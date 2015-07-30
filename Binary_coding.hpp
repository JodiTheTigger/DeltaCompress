// Copyright 2015 Richard Maxwell, all rights reserved.

// Interesting reading:
// http://www.arturocampos.com/ac_range.html
// http://cbloomrants.blogspot.co.nz/2010/09/09-16-10-modern-arithmetic-coding-from.html
//
// But you really should be spending all your time here:
// http://mattmahoney.net/dc/text.html
#pragma once

#include <vector>
#include <array>
#include <cstdint>
#include <cassert>
#include <cmath>

// //////////////////////////////////////////////
// Types
// //////////////////////////////////////////////

namespace Binary_types
{
    // Maths means for 32 bits, max range is 15. Can't remember why
    // Maybe test to see if I can get away with > 15?
    // Maybe that's only for arth coders? and not range coders?
    static const uint32_t PROBABILITY_RANGE_BITS = 13;
    static const uint32_t PROBABILITY_RANGE      = 1 << PROBABILITY_RANGE_BITS;

    using Bytes = std::vector<uint8_t>;
}

// //////////////////////////////////////////////
// Coders
// //////////////////////////////////////////////

namespace Binary_coders
{
    using namespace Binary_types;

    // RAM: TODO: clean up this class to share common things.
    class fpaq0p_encoder
    {
    public:
        fpaq0p_encoder(Bytes& bytes)
            : m_bytes(&bytes)
        {}

        void Encode(unsigned bit, uint16_t one_probability)
        {
            assert(one_probability < (1 << PROBABILITY_RANGE_BITS));

            // cast to 64 bits for more precision.
            const uint32_t middle =
                m_low
                +
                (
                    (
                        static_cast<uint64_t>(m_high - m_low)
                        * one_probability
                    )
                    >> PROBABILITY_RANGE_BITS
                );

            if (bit)
            {
                m_high = middle;
            }
            else
            {
                // +1 has something to do with squashing carrry bits and
                // underflow. I'm fuzzy with the maths.
                m_low = middle + 1;
            }

            // flush while the top byte is all set for both low and high
            // i.e. normalise.
            while ((m_low ^ m_high) < (1u << 24))
            {
                m_bytes->push_back(m_low >> 24);
                m_low <<= 8;

                // Some people ignore the 0xff at the end, can't remember
                // why it was there. *shrug*
                m_high = (m_high << 8) | 0xff;
            }
        }

        ~fpaq0p_encoder()
        {
            // Just copy Fabian's flush code.

            // Find shortest encoding that still decodes to the
            // right symbols. This is assume the decoder implicitly
            // zero-pads out of bounds array accesses.
            uint32_t round_up = 0xffffffu;
            while (round_up)
            {
                if ((m_low | round_up) != ~0u)
                {
                    uint32_t rounded = (m_low + round_up) & ~round_up;

                    if (rounded <= m_high)
                    {
                        // inside interval, we're good!
                        m_low = rounded;
                        break;
                    }
                }

                round_up >>= 8;
            }

            while (m_low)
            {
                m_bytes->push_back(m_low >> 24);
                m_low <<= 8;
            }
        }

    private:
        Bytes*      m_bytes = 0;
        uint32_t    m_low   = 0;
        uint32_t    m_high  = 0xffffffff;
    };

    class fpaq0p_decoder
    {
    public:
        fpaq0p_decoder(const Bytes& bytes)
            : m_bytes(&bytes)
            , m_byteSize(bytes.size())
        {
            m_value = (m_value << 8) + Read();
            m_value = (m_value << 8) + Read();
            m_value = (m_value << 8) + Read();
            m_value = (m_value << 8) + Read();
        }

        unsigned Decode(uint16_t one_probability)
        {
            assert(one_probability < (1 << PROBABILITY_RANGE_BITS));

            const uint32_t middle =
                m_low
                +
                (
                    (
                        static_cast<uint64_t>(m_high - m_low)
                        * one_probability
                    )
                    >> PROBABILITY_RANGE_BITS
                );

            auto bit = 0u;

            if (m_value <= middle)
            {
                bit = 1;
                m_high = middle;
            }
            else
            {
                m_low = middle + 1;
            }

            // Renormaise
            while ((m_low ^ m_high) < (1u << 24))
            {
                m_low <<= 8;
                m_high = (m_high << 8) | 0xff;

                m_value <<= 8;
                m_value += Read();
            }

            return bit;
        }

        uint32_t FlushAndGetBytesRead()
        {
            // Do I need to renormalise?
            return m_read_index;
        }

    private:
        const Bytes*    m_bytes         = 0;
        uint32_t        m_value         = 0;
        uint32_t        m_low           = 0;
        uint32_t        m_high          = 0xffffffff;
        uint32_t        m_read_index    = 0;
        uint32_t        m_byteSize      = 0;

        uint8_t Read()
        {
            if (m_read_index < m_byteSize)
            {
                return (*m_bytes)[m_read_index++];
            }

            return 0;
        }
    };

}

// //////////////////////////////////////////////
// Models
// //////////////////////////////////////////////

namespace Models
{
    using namespace Binary_types;
    using namespace Binary_coders;

    // Ideas from https://github.com/rygorous/gaffer_net/blob/master/main.cpp    

    template<typename ENCODER, typename DECODER>
    class Binary_two_speed
    {
    public:
        typedef ENCODER Encoder;
        typedef DECODER Decoder;

        Binary_two_speed(
                unsigned inertia_1 = 4,
                unsigned inertia_2 = 4,
                unsigned initial_probability_1 = QUARTER_RANGE,
                unsigned initial_probability_2 = QUARTER_RANGE)
            : m_inertia_1(inertia_1)
            , m_inertia_2(inertia_2)
            , m_probabilities_1(initial_probability_1)
            , m_probabilities_2(initial_probability_2)
        {}

        void Encode(ENCODER& coder, unsigned value)
        {
            coder.Encode(value, m_probabilities_1 + m_probabilities_2);
            Adapt(value);
        }

        unsigned Decode(DECODER& coder)
        {
            auto result = coder.Decode(m_probabilities_1 + m_probabilities_2);
            Adapt(result);

            return result;
        }

    private:
        static const unsigned HALF_RANGE    = PROBABILITY_RANGE / 2;
        static const unsigned QUARTER_RANGE = HALF_RANGE / 2;

        unsigned m_inertia_1;
        unsigned m_inertia_2;
        uint16_t m_probabilities_1;
        uint16_t m_probabilities_2;

        void Adapt(unsigned value)
        {
            if (value)
            {
                m_probabilities_1 += (HALF_RANGE - m_probabilities_1) >> m_inertia_1;
                m_probabilities_2 += (HALF_RANGE - m_probabilities_2) >> m_inertia_2;
            }
            else
            {
                m_probabilities_1 -= m_probabilities_1 >> m_inertia_1;
                m_probabilities_2 -= m_probabilities_2 >> m_inertia_2;
            }
        }
    };

    template<typename ENCODER, typename DECODER>
    class Binary_bitstream_general
    {
    public:
        static const constexpr unsigned HALF_RANGE = PROBABILITY_RANGE / 2;

        static void Encode(ENCODER& coder, unsigned value)
        {
            coder.Encode(value, HALF_RANGE);
        }

        static unsigned Decode(DECODER& coder)
        {
            auto result = coder.Decode(HALF_RANGE);
            return result;
        }
    };

    template<class BINARY_MODEL, unsigned BITS>
    class Binary_tree
    {
    public:
        Binary_tree() = default;

        // This is getting too c++'y. find an easier way..
        template<typename... Args>
        Binary_tree(Args&&... args)
        {
            for(auto& model : m_models)
            {
                model = BINARY_MODEL(std::forward<Args>(args)...);
            }
        }

        void Encode(typename BINARY_MODEL::Encoder& coder, unsigned value, unsigned bits = BITS)
        {
            assert(value < MODEL_COUNT);
            assert(bits <= BITS);

            // Model the MSB first, then work our way down.
            // Seems adds are better than << 1.
            unsigned rebuilt = 1;

            unsigned skips = BITS - bits;
            value <<= skips;
            rebuilt <<= skips;

            while (rebuilt < MODEL_COUNT)
            {
                unsigned bit = ((value & TOP_BIT) != 0);
                value += value;
                m_models[rebuilt - 1].Encode(coder, bit);
                rebuilt += rebuilt + bit;
            }
        }

        unsigned Decode(typename BINARY_MODEL::Decoder& coder, unsigned bits = BITS)
        {
            assert(bits <= BITS);

            unsigned rebuilt = 1 << (BITS - bits);
            unsigned count = MODEL_COUNT;

            while (rebuilt < count)
            {
                rebuilt += rebuilt + m_models[rebuilt - 1].Decode(coder);
            }

            // Clear the top bit due to starting rebuilt with 1.
            return rebuilt - count;
        }

    private:
        static const unsigned MODEL_COUNT   = 1 << BITS;
        static const unsigned TOP_BIT       = MODEL_COUNT / 2;

        std::array<BINARY_MODEL, MODEL_COUNT - 1> m_models;
    };

    template<class BINARY_MODEL, unsigned BITS, unsigned MAX_VALUE>
    class Binary_tree_max_value
    {
    public:
        Binary_tree_max_value() = default;

        // This is getting too c++'y. find an easier way..
        template<typename... Args>
        Binary_tree_max_value(Args&&... args)
        {
            for(auto& model : m_models)
            {
                model = BINARY_MODEL(std::forward<Args>(args)...);
            }
        }

        void Encode(typename BINARY_MODEL::Encoder& coder, unsigned value)
        {
            assert(value < MODEL_COUNT);

            // Model the MSB first, then work our way down.
            // Seems adds are better than << 1.
            unsigned rebuilt = 1;
            unsigned mask = MAX_VALUE;

            while (rebuilt < MODEL_COUNT)
            {
                unsigned bit = ((value & TOP_BIT) != 0);
                value += value;

                if (mask & TOP_BIT)
                {
                    m_models[rebuilt - 1].Encode(coder, bit);

                    // At any point we get a zero, then that means
                    // we are not constrained by max value anymore.
                    // RAM: TODO: algorithm is clunky - fix.
                    if (!bit)
                    {
                        mask = MODEL_COUNT - 1;
                    }
                }

                rebuilt += rebuilt + bit;
                mask += mask;
            }
        }

        unsigned Decode(typename BINARY_MODEL::Decoder& coder)
        {
            unsigned rebuilt = 1;
            unsigned count = MODEL_COUNT;
            unsigned mask = MAX_VALUE;

            while (rebuilt < count)
            {
                auto new_bit = 0u;

                if (mask & TOP_BIT)
                {
                    new_bit = m_models[rebuilt - 1].Decode(coder);

                    if (!new_bit)
                    {
                        mask = MODEL_COUNT - 1;
                    }
                }

                rebuilt += rebuilt + new_bit;
                mask += mask;
            }

            // Clear the top bit due to starting rebuilt with 1.
            return rebuilt - count;
        }

    private:
        static const unsigned MODEL_COUNT   = 1 << BITS;
        static const unsigned TOP_BIT       = MODEL_COUNT / 2;

        std::array<BINARY_MODEL, MODEL_COUNT - 1> m_models;
    };

    // RAM: TODO: Pre init the probabilities for values
    // larger than some MAX value to 0.
    template
    <
        class BINARY_MODEL,
        unsigned BITS_FOR_BITS,
        unsigned MAX_VALUE=((1 << BITS_FOR_BITS) - 1)
    >
    class Unsigned_golomb_binary
    {
    public:
        Unsigned_golomb_binary() = default;

        // This is getting too c++'y. find an easier way..
        template<typename... Args>
        Unsigned_golomb_binary(Args&&... args)
        {
            m_bits =
                Binary_tree_max_value<BINARY_MODEL, BITS_FOR_BITS, MAX_VALUE>
                (
                    std::forward<Args>(args)...
                );
        }

        void Encode(typename BINARY_MODEL::Encoder& coder, unsigned value)
        {
            // Enocde how many bits we are sending
            unsigned min_bits = 0;

            while ((1u << min_bits) <= value)
            {
                ++min_bits;
            }

            assert(min_bits <= MAX_VALUE);

            {
                unsigned min_bits_2 = 0;

                while ((1u << min_bits_2) <= min_bits)
                {
                    ++min_bits_2;
                }

                assert(min_bits_2 <= BITS_FOR_BITS);
            }

            m_bits.Encode(coder, min_bits);

            if (min_bits)
            {
                // No need to send the top bit
                // as we can assume it's set due
                // to min_bits.
                --min_bits;

                // Send the bits, no probabilities.
                while (min_bits--)
                {
                    const auto mask = (1 << min_bits);

                    Binary_bitstream_general
                    <
                        typename BINARY_MODEL::Encoder,
                        typename BINARY_MODEL::Decoder
                    >
                    ::Encode
                    (
                        coder,
                        value & mask
                    );
                }
            }
        }

        unsigned Decode(typename BINARY_MODEL::Decoder& coder)
        {
            auto min_bits = m_bits.Decode(coder);

            unsigned result = 0;

            if (min_bits)
            {
                // top bit is always set.
                result |= 1;

                --min_bits;

                while (min_bits)
                {
                    result <<= 1;

                    result |=
                        Binary_bitstream_general
                        <
                            typename BINARY_MODEL::Encoder,
                            typename BINARY_MODEL::Decoder
                        >
                        ::Decode(coder);

                    --min_bits;
                }
            }

            return result;
        }

    private:
        Binary_tree_max_value<BINARY_MODEL, BITS_FOR_BITS, MAX_VALUE> m_bits;
    };

}

void binary_tests()
{
    using namespace Binary_types;
    using namespace Binary_coders;
    using namespace Models;

    {
        Bytes data;

        auto tests = {0,1,0,1,0,1,1,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0};

        {
            fpaq0p_encoder encoder(data);

            for (const unsigned t : tests)
            {
                encoder.Encode(t, (PROBABILITY_RANGE * 2) / 3);
            }
        }

        {
            fpaq0p_decoder decoder(data);

            for (const unsigned t : tests)
            {
                auto value = decoder.Decode((PROBABILITY_RANGE * 2) / 3);
                assert(value == t);
            }

            auto read = decoder.FlushAndGetBytesRead();
            assert(read == data.size());
        }
    }

    // Binary Model Tests
    {
        auto tests = {0,1,0,1,0,1,1,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0};

        auto Binary_test = [](auto tests, auto model_in, auto model_out)
        {
            Bytes data;

            {
                fpaq0p_encoder test_encoder(data);

                for (auto t : tests)
                {
                    model_in.Encode(test_encoder, t);
                }
            }
            {
                fpaq0p_decoder test_decoder(data);

                for (unsigned t : tests)
                {
                    auto value = model_out.Decode(test_decoder);
                    assert(value == t);
                }

                auto read = test_decoder.FlushAndGetBytesRead();
                assert(read == data.size());
            }
        };

        Binary_test
        (
            tests,
            Binary_bitstream_general<fpaq0p_encoder, fpaq0p_decoder>(),
            Binary_bitstream_general<fpaq0p_encoder, fpaq0p_decoder>()
        );

        Binary_test(
            tests,
            Binary_two_speed<fpaq0p_encoder, fpaq0p_decoder>(5,2),
            Binary_two_speed<fpaq0p_encoder, fpaq0p_decoder>(5,2));
        Binary_test(
            tests,
            Binary_two_speed<fpaq0p_encoder, fpaq0p_decoder>(6,1),
            Binary_two_speed<fpaq0p_encoder, fpaq0p_decoder>(6,1));
        Binary_test(
            tests,
            Binary_two_speed<fpaq0p_encoder, fpaq0p_decoder>(3,4),
            Binary_two_speed<fpaq0p_encoder, fpaq0p_decoder>(3,4));
    }

    // Binary tree tests.
    {
        Bytes data;
        Bytes tests{0,1,2,3,4,5,6,7, 0,1,2,3,4,5,6,7, 0,1,2,3,4,5,6,7, 0,1,2};

        {
            fpaq0p_encoder test_encoder(data);
            Binary_tree_max_value
            <
                Binary_two_speed<fpaq0p_encoder, fpaq0p_decoder>,
                3,
                5
            > tree;

            for (auto t : tests)
            {
                tree.Encode(test_encoder, t % 6);
            }
        }
        {
            fpaq0p_decoder test_decoder(data);
            Binary_tree_max_value
            <
                Binary_two_speed<fpaq0p_encoder, fpaq0p_decoder>,
                3,
                5
            > tree;

            for (unsigned t : tests)
            {
                auto value = tree.Decode(test_decoder);
                assert(value == (t % 6));
            }

            auto read = test_decoder.FlushAndGetBytesRead();
            assert(read == data.size());
        }
    }

    // Binary coder.
    {
        Bytes data;

        auto tests = {0,1,0,1,0,1,1,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0};

        {
            fpaq0p_encoder encoder(data);

            for (const unsigned t : tests)
            {
                encoder.Encode(t, PROBABILITY_RANGE / 10);
            }
        }

        {
            fpaq0p_decoder decoder(data);

            for (const unsigned t : tests)
            {
                auto value =
                    decoder.Decode(PROBABILITY_RANGE / 10);

                assert(value == t);
            }

            auto read = decoder.FlushAndGetBytesRead();
            assert(read == data.size());
        }
    }
}
