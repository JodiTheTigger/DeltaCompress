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

namespace Range_types
{
    static const uint32_t PROBABILITY_RANGE_BITS = 15;
    static const uint32_t PROBABILITY_RANGE      = 1 << PROBABILITY_RANGE_BITS;

    using Bytes = std::vector<uint8_t>;

    struct Range
    {
        uint32_t min;
        uint32_t count;
    };
}

// //////////////////////////////////////////////
// Coders
// //////////////////////////////////////////////

// NOTE: don't encode more than 2 << 24 values or you will
// have a bad time.
namespace Range_coders
{
    using namespace Range_types;

    // Uses 31 instead of 32 bits to reserve space for a possible carry.
    static const uint32_t CODE_BITS     = 32;
    static const uint32_t TOP_VALUE     = (1ul << (CODE_BITS - 1));
    static const uint32_t SHIFT_BITS    = (CODE_BITS - 9);
    static const uint32_t EXTRA_BITS    = ((CODE_BITS - 2) % 8 + 1);
    static const uint32_t BOTTOM_VALUE  = (TOP_VALUE >> 8);

    class Encoder
    {
    public:
        Encoder(Bytes& bytes)
            : m_bytes(&bytes)
        {}

        ~Encoder()
        {
            Renormalise();

            // If we still haven't converged to a top bit, we need to
            // force it to converge and flush the underflows.
            auto temp = m_min >> SHIFT_BITS;

            // Still NFI why the comparison with bytes written is made.
            if  (
                    (m_min & (BOTTOM_VALUE - 1)) >=
                    ((m_bytes->size() & 0xffffffL) >> 1)
                )
            {
                ++temp;
            }

            if (temp > 0xFF)
            {
                Flush_underflow_with_carry(1);
            }
            else
            {
                Flush_underflow_with_carry(0);
            }

            // Ok, so now we should be in the state with no underflows left.
            // Change the range so that we write the minimum amount of tail
            // bytes.
            auto min        = m_min;
            auto max        = m_min + m_range;
            auto round_up   = BOTTOM_VALUE - 1;

            while (round_up)
            {
                auto rounded = (min + round_up) & ~round_up;
                if (rounded <= max)
                {
                    min = rounded;
                    break;
                }

                round_up >>= 8;
            }

            while (min)
            {
                auto to_write = (min >> SHIFT_BITS) & 0xFF;
                Write(to_write);

                min <<= 8;
                min &= (TOP_VALUE - 1);
            }
        }

        void Encode(Range range)
        {
            Renormalise();

            auto new_range          = m_range >> PROBABILITY_RANGE_BITS;
            auto new_range_start    = range.min * new_range;

            m_min += new_range_start;

            if (range.min + range.count < PROBABILITY_RANGE)
            {
                m_range = new_range * range.count;
            }
            else
            {
                m_range -= new_range_start;
            }
        }

    private:
        Bytes*      m_bytes                 = 0;
        uint32_t    m_range                 = TOP_VALUE;
        uint32_t    m_min                   = 0;
        uint8_t     m_underflow_byte_count  = 0;
        uint8_t     m_buffer                = 0;
        bool        m_first                 = true;

        void Write(uint8_t value)
        {
            if (!m_first)
            {
                m_bytes->push_back(value);
            }
            else
            {
                m_first = false;
            }
        }

        void Flush_underflow_with_carry(uint8_t carry)
        {
            // Silly way of avoiding a branch
            const uint8_t buffered_bits = 0xFF * (1 - carry);

            Write(m_buffer + carry);

            auto underfloat_count = m_underflow_byte_count;
            while (underfloat_count)
            {
                Write(buffered_bits);
                --underfloat_count;
            }

            m_buffer                = m_min >> SHIFT_BITS;
            m_underflow_byte_count  = 0;
        }

        void Renormalise()
        {
            while (m_range <= BOTTOM_VALUE)
            {
                if (m_min < (0xFF << SHIFT_BITS))
                {
                    Flush_underflow_with_carry(0);
                }
                else
                {
                    if (m_min & TOP_VALUE)
                    {
                        Flush_underflow_with_carry(1);
                    }
                    else
                    {
                        ++m_underflow_byte_count;
                    }
                }

                m_range    <<= 8;
                m_min      <<= 8;
                m_min      &= (TOP_VALUE - 1);
            }
        }
    };

    class Decoder
    {
    public:
        Decoder(const Bytes& bytes)
            : m_bytes(&bytes)
            , m_byteSize(bytes.size())
        {
            m_buffer    = Read();
            m_min       = m_buffer >> (8 - EXTRA_BITS);
        }

        uint32_t Decode()
        {
            Renormalise();

            m_next_range = m_range >> PROBABILITY_RANGE_BITS;
            auto symbol_range = m_min / m_next_range;

            if (symbol_range >= PROBABILITY_RANGE)
            {
                return PROBABILITY_RANGE - 1;
            }
            else
            {
                return symbol_range;
            }
        }

        void Update(Range range)
        {
            auto temp = m_next_range * range.min;
            m_min -= temp;
            if ((range.min + range.count) < PROBABILITY_RANGE)
            {
                m_range = m_next_range * range.count;
            }
            else
            {
                m_range -= temp;
            }
        }

        uint32_t FlushAndGetBytesRead()
        {
            Renormalise();
            return m_read_index;
        }

    private:
        const Bytes*    m_bytes;
        uint32_t        m_range         = 1 << EXTRA_BITS;
        uint32_t        m_min           = 0;
        uint32_t        m_next_range    = 0;
        uint32_t        m_read_index    = 0;
        uint32_t        m_byteSize      = 0;
        uint8_t         m_buffer        = 0;

        uint8_t Read()
        {
            if (m_read_index < m_byteSize)
            {
                return (*m_bytes)[m_read_index++];
            }

            return 0;
        }

        void Renormalise()
        {
            while (m_range <= BOTTOM_VALUE)
            {
                m_min =
                    (m_min << 8) |
                    ((m_buffer << EXTRA_BITS) & 0xFF);

                m_buffer    =   Read();
                m_min       |=  (m_buffer >> (8 - EXTRA_BITS));
                m_range     <<= 8;
            }
        }
    };

    class Binary_encoder
    {
    public:
        Binary_encoder(Encoder& coder)
            : m_encoder(&coder)
        {}

        void Encode(unsigned value, uint16_t one_probability)
        {
            // Note: First range is for one_probability, not zero.
            m_encoder->Encode(
                value ?
                    Range{0, one_probability} :
                    Range{one_probability, PROBABILITY_RANGE - one_probability});
        }

    private:
        Encoder* m_encoder;
    };

    class Binary_decoder
    {
    public:
        Binary_decoder(Decoder& coder)
            : m_decoder(&coder)
        {}

        unsigned Decode(unsigned one_probability)
        {
            auto symbol = m_decoder->Decode();
            auto result = (symbol < one_probability);

            m_decoder->Update(
                result ?
                    Range{0, one_probability} :
                    Range{one_probability, PROBABILITY_RANGE - one_probability});

            return result;
        }

        uint32_t FlushAndGetBytesRead()
        {
            return m_decoder->FlushAndGetBytesRead();
        }

    private:
        Decoder* m_decoder;
    };

    class fpaq0p_encoder
    {
    public:
        // Maths means for 32 bits, max range is 15. Can't remember why
        // Maybe test to see if I can get away with > 15?
        const unsigned PROBABILITY_BITS = 15;

        fpaq0p_encoder(Bytes& bytes)
            : m_bytes(&bytes)
        {}

        void Encode(unsigned bit, uint16_t one_probability)
        {
            assert(one_probability < (1 << PROBABILITY_BITS));

            // cast to 64 bits for more precision.
            const uint32_t middle =
                m_low
                +
                (
                    (
                        static_cast<uint64_t>(m_high - m_low)
                        * one_probability
                    )
                    >> PROBABILITY_BITS
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
        // Maths means for 32 bits, max range is 15. Can't remember why
        // Maybe test to see if I can get away with > 15?
        const unsigned PROBABILITY_BITS = 15;

        fpaq0p_decoder(Bytes& bytes)
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
            assert(one_probability < (1 << PROBABILITY_BITS));

            const uint32_t middle =
                m_low
                +
                (
                    (
                        static_cast<uint64_t>(m_high - m_low)
                        * one_probability
                    )
                    >> PROBABILITY_BITS
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
        Bytes*      m_bytes         = 0;
        uint32_t    m_value         = 0;
        uint32_t    m_low           = 0;
        uint32_t    m_high          = 0xffffffff;
        uint32_t    m_read_index    = 0;
        uint32_t    m_byteSize      = 0;

        uint8_t Read()
        {
            if (m_read_index < m_byteSize)
            {
                return (*m_bytes)[m_read_index++];
            }

            return 0;
        }
    };

} // namespace Range_coders

// //////////////////////////////////////////////
// Models
// //////////////////////////////////////////////

namespace Range_models
{
    using namespace Range_types;
    using namespace Range_coders;

    // Ideas from https://github.com/rygorous/gaffer_net/blob/master/main.cpp

    class Binary
    {
    public:
        Binary(
                unsigned inertia = 4,
                unsigned initial_probability = (PROBABILITY_RANGE / 2))
            : m_inertia(inertia)
            , m_one_probability(initial_probability)
        {
            assert(inertia > 0);
        }

        void Encode(Binary_encoder& coder, unsigned value)
        {
            coder.Encode(value, m_one_probability);
            Adapt(value);
        }

        unsigned Decode(Binary_decoder& coder)
        {
            auto result = coder.Decode(m_one_probability);
            Adapt(result);

            return result;
        }

    private:
        uint16_t m_inertia;
        uint16_t m_one_probability;

        void Adapt(unsigned value)
        {
            if (value)
            {
                m_one_probability +=
                    (PROBABILITY_RANGE - m_one_probability) >> m_inertia;
            }
            else
            {
                m_one_probability -= m_one_probability >> m_inertia;
            }
        }
    };

    class Binary_two_speed
    {
    public:
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

        void Encode(Binary_encoder& coder, unsigned value)
        {
            coder.Encode(value, m_probabilities_1 + m_probabilities_2);
            Adapt(value);
        }

        unsigned Decode(Binary_decoder& coder)
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

    class Binary_bitstream
    {
    public:
        static const constexpr unsigned HALF_RANGE = PROBABILITY_RANGE / 2;

        static void Encode(Binary_encoder& coder, unsigned value)
        {
            coder.Encode(value, HALF_RANGE);
        }

        static void Encode(Encoder& coder, unsigned value)
        {
            coder.Encode
            (
                value ?
                    Range{0, HALF_RANGE} :
                    Range{HALF_RANGE, HALF_RANGE}
            );
        }

        static unsigned Decode(Binary_decoder& coder)
        {
            auto result = coder.Decode(HALF_RANGE);
            return result;
        }

        static unsigned Decode(Decoder& coder)
        {
            auto symbol = coder.Decode();
            auto result = (symbol < HALF_RANGE);

            coder.Update
            (
                result ?
                    Range{0, HALF_RANGE} :
                    Range{HALF_RANGE, HALF_RANGE}
            );

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

        void Encode(Binary_encoder& coder, unsigned value, unsigned bits = BITS)
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

        unsigned Decode(Binary_decoder& coder, unsigned bits = BITS)
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

        void Encode(Binary_encoder& coder, unsigned value)
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

        unsigned Decode(Binary_decoder& coder)
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

        void Encode(Binary_encoder& coder, unsigned value)
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

                    Binary_bitstream::Encode
                    (
                        coder,
                        value & mask
                    );
                }
            }
        }

        unsigned Decode(Binary_decoder& coder)
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

                    result |= Binary_bitstream::Decode(coder);

                    --min_bits;
                }
            }

            return result;
        }

    private:
        Binary_tree_max_value<BINARY_MODEL, BITS_FOR_BITS, MAX_VALUE> m_bits;
    };

    class Periodic_update
    {
    public:
        typedef std::vector<unsigned>  Freqencies;
        typedef std::vector<Range>     Ranges;

        Periodic_update
        (
            unsigned size,
            unsigned slowest_update_rate
        )
            : m_f()
            , m_r(size)
            , m_size(size)
            , m_slowest_update_rate(slowest_update_rate)
        {
            m_f.reserve(size);

            for (unsigned i = 0; i < size; ++i)
            {
                m_f.push_back(1);
                ++m_last_total;
            }

            Recalculate_ranges();
        }

        Periodic_update
        (
            Freqencies frequencies,
            unsigned slowest_update_rate
        )
            : m_f(frequencies)
            , m_r(frequencies.size())
            , m_size(frequencies.size())
            , m_slowest_update_rate(slowest_update_rate)
        {
            for (const auto f : m_f)
            {
                assert(f);

                m_last_total += f;
            }

            Recalculate_ranges();
        }

        void Encode(Encoder& coder, unsigned value)
        {
            assert(value < m_size);

            coder.Encode(m_r[value]);

            Adapt(value);
        }

        void Encode(Encoder& coder, unsigned value, unsigned max_value)
        {
            assert(value < m_size);
            assert(max_value < m_size);

            // Hmm, encoding one bit bad give me issues.
            max_value = std::max(max_value, 2u);

            // Multiply the ranges instead of the value
            // so the frequencies still get updated correctly.
            auto old_range  = m_r[value];
            auto last_range = m_r[max_value];
            float max_cf = last_range.min + last_range.count;
            float multiplier = PROBABILITY_RANGE / max_cf;

            // Man, rounding errors and off my one are my bane.
            unsigned next_min = 0;
            if (value < max_value)
            {
                next_min =
                    static_cast<uint32_t>
                    (
                        std::round(m_r[value + 1].min * multiplier)
                    );
            }
            else
            {
                next_min = PROBABILITY_RANGE;
            }

            auto new_min =
                static_cast<uint32_t>(std::round(old_range.min * multiplier));

            auto new_range = Range
            {
                new_min,
                next_min - new_min,
            };

            coder.Encode(new_range);

            Adapt(value);
        }

        unsigned Decode(Decoder& coder)
        {
            auto range = coder.Decode();
            unsigned result = 0;

            const auto size = m_size;
            while ((m_r[result].min <= range) && (result < size))
            {
                ++result;
            }

            --result;

            coder.Update(m_r[result]);
            Adapt(result);

            return result;
        }

        unsigned Decode(Decoder& coder, unsigned max_value)
        {
            assert(max_value < m_size);

            // Hmm, encoding one bit bad give me issues.
            max_value = std::max(max_value, 2u);

            auto range = coder.Decode();
            unsigned result = 0;

            auto last_range = m_r[max_value];
            float max_cf = last_range.min + last_range.count;
            float multiplier = PROBABILITY_RANGE / max_cf;

            const auto size = m_size;

            while (result < size)
            {
                auto new_min =
                    static_cast<uint32_t>
                    (
                        std::round(m_r[result].min * multiplier)
                    );

                if (new_min > range)
                {
                    break;
                }

                ++result;
            }

            --result;

            unsigned next_min = 0;
            if (result < max_value)
            {
                next_min =
                    static_cast<uint32_t>
                    (
                        std::round(m_r[result + 1].min * multiplier)
                    );
            }
            else
            {
                next_min = PROBABILITY_RANGE;
            }

            auto new_min =
                static_cast<uint32_t>(std::round(m_r[result].min * multiplier));

            auto new_range = Range
            {
                new_min,
                next_min - new_min,
            };

            coder.Update(new_range);
            Adapt(result);

            return result;
        }

    private:
        Freqencies  m_f;
        Ranges      m_r;

        unsigned m_updates          = 0;
        unsigned m_update_trigger   = 1;
        unsigned m_last_total       = 0;

        const unsigned m_size                   = 0;
        const unsigned m_slowest_update_rate    = 0;

        void Adapt(unsigned value)
        {
            m_f[value] += 2;
            m_updates += 2;

            if (m_updates >= m_update_trigger)
            {
                m_update_trigger += m_update_trigger;
                if (m_update_trigger > m_slowest_update_rate)
                {
                    m_update_trigger = m_slowest_update_rate;
                }

                Recalculate_ranges();
                m_updates = 0;
            }
        }

        void Recalculate_ranges()
        {
            unsigned total = m_last_total + m_updates;

            if (total > PROBABILITY_RANGE)
            {
                // Dunno how to do this nice for now. Brute it.
                total = 0;
                for (auto& f : m_f)
                {
                    f >>= 1;
                    ++f;

                    total += f;
                }
            }

            const auto size     = m_size;
            auto multiple       = PROBABILITY_RANGE / total;
            auto reminder       = PROBABILITY_RANGE % total;
            auto global_adjust  = reminder / size;
            auto reminder_count = reminder % size;
            unsigned last_min   = 0;

            for (unsigned i = 0; i < reminder_count; ++i)
            {
                m_r[i].min      = last_min;
                m_r[i].count    = 1 + (global_adjust + m_f[i] * multiple);

                last_min += m_r[i].count;
            }

            for (unsigned i = reminder_count; i < size; ++i)
            {
                m_r[i].min      = last_min;
                m_r[i].count    = 0 + (global_adjust + m_f[i] * multiple);

                last_min += m_r[i].count;
            }

            assert(last_min == PROBABILITY_RANGE);

            m_last_total = total;
        }
    };

    class No_update
    {
    public:
        typedef std::vector<unsigned>  Freqencies;
        typedef std::vector<Range>     Ranges;

        No_update
        (
            unsigned size
        )
            : m_r(size)
            , m_size(size)
        {
            Freqencies f;
            f.reserve(size);

            unsigned total = 0;
            for (unsigned i = 0; i < size; ++i)
            {
                f.push_back(1);
                ++total;
            }

            Recalculate_ranges(f, total);
        }

        No_update
        (
            Freqencies frequencies
        )
            : m_r(frequencies.size())
            , m_size(frequencies.size())
        {
            unsigned total = 0;
            for (const auto f : frequencies)
            {
                assert(f);

                total += f;
            }

            Recalculate_ranges(frequencies, total);
        }

        void Encode(Encoder& coder, unsigned value)
        {
            assert(value < m_size);

            coder.Encode(m_r[value]);
        }

        void Encode(Encoder& coder, unsigned value, unsigned max_value)
        {
            assert(value < m_size);
            assert(max_value < m_size);

            // Hmm, encoding one bit bad give me issues.
            max_value = std::max(max_value, 2u);

            // Multiply the ranges instead of the value
            // so the frequencies still get updated correctly.
            auto old_range  = m_r[value];
            auto last_range = m_r[max_value];
            float max_cf = last_range.min + last_range.count;
            float multiplier = PROBABILITY_RANGE / max_cf;

            // Man, rounding errors and off my one are my bane.
            unsigned next_min = 0;
            if (value < max_value)
            {
                next_min =
                    static_cast<uint32_t>
                    (
                        std::round(m_r[value + 1].min * multiplier)
                    );
            }
            else
            {
                next_min = PROBABILITY_RANGE;
            }

            auto new_min =
                static_cast<uint32_t>(std::round(old_range.min * multiplier));

            auto new_range = Range
            {
                new_min,
                next_min - new_min,
            };

            coder.Encode(new_range);
        }

        unsigned Decode(Decoder& coder)
        {
            auto range = coder.Decode();
            unsigned result = 0;

            const auto size = m_size;
            while ((m_r[result].min <= range) && (result < size))
            {
                ++result;
            }

            --result;

            coder.Update(m_r[result]);

            return result;
        }

        unsigned Decode(Decoder& coder, unsigned max_value)
        {
            assert(max_value < m_size);

            // Hmm, encoding one bit bad give me issues.
            max_value = std::max(max_value, 2u);

            auto range = coder.Decode();
            unsigned result = 0;

            auto last_range = m_r[max_value];
            float max_cf = last_range.min + last_range.count;
            float multiplier = PROBABILITY_RANGE / max_cf;

            const auto size = m_size;

            while (result < size)
            {
                auto new_min =
                    static_cast<uint32_t>
                    (
                        std::round(m_r[result].min * multiplier)
                    );

                if (new_min > range)
                {
                    break;
                }

                ++result;
            }

            --result;

            unsigned next_min = 0;
            if (result < max_value)
            {
                next_min =
                    static_cast<uint32_t>
                    (
                        std::round(m_r[result + 1].min * multiplier)
                    );
            }
            else
            {
                next_min = PROBABILITY_RANGE;
            }

            auto new_min =
                static_cast<uint32_t>(std::round(m_r[result].min * multiplier));

            auto new_range = Range
            {
                new_min,
                next_min - new_min,
            };

            coder.Update(new_range);

            return result;
        }

    private:
        Ranges      m_r;

        const unsigned m_size = 0;

        void Recalculate_ranges(Freqencies frequencies, unsigned total)
        {
            if (total > PROBABILITY_RANGE)
            {
                // Dunno how to do this nice for now. Brute it.
                total = 0;
                for (auto& f : frequencies)
                {
                    f >>= 1;
                    ++f;

                    total += f;
                }
            }

            const auto size     = m_size;
            auto multiple       = PROBABILITY_RANGE / total;
            auto reminder       = PROBABILITY_RANGE % total;
            auto global_adjust  = reminder / size;
            auto reminder_count = reminder % size;
            unsigned last_min   = 0;

            for (unsigned i = 0; i < reminder_count; ++i)
            {
                m_r[i].min   = last_min;
                m_r[i].count = 1 + (global_adjust + frequencies[i] * multiple);

                last_min += m_r[i].count;
            }

            for (unsigned i = reminder_count; i < size; ++i)
            {
                m_r[i].min   = last_min;
                m_r[i].count = 0 + (global_adjust + frequencies[i] * multiple);

                last_min += m_r[i].count;
            }

            assert(last_min == PROBABILITY_RANGE);
        }
    };

    class Periodic_update_with_kernel
    {
    public:
        typedef std::vector<unsigned>  Freqencies;
        typedef std::vector<Range>     Ranges;
        typedef std::vector<unsigned>  Kernel;

        Periodic_update_with_kernel
        (
            unsigned size,
            unsigned slowest_update_rate,
            Kernel kernel = {1, 1, 2, 8, 2, 1, 1}
        )
            : m_f()
            , m_r(size)
            , m_k(kernel)
            , m_size(size)
            , m_slowest_update_rate(slowest_update_rate)
            , m_kernal_offset(kernel.size() >> 1)
        {
            m_f.reserve(size);

            for (unsigned i = 0; i < size; ++i)
            {
                m_f.push_back(1);
                ++m_last_total;
            }

            Recalculate_ranges();
        }

        Periodic_update_with_kernel
        (
            Freqencies frequencies,
            unsigned slowest_update_rate,
            Kernel kernel = {1, 1, 2, 8, 2, 1, 1}
        )
            : m_f(frequencies)
            , m_r(frequencies.size())
            , m_k(kernel)
            , m_size(frequencies.size())
            , m_slowest_update_rate(slowest_update_rate)
            , m_kernal_offset(kernel.size() >> 1)
        {
            for (const auto f : m_f)
            {
                assert(f);

                m_last_total += f;
            }

            Recalculate_ranges();
        }

        void Encode(Encoder& coder, unsigned value)
        {
            assert(value < m_size);

            coder.Encode(m_r[value]);

            Adapt(value);
        }

        void Encode(Encoder& coder, unsigned value, unsigned max_value)
        {
            assert(value < m_size);
            assert(max_value < m_size);

            // Hmm, encoding one bit bad give me issues.
            max_value = std::max(max_value, 2u);

            // Multiply the ranges instead of the value
            // so the frequencies still get updated correctly.
            auto old_range  = m_r[value];
            auto last_range = m_r[max_value];
            float max_cf = last_range.min + last_range.count;
            float multiplier = PROBABILITY_RANGE / max_cf;

            // Man, rounding errors and off my one are my bane.
            unsigned next_min = 0;
            if (value < max_value)
            {
                next_min =
                    static_cast<uint32_t>
                    (
                        std::round(m_r[value + 1].min * multiplier)
                    );
            }
            else
            {
                next_min = PROBABILITY_RANGE;
            }

            auto new_min =
                static_cast<uint32_t>(std::round(old_range.min * multiplier));

            auto new_range = Range
            {
                new_min,
                next_min - new_min,
            };

            coder.Encode(new_range);

            Adapt(value);
        }

        unsigned Decode(Decoder& coder)
        {
            auto range = coder.Decode();
            unsigned result = 0;

            const auto size = m_size;
            while ((m_r[result].min <= range) && (result < size))
            {
                ++result;
            }

            --result;

            coder.Update(m_r[result]);
            Adapt(result);

            return result;
        }

        unsigned Decode(Decoder& coder, unsigned max_value)
        {
            assert(max_value < m_size);

            // Hmm, encoding one bit bad give me issues.
            max_value = std::max(max_value, 2u);

            auto range = coder.Decode();
            unsigned result = 0;

            auto last_range = m_r[max_value];
            float max_cf = last_range.min + last_range.count;
            float multiplier = PROBABILITY_RANGE / max_cf;

            const auto size = m_size;

            while (result < size)
            {
                auto new_min =
                    static_cast<uint32_t>
                    (
                        std::round(m_r[result].min * multiplier)
                    );

                if (new_min > range)
                {
                    break;
                }

                ++result;
            }

            --result;

            unsigned next_min = 0;
            if (result < max_value)
            {
                next_min =
                    static_cast<uint32_t>
                    (
                        std::round(m_r[result + 1].min * multiplier)
                    );
            }
            else
            {
                next_min = PROBABILITY_RANGE;
            }

            auto new_min =
                static_cast<uint32_t>(std::round(m_r[result].min * multiplier));

            auto new_range = Range
            {
                new_min,
                next_min - new_min,
            };

            coder.Update(new_range);
            Adapt(result);

            return result;
        }

    private:
        Freqencies  m_f;
        Ranges      m_r;
        Kernel      m_k;

        unsigned m_updates          = 0;
        unsigned m_update_total     = 0;
        unsigned m_update_trigger   = 1;
        unsigned m_last_total       = 0;

        const unsigned m_size                   = 0;
        const unsigned m_slowest_update_rate    = 0;
        const unsigned m_kernal_offset          = 0;

        void Adapt(unsigned value)
        {
            const auto kernal_offset = m_kernal_offset;
            const auto kernal_size = m_k.size();

            // Maybe add kernal_offset padding to both
            // sides of m_f to remove the ifs.
            const auto size = m_f.size();
            auto start =
                (value >= kernal_offset) ?
                    0 :
                    kernal_offset - value;

            auto tail_offset = (size - value) - 1;
            auto end =
                (tail_offset >= kernal_offset) ?
                    kernal_size :
                    kernal_size - (kernal_offset - tail_offset);

            const auto& k = m_k;
            for (auto i = start; i < end; ++i)
            {
                m_f[(value - kernal_offset) + i] += k[i];
                m_updates += k[i];
            }

            auto updates = m_update_total + 1;
            if (updates >= m_update_trigger)
            {
                m_update_trigger += m_update_trigger;
                if (m_update_trigger > m_slowest_update_rate)
                {
                    m_update_trigger = m_slowest_update_rate;
                }

                Recalculate_ranges();
                updates = 0;
            }
            m_update_total = updates;
        }

        void Recalculate_ranges()
        {
            unsigned total = m_last_total + m_updates;

            if (total > PROBABILITY_RANGE)
            {
                // Dunno how to do this nice for now. Brute it.
                total = 0;
                for (auto& f : m_f)
                {
                    f >>= 1;
                    ++f;

                    total += f;
                }
            }

            const auto size     = m_size;
            auto multiple       = PROBABILITY_RANGE / total;
            auto reminder       = PROBABILITY_RANGE % total;
            auto global_adjust  = reminder / size;
            auto reminder_count = reminder % size;
            unsigned last_min   = 0;

            for (unsigned i = 0; i < reminder_count; ++i)
            {
                m_r[i].min      = last_min;
                m_r[i].count    = 1 + (global_adjust + m_f[i] * multiple);

                last_min += m_r[i].count;
            }

            for (unsigned i = reminder_count; i < size; ++i)
            {
                m_r[i].min      = last_min;
                m_r[i].count    = 0 + (global_adjust + m_f[i] * multiple);

                last_min += m_r[i].count;
            }

            assert(last_min == PROBABILITY_RANGE);

            m_last_total = total;
            m_updates = 0;
        }
    };

    class Exp_update
    {
    public:
        static const constexpr unsigned DEFAULT_SPEED = 3;

        typedef std::vector<unsigned>  Freqencies;
        typedef std::vector<Range>     Ranges;

        Exp_update
        (
            unsigned size,
            unsigned speed = DEFAULT_SPEED
        )
            : m_r()
            , m_size(size)
            , m_speed(speed)
        {
            m_r.reserve(size);

            const unsigned number = PROBABILITY_RANGE / size;

            unsigned total = 0;
            for (unsigned i = 0; i < size; ++i)
            {
                m_r.push_back
                ({
                    total,
                    number
                });

                total += number;
            }

            Recalculate_ranges(total);
        }

        Exp_update
        (
            Freqencies frequencies,
            unsigned speed = DEFAULT_SPEED
        )
            : m_r()
            , m_size(frequencies.size())
            , m_speed(speed)
        {
            m_r.reserve(m_size);

            unsigned total = 0;
            for (const auto f : frequencies)
            {
                assert(f);

                m_r.push_back
                ({
                    total,
                    f
                });

                total += f;
            }

            Recalculate_ranges(total);
        }

        void Encode(Encoder& coder, unsigned value)
        {
            assert(value < m_size);

            coder.Encode(m_r[value]);

            Adapt(value);
        }

        void Encode(Encoder& coder, unsigned value, unsigned max_value)
        {
            assert(value < m_size);
            assert(max_value < m_size);

            // Hmm, encoding one bit bad give me issues.
            max_value = std::max(max_value, 2u);

            // Multiply the ranges instead of the value
            // so the frequencies still get updated correctly.
            auto old_range  = m_r[value];
            auto last_range = m_r[max_value];
            float max_cf = last_range.min + last_range.count;
            float multiplier = PROBABILITY_RANGE / max_cf;

            // Man, rounding errors and off my one are my bane.
            unsigned next_min = 0;
            if (value < max_value)
            {
                next_min =
                    static_cast<uint32_t>
                    (
                        std::round(m_r[value + 1].min * multiplier)
                    );
            }
            else
            {
                next_min = PROBABILITY_RANGE;
            }

            auto new_min =
                static_cast<uint32_t>(std::round(old_range.min * multiplier));

            auto new_range = Range
            {
                new_min,
                next_min - new_min,
            };

            coder.Encode(new_range);

            Adapt(value);
        }

        unsigned Decode(Decoder& coder)
        {
            auto range = coder.Decode();
            unsigned result = 0;

            const auto size = m_size;
            while ((m_r[result].min <= range) && (result < size))
            {
                ++result;
            }

            --result;

            coder.Update(m_r[result]);
            Adapt(result);

            return result;
        }

        unsigned Decode(Decoder& coder, unsigned max_value)
        {
            assert(max_value < m_size);

            // Hmm, encoding one bit bad give me issues.
            max_value = std::max(max_value, 2u);

            auto range = coder.Decode();
            unsigned result = 0;

            auto last_range = m_r[max_value];
            float max_cf = last_range.min + last_range.count;
            float multiplier = PROBABILITY_RANGE / max_cf;

            const auto size = m_size;

            while (result < size)
            {
                auto new_min =
                    static_cast<uint32_t>
                    (
                        std::round(m_r[result].min * multiplier)
                    );

                if (new_min > range)
                {
                    break;
                }

                ++result;
            }

            --result;

            unsigned next_min = 0;
            if (result < max_value)
            {
                next_min =
                    static_cast<uint32_t>
                    (
                        std::round(m_r[result + 1].min * multiplier)
                    );
            }
            else
            {
                next_min = PROBABILITY_RANGE;
            }

            auto new_min =
                static_cast<uint32_t>(std::round(m_r[result].min * multiplier));

            auto new_range = Range
            {
                new_min,
                next_min - new_min,
            };

            coder.Update(new_range);
            Adapt(result);

            return result;
        }

    private:
        Ranges      m_r;

        const unsigned m_size       = 0;
        const unsigned m_speed      = 3;

        void Adapt(unsigned value)
        {
            unsigned total = 0;
            const unsigned spread = m_size - 1;
            for (unsigned i = 0; i < m_size; ++i)
            {
                if (value == i)
                {
                    m_r[i].count +=
                    (
                        ((PROBABILITY_RANGE - m_r[i].count) >> m_speed)
                        / spread
                    )
                    * spread;
                }
                else
                {
                    m_r[i].count -=
                    (
                        (m_r[i].count >> m_speed)
                        / spread
                    )
                    * spread;
                }

                assert(m_r[i].min < PROBABILITY_RANGE);

                m_r[i].min = total;
                total += m_r[i].count;
            }

            Recalculate_ranges(total);
        }

        void Recalculate_ranges(unsigned total)
        {
            while (total > PROBABILITY_RANGE)
            {
                unsigned sub = 1 + ((total - PROBABILITY_RANGE) / m_size);

                total = 0;
                for (auto& r : m_r)
                {
                    r.count -= std::min(sub, r.count - 1);
                    r.min = total;
                    total += r.count;
                }
            }

            if (total == PROBABILITY_RANGE)
            {
                // yay, nothing to do!
                return;
            }

            const auto size     = m_size;
            auto multiple       = PROBABILITY_RANGE / total;
            auto reminder       = PROBABILITY_RANGE % total;
            auto global_adjust  = reminder / size;
            auto reminder_count = reminder % size;
            unsigned last_min   = 0;

            for (unsigned i = 0; i < reminder_count; ++i)
            {
                m_r[i].min      = last_min;
                m_r[i].count    = 1 + (global_adjust + m_r[i].count * multiple);

                assert(m_r[i].min < PROBABILITY_RANGE);

                last_min += m_r[i].count;
            }

            for (unsigned i = reminder_count; i < size; ++i)
            {
                m_r[i].min      = last_min;
                m_r[i].count    = 0 + (global_adjust + m_r[i].count * multiple);

                assert(m_r[i].min < PROBABILITY_RANGE);

                last_min += m_r[i].count;
            }

            assert(last_min == PROBABILITY_RANGE);
        }
    };

    class Unsigned_golomb_range
    {
    public:
        Unsigned_golomb_range(unsigned max_value_bits, unsigned speed = 5)
            : m_bits(max_value_bits, speed)
            , m_max_value_bits(max_value_bits)
        {
        }

        void Encode(Encoder& coder, unsigned value)
        {
            // Enocde how many bits we are sending
            unsigned min_bits = 0;

            while ((1u << min_bits) <= value)
            {
                ++min_bits;
            }

            assert(min_bits < m_max_value_bits);

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

                    Binary_bitstream::Encode
                    (
                        coder,
                        value & mask
                    );
                }
            }
        }

        unsigned Decode(Decoder& coder)
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

                    result |= Binary_bitstream::Decode(coder);

                    --min_bits;
                }
            }

            return result;
        }

    private:
        Exp_update m_bits = {4};
        unsigned m_max_value_bits;
    };

} // Namespace models

void range_tests()
{
    using namespace Range_types;
    using namespace Range_coders;
    using namespace Range_models;

    {
        Bytes data;

        // just test 3 ranges
        std::vector<Range> ranges
        {
            {0, 2000},
            {2000, 10000},
            {12000, (65536 - 12000)},
        };

        auto tests = {1,2,2,2,2,0,1,1,2,2,2,2,2,2,2,1,2,2,1};

        {
            Encoder encoder(data);

            for (const auto t : tests)
            {
                encoder.Encode(ranges[t]);
            }
        }

        {
            Decoder decoder(data);

            auto k =0;
            for (const auto t : tests)
            {
                auto value = decoder.Decode();
                assert(value >= ranges[t].min);
                assert(value < (ranges[t].min + ranges[t].count));
                decoder.Update(ranges[t]);
                k++;
            }

            auto read = decoder.FlushAndGetBytesRead();
            assert(read == data.size());
        }
    }

    {
        Bytes data;

        auto tests = {0,1,0,1,0,1,1,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0};

        {
            Encoder range_encoder(data);
            Binary_encoder encoder(range_encoder);

            for (const unsigned t : tests)
            {
                encoder.Encode(t, (PROBABILITY_RANGE * 2) / 3);
            }
        }

        {
            Decoder range_decoder(data);
            Binary_decoder decoder(range_decoder);

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
                Encoder range_encoder(data);
                Binary_encoder test_encoder(range_encoder);

                for (auto t : tests)
                {
                    model_in.Encode(test_encoder, t);
                }
            }
            {
                Decoder range_decoder(data);
                Binary_decoder test_decoder(range_decoder);

                for (unsigned t : tests)
                {
                    auto value = model_out.Decode(test_decoder);
                    assert(value == t);
                }

                auto read = test_decoder.FlushAndGetBytesRead();
                assert(read == data.size());
            }
        };

        Binary_test(tests, Binary_bitstream(), Binary_bitstream());
        Binary_test(tests, Binary(5), Binary(5));
        Binary_test(tests, Binary(1), Binary(1));
        Binary_test(tests, Binary(2), Binary(2));


        Binary_test(
            tests,
            Binary_two_speed(5,2),
            Binary_two_speed(5,2));
        Binary_test(
            tests,
            Binary_two_speed(6,1),
            Binary_two_speed(6,1));
        Binary_test(
            tests,
            Binary_two_speed(3,4),
            Binary_two_speed(3,4));
    }

    // Range Model Tests
    {
        auto range_data =
        {
            Bytes{0,1,2,6,4,5,3,7, 4,3,4,3,3,3,3,0, 1,2,0,0,0,3,3,3, 3,3,2},
            Bytes{0,0,0,0,0,0,0,0, 0,4,4,4,4,4,4,4, 2,2,4,4,4,4,4,0, 0,0,0},
            Bytes{0,1,2,3,4,5,6,7, 0,1,2,3,4,5,6,7, 0,1,2,3,4,5,6,7, 0,1,2},
        };

        auto Range_test = []
        (
            auto mod,
            auto tests,
            auto model_in,
            auto model_out,
            uint32_t max_value = 0
        )
        {
            Bytes data;

            if (!max_value)
            {
                Encoder test_encoder(data);

                for (auto t : tests)
                {
                    model_in.Encode(test_encoder, t % mod);
                }
            }
            else
            {
                Encoder test_encoder(data);

                for (auto t : tests)
                {
                    model_in.Encode(test_encoder, t % max_value, max_value);
                }
            }

            if (!max_value)
            {
                Decoder test_decoder(data);

                for (unsigned t : tests)
                {
                    auto value = model_out.Decode(test_decoder);
                    assert(value == (t % mod));
                }

                auto read = test_decoder.FlushAndGetBytesRead();
                assert(read == data.size());
            }
            else
            {
                Decoder test_decoder(data);

                for (unsigned t : tests)
                {
                    auto value = model_out.Decode(test_decoder, max_value);
                    assert(value == (t % max_value));
                }

                auto read = test_decoder.FlushAndGetBytesRead();
                assert(read == data.size());
            }
        };

        // RAM: Why the copy paste?

        for (const auto& range_set : range_data)
        {
            Exp_update::Freqencies
                    frequencies{1,1,1,1};
            Exp_update::Freqencies
                    overflow{65536,44,100000,34567};

            Range_test(
                4,
                range_set,
                Exp_update(frequencies),
                Exp_update(frequencies));

            Range_test(
                4,
                range_set,
                Exp_update(overflow),
                Exp_update(overflow));

            Range_test(
                8,
                range_set,
                Exp_update(8),
                Exp_update(8),
                3);
        }

        for (const auto& range_set : range_data)
        {
            Periodic_update_with_kernel::Freqencies
                    frequencies{1,1,1,1};
            Periodic_update_with_kernel::Freqencies
                    overflow{65536,44,100000,34567};

            Range_test(
                4,
                range_set,
                Periodic_update_with_kernel(frequencies, 8),
                Periodic_update_with_kernel(frequencies, 8));

            Range_test(
                4,
                range_set,
                Periodic_update_with_kernel(overflow, 8),
                Periodic_update_with_kernel(overflow, 8));

            Range_test(
                8,
                range_set,
                Periodic_update_with_kernel(8,8),
                Periodic_update_with_kernel(8,8),
                3);
        }

        for (const auto& range_set : range_data)
        {
            Periodic_update::Freqencies
                    frequencies{1,1,1,1};
            Periodic_update::Freqencies
                    overflow{65536,44,100000,34567};

            Range_test(
                4,
                range_set,
                Periodic_update(frequencies, 8),
                Periodic_update(frequencies, 8));

            Range_test(
                4,
                range_set,
                Periodic_update(overflow, 8),
                Periodic_update(overflow, 8));

            Range_test(
                8,
                range_set,
                Periodic_update(8,8),
                Periodic_update(8,8),
                3);
        }


        auto Range_test2 = []
        (
            auto mod,
            auto tests,
            auto model_in,
            auto model_out
        )
        {
            Bytes data;

            {
                Encoder test_encoder(data);

                for (auto t : tests)
                {
                    model_in.Encode(test_encoder, t % mod);
                }
            }

            Decoder test_decoder(data);

            for (unsigned t : tests)
            {
                auto value = model_out.Decode(test_decoder);
                assert(value == (t % mod));
            }

            auto read = test_decoder.FlushAndGetBytesRead();
            assert(read == data.size());
        };

        for (const auto& range_set : range_data)
        {
            Range_test2
            (
                1,
                range_set,
                Binary_bitstream(),
                Binary_bitstream()
            );

            Range_test2
            (
                8,
                range_set,
                Unsigned_golomb_range(8),
                Unsigned_golomb_range(8)
            );
        }
    }    

    // Binary tree tests.
    {
        Bytes data;
        Bytes tests{0,1,2,3,4,5,6,7, 0,1,2,3,4,5,6,7, 0,1,2,3,4,5,6,7, 0,1,2};

        {
            Encoder range_encoder(data);
            Binary_encoder test_encoder(range_encoder);
            Binary_tree_max_value<Binary_two_speed, 3, 5> tree;

            for (auto t : tests)
            {
                tree.Encode(test_encoder, t % 6);
            }
        }
        {
            Decoder range_decoder(data);
            Binary_decoder test_decoder(range_decoder);
            Binary_tree_max_value<Binary_two_speed, 3, 5> tree;

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
                encoder.Encode(t, 4000);
            }
        }

        {
            fpaq0p_decoder decoder(data);

            for (const unsigned t : tests)
            {
                auto value = decoder.Decode(4000);
                assert(value == t);
            }

            auto read = decoder.FlushAndGetBytesRead();
            assert(read == data.size());
        }
    }
}
