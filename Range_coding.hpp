// Copyright 2015 Richard Maxwell, all rights reserved.

// Try doing http://www.arturocampos.com/ac_range.html
// NOTE: don't encode more than 2 << 24 values or you will
// have a bad time.

#pragma once

#include <vector>
#include <array>
#include <cstdint>
#include <cassert>

// //////////////////////////////////////////////
// Types
// //////////////////////////////////////////////

namespace Range_types
{
    static const uint32_t TOTAL_RANGE_BITS  = 16;
    static const uint32_t TOTAL_RANGE       = 1 << TOTAL_RANGE_BITS;

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

namespace Range_coders
{
    using namespace Range_types;

    static const uint32_t TOP_VALUE         = 0x80000000;
    static const uint32_t BOTTOM_VALUE      = 0x00800000;
    static const uint32_t SHIFT_BITS        = 23;
    static const uint32_t EXTRA_BITS        = 7;

    class Encoder
    {
    public:
        Encoder(Bytes& bytes)
            : m_bytes(&bytes)
        {}

        ~Encoder()
        {
            Renormalise();

            auto temp = m_min >> SHIFT_BITS;

            if  (
                    (m_min & (BOTTOM_VALUE - 1)) >=
                    ((m_bytes->size() & 0xffffffL) >> 1)
                )
            {
                ++temp;
            }

            if (temp > 0xFF)
            {
                Flush_underflow(0, 1);
            }
            else
            {
                Flush_underflow(0xFF, 0);
            }

            Write(temp & 0xFF);
            Write(((temp = m_min) >> (SHIFT_BITS - 8)) & 0xFF);

            // I think we can get away with not needing these.
            // but whenever I remove them, it goes wrong.
            // I don't know why.
            Write(0);
            Write(0);
        }

        void Encode(Range range)
        {
            Renormalise();

            auto new_range = m_range >> TOTAL_RANGE_BITS;
            auto new_range_start = new_range * range.min;

            m_min += new_range_start;

            if (range.min + range.count < TOTAL_RANGE)
            {
                m_range = new_range * range.count;
            }
            else
            {
                m_range -= new_range_start;
            }
        }

    private:
        Bytes*      m_bytes         = 0;
        uint32_t    m_range         = TOP_VALUE;
        uint32_t    m_min           = 0;
        uint8_t     m_underflow     = 0;
        uint8_t     m_buffer        = 0;
        bool        m_first         = true;

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

        void Flush_underflow(uint8_t value, uint8_t overflow)
        {
            Write(m_buffer + overflow);

            while (m_underflow)
            {
                Write(value);
                --m_underflow;
            }

            m_buffer = m_min >> SHIFT_BITS;
        }

        void Renormalise()
        {
            while (m_range <= BOTTOM_VALUE)
            {
                if (m_min < (0xFF << SHIFT_BITS))
                {
                    Flush_underflow(0xFF, 0);
                }
                else
                {
                    if (m_min & TOP_VALUE)
                    {
                        Flush_underflow(0, 1);
                    }
                    else
                    {
                        ++m_underflow;
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

            m_next_range = m_range >> TOTAL_RANGE_BITS;
            auto symbol_range = m_min / m_next_range;

            if (symbol_range >= TOTAL_RANGE)
            {
                return TOTAL_RANGE - 1;
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
            if ((range.min + range.count) < TOTAL_RANGE)
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
                    Range{one_probability, TOTAL_RANGE - one_probability});
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
                    Range{one_probability, TOTAL_RANGE - one_probability});

            return result;
        }

        uint32_t FlushAndGetBytesRead()
        {
            return m_decoder->FlushAndGetBytesRead();
        }

    private:
        Decoder* m_decoder;
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
                unsigned initial_probability = (TOTAL_RANGE / 2))
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
                    (TOTAL_RANGE - m_one_probability) >> m_inertia;
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
        static const unsigned HALF_RANGE    = TOTAL_RANGE / 2;
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

    template<class BINARY_MODEL, unsigned BITS>
    class Binary_tree
    {
    public:
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

    template<class BINARY_MODEL_ZERO, class BINARY_MODEL_ONE>
    class Binary_history
    {
    public:
        Binary_history(
                unsigned last_value = 0,
                BINARY_MODEL_ZERO zero = BINARY_MODEL_ZERO(),
                BINARY_MODEL_ONE one = BINARY_MODEL_ONE())
            : m_last_value(last_value)
            , m_model_zero(zero)
            , m_model_one(one)
        {}

        void Encode(Binary_encoder& coder, unsigned value)
        {
            if (!m_last_value)
            {
                m_model_zero.Encode(coder, value);
            }
            else
            {
                m_model_one.Encode(coder, value);
            }

            m_last_value = value;
        }

        unsigned Decode(Binary_decoder& coder)
        {
            auto value =
                (!m_last_value) ?
                    m_model_zero.Decode(coder) :
                    m_model_one.Decode(coder);

            m_last_value = value;

            return value;
        }

    private:
        unsigned            m_last_value;
        BINARY_MODEL_ZERO   m_model_zero;
        BINARY_MODEL_ONE    m_model_one;
    };

    template<unsigned SIZE, unsigned SLOWEST_UPDATE_RATE>
    class Perodic_renomalisation
    {
    public:
        typedef std::array<unsigned, SIZE>  Freqencies;
        typedef std::array<Range, SIZE>     Ranges;

        Perodic_renomalisation()
            : m_f()
        {
            for (auto& f : m_f)
            {
                f = 1;
                m_last_total += f;
            }

            Recalculate_ranges();
        }

        Perodic_renomalisation(Freqencies frequencies)
            : m_f(frequencies)
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
            assert(value < SIZE);

            coder.Encode(m_r[value]);

            Adapt(value);
        }

        unsigned Decode(Decoder& coder)
        {
            auto range = coder.Decode();
            unsigned result = 0;

            while ((m_r[result].min <= range) && (result < SIZE))
            {
                ++result;
            }

            --result;

            coder.Update(m_r[result]);
            Adapt(result);

            return result;
        }

    private:
        Ranges      m_r;
        Freqencies  m_f;
        unsigned m_updates          = 0;
        unsigned m_update_trigger   = 1;
        unsigned m_last_total       = 0;

        void Adapt(unsigned value)
        {
            ++m_f[value];
            ++m_updates;

            if (m_updates >= m_update_trigger)
            {
                m_update_trigger += m_update_trigger;
                if (m_update_trigger > SLOWEST_UPDATE_RATE)
                {
                    m_update_trigger = SLOWEST_UPDATE_RATE;
                }

                Recalculate_ranges();
                m_updates = 0;
            }
        }

        void Recalculate_ranges()
        {
            unsigned total      = m_last_total + m_updates;

            if (total > TOTAL_RANGE)
            {
                // Dunno how to do this nice for now. Brute it.
                total = 0;
                for (auto& f : m_f)
                {
                    f >>= 8;
                    ++f;

                    total += f;
                }
            }

            auto multiple       = TOTAL_RANGE / total;
            auto reminder       = TOTAL_RANGE % total;
            auto global_adjust  = reminder / SIZE;
            auto reminder_count = reminder % SIZE;
            unsigned last_min   = 0;

            for (unsigned i = 0; i < SIZE; ++i)
            {
                m_r[i].min      = last_min;
                m_r[i].count    = global_adjust + m_f[i] * multiple;

                if (reminder_count)
                {
                    reminder_count--;
                    ++m_r[i].count;
                }

                last_min += m_r[i].count;
            }

            assert(last_min == TOTAL_RANGE);

            m_last_total = total;
        }
    };
} // Namespace models

void Range_tests()
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

            for (const auto t : tests)
            {
                auto value = decoder.Decode();
                assert(value >= ranges[t].min);
                assert(value < (ranges[t].min + ranges[t].count));
                decoder.Update(ranges[t]);
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
                encoder.Encode(t, 40000);
            }
        }

        {
            Decoder range_decoder(data);
            Binary_decoder decoder(range_decoder);

            for (const unsigned t : tests)
            {
                auto value = decoder.Decode(40000);
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

        Binary_test(
            tests,
            Binary_history<Binary, Binary>(0, Binary(3), Binary(3)),
            Binary_history<Binary, Binary>(0, Binary(3), Binary(3)));
    }

    // Range Model Tests
    {
        auto range_data =
        {
            Bytes{0,1,2,6,4,5,3,7,4,3,4,3,3,3,3,0,1,2,0,0,0,3,3,3,3,3,2},
            Bytes{0,0,0,0,0,0,0,0,0,4,4,4,4,4,4,4,2,2,4,4,4,4,4,0,0,0,0},
            Bytes{0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7,0,1,2},
        };

        auto Range_test = [](auto mod, auto tests, auto model_in, auto model_out)
        {
            Bytes data;

            {
                Encoder test_encoder(data);

                for (auto t : tests)
                {
                    model_in.Encode(test_encoder, t % mod);
                }
            }
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
        };

        for (const auto& range_set : range_data)
        {
            Perodic_renomalisation<4,8>::Freqencies
                    frequencies{1,1,1,1};
            Perodic_renomalisation<4,8>::Freqencies
                    overflow{65536,44,100000,34567};

            Range_test(
                4,
                range_set,
                Perodic_renomalisation<4,8>(frequencies),
                Perodic_renomalisation<4,8>(frequencies));

            Range_test(
                4,
                range_set,
                Perodic_renomalisation<4,8>(overflow),
                Perodic_renomalisation<4,8>(overflow));
        }
    }
}