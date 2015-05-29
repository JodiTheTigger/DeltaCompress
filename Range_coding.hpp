// Copyright 2015 Richard Maxwell, all rights reserved.

// Try doing http://www.arturocampos.com/ac_range.html
// NOTE: don't encode more than 2 << 24 values or you will
// have a bad time.

// RAM: TODO: Check out http://ezcodesample.com/reanatomy.html

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
            // Deal with overflow and underflow
            unsigned overflow = 0;
            {
                auto temp = m_min >> SHIFT_BITS;

                auto round_up = BOTTOM_VALUE - 1;
                auto a = m_min & round_up;
                auto s = m_bytes->size();
                auto sm = s & 0x00ffffff;
                auto b = sm >> 1;

                if  (
                        a >= b
                    )
                {
                    ++temp;
                }

                auto underflow_value = 0xFF;
                if (temp > 0xFF)
                {
                    overflow++;
                    underflow_value = 0;
                }

                Write(m_buffer + overflow);

                auto underflow = m_underflow;
                while (underflow)
                {
                    Write(underflow_value);
                    underflow--;
                }
            }

            // RAM: Let's try Fabian's trick.
            auto round_up = BOTTOM_VALUE - 1;
            auto min = m_min;
            auto hi = m_min + m_range;
            while (round_up)
            {
                auto rounded = (min + round_up) & ~round_up;

                if (rounded <= hi)
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

                min      <<= 8;
                min      &= (TOP_VALUE - 1);
            }
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

            if (symbol_range < TOTAL_RANGE)
            {
                return symbol_range;
            }
            else
            {
                return TOTAL_RANGE - 1;
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
                unsigned inertia_1 = 1,
                unsigned inertia_2 = 6,
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
            float multiplier = TOTAL_RANGE / max_cf;

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
                next_min = TOTAL_RANGE;
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
            float multiplier = TOTAL_RANGE / max_cf;

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
                next_min = TOTAL_RANGE;
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

            if (total > TOTAL_RANGE)
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
            auto multiple       = TOTAL_RANGE / total;
            auto reminder       = TOTAL_RANGE % total;
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

            assert(last_min == TOTAL_RANGE);

            m_last_total = total;
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
            float multiplier = TOTAL_RANGE / max_cf;

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
                next_min = TOTAL_RANGE;
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
            float multiplier = TOTAL_RANGE / max_cf;

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
                next_min = TOTAL_RANGE;
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

            if (total > TOTAL_RANGE)
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
            auto multiple       = TOTAL_RANGE / total;
            auto reminder       = TOTAL_RANGE % total;
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

            assert(last_min == TOTAL_RANGE);

            m_last_total = total;
            m_updates = 0;
        }
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

            auto k = 0;
            for (const auto t : tests)
            {
                auto value = decoder.Decode();
                assert(value >= ranges[t].min);
                assert(value < (ranges[t].min + ranges[t].count));
                decoder.Update(ranges[t]);
                ++k;
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
    }

    // Range Model Tests
    {
        auto range_data =
        {
            Bytes{0,1,2,6,4,5,3,7,4,3,4,3,3,3,3,0,1,2,0,0,0,3,3,3,3,3,2},
            Bytes{0,0,0,0,0,0,0,0,0,4,4,4,4,4,4,4,2,2,4,4,4,4,4,0,0,0,0},
            Bytes{0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7,0,1,2},
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

                auto t2 = 0;
                for (unsigned t : tests)
                {
                    auto value = model_out.Decode(test_decoder, max_value);
                    assert(value == (t % max_value));
                    ++t2;
                }

                auto read = test_decoder.FlushAndGetBytesRead();
                assert(read == data.size());
            }
        };

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
    }
}
