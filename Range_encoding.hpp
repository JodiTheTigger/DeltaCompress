// Copyright 2015 Richard Maxwell, all rights reserved.

// Try doing http://www.arturocampos.com/ac_range.html
// NOTE: don't encode more than 2 << 24 values or you will
// have a bad time.

#pragma once

#include <vector>
#include <array>
#include <cstdint>
#include <cassert>

namespace Range_encoding
{
static const uint32_t TOTAL_RANGE_BITS  = 16;

static const uint32_t TOP_VALUE         = 0x80000000;
static const uint32_t BOTTOM_VALUE      = 0x00800000;
static const uint32_t SHIFT_BITS        = 23;
static const uint32_t EXTRA_BITS        = 7;
static const uint32_t TOTAL_RANGE       = 1 << TOTAL_RANGE_BITS;

using Bytes = std::vector<uint8_t>;

struct Range
{
    uint32_t min;
    uint32_t count;
};

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
        // but whenever I remove them both, it goes wrong.
        // I don't know why.
        Write(0);
        //Write(0);
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

        // -1 since I don't send the last byte
        // when I encode.
        return m_read_index - 1;
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

class BinaryEncoder
{
public:
    BinaryEncoder(Bytes& bytes)
        : m_encoder(bytes)
    {}

    void Encode(unsigned value, uint16_t probability)
    {
        m_encoder.Encode(
            !value ?
                Range{0, probability} :
                Range{probability, TOTAL_RANGE - probability});
    }

private:
    Encoder m_encoder;
};

class BinaryDecoder
{
public:
    BinaryDecoder(const Bytes& bytes)
        : m_decoder(bytes)
    {}

    unsigned Decode(unsigned probability)
    {
        auto symbol = m_decoder.Decode();
        auto result = (symbol >= probability);

        m_decoder.Update(
            !result ?
                Range{0, probability} :
                Range{probability, TOTAL_RANGE - probability});

        return result;
    }

    uint32_t FlushAndGetBytesRead()
    {
        return m_decoder.FlushAndGetBytesRead();
    }

private:
    Decoder m_decoder;
};

// //////////////////////////////////////////////
// Models
// //////////////////////////////////////////////

namespace Models
{

// Ideas from https://github.com/rygorous/gaffer_net/blob/master/main.cpp

template<unsigned INERTIA>
class Binary
{
public:
    Binary(unsigned initial_probability = (TOTAL_RANGE / 2))
        : m_probability(initial_probability)
    {}

    void Encode(BinaryEncoder& coder, unsigned value)
    {
        coder.Encode(value, m_probability);
        Adapt(value);
    }

    unsigned Decode(BinaryDecoder& coder)
    {
        auto result = coder.Decode(m_probability);
        Adapt(result);

        return result;
    }

private:
    uint16_t  m_probability;

    void Adapt(unsigned value)
    {
        if (value)
        {
            m_probability += (TOTAL_RANGE - m_probability) >> INERTIA;
        }
        else
        {
            m_probability -= m_probability >> INERTIA;
        }
    }
};

template<unsigned INERTIA_1, unsigned INERTIA_2>
class Binary_two_speed
{
public:
    Binary_two_speed(
            unsigned initial_probability_1 = (TOTAL_RANGE / 4),
            unsigned initial_probability_2 = (TOTAL_RANGE / 4))
        : m_probability_1(initial_probability_1)
        , m_probability_2(initial_probability_2)
    {}

    void Encode(BinaryEncoder& coder, unsigned value)
    {
        coder.Encode(value, m_probability_1 + m_probability_2);
        Adapt(value);

    }

    unsigned Decode(BinaryDecoder& coder)
    {
        auto result = coder.Decode(m_probability_1 + m_probability_2);
        Adapt(result);

        return result;
    }

private:
    uint16_t  m_probability_1;
    uint16_t  m_probability_2;

    void Adapt(unsigned value)
    {
        if (value)
        {
            m_probability_1 += (TOTAL_RANGE - m_probability_1) >> INERTIA_1;
            m_probability_2 += (TOTAL_RANGE - m_probability_2) >> INERTIA_2;
        }
        else
        {
            m_probability_1 -= m_probability_1 >> INERTIA_1;
            m_probability_2 -= m_probability_2 >> INERTIA_2;
        }
    }
};

template<class BINARY_MODEL, unsigned BITS>
class Binary_tree
{
    void Encode(BinaryEncoder& coder, unsigned value)
    {
        assert(value < MODEL_COUNT);

        // Model the MSB first, then work our way down.
        // Seems adds are better than << 1.
        unsigned rebuilt = 1;

        while (rebuilt < MODEL_COUNT)
        {
            unsigned bit = ((value & TOP_BIT) != 0);
            value += value;
            m_models[rebuilt - 1].Encode(coder, bit);
            rebuilt += rebuilt + bit;
        }
    }

    unsigned Decode(BinaryDecoder& coder)
    {
        unsigned rebuilt = 1;
        while (rebuilt < MODEL_COUNT)
        {
            rebuilt += rebuilt + m_models[rebuilt - 1].decode(coder);
        }

        // Clear the top bit due to starting rebuilt with 1.
        return rebuilt - MODEL_COUNT;
    }

private:
    static const unsigned MODEL_COUNT   = 1 << BITS;
    static const unsigned TOP_BIT       = MODEL_COUNT / 2;

    std::array<BINARY_MODEL, MODEL_COUNT - 1> m_models;
};

template<unsigned SIZE, unsigned SLOWEST_UPDATE_RATE>
class PerodicRenomalisation
{
    typedef std::array<unsigned, SIZE>  Freqencies;
    typedef std::array<Range, SIZE>     Ranges;

    // How to init the entire array to 1?
    PerodicRenomalisation(Freqencies frequencies = {1})
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
        auto result = 0;

        while (m_r[result].min <= range)
        {
            ++result;
        }

        --result;

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
        auto multiple       = TOTAL_RANGE / total;
        auto reminder       = TOTAL_RANGE % total;
        unsigned last_min   = 0;

        for (unsigned i = 0; i < SIZE; ++i)
        {
            m_r[i].min      = last_min;
            m_r[i].count    = m_f[i] * multiple;

            if (reminder--)
            {
                ++m_r[i].count;
            }

            last_min += m_r[i].count;
        }

        assert(last_min == TOTAL_RANGE);
    }
};

} // Namespace models

void Tests()
{

}

} // namespace
