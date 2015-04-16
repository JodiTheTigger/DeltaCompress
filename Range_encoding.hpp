// Copyright 2015 Richard Maxwell, all rights reserved.

// Try doing http://www.arturocampos.com/ac_range.html
// NOTE: don't encode more than 2 << 24 values or you will
// have a bad time.

#pragma once

#include <vector>
#include <cstdint>

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
    }

    void Renormalise()
    {
        while (m_range <= BOTTOM_VALUE)
        {
            if (m_min < (0xFF << SHIFT_BITS))
            {
                Flush_underflow(0xFF, 0);

                m_buffer = m_min >> SHIFT_BITS;
            }
            else
            {
                if (m_min & TOP_VALUE)
                {
                    Flush_underflow(0, 1);

                    m_buffer = m_min >> SHIFT_BITS;
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
    uint8_t         m_buffer        = 0;

    uint8_t Read()
    {
        if (m_read_index < m_bytes->size())
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


} // namespace
