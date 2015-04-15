// Copyright 2015 Richard Maxwell, all rights reserved.

// Try doing http://www.arturocampos.com/ac_range.html
// NOTE: don't encode more than 2 << 24 values or you will
// have a bad time.

#pragma once

#include <vector>
#include <cstdint>

namespace Range_encoding
{
static const uint32_t TOP_VALUE     = 0x80000000;
static const uint32_t BOTTOM_VALUE  = 0x00800000;
static const uint32_t SHIFT_BITS    = 23;
static const uint32_t EXTRA_BITS    = 7;
static const uint32_t TOTAL_RANGE   = 65536;

using Bytes = std::vector<uint8_t>;

class Encoder
{
private:
    Bytes* m_bytes;
};

class Decoder
{
private:
    const Bytes* m_bytes;
};


} // namespace
