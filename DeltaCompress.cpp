// Copyright 2015 Richard Maxwell, all rights reserved.

// http://gafferongames.com/2015/03/14/the-networked-physics-data-compression-challenge/

// g++ -std=c++14 DeltaCompress.cpp -Wall -Wextra -Werror -g -o DeltaCompress

// //////////////////////////////////////////////////////

#define _CRT_SECURE_NO_WARNINGS

#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <cstdlib>

#include <array>
#include <vector>
#include <functional>
#include <cmath>

// //////////////////////////////////////////////////////

// http://the-witness.net/news/2012/11/scopeexit-in-c11/

template <typename F>
struct ScopeExit {
    ScopeExit(F f) : f(f) {}
    ~ScopeExit() { f(); }
    F f;
};

template <typename F>
ScopeExit<F> MakeScopeExit(F f) {
    return ScopeExit<F>(f);
};

#define STRING_JOIN2(arg1, arg2) DO_STRING_JOIN2(arg1, arg2)
#define DO_STRING_JOIN2(arg1, arg2) arg1 ## arg2
#define SCOPED_EXIT(code) \
    auto STRING_JOIN2(scope_exit_, __LINE__) = MakeScopeExit([=](){code;})

// //////////////////////////////////////////////////////

struct DeltaData
{
    int orientation_largest;
    int orientation_a;
    int orientation_b;
    int orientation_c;
    int position_x;
    int position_y;
    int position_z;
    int interacting;
};

inline bool operator==(const DeltaData& lhs, const DeltaData& rhs)
{
    return
    (
            (lhs.interacting    == rhs.interacting)
        &&  (lhs.orientation_a  == rhs.orientation_a)
        &&  (lhs.orientation_b  == rhs.orientation_b)
        &&  (lhs.orientation_c  == rhs.orientation_c)
        &&  (lhs.position_x     == rhs.position_x)
        &&  (lhs.position_y     == rhs.position_y)
        &&  (lhs.position_z     == rhs.position_z)
    );
}

inline bool operator!=(const DeltaData& lhs, const DeltaData& rhs){return !operator==(lhs,rhs);}

// //////////////////////////////////////////////////////

using ByteVector = std::vector<uint8_t>;

// //////////////////////////////////////////////////////

static const size_t Cubes = 901;
static const size_t CubeBits = 10;

static_assert(Cubes < ((1 << CubeBits) - 1), "CubeBits is too small.");

// Using contraints from
// http://gafferongames.com/networked-physics/snapshot-compression/

// With this technique I’ve found that minimum sufficient precision for
// my simulation is 9 bits per-smallest component. This gives a result
// of 2 + 9 + 9 + 9 = 29 bits per-orientation (originally 128!).

static const unsigned RotationMaxBits = 9;
static const unsigned RotationIndexMaxBits = 2;

// I found a maximum speed of 32 meters per-second is a nice power of two
// and doesn’t negatively affect the player experience in the cube simulation.
static const unsigned MaxSpeedMetersPerSecond = 32;

// Now we have only position to compress. We’ll use the same trick that we used
// for linear velocity: bound and quantize. Most game worlds are reasonably big
// so I chose a position bound of [-256,255] meters in the horizontal plane (xy)
// and since in my cube simulation the floor is at z=0, I chose for z a range of
// [0,32] meters.

// Now we need to work out how much precision is required. With some
// experimentation I found that 512 values per-meter (roughly 0.5mm precision)
// provides sufficient precision. This gives position x and y components
// in [-131072,+131071] and z components in range [0,16383].

// RAM: NOTE: ambiguous statement. So I'll use [-131072,+131071]

static const unsigned   ValuesPerMeter    = 512;
static const unsigned   MaxZInclusize     = (32 * ValuesPerMeter) - 1;
static const unsigned   MinZInclusize     = 0;
static const int        MaxXYInclusize    = 131071;  // (256 * 512) - 1
static const int        MinXYInclusize    = -131072; // (-256 * 512)
static const int        XYRange           = MaxXYInclusize - MinXYInclusize;

static const unsigned   MaxBitsZ          = 14;
static const unsigned   MaxBitsXY         = 18;

static const unsigned   MaxSnapshotsPerSecond = 60;

static const unsigned   FirstBase           = 6;
static const unsigned   PacketDelta         = 6;

// This one is important
static const float      MaxPositionChangePerSnapshot =
        MaxSpeedMetersPerSecond * ValuesPerMeter /
        (float) MaxSnapshotsPerSecond;

static_assert(
    ((1 << MaxBitsZ) - 1) == MaxZInclusize,
    "MaxBitZ doesn't match MaxZInclusize");

static_assert(
    ((1 << MaxBitsXY) - 1) == XYRange,
    "MaxBitsXY doesn't match XYRange");

typedef std::array<DeltaData, Cubes> Frame;

inline bool operator==(const Frame& lhs, const Frame& rhs)
{
    auto size = lhs.size();

    if (size != rhs.size())
    {
        return false;
    }

    for (size_t i = 0; i < size; ++i)
    {
        if (lhs[i] != rhs[i])
        {
            return false;
        }
    }

    return true;
}

inline bool operator!=(const Frame& lhs, const Frame& rhs){return !operator==(lhs,rhs);}

// //////////////////////////////////////////////////////

struct MinMaxSum
{
    unsigned min;
    unsigned max;
    unsigned sum;
    unsigned count;

    void Update(unsigned value)
    {
        min = std::min(min, value);
        max = std::max(max, value);
        sum += value;
        ++count;
    }

    float Average()
    {
        if (!count)
        {
            return 0;
        }

        return sum / count;
    }
};

struct Stats
{
    unsigned itemCount;

    unsigned notChangedCount;

    unsigned interactingTotal;
    unsigned interactingNotChanged;
    unsigned notInteractingNotChanged;

    std::array<unsigned, Cubes> changed0DistanceRunHistogram;
    std::array<unsigned, Cubes> changed1DistanceRunHistogram;

    MinMaxSum rle;
    MinMaxSum bitpack;
    MinMaxSum bitexprle;
    MinMaxSum bitbitpack;
    MinMaxSum bitbitfullpack;

    unsigned changed;
    unsigned quatChangedNoPos;
    unsigned PosChangedNotQuat;
    unsigned wtfBothSame;
    MinMaxSum deltaX;
    MinMaxSum deltaY;
    MinMaxSum deltaZ;
    MinMaxSum deltaTotal;

    unsigned posDeltaUnpackedBitCount;
    MinMaxSum posDeltaPackedBitCount;
};

enum class ChangedArrayEncoding
{
    None,
    Rle,
    BitPack,
    Exp,
    BitBitPack,
    BitBitFullPack
};

enum class PosVector3Packer
{
    None,
    BitVector3,
    BitVector3_2BitExpPrefix,
    BitVector3Truncated
};

void DoSomeStats(const Frame& base, const Frame& target, Stats& stats)
{
    assert(base.size() == target.size());
    assert(!base.empty());

    auto size = base.size();
    for (size_t i = 0; i < size; ++i)
    {
        stats.itemCount++;

        if (base[i] == target[i])
        {
            stats.notChangedCount++;
        }

        if (base[i].interacting)
        {
            stats.interactingTotal++;

            if (base[i].interacting == target[i].interacting)
            {
                stats.interactingNotChanged++;
            }
        }

        if (!base[i].interacting)
        {
            if (base[i].interacting == target[i].interacting)
            {
                stats.notInteractingNotChanged++;
            }
        }
    }
}

// //////////////////////////////////////////////////////

inline constexpr uint32_t ZigZag(int32_t n)
{
    return (n << 1) ^ (n >> 31);
}

inline constexpr int32_t ZigZag(uint32_t n)
{
    return (n >> 1) ^ (-(static_cast<int>(n) & 1));
}

void ZigZagTest()
{
    assert(42 == ZigZag(ZigZag(42)));
    assert(-42 == ZigZag(ZigZag(-42)));
    assert(0 == ZigZag(ZigZag(0)));
    assert(-12345 == ZigZag(ZigZag(-12345)));
    assert(30654 == ZigZag(ZigZag(30654)));
    assert(-31654 == ZigZag(ZigZag(-31654)));
}

// //////////////////////////////////////////////////////

class BitStream
{
public:
    BitStream()
    {
    }

    BitStream(std::vector<uint8_t> newData)
        : m_data{newData}
    {
    }

    BitStream(std::vector<uint8_t> newData, size_t bits)
        : m_bitOffset(bits)
        , m_data{newData}
    {
    }

    size_t Bits() const
    {
        // Next place to write, so in affect == existing size.
        return m_bitOffset;
    }

    void SetOffset(unsigned offset)
    {
        m_bitOffset = offset;
    }

    void Reset()
    {
        SetOffset(0);
    }

    void Write(const BitStream& other)
    {
        const auto newOffset = m_bitOffset + other.m_bitOffset;
        const auto newSize = ((newOffset + 7) / 8);

        for (const auto& b : other.m_data)
        {
            Write(b, 8);
        }

        while (m_data.size() > newSize)
        {
            m_data.pop_back();
        }

        m_bitOffset = newOffset;
    }

    void Write(unsigned value, unsigned bitsToWrite)
    {
        assert(value <= ((1u << bitsToWrite) - 1));

        auto bitOffset = m_bitOffset;
        auto byteOffset = bitOffset / 8u;
        auto byteOffsetAfter = (bitOffset + bitsToWrite + -1) / 8u;
        auto index = bitOffset - (byteOffset * 8);
        auto mask = ((1 << index) - 1);

        while (m_data.size() <= byteOffsetAfter)
        {
            m_data.push_back(0);
        }

        while (bitsToWrite)
        {
            m_data[byteOffset] |= (value & 0xFF) << index;

            if ((index + bitsToWrite) > 8)
            {
                m_data[byteOffset + 1] |= (value >> (8 - index)) & mask;
            }

            if (bitsToWrite >= 8)
            {
                value >>= 8;
                byteOffset++;

                bitOffset+= 8;
                bitsToWrite -= 8;
            }
            else
            {
                bitOffset+= bitsToWrite;
                bitsToWrite = 0;
            }
        }

        m_bitOffset = bitOffset;
    }

    unsigned Read(unsigned bitsToRead)
    {
        auto bitOffset = m_bitOffset;
        auto byteOffset = bitOffset / 8u;
        auto index = bitOffset - (byteOffset * 8);
        auto size = m_data.size();
        auto mask = ((1 << bitsToRead) - 1);
        auto shift = 0;

        unsigned value = 0;

        while (bitsToRead)
        {
            if (byteOffset < size)
            {
                value |= (m_data[byteOffset] >> index) << shift;
            }

            if (((index + bitsToRead) > 8) && ((byteOffset + 1) < size))
            {
                value |= ((m_data[byteOffset + 1] << (8 - index) & 0xFF)) << shift;
            }

            if (bitsToRead >= 8)
            {
                byteOffset++;
                shift+=8;

                bitOffset+= 8;
                bitsToRead -= 8;
            }
            else
            {
                bitOffset += bitsToRead;
                bitsToRead = 0;
            }
        }

        value &= mask;
        m_bitOffset = bitOffset;

        return value;
    }

    BitStream ReadArray(unsigned bitsToRead)
    {
        BitStream result;

        while (bitsToRead > 7)
        {
            result.Write(Read(8), 8);
            bitsToRead -= 8;
        }

        if (bitsToRead)
        {
            result.Write(Read(bitsToRead), bitsToRead);
        }

        return result;
    }

    void TrimZerosFromBack()
    {
        while (m_data.back() == 0)
        {
            m_data.pop_back();
        }
    }

    std::vector<uint8_t> Data() { return m_data; }

    inline bool operator==(const BitStream& rhs)
    {
        if (this->m_bitOffset != rhs.m_bitOffset)
        {
            return false;
        }

        return m_data == rhs.m_data;
    }

    inline bool operator!=(const BitStream& rhs)
    {
        return !operator==(rhs);
    }

private:
    size_t                  m_bitOffset   = 0;
    std::vector<uint8_t>    m_data        = {};
};

void BitStreamTest()
{
    BitStream bitsIn;

    bitsIn.Write(46, 6);
    bitsIn.Write(666, 10);
    bitsIn.Write(169, 8);

    BitStream bitsOut(bitsIn.Data());

    auto a = bitsOut.Read(6);
    auto b = bitsOut.Read(10);
    auto c = bitsOut.Read(8);

    assert(a == 46);
    assert(b == 666);
    assert(c == 169);

    // ////////////////////////////////

    BitStream bitsInB;

    bitsInB.Write(3, 2);
    bitsInB.Write(256, 9);
    bitsInB.Write(256, 9);

    BitStream bitsOutB(bitsInB.Data());

    auto ab = bitsOutB.Read(2);
    auto bb = bitsOutB.Read(9);
    auto cb = bitsOutB.Read(9);

    assert(ab == 3);
    assert(bb == 256);
    assert(cb == 256);

    // ////////////////////////////////

    BitStream bitsOutToSplit(bitsIn.Data());
    BitStream bitsSplit = bitsOutToSplit.ReadArray(6);
    bitsSplit.Reset();

    auto as = bitsSplit.Read(6);
    auto bs = bitsOutToSplit.Read(10);
    auto cs = bitsOutToSplit.Read(8);

    assert(as == 46);
    assert(bs == 666);
    assert(cs == 169);

    // ////////////////////////////////

    BitStream bitsIn1;
    BitStream bitsIn2;

    bitsIn1.Write(46, 6);
    bitsIn2.Write(666, 10);
    bitsIn2.Write(169, 8);

    BitStream bitJoin;
    bitJoin.Write(bitsIn1);
    bitJoin.Write(bitsIn2);

    bitJoin.Reset();

    auto aj = bitJoin.Read(6);
    auto bj = bitJoin.Read(10);
    auto cj = bitJoin.Read(8);

    assert(aj == 46);
    assert(bj == 666);
    assert(cj == 169);
}

// //////////////////////////////////////////////////////

unsigned MinBits(unsigned value)
{
    unsigned result = 0;

    while ((1u << result) <= value)
    {
        ++result;
    }

    return result;
}

// //////////////////////////////////////////////////////

// Unneeded if we change the endianess of our bitvector.
unsigned Flip(unsigned value, unsigned bits)
{
    unsigned result = 0;

    while (bits--)
    {
        if (value & 1)
        {
            result |= 1 << bits;
        }

        value >>= 1;
    }

    return result;
}

unsigned TruncateEncode(unsigned value, unsigned maxValue, BitStream& target)
{
    // assume 0-maxValue inclusive.
    unsigned count = maxValue + 1;

    // http://en.wikipedia.org/wiki/Truncated_binary_encoding
    auto bits = MinBits(count);
    auto k = bits - 1;
    auto u = (1u << bits) - count;

    if (value < u)
    {
        target.Write(Flip(value, k), k);
        return k;
    }
    else
    {
        value += u;
        target.Write(Flip(value, bits), bits);
        return bits;
    }
}

unsigned TruncateDecode(unsigned maxValue, BitStream& source)
{
    unsigned count = maxValue + 1;
    auto bits = MinBits(count);
    auto k = bits - 1;
    auto u = (1u << bits) - count;

    // Bah, by bitstream is the wrong endianess, so I have to
    // do this a hacky way.
    auto value = Flip(source.Read(k), k);

    if (value < u)
    {
        return value;
    }

    value <<= 1;
    value |= source.Read(1);
    return value - u;
}

void TestTruncate()
{
    const unsigned count = 10;
    const unsigned max = count - 1;

    for (unsigned i = 0; i < count; ++i)
    {
        BitStream coded;
        TruncateEncode(i, max, coded);
        coded.Reset();
        auto result = TruncateDecode(max, coded);

        assert(result == i);
    }
}

// //////////////////////////////////////////////////////

struct IntVec3
{
    int x;
    int y;
    int z;
};

unsigned BitVector3Encode(
        IntVec3 vec,
        unsigned maxMagnitude,
        BitStream& target)
{
    unsigned bitsUsed = 0;

    // +1 for the sign bit.
    unsigned maxBitsRequired = 1 + MinBits(maxMagnitude);
    unsigned maxPrefixSize = MinBits(maxBitsRequired);

    assert(abs(vec.x) <= static_cast<int>(maxMagnitude));
    auto zx = ZigZag(vec.x);

    unsigned bitsForzx = MinBits(zx);

    // TODO: Use truncation binary encoding to shave off a few bits.
    target.Write(bitsForzx, maxPrefixSize);
    bitsUsed += maxPrefixSize;

    // Don't need to send the top bit, as we know it's set since
    // otherwise bitsForzx would be smaller.
    if (bitsForzx)
    {
        zx &= (1 << (bitsForzx - 1)) - 1;
        target.Write(zx, bitsForzx - 1);
        bitsUsed += bitsForzx - 1;
    }

    float next = maxMagnitude * maxMagnitude;
    next -= vec.x * vec.x;
    assert(next >= 0);
    maxMagnitude = static_cast<unsigned>(sqrt(next) + 1);

    // //////////////////////////////////////////

    maxBitsRequired = 1 + MinBits(maxMagnitude);
    maxPrefixSize = MinBits(maxBitsRequired);

    assert(abs(vec.y) <= static_cast<int>(maxMagnitude));
    auto zy = ZigZag(vec.y);

    unsigned bitsForzy = MinBits(zy);

    target.Write(bitsForzy, maxPrefixSize);
    bitsUsed += maxPrefixSize;

    if (bitsForzy)
    {
        zy &= (1 << (bitsForzy - 1)) - 1;
        target.Write(zy, bitsForzy - 1);
        bitsUsed += bitsForzy - 1;
    }

    next = maxMagnitude * maxMagnitude;
    next -= vec.y * vec.y;
    assert(next >= 0);
    maxMagnitude = static_cast<unsigned>(sqrt(next) + 1);

    // //////////////////////////////////////////

    maxBitsRequired = 1 + MinBits(maxMagnitude);
    //maxPrefixSize = MinBits(maxBitsRequired);

    assert(abs(vec.z) <= static_cast<int>(maxMagnitude));
    auto zz = ZigZag(vec.z);

    //unsigned bitsForzz = MinBits(zz);

    //target.Write(bitsForzz, maxPrefixSize);
    target.Write(zz, maxBitsRequired);
    bitsUsed += maxBitsRequired;

    return bitsUsed;
}

IntVec3 BitVector3Decode(
        unsigned maxMagnitude,
        BitStream& source)
{
    IntVec3 result = {0,0,0};

    // +1 for the sign bit.
    unsigned maxBitsRequired = 1 + MinBits(maxMagnitude);
    unsigned maxPrefixSize = MinBits(maxBitsRequired);

    auto bitsX = source.Read(maxPrefixSize);

    if (bitsX)
    {
        auto zx = source.Read(bitsX - 1);
        // Don't need to send the top bit, as we know it's set since
        // otherwise bitsForzx would be smaller.
        zx |= 1 << (bitsX - 1);
        result.x = ZigZag(zx);
    }

    float next = maxMagnitude * maxMagnitude;
    next -= result.x * result.x;
    maxMagnitude = static_cast<unsigned>(sqrt(next) + 1);

    // //////////////////////////////////////////

    maxBitsRequired = 1 + MinBits(maxMagnitude);
    maxPrefixSize = MinBits(maxBitsRequired);

    auto bitsY = source.Read(maxPrefixSize);
    if (bitsY)
    {
        auto zy = source.Read(bitsY - 1);
        zy |= 1 << (bitsY - 1);
        result.y = ZigZag(zy);
    }

    next = maxMagnitude * maxMagnitude;
    next -= result.y * result.y;
    maxMagnitude = static_cast<unsigned>(sqrt(next) + 1);

    // //////////////////////////////////////////

    maxBitsRequired = 1 + MinBits(maxMagnitude);
    //maxPrefixSize = MinBits(maxBitsRequired);

    //auto bitsZ = source.Read(maxPrefixSize);
    auto zz = source.Read(maxBitsRequired);
    result.z = ZigZag(zz);

    return result;
}

unsigned BitVector3TruncatedEncode(
        IntVec3 vec,
        unsigned maxMagnitude,
        BitStream& target)
{
    unsigned bitsUsed = 0;

    // +1 for the sign bit.
    unsigned maxBitsRequired = 1 + MinBits(maxMagnitude);

    assert(abs(vec.x) <= static_cast<int>(maxMagnitude));
    auto zx = ZigZag(vec.x);

    unsigned bitsForzx = MinBits(zx);
    bitsUsed += TruncateEncode(bitsForzx, maxBitsRequired, target);

    // Don't need to send the top bit, as we know it's set since
    // otherwise bitsForzx would be smaller.
    if (bitsForzx)
    {
        zx &= (1 << (bitsForzx - 1)) - 1;
        target.Write(zx, bitsForzx - 1);
        bitsUsed += bitsForzx - 1;
    }

    float next = maxMagnitude * maxMagnitude;
    next -= vec.x * vec.x;
    assert(next >= 0);
    maxMagnitude = static_cast<unsigned>(sqrt(next) + 1);

    // //////////////////////////////////////////

    maxBitsRequired = 1 + MinBits(maxMagnitude);

    assert(abs(vec.y) <= static_cast<int>(maxMagnitude));
    auto zy = ZigZag(vec.y);

    unsigned bitsForzy = MinBits(zy);
    bitsUsed += TruncateEncode(bitsForzy, maxBitsRequired, target);

    if (bitsForzy)
    {
        zy &= (1 << (bitsForzy - 1)) - 1;
        target.Write(zy, bitsForzy - 1);
        bitsUsed += bitsForzy - 1;
    }

    next = maxMagnitude * maxMagnitude;
    next -= vec.y * vec.y;
    assert(next >= 0);
    maxMagnitude = static_cast<unsigned>(sqrt(next) + 1);

    // //////////////////////////////////////////

    maxBitsRequired = 1 + MinBits(maxMagnitude);
    //maxPrefixSize = MinBits(maxBitsRequired);

    assert(abs(vec.z) <= static_cast<int>(maxMagnitude));
    auto zz = ZigZag(vec.z);

    //unsigned bitsForzz = MinBits(zz);

    //target.Write(bitsForzz, maxPrefixSize);
    target.Write(zz, maxBitsRequired);
    bitsUsed += maxBitsRequired;

    return bitsUsed;
}

IntVec3 BitVector3TruncatedDecode(
        unsigned maxMagnitude,
        BitStream& source)
{
    IntVec3 result = {0,0,0};

    // +1 for the sign bit.
    unsigned maxBitsRequired = 1 + MinBits(maxMagnitude);

    auto bitsX = TruncateDecode(maxBitsRequired, source);

    if (bitsX)
    {
        auto zx = source.Read(bitsX - 1);
        // Don't need to send the top bit, as we know it's set since
        // otherwise bitsForzx would be smaller.
        zx |= 1 << (bitsX - 1);
        result.x = ZigZag(zx);
    }

    float next = maxMagnitude * maxMagnitude;
    next -= result.x * result.x;
    maxMagnitude = static_cast<unsigned>(sqrt(next) + 1);

    // //////////////////////////////////////////

    maxBitsRequired = 1 + MinBits(maxMagnitude);

    auto bitsY = TruncateDecode(maxBitsRequired, source);
    if (bitsY)
    {
        auto zy = source.Read(bitsY - 1);
        zy |= 1 << (bitsY - 1);
        result.y = ZigZag(zy);
    }

    next = maxMagnitude * maxMagnitude;
    next -= result.y * result.y;
    maxMagnitude = static_cast<unsigned>(sqrt(next) + 1);

    // //////////////////////////////////////////

    maxBitsRequired = 1 + MinBits(maxMagnitude);
    //maxPrefixSize = MinBits(maxBitsRequired);

    //auto bitsZ = source.Read(maxPrefixSize);
    auto zz = source.Read(maxBitsRequired);
    result.z = ZigZag(zz);

    return result;
}

static const unsigned prefixBits = 2;

unsigned ShiftsRequired(unsigned value, unsigned maxBits)
{
    assert(value < (1u << maxBits));

    unsigned maxShift = ((1u << prefixBits) - 1u);
    unsigned shifts = 1;

    if (!value)
    {
        return maxShift;
    }

    while (value < (1u << (maxBits >> shifts)))
    {
        shifts++;
    }

    shifts--;

    if (shifts > maxShift)
    {
        shifts = maxShift;
    }

    return shifts;
}

unsigned BitVector3Encode2BitExpPrefix(
        IntVec3 vec,
        unsigned maxMagnitude,
        BitStream& target)
{
    unsigned bitsUsed = 0;
    unsigned maxBitsRequired = 1 + MinBits(maxMagnitude);

    assert(abs(vec.x) <= static_cast<int>(maxMagnitude));
    auto zx = ZigZag(vec.x);

    unsigned shiftsForzx = ShiftsRequired(zx, maxBitsRequired);

    target.Write(shiftsForzx, prefixBits);
    bitsUsed += prefixBits;

    target.Write(zx, maxBitsRequired >> shiftsForzx);
    bitsUsed += maxBitsRequired >> shiftsForzx;

    float next = maxMagnitude * maxMagnitude;
    next -= vec.x * vec.x;
    assert(next >= 0);
    maxMagnitude = static_cast<unsigned>(sqrt(next) + 1);

    // //////////////////////////////////////////

    maxBitsRequired = 1 + MinBits(maxMagnitude);

    assert(abs(vec.y) <= static_cast<int>(maxMagnitude));
    auto zy = ZigZag(vec.y);

    unsigned shiftsForzy = ShiftsRequired(zy, maxBitsRequired);

    target.Write(shiftsForzy, prefixBits);
    bitsUsed += prefixBits;

    target.Write(zy, maxBitsRequired >> shiftsForzy);
    bitsUsed += maxBitsRequired >> shiftsForzy;

    next = maxMagnitude * maxMagnitude;
    next -= vec.y * vec.y;
    assert(next >= 0);
    maxMagnitude = static_cast<unsigned>(sqrt(next) + 1);

    // //////////////////////////////////////////

    maxBitsRequired = 1 + MinBits(maxMagnitude);
    //maxPrefixSize = MinBits(maxBitsRequired);

    assert(abs(vec.z) <= static_cast<int>(maxMagnitude));
    auto zz = ZigZag(vec.z);

    //unsigned bitsForzz = MinBits(zz);

    //target.Write(bitsForzz, maxPrefixSize);
    target.Write(zz, maxBitsRequired);
    bitsUsed += maxBitsRequired;

    return bitsUsed;
}

IntVec3 BitVector3Decode2BitExpPrefix(
        unsigned maxMagnitude,
        BitStream& source)
{
    IntVec3 result = {0,0,0};

    unsigned maxBitsRequired = 1 + MinBits(maxMagnitude);

    auto shiftX = source.Read(prefixBits);

    auto zx = source.Read(maxBitsRequired >> shiftX);
    result.x = ZigZag(zx);

    float next = maxMagnitude * maxMagnitude;
    next -= result.x * result.x;
    maxMagnitude = static_cast<unsigned>(sqrt(next) + 1);

    // //////////////////////////////////////////

    maxBitsRequired = 1 + MinBits(maxMagnitude);

    auto shiftY = source.Read(prefixBits);

    auto zy = source.Read(maxBitsRequired >> shiftY);
    result.y = ZigZag(zy);

    next = maxMagnitude * maxMagnitude;
    next -= result.y * result.y;
    maxMagnitude = static_cast<unsigned>(sqrt(next) + 1);

    // //////////////////////////////////////////

    maxBitsRequired = 1 + MinBits(maxMagnitude);
    //maxPrefixSize = MinBits(maxBitsRequired);

    //auto bitsZ = source.Read(maxPrefixSize);
    auto zz = source.Read(maxBitsRequired);
    result.z = ZigZag(zz);

    return result;
}

void BitVector3Tests()
{
    struct Test
    {
        std::function<unsigned(IntVec3, unsigned, BitStream&)> encode;
        std::function<IntVec3(unsigned, BitStream&)> decode;
    };

    std::vector<Test> tests;

    tests.push_back(
    {
        BitVector3TruncatedEncode,
        BitVector3TruncatedDecode,
    });

    tests.push_back(
    {
        BitVector3Encode2BitExpPrefix,
        BitVector3Decode2BitExpPrefix,
    });

    tests.push_back(
    {
        BitVector3Encode,
        BitVector3Decode,
    });

    int const max = (MaxPositionChangePerSnapshot) * 6 + 1;
    for (const auto& test : tests)
    {
        {
            IntVec3 data =
            {
                -1,
                -1586,
                0,
            };

            BitStream encoded;

            test.encode(
                data,
                max,
                encoded);

            encoded.Reset();

            auto decoded = test.decode(max, encoded);

            assert(data.x == decoded.x);
            assert(data.y == decoded.y);
            assert(data.z == decoded.z);
        }

        {
            IntVec3 data =
            {
                0,
                0,
                1,
            };

            BitStream encoded;

            test.encode(
                data,
                max,
                encoded);

            encoded.Reset();

            auto decoded = test.decode(max, encoded);

            assert(data.x == decoded.x);
            assert(data.y == decoded.y);
            assert(data.z == decoded.z);
        }

        for (auto i = 5 - max; i < max; i+=23)
        {
            for (auto j = 53 - max; j < max; j+=37)
            {
                for (auto k = 0; k < max; k+=17)
                {
                    auto mag = sqrt(i*i + j*j + k*k);

                    if (mag > max)
                    {
                        continue;
                    }

                    IntVec3 data =
                    {
                        i,
                        j,
                        k,
                    };

                    BitStream encoded;

                    test.encode(
                        data,
                        max,
                        encoded);

                    encoded.Reset();

                    auto decoded = test.decode(max, encoded);

                    assert(data.x == decoded.x);
                    assert(data.y == decoded.y);
                    assert(data.z == decoded.z);
                }
            }
        }
    }
}

// //////////////////////////////////////////////////////

ByteVector RunLengthEncode(const ByteVector& data)
{
    auto size = data.size();

    if (size < 3)
    {
        return data;
    }

    ByteVector  result;
    auto        previous    = data[0];
    unsigned    index       = 1;

    result.push_back(previous);

    while (index < size)
    {
        auto current = data[index++];

        if (previous == current)
        {
            unsigned run = 0;

            while (index < size)
            {
                current = data[index++];

                if (current != previous)
                {
                    break;
                }

                if (run == 255)
                {
                    break;
                }

                ++run;
            }

            result.push_back(previous);
            result.push_back(run);

            if (run == 255)
            {
                result.push_back(current);
            }
        }

        if (current != previous)
        {
            result.push_back(current);
        }

        previous = current;
    }

    return result;
}

ByteVector RunLengthDecode(
    const ByteVector& data,
    unsigned& bytesConsumed,
    unsigned maxResultBytes = 0)
{
    auto size = data.size();

    if (size < 3)
    {
        return data;
    }

    ByteVector  result;
    auto        previous    = data[0];
    unsigned    index       = 1;

    result.push_back(previous);

    while (index < size)
    {
        if (result.size() == maxResultBytes)
        {
            break;
        }

        auto current = data[index++];

        if (previous == current)
        {
            unsigned run = 0;

            if (index < size)
            {
                run = data[index++];
            }

            bool max = (run == 255);

            result.push_back(current);
            while (run--)
            {
                result.push_back(current);
            }

            if (result.size() == maxResultBytes)
            {
                continue;
            }

            if (index < size)
            {
                current = data[index++];
            }

            if (max)
            {
                result.push_back(current);
            }
        }

        if (current != previous)
        {
            result.push_back(current);
        }

        previous = current;
    }

    bytesConsumed = index;
    return result;
}

ByteVector RunLengthDecode(const ByteVector& data)
{
    unsigned unused;
    return RunLengthDecode(data, unused);
}

// //////////////////////////////////////////////////////

// http://michael.dipperstein.com/rle/

ByteVector BitPackEncode(const ByteVector& data)
{
    auto size = data.size();

    if (size < 2)
    {
        return data;
    }

    ByteVector  result;
    unsigned startIndex = 0;

    while (startIndex < size)
    {
        unsigned i = 1;
        unsigned run = 1;

        while (i != 128)
        {
            auto index = startIndex + i;
            if (data[index] == data[index - 1])
            {
                run++;
            }
            else
            {
                run = 1;
            }

            if (run == 3)
            {
                break;
            }

            ++i;

            if ((startIndex + i) == size)
            {
                break;
            }
        }

        unsigned unRunCount = 1 + (i - run);

        if (unRunCount)
        {
            auto header = ZigZag(static_cast<int32_t>(unRunCount) - 1);

            result.push_back(header);

            unsigned index = 0;
            while (index != unRunCount)
            {
                result.push_back(data[startIndex + index]);

                index++;
            }
        }

        startIndex += unRunCount;

        if (run == 3)
        {
            while (run < 130)
            {
                auto index = startIndex + run;

				if ((startIndex + run) == size)
				{
					break;
				}

                if (data[index] != data[index - 1])
                {
                    break;
                }

                ++run;
            }

            auto header = ZigZag(-(static_cast<int32_t>(run) - 2));

            result.push_back(header);
            result.push_back(data[startIndex]);

            startIndex += run;
        }

        if ((startIndex + 1) == size)
        {
            auto header = ZigZag(static_cast<int32_t>(0));

            result.push_back(header);
            result.push_back(data[startIndex]);
            startIndex++;
        }
    }

    return result;
}

ByteVector BitPackDecode(
    const ByteVector& data,
    unsigned& bytesConsumed,
    unsigned maxResultBytes = 0)
{
    auto size = data.size();

    if (size < 2)
    {
        return data;
    }

    ByteVector  result;
    unsigned index = 0;
    while (index < size)
    {
        auto header = ZigZag(static_cast<uint32_t>(data[index++]));

        if (header >= 0)
        {
            auto count = header + 1;

            while (count--)
            {
                result.push_back(data[index++]);
            }
        }
        else
        {
            auto count = -header + 2u;

            while (count--)
            {
                result.push_back(data[index]);
            }

            index++;
        }

        if (maxResultBytes > 0)
        {
            if (result.size() == maxResultBytes)
            {
                break;
            }
        }
    }

    bytesConsumed = index;

    return result;
}

ByteVector BitPackDecode(const ByteVector& data)
{
    unsigned unused;
    return BitPackDecode(data, unused);
}

// //////////////////////////////////////////////////////

static const std::array<unsigned, 8> exponentialBitLevelRunLengthLookup
{
    0, 2, 10, 27, 57, 90, 256, 898,
};

BitStream ExponentialBitLevelRunLengthEncode(BitStream data)
{
    auto size = data.Bits();

    if (size < 2)
    {
        return data;
    }

    unsigned bitsReadCount = 1;
    data.Reset();
    auto previous = data.Read(1);
    BitStream result;

    while (bitsReadCount < size)
    {
        auto current = data.Read(1);
        bitsReadCount++;

        unsigned run = 1;

        while (current == previous)
        {
            run++;

            if (bitsReadCount == size)
            {
                break;
            }

            bitsReadCount++;
            current = data.Read(1);
        }

        result.Write(previous, 1);

        if (run > 1)
        {
            while (run > 1)
            {
                result.Write(previous, 1);

                run -= 2;
                unsigned expRun = 0;
                unsigned ordinal = 0;
                for (const auto& exp : exponentialBitLevelRunLengthLookup)
                {
                    if (run >= exp)
                    {
                        expRun = ordinal;
                        ordinal++;
                    }
                }

                for (unsigned i = 0; i < expRun; ++i)
                {
                    result.Write(previous, 1);
                }

                result.Write(1 - previous, 1);

                run -= exponentialBitLevelRunLengthLookup[expRun];

                if (run > 0)
                {
                    result.Write(previous, 1);
                }

                if ((run < 2) && (bitsReadCount == size) && (current != previous))
                {
                    result.Write(current, 1);
                }
            }
        }
        else
        {
            if (bitsReadCount == size)
            {
                result.Write(current, 1);
            }
        }

        if ((bitsReadCount - 1) == size)
        {
            break;
        }

        previous = current;
    }

    return result;
}

BitStream ExponentialBitLevelRunLengthDecode(BitStream& data, unsigned targetBits = 0)
{
    auto size = data.Bits();

    if (size < 2)
    {
        return data;
    }

    unsigned bitsReadCount = 1;
    data.Reset();
    auto previous = data.Read(1);
    BitStream result;

    while (bitsReadCount != size)
    {
        if (targetBits > 0)
        {
            if (result.Bits() == targetBits)
            {
                break;
            }

            if (result.Bits() == (targetBits - 1))
            {
                result.Write(previous, 1);
                break;
            }
        }

        auto current = data.Read(1);
        bitsReadCount++;

        unsigned run = 1;

        while (current == previous)
        {
            run++;

            if (bitsReadCount == size)
            {
                break;
            }

            bitsReadCount++;
            current = data.Read(1);
        }

        if (run > 1)
        {
            result.Write(previous, 1);
            result.Write(previous, 1);

            run -= 2;
            assert(run < exponentialBitLevelRunLengthLookup.size());
            auto runCount = exponentialBitLevelRunLengthLookup[run];

            while (runCount--)
            {
                result.Write(previous, 1);
            }

            if (targetBits > 0)
            {
                if (result.Bits() == targetBits)
                {
                    continue;
                }
            }

            if (bitsReadCount < size)
            {
                current = data.Read(1);
                bitsReadCount++;

                if (bitsReadCount == size)
                {
                    result.Write(current, 1);
                }

                if (targetBits > 0)
                {
                    if (result.Bits() == (targetBits - 1))
                    {
                        result.Write(current, 1);
                    }
                }
            }
        }
        else
        {
            result.Write(previous, 1);

            if (bitsReadCount == size)
            {
                result.Write(current, 1);
            }
        }

        previous = current;
    }

    return result;
}

// //////////////////////////////////////////////////////

static const unsigned minRun = 8;
static const unsigned runCountBits = 7;
static const unsigned maxUnrun = 1 + ((1 << (runCountBits - 1)) - 1);
static const unsigned maxRun = minRun + ((1 << (runCountBits - 1)) - 1);

BitStream BitBitPackFullEncode(BitStream data)
{
    auto size = data.Bits();

    if (size < 2)
    {
        return data;
    }

    BitStream result;
    unsigned startIndex = 0;

    auto BitRead = [&data](unsigned index) -> auto
    {
        data.SetOffset(index);
        return data.Read(1);
    };

    while (startIndex < size)
    {
        unsigned i = 1;
        unsigned run = 1;

        while (i != maxUnrun)
        {
            auto index = startIndex + i;
            if (BitRead(index) == BitRead(index - 1))
            {
                run++;
            }
            else
            {
                run = 1;
            }

            if (run == (minRun + 1))
            {
                break;
            }

            ++i;

            if ((startIndex + i) == size)
            {
                break;
            }
        }

        unsigned unRunCount = 1 + (i - run);

        if (unRunCount)
        {
            auto header = ZigZag(static_cast<int32_t>(unRunCount) - 1);

            result.Write(header, runCountBits);

            unsigned index = 0;
            while (index != unRunCount)
            {
                result.Write(BitRead(startIndex + index),1);

                index++;
            }
        }

        startIndex += unRunCount;

        if (run == (minRun + 1))
        {
            while (run < (maxRun - minRun))
            {
                auto index = startIndex + run;

                if ((startIndex + run) == size)
                {
                    break;
                }

                if (BitRead(index) != BitRead(index - 1))
                {
                    break;
                }

                ++run;
            }

            int h = -(static_cast<int32_t>(run - minRun));
            auto header = ZigZag(h);

            result.Write(header, runCountBits);
            result.Write(BitRead(startIndex), 1);

            startIndex += run;
        }

        if ((startIndex + 1) == size)
        {
            auto header = ZigZag(static_cast<int32_t>(0));

            result.Write(header, runCountBits);
            result.Write(BitRead(startIndex), 1);
            startIndex++;
        }
    }

    return result;
}


BitStream BitBitPackFullDecode(
    BitStream data,
    unsigned& bytesConsumed,
    unsigned maxResultBytes = 0)
{
    auto size = data.Bits();

    if (size < 2)
    {
        return data;
    }

    auto BitRead = [&data](unsigned index, unsigned count = 1) -> auto
    {
        data.SetOffset(index);
        return data.Read(count);
    };

    BitStream result;
    unsigned index = 0;
    while (index < size)
    {
        auto rawHeader = BitRead(index, runCountBits);
        auto header = ZigZag(static_cast<uint32_t>(rawHeader));
        index += runCountBits;

        if (header >= 0)
        {
            auto count = header + 1;

            while (count--)
            {
                result.Write(BitRead(index++),1);
            }
        }
        else
        {
            unsigned count = -header + minRun;

            while (count--)
            {
                result.Write(BitRead(index),1);
            }

            index++;
        }

        if (maxResultBytes > 0)
        {
            if (result.Bits() == maxResultBytes)
            {
                break;
            }
        }
    }

    bytesConsumed = index;

    return result;
}

BitStream BitBitPackFullDecode(BitStream data)
{
    unsigned unused;
    return BitBitPackFullDecode(data, unused);
}

// //////////////////////////////////////////////////////

static const unsigned long expMaxUnrun = 5;
static const unsigned allSet = (1 << expMaxUnrun) - 1;

BitStream BitBitPackEncode(BitStream data)
{
    auto size = data.Bits();

    if (size < 2)
    {
        return data;
    }

    data.Reset();
    BitStream result;

    unsigned count = size;
    unsigned read = 0;

    while (count > (expMaxUnrun - 1))
    {
        auto chunk = data.Read(expMaxUnrun);
        count -= expMaxUnrun;
        read += expMaxUnrun;

        if ((chunk == allSet) || (chunk == 0))
        {
            auto previous = chunk & 1;
            auto current = data.Read(1);
            auto run = expMaxUnrun;

            while (((previous == current) && (read != size)) && (run < 255))
            {
                run++;
                current = data.Read(1);
                count--;
                read++;
            }

            if (read != size)
            {
                data.SetOffset(data.Bits() - 1);
            }

            result.Write(1,1);
            result.Write(previous, 1);
            result.Write(run, 8);
        }
        else
        {
            result.Write(0, 1);
            result.Write(chunk, expMaxUnrun);
        }
    }

    if (count)
    {
        auto theRest = data.Read(count);
        result.Write(0, 1);
        result.Write(theRest, count);
        result.Write(0, expMaxUnrun - count);
    }

    return result;
}

BitStream BitBitPackDecode(BitStream& data, unsigned targetBits = 0)
{
    auto size = data.Bits();

    if (size < 5)
    {
        return data;
    }

    unsigned bitsReadCount = 0;
    data.Reset();
    BitStream result;

    while (bitsReadCount < size)
    {
        auto tag = data.Read(1);
        bitsReadCount++;

        if (tag)
        {
            auto repeat = data.Read(1);
            bitsReadCount++;

            auto run = data.Read(8);
            bitsReadCount += 8;

            while(run--)
            {
                result.Write(repeat, 1);
            }
        }
        else
        {
            auto toWrite = targetBits ? targetBits - result.Bits() : expMaxUnrun;
            toWrite = std::min(toWrite, expMaxUnrun);

            auto unRun = data.Read(expMaxUnrun);
            bitsReadCount += expMaxUnrun;
            result.Write(unRun, toWrite);
        }

        if (targetBits)
        {
            if (targetBits == result.Bits())
            {
                break;
            }
        }
    }

    return result;
}

// //////////////////////////////////////////////////////

void RunLengthTests()
{
    std::vector<std::function<void(BitStream)>> tests;

    tests.push_back([](BitStream data)
    {
        auto encoded = RunLengthEncode(data.Data());
        auto decoded = RunLengthDecode(encoded);

        assert(data.Data() == decoded);
    });

    tests.push_back([](BitStream data)
    {
        auto encoded = BitPackEncode(data.Data());
        auto decoded = BitPackDecode(encoded);

        assert(data.Data() == decoded);
    });

    tests.push_back([](BitStream data)
    {
        auto encoded = ExponentialBitLevelRunLengthEncode(data);
        auto decoded = ExponentialBitLevelRunLengthDecode(encoded);

        assert(data == decoded);
    });

    tests.push_back([](BitStream data)
    {
        auto bits = data.Bits();
        auto encoded = BitBitPackEncode(data);
        auto decoded = BitBitPackDecode(encoded, bits);

        assert(data == decoded);
    });

    tests.push_back([](BitStream data)
    {
        unsigned unused;
        auto bits = data.Bits();
        auto encoded = BitBitPackFullEncode(data);
        auto decoded = BitBitPackFullDecode(encoded, unused, bits);

        assert(data == decoded);
    });

    std::vector<BitStream> testData;

    for (uint8_t i = 0; i < 32; ++i)
    {
        testData.push_back(BitStream(ByteVector{i}, 5));
    }

    testData.push_back(BitStream(ByteVector{0, 248, 11}, (2 * 8) + 5));

    testData.push_back(BitStream(
                        ByteVector{1,2,3,4,5,0,0},
                        7 * 8));

    testData.push_back(BitStream(
                        ByteVector{0,1,3,3,0,0,0,0,0,5,6,6,6,5,4,3,3,3,3,4,},
                        20 * 8));

    testData.push_back(BitStream(
                        ByteVector{0,1,3,3,0,0,0,0,0,5,6,6,6,5,4,3,3,3,3,},
                        19 * 8));

    testData.push_back(BitStream(
                        ByteVector(100, 0),
                        100 * 8));

    testData.push_back(BitStream(
                        ByteVector(300, 5),
                        300 * 8));

    for (auto tester : tests)
    {
        for (auto data : testData)
        {
            tester(data);
        }
    }
}

// //////////////////////////////////////////////////////

struct Config
{
    ChangedArrayEncoding    rle;
    PosVector3Packer        posPacker;
};

std::vector<uint8_t> Encode(
    const Frame& base,
    const Frame& target,
    unsigned frameDelta,
    Stats& stats,
    Config config = {ChangedArrayEncoding::None, PosVector3Packer::None})
{
    const auto count = base.size();

    assert(count > 0);
    assert(count == target.size());

    // ////////////////////////////////
    // A. If nothing changed, don't send anything at all.

    bool same = true;
    size_t firstChanged = 0;
    for (;firstChanged < count; ++firstChanged)
    {
        if (base[firstChanged] != target[firstChanged])
        {
            same = false;
            break;
        }
    }

    if (same)
    {
        return {};
    }

    // ////////////////////////////////

    BitStream result;

    BitStream changed;
    BitStream deltas;

    unsigned runLength0 = 0;
    unsigned runLength1 = 0;

    for (size_t i = 0; i < count; ++i)
    {
        if (base[i] == target[i])
        {
            changed.Write(0, 1);
            ++runLength0;

            ++(stats.changed1DistanceRunHistogram[runLength1]);
            runLength1 = 0;
        }
        else
        {
            changed.Write(1, 1);
            ++runLength1;

            ++(stats.changed0DistanceRunHistogram[runLength0]);
            runLength0 = 0;

            assert(target[i].orientation_a >= 0);
            assert(target[i].orientation_b >= 0);
            assert(target[i].orientation_c >= 0);

            // Get stats on:
            // If the quat changes but the position doesn't
            // if the position changes but the quat doesn't
            // largest distance change for x,y,z,total
            bool posSame =
                    (base[i].position_x == target[i].position_x) &&
                    (base[i].position_y == target[i].position_y) &&
                    (base[i].position_z == target[i].position_z);

            bool quatSame =
                    (base[i].orientation_largest == target[i].orientation_largest) &&
                    (base[i].orientation_a == target[i].orientation_a) &&
                    (base[i].orientation_b == target[i].orientation_b) &&
                    (base[i].orientation_c == target[i].orientation_c);

            if (posSame && quatSame)
            {
                ++stats.wtfBothSame;
            }

            if (!posSame && quatSame)
            {
                ++stats.quatChangedNoPos;
            }

            if (posSame && !quatSame)
            {
                ++stats.PosChangedNotQuat;
            }

            ++stats.changed;

            if (!posSame)
            {
                auto dx = abs(base[i].position_x - target[i].position_x);
                auto dy = abs(base[i].position_y - target[i].position_y);
                auto dz = abs(base[i].position_z - target[i].position_z);

                stats.deltaX.Update(dx);
                stats.deltaY.Update(dy);
                stats.deltaZ.Update(dz);

                stats.deltaTotal.Update(dx + dy + dz);
            }

            // Only write what's changed
            if (quatSame)
            {
                deltas.Write(0, 1);
            }
            else
            {
                deltas.Write(1, 1);
                deltas.Write(target[i].orientation_largest, RotationIndexMaxBits);
                deltas.Write(target[i].orientation_a, RotationMaxBits);
                deltas.Write(target[i].orientation_b, RotationMaxBits);
                deltas.Write(target[i].orientation_c, RotationMaxBits);
            }

            if (posSame)
            {
                deltas.Write(0, 1);
            }
            else
            {
                deltas.Write(1, 1);
                stats.posDeltaUnpackedBitCount += ((MaxBitsXY * 2) + MaxBitsZ);

                switch (config.posPacker)
                {
                    default:
                    case PosVector3Packer::None:
                    {
                        deltas.Write(ZigZag(target[i].position_x), MaxBitsXY);
                        deltas.Write(ZigZag(target[i].position_y), MaxBitsXY);
                        deltas.Write(target[i].position_z, MaxBitsZ);
                        break;
                    }

                    case PosVector3Packer::BitVector3:
                    case PosVector3Packer::BitVector3_2BitExpPrefix:
                    case PosVector3Packer::BitVector3Truncated:
                    {
                        // GLEE, actually writing a delta this time!
                        auto vec = IntVec3
                        {
                            target[i].position_x - base[i].position_x,
                            target[i].position_y - base[i].position_y,
                            target[i].position_z - base[i].position_z,
                        };

                        auto maxMagnitude = MaxPositionChangePerSnapshot * frameDelta;
                        auto mag = sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
                        assert(mag < maxMagnitude);

                        unsigned bitsWritten;
                        switch (config.posPacker)
                        {
                            default:
                            case PosVector3Packer::BitVector3:
                            {
                                bitsWritten = BitVector3Encode(
                                    vec,
                                    maxMagnitude,
                                    deltas);

                                break;
                            }

                            case PosVector3Packer::BitVector3_2BitExpPrefix:
                            {
                                bitsWritten = BitVector3Encode2BitExpPrefix(
                                    vec,
                                    maxMagnitude,
                                    deltas);

                                break;
                            }

                            case PosVector3Packer::BitVector3Truncated:
                            {
                                bitsWritten = BitVector3TruncatedEncode(
                                    vec,
                                    maxMagnitude,
                                    deltas);

                                break;
                            }
                        }

                        assert(bitsWritten < ((MaxBitsXY * 2) + MaxBitsZ));
                        stats.posDeltaPackedBitCount.Update(bitsWritten);
                        break;
                    }
                }
            }

            deltas.Write(target[i].interacting, 1);
        }
    }

    ++(stats.changed0DistanceRunHistogram[runLength0]);
    ++(stats.changed1DistanceRunHistogram[runLength1]);

    // Oooh, finally, some research!
    auto rle = RunLengthEncode(changed.Data());
    auto bitpack = BitPackEncode(changed.Data());
    auto bitexprle = ExponentialBitLevelRunLengthEncode(changed);
    auto expbitpack = BitBitPackEncode(changed);
    auto bitbitpackfull = BitBitPackFullEncode(changed);

    assert((rle.size() * 8) < Cubes);
    assert((bitpack.size() * 8) < Cubes);
    assert((bitexprle.Bits()) < Cubes);
    assert((expbitpack.Bits()) < Cubes);
    assert((bitbitpackfull.Bits()) < Cubes);

    stats.rle.Update(rle.size() * 8);
    stats.bitpack.Update(bitpack.size() * 8);
    stats.bitexprle.Update(bitexprle.Bits());
    stats.bitbitpack.Update(expbitpack.Bits());
    stats.bitbitfullpack.Update(bitbitpackfull.Bits());

    auto changedCompressed = [&changed, &config]()
    {
        switch (config.rle)
        {
            default:
            case ChangedArrayEncoding::None:
            {
                return changed;
            }

            case ChangedArrayEncoding::Rle:
            {
                auto result = RunLengthEncode(changed.Data());
                return BitStream(result, result.size() * 8);
            }

            case ChangedArrayEncoding::BitPack:
            {
                auto result = BitPackEncode(changed.Data());
                return BitStream(result, result.size() * 8);
            }

            case ChangedArrayEncoding::Exp:
            {
                return ExponentialBitLevelRunLengthEncode(changed);
            }

            case ChangedArrayEncoding::BitBitPack:
            {
                return BitBitPackEncode(changed);
            }

            case ChangedArrayEncoding::BitBitFullPack:
            {
                return BitBitPackFullEncode(changed);
            }
        }
    }();

    assert(changedCompressed.Bits() <= 901);

    result.Write(changedCompressed);
    result.Write(deltas);

    result.TrimZerosFromBack();

    return result.Data();
}

Frame Decode(
    const Frame& base,
    std::vector<uint8_t>& buffer,
    unsigned frameDelta,
    Config config = {ChangedArrayEncoding::None, PosVector3Packer::None})
{
    // A.
    if (buffer.empty())
    {
        return base;
    }

    BitStream bits(buffer, buffer.size() * 8);
    Frame result;

    auto changed = [&bits, &config]()
    {
        switch (config.rle)
        {
            default:
            case ChangedArrayEncoding::None:
            {
                bits.Reset();
                return bits.ReadArray(Cubes);
            }

            case ChangedArrayEncoding::Rle:
            {
                unsigned bytesUsed = 0;
                auto decode = RunLengthDecode(bits.Data(), bytesUsed, 113);

                bits.SetOffset(bytesUsed * 8);

                return BitStream(decode);
            }

            case ChangedArrayEncoding::BitPack:
            {
                unsigned bytesUsed = 0;
                auto decode = BitPackDecode(bits.Data(), bytesUsed, 113);

                bits.SetOffset(bytesUsed * 8);

                return BitStream(decode);
            }

            case ChangedArrayEncoding::Exp:
            {
                return ExponentialBitLevelRunLengthDecode(bits, Cubes);
            }

            case ChangedArrayEncoding::BitBitPack:
            {
                return BitBitPackDecode(bits, Cubes);
            }

            case ChangedArrayEncoding::BitBitFullPack:
            {
                unsigned bitsUsed = 0;
                auto decode = BitBitPackFullDecode(bits, bitsUsed, Cubes);

                bits.SetOffset(bitsUsed);

                return decode;
            }
        }
    }();

    changed.Reset();

    for (size_t i = 0; i < Cubes; ++i)
    {
        if (changed.Read(1))
        {
            bool quadChanged = bits.Read(1);

            if (quadChanged)
            {
                result[i].orientation_largest   = bits.Read(RotationIndexMaxBits);
                result[i].orientation_a         = bits.Read(RotationMaxBits);
                result[i].orientation_b         = bits.Read(RotationMaxBits);
                result[i].orientation_c         = bits.Read(RotationMaxBits);
            }
            else
            {
                result[i].orientation_largest   = base[i].orientation_largest;
                result[i].orientation_a         = base[i].orientation_a;
                result[i].orientation_b         = base[i].orientation_b;
                result[i].orientation_c         = base[i].orientation_c;
            }

            bool posChanged = bits.Read(1);

            if (posChanged)
            {
                switch (config.posPacker)
                {
                    default:
                    case PosVector3Packer::None:
                    {
                        result[i].position_x = ZigZag(bits.Read(MaxBitsXY));
                        result[i].position_y = ZigZag(bits.Read(MaxBitsXY));
                        result[i].position_z = bits.Read(MaxBitsZ);
                        break;
                    }

                    case PosVector3Packer::BitVector3:
                    case PosVector3Packer::BitVector3_2BitExpPrefix:
                    case PosVector3Packer::BitVector3Truncated:
                    {
                        auto maxMagnitude = MaxPositionChangePerSnapshot * frameDelta;

                        IntVec3 vec;

                        switch (config.posPacker)
                        {
                            default:
                            case PosVector3Packer::BitVector3:
                            {
                                vec = BitVector3Decode(
                                    maxMagnitude,
                                    bits);

                                break;
                            }

                            case PosVector3Packer::BitVector3_2BitExpPrefix:
                            {
                                vec = BitVector3Decode2BitExpPrefix(
                                    maxMagnitude,
                                    bits);

                                break;
                            }

                            case PosVector3Packer::BitVector3Truncated:
                            {
                                vec = BitVector3TruncatedDecode(
                                    maxMagnitude,
                                    bits);

                                break;
                            }
                        }

                        auto mag = sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
                        assert(mag < maxMagnitude);

                        result[i].position_x = vec.x + base[i].position_x;
                        result[i].position_y = vec.y + base[i].position_y;
                        result[i].position_z = vec.z + base[i].position_z;

                        break;
                    }
                }
            }
            else
            {
                result[i].position_x            = base[i].position_x;
                result[i].position_y            = base[i].position_y;
                result[i].position_z            = base[i].position_z;
            }

            result[i].interacting           = bits.Read(1);
        }
        else
        {
            result[i] = base[i];
        }
    }

    return result;
}

// //////////////////////////////////////////////////////

int main(int, char**)
{
    TestTruncate();
    BitVector3Tests();
    RunLengthTests();
    BitStreamTest();
    ZigZagTest();

    // //////////////////////////////////////////////////////

    auto frames = []() -> std::vector<Frame>
    {
        const auto fileHandle = fopen("delta_data.bin", "rb");
        if (!fileHandle)
        {
            printf("ERROR: Cannot open 'delta_data.bin' for reading.\n");
            return {};
        }
        SCOPED_EXIT(fclose(fileHandle));

        fseek(fileHandle, 0, SEEK_END);
        const auto size = ftell(fileHandle);
        fseek(fileHandle, 0, SEEK_SET);

        const auto frameCount = size / sizeof(Frame);

        // I really hate how I cannot make an array without initilising it first.
        std::vector<Frame> frames(frameCount);

        fread(
            frames.data(),
            sizeof(Frame),
            frameCount,
            fileHandle);

        return frames;
    }();

    if (frames.empty())
    {
        return 1;
    }

    // Lets do some research
    auto size = frames.size();
    Stats stats
    {
        0,
        0,
        0,
        0,
        0,
        {0},
        {0},
        {Cubes,0,0,0},
        {Cubes,0,0,0},
        {Cubes,0,0,0},
        {Cubes,0,0,0},
        {Cubes,0,0,0},
        0,
        0,
        0,
        0,
        {static_cast<unsigned>(MaxPositionChangePerSnapshot * 20),0,0,0},
        {static_cast<unsigned>(MaxPositionChangePerSnapshot * 20),0,0,0},
        {static_cast<unsigned>(MaxPositionChangePerSnapshot * 20),0,0,0},
        {static_cast<unsigned>(MaxPositionChangePerSnapshot * 20),0,0,0},
        0,
        {static_cast<unsigned>(MaxPositionChangePerSnapshot * 20),0,0,0},
    };

    for (size_t i = FirstBase; i < size; ++i)
    {
         DoSomeStats(frames[i-FirstBase], frames[i], stats);
    }

    #define PRINT_INT(x) printf("%-32s\t%d\n", #x, x);

    PRINT_INT(stats.itemCount)
    PRINT_INT(stats.notChangedCount)
    PRINT_INT(stats.interactingTotal)
    PRINT_INT(stats.interactingNotChanged)
    PRINT_INT(stats.notInteractingNotChanged)

    auto percentUnchanged = (100.0f * stats.notChangedCount) / stats.itemCount;
    auto percentIUnchanged = (100.0f * stats.interactingNotChanged) / stats.interactingTotal;
    auto percentSUnchanged = (100.0f * stats.notInteractingNotChanged) / (stats.itemCount - stats.interactingTotal);

    #define PRINT_FLOAT(x) printf("%-32s\t%f\n", #x, x);

    PRINT_FLOAT(percentUnchanged)
    PRINT_FLOAT(percentIUnchanged)
    PRINT_FLOAT(percentSUnchanged)

    // Lets actually do the stuff.
    unsigned bytes = 0;
    unsigned packetsCoded = 0;
    Config config
    {
        ChangedArrayEncoding::Exp,
        PosVector3Packer::BitVector3Truncated,
    };

    for (size_t i = FirstBase; i < size; ++i)
    {
        auto buffer = Encode(
            frames[i-FirstBase],
            frames[i],
            FirstBase,
            stats,
            config);

        bytes += buffer.size();

        auto back = Decode(
            frames[i-FirstBase],
            buffer,
            FirstBase,
            config);

        assert(back == frames[i]);

        packetsCoded++;
    }

    float packetSizeAverge = ((float) bytes) / (size - FirstBase);
    float bytesPerSecondAverage = packetSizeAverge * 60.0f;
    float kbps = bytesPerSecondAverage * 8 / 1000.0f;

    PRINT_INT(bytes)
    PRINT_INT(packetsCoded)
    PRINT_FLOAT(packetSizeAverge)
    PRINT_FLOAT(bytesPerSecondAverage)
    PRINT_FLOAT(kbps)

    float changedBitsTotal  = (901.0f * packetsCoded);
    float rle               = 100 * (stats.rle.sum / changedBitsTotal);
    float bitpack           = 100 * (stats.bitpack.sum / changedBitsTotal);
    float bitexprle         = 100 * (stats.bitexprle.sum / changedBitsTotal);
    float expbitpack        = 100 * (stats.bitbitpack.sum / changedBitsTotal);
    float bitbitpackfull    = 100 * (stats.bitbitfullpack.sum / changedBitsTotal);

    PRINT_FLOAT(rle)
    PRINT_INT(stats.rle.min)
    PRINT_INT(stats.rle.max)
    PRINT_FLOAT(bitpack)
    PRINT_INT(stats.bitpack.min)
    PRINT_INT(stats.bitpack.max)
    PRINT_FLOAT(bitexprle)
    PRINT_INT(stats.bitexprle.min)
    PRINT_INT(stats.bitexprle.max)
    PRINT_FLOAT(expbitpack)
    PRINT_INT(stats.bitbitpack.min)
    PRINT_INT(stats.bitbitpack.max)
    PRINT_FLOAT(bitbitpackfull)
    PRINT_INT(stats.bitbitfullpack.min)
    PRINT_INT(stats.bitbitfullpack.max)

    for (const auto h : stats.changed0DistanceRunHistogram)
    {
        printf("%d, ", h);
    }
    printf("\n");

    for (const auto h : stats.changed1DistanceRunHistogram)
    {
        printf("%d, ", h);
    }
    printf("\n");
    printf("\n");
    auto dxa =  stats.deltaX.Average();
    auto dya =  stats.deltaY.Average();
    auto dza =  stats.deltaZ.Average();
    auto dta =  stats.deltaTotal.Average();
    PRINT_FLOAT(dxa)
    PRINT_INT(stats.deltaX.min)
    PRINT_INT(stats.deltaX.max)
    PRINT_FLOAT(dya)
    PRINT_INT(stats.deltaY.min)
    PRINT_INT(stats.deltaY.max)
    PRINT_FLOAT(dza)
    PRINT_INT(stats.deltaZ.min)
    PRINT_INT(stats.deltaZ.max)
    PRINT_FLOAT(dta)
    PRINT_INT(stats.deltaTotal.min)
    PRINT_INT(stats.deltaTotal.max)
    PRINT_FLOAT(MaxPositionChangePerSnapshot*6)
    PRINT_INT(stats.quatChangedNoPos)
    PRINT_INT(stats.PosChangedNotQuat)
    PRINT_INT(stats.wtfBothSame)
    PRINT_INT(stats.changed)

    auto posPackedAvg = stats.posDeltaPackedBitCount.Average();
    PRINT_FLOAT(posPackedAvg)
    PRINT_INT(stats.posDeltaPackedBitCount.min)
    PRINT_INT(stats.posDeltaPackedBitCount.max)

    return 0;
}
