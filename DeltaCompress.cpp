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
#include <algorithm>
#include <functional>
#include <cmath>

#include "Range_coding.hpp"

// //////////////////////////////////////////////////////

bool doTests            = false;
bool doStats            = true;
bool doCompression      = true;
bool doRangeCompression = false;
bool doRangeSearch      = false;

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
            (lhs.interacting            == rhs.interacting)
        &&  (lhs.orientation_a          == rhs.orientation_a)
        &&  (lhs.orientation_b          == rhs.orientation_b)
        &&  (lhs.orientation_c          == rhs.orientation_c)
        &&  (lhs.orientation_largest    == rhs.orientation_largest)
        &&  (lhs.position_x             == rhs.position_x)
        &&  (lhs.position_y             == rhs.position_y)
        &&  (lhs.position_z             == rhs.position_z)
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

// With this technique I've found that minimum sufficient precision for
// my simulation is 9 bits per-smallest component. This gives a result
// of 2 + 9 + 9 + 9 = 29 bits per-orientation (originally 128!).

// NOTE: All quanternion values are positive even though they are stored
//       as ints.

static const unsigned RotationMaxBits = 9;
static const unsigned RotationIndexMaxBits = 2;

// I found a maximum speed of 32 meters per-second is a nice power of two
// and doesn't negatively affect the player experience in the cube simulation.
static const unsigned MaxSpeedMetersPerSecond = 32;

// Now we have only position to compress. We'll use the same trick that we used
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
//static const unsigned   MinZInclusize     = 0;
static const int        MaxXYInclusize    = 131071;  // (256 * 512) - 1
static const int        MinXYInclusize    = -131072; // (-256 * 512)
static const int        XYRange           = MaxXYInclusize - MinXYInclusize;

static const unsigned   MaxBitsZ          = 14;
static const unsigned   MaxBitsXY         = 18;

static const unsigned   MaxSnapshotsPerSecond = 60;

static const unsigned   PacketDelta           = 6;

static const unsigned   QuanternionUncompressedBitThreshold = 8;

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
};

float Average(const MinMaxSum& mms)
{
    auto count = mms.count;

    if (!count)
    {
        return 0;
    }

    return mms.sum / static_cast<float>(count);
}

struct Stats
{
    MinMaxSum bytesPerPacket;

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
    MinMaxSum range_simple;
    MinMaxSum range_smarter;
    MinMaxSum range_simple_adaptive;
    MinMaxSum range_smarter_adaptive;

    unsigned changed;
    unsigned quatChangedNoPos;
    unsigned PosChangedNotQuat;
    unsigned wtfBothSame;
    MinMaxSum deltaX;
    MinMaxSum deltaY;
    MinMaxSum deltaZ;
    MinMaxSum deltaTotal;

    unsigned quatChanged;
    unsigned largestDifferent;
    MinMaxSum deltaA;
    MinMaxSum deltaB;
    MinMaxSum deltaC;
    MinMaxSum deltaABC;

    unsigned posDeltaUnpackedBitCount;
    MinMaxSum posDeltaPackedBitCount;

    unsigned quatDeltaUnpackedBitCount;
    MinMaxSum quatDeltaPackedBitCount;

    unsigned bitStream0;
    unsigned bitStream1;
    unsigned bitStream00;
    unsigned bitStream01;
    unsigned bitStream10;
    unsigned bitStream11;

    unsigned zeros;
    unsigned ones;
    unsigned one_one;
    unsigned zero_zero;
    unsigned one_zero;
    unsigned zero_one;

    MinMaxSum rotor;
    MinMaxSum rotor_bits;
    unsigned rotor_bits_8;
    unsigned rotor_bits_9;
    unsigned rotor_bits_10;
    unsigned rotor_bits_11;
    unsigned rotor_bits_12;
    unsigned rotor_bits_13;
    unsigned rotor_bits_14;
    unsigned rotor_bits_15;
    unsigned rotor_bits_wtf;
};

enum class ChangedArrayEncoding
{
    None,
    Rle,
    BitPack,
    Exp,
    BitBitPack,
    BitBitFullPack,
    RangeSimple,
    RangeSmarter,
    RangeSimpleAdaptive,
    RangeSmarterAdaptive,
};

enum class PosVector3Packer
{
    None,
    BitVector3,
    BitVector3_2BitExpPrefix,
    BitVector3Truncated,
    BitVector3Sorted,
    Sorted_no_bit_count,
};

enum class QuatPacker
{
    None,
    BitVector3Unrelated,
    BitVector3BitCount,
    BitVector3ModifiedZigZag,
    BitVector3ModifiedGaffer,
    Gaffer,
    BitVector3Sorted,
    Sorted_no_bit_count,
};

// //////////////////////////////////////////////////////

struct Quat2
{
    std::array<float, 4> q;

    float constexpr operator[](unsigned index) const
    {
        return q[index];
    }

    float& operator[](unsigned index)
    {
        return q[index];
    }
};

float* begin(Quat2& q)
{
    return &(q.q[0]);
}

float* end(Quat2& q)
{
    return &(q.q[4]);
}

struct Rotor
{
    std::array<float, 3> r;

    float constexpr operator[](unsigned index) const
    {
        return r[index];
    }
};

float* begin(Rotor& r)
{
    return &(r.r[0]);
}

float* end(Rotor& r)
{
    return &(r.r[4]);
}

float constexpr Magnitude_squared(const Quat2& q)
{
    return
        q[0] * q[0] +
        q[1] * q[1] +
        q[2] * q[2] +
        q[3] * q[3];
}

float constexpr Magnitude_squared(const Rotor& r)
{
    return
        r[0] * r[0] +
        r[1] * r[1] +
        r[2] * r[2];
}

Quat2 constexpr Mul(const Quat2& lhs, float rhs)
{
    return
    {
        lhs[0]*rhs,
        lhs[1]*rhs,
        lhs[2]*rhs,
        lhs[3]*rhs
    };
}

Rotor constexpr Mul(const Rotor& lhs, float rhs)
{
    return
    {
        lhs[0]*rhs,
        lhs[1]*rhs,
        lhs[2]*rhs
    };
}

Quat2 constexpr Mul(const Quat2& lhs, const Quat2& rhs)
{
    return
    {
        (lhs[0]*rhs[0] - lhs[1]*rhs[1] - lhs[2]*rhs[2] - lhs[3]*rhs[3]),
        (lhs[0]*rhs[1] + lhs[1]*rhs[0] + lhs[2]*rhs[3] - lhs[3]*rhs[2]),
        (lhs[0]*rhs[2] - lhs[1]*rhs[3] + lhs[2]*rhs[0] + lhs[3]*rhs[1]),
        (lhs[0]*rhs[3] + lhs[1]*rhs[2] - lhs[2]*rhs[1] + lhs[3]*rhs[0])
    };
}

Quat2 Normalise(const Quat2& q)
{
    return Mul(q, 1.0f / sqrt(Magnitude_squared(q)));
}

const float EPISLON_U16 = 1.0f / 65536.5f;

Quat2 Make_largest_positive(const Quat2& q)
{
    float largest = 0;
    for (const auto qi : q.q)
    {
        if (std::abs(qi) > (largest + EPISLON_U16))
        {
            largest = qi;
        }
    }

    if (largest < 0)
    {
        return
        {
            -q[0],
            -q[1],
            -q[2],
            -q[3],
        };
    }

    return q;
}

Quat2 constexpr q_star(const Quat2& q)
{
    return
    {
        q[0],
        -q[1],
        -q[2],
        -q[3]
    };
}

// http://www.geomerics.com/blogs/quaternions-rotations-and-compression/
Quat2 constexpr R(const Quat2& base, const Quat2& target)
{
    return Mul(target, q_star(base));
}

Rotor constexpr to_rotor(const Quat2& r)
{
    return
    {
        r[1] / (1.0f + r[0]),
        r[2] / (1.0f + r[0]),
        r[3] / (1.0f + r[0]),
    };
}

Quat2 constexpr to_quat(const Rotor& b)
{
    return
    {
        (1.0f - Magnitude_squared(b)) / (1 + Magnitude_squared(b)),
        (b[0] * 2.0f) / (1 + Magnitude_squared(b)),
        (b[1] * 2.0f) / (1 + Magnitude_squared(b)),
        (b[2] * 2.0f) / (1 + Magnitude_squared(b)),
    };
}

bool Equal(const Quat2& lhs, const Quat2& rhs, float epislon = EPISLON_U16)
{
    const auto size = lhs.q.size();
    for (unsigned i = 0; i < size; ++i)
    {
        if ((lhs[i] - rhs[i]) >= epislon)
        {
            return false;
        }
        if ((rhs[i] - lhs[i]) >= epislon)
        {
            return false;
        }
    }

    return true;
}

void assert_float_eq(float a, float b, float epislon = EPISLON_U16)
{
    assert((a - b) < epislon);
    assert((b - a) < epislon);
}

void Quat_tests()
{
    for (float i = -10; i < 20; ++i)
    {
        for (float j  = -10; j < 15; ++j)
        {
            for (float k = -5; k < 16; ++k)
            {
                for (float l = -55; l < 10; ++l)
                {                    // Skip malformed quats (vector == 0)
                    if (!j)
                    {
                        continue;
                    }

                    auto base = Normalise(Quat2{i,j,k,l});
                    auto mag_squared = Magnitude_squared(base);

                    assert_float_eq(mag_squared, 1.0f);

                    auto target =
                        Make_largest_positive(
                            Normalise(
                                Quat2
                                {
                                    i,
                                    k,
                                    j,
                                    j
                                }));

                    // Skip if the same
                    if (Equal(target, Mul(base, -1.0f)))
                    {
                        continue;
                    }
                    if (Equal(target, base))
                    {
                        continue;
                    }

                    auto difference = R(base, target);
                    auto b = to_rotor(difference);

                    auto d2 = to_quat(b);
                    auto result = Mul(d2, base);

                    auto should_be_identity = R(result, target);
                    // get rid of sign bit.
                    should_be_identity.q[0] *= should_be_identity.q[0];

                    assert_float_eq(should_be_identity[0], 1.0f);
                    assert_float_eq(should_be_identity[1], 0.0f);
                    assert_float_eq(should_be_identity[2], 0.0f);
                    assert_float_eq(should_be_identity[3], 0.0f);

                    result = Make_largest_positive(result);

                    assert_float_eq(result[0], target[0]);
                    assert_float_eq(result[1], target[1]);
                    assert_float_eq(result[2], target[2]);
                    assert_float_eq(result[3], target[3]);
                }
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

uint32_t ZigZagEncode(unsigned target, unsigned base, unsigned maxBits)
{
    {
        unsigned max = 1 << maxBits;
        unsigned half = max >> 1;
        int diff = static_cast<int>(target) - base;
        int iBase = static_cast<int>(base);

        if (base == half)
        {
            return ZigZag(diff);
        }

        if (base < half)
        {
            if ((diff < iBase) && (diff >= -iBase))
            {
                return ZigZag(diff);
            }

            return diff + base;
        }

        int antiBase = (max - base);

        if ((diff < antiBase) && (diff >= -antiBase))
        {
            return ZigZag(diff);
        }

        return (antiBase - 1) + -diff;
    }
}

uint32_t ZigZagDecode(unsigned value, unsigned base, unsigned maxBits)
{
    {
        unsigned max = 1 << maxBits;
        unsigned half = max >> 1;
        int iBase = static_cast<int>(base);
        int wtf = ((iBase - 1) * 2) + 1;

        if (half == base)
        {
            return iBase + ZigZag(value);
        }

        if (base < half)
        {
            if (static_cast<int>(value) <= wtf)
            {
                return iBase + ZigZag(value);
            }

            return value;
        }

        unsigned antiBase = (((max - 1) - base) * 2) + 1;

        if (value < antiBase)
        {
            return iBase + ZigZag(value);
        }

        return (max - 1) - value;
    }
}

void ZigZagTest()
{
    assert(42 == ZigZag(ZigZag(42)));
    assert(-42 == ZigZag(ZigZag(-42)));
    assert(0 == ZigZag(ZigZag(0)));
    assert(-12345 == ZigZag(ZigZag(-12345)));
    assert(30654 == ZigZag(ZigZag(30654)));
    assert(-31654 == ZigZag(ZigZag(-31654)));

    {
        auto encoded = ZigZagEncode(9, 7, 4);
        auto decoded = ZigZagDecode(encoded, 7, 4);

        assert(decoded == 9);
    }

    {
        auto encoded = ZigZagEncode(8, 15, 4);
        auto decoded = ZigZagDecode(encoded, 15, 4);

        assert(decoded == 8);
    }

    for (unsigned i = 0; i < 16; ++i)
    {
        for (unsigned j = 0; j < 16; ++j)
        {
            auto encoded = ZigZagEncode(i, j, 4);
            auto decoded = ZigZagDecode(encoded, j, 4);

            assert(decoded == i);
        }
    }
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

void TruncateTest()
{
    const unsigned count = 10;
    const unsigned max = count - 1;

    struct Set
    {
        unsigned i;
        unsigned max;
    };

    auto Test = [](const std::vector<Set>& sets)
    {
        BitStream coded2;
        BitStream coded;

        for (const auto& set : sets)
        {
            TruncateEncode(set.i, set.max, coded);
        }

        coded.Reset();

        for (const auto& set : sets)
        {
            auto result = TruncateDecode(set.max, coded);

            assert(result == set.i);
        }
    };

    {
        std::vector<Set> set
        {
            {446, 1021},
            {240, 359},
            {62, 240},
        };

        Test(set);
    }

    for (unsigned i = 0; i < count; ++i)
    {
        std::vector<Set> set = {{i, max}};
        Test(set);
    }
}
// //////////////////////////////////////////////////////

namespace {
    using namespace Range_models;
    struct Everything_model
    {
        // Note: I == investigate if dependency exists.
        // code quat changed as adaptive history binary
        // {
        //   assume no history on largest quat changed, simple binary adaptive
        //   3 models dependent on largest changed (I)
        //   {
        //     dependent model on largest changed to find largest vector index (I)
        //     dependent model on largest, to get next largest. (I)
        //     3 models per index
        //     {
        //       model on number of bits for first item (per max magnitude value)
        //       tree model on bits left, but with force 0 so we dictate how many
        //       bits are actually used. To avoid coding them in the first place.
        //       don't need to code last bits needed as we know already.
        //     }
        //   }
        //
        // code position changed dependent on quat changed. Hmm, rygorous does
        // 8 models, 2x dependent on quat changed, then 2x per quat different
        // component. weird. (I)
        // Code position just like quat, but with max magnitude specifying the bits
        // as well! yay! force zero.
        // interacting based on current interacting + quat/position changed.
        using Simple = Binary;

        Binary_history<Simple, Simple> quant_changed =
        {
            1,
            {4, 31}
        };

//        Binary_history<
//            Binary_history<Simple, Simple>,
//            Binary_history<Simple, Simple>> quant_changed =
//        {
//            0,
//            {0,6,7},
//            {1,7,7}
//        };

        // Test range vs tree
        // test if needed multiple models depending on previous model.
        Simple largest_index_quant_changed;
        Perodic_renomalisation<4, 8> largest_index_quant;

        struct Vector_model
        {
            Perodic_renomalisation<3, 8> largest_index;
            std::array<Simple, 3> next_largest_index;

            // Note: even though I hard code 12, it shouldn't be hard coded.
            // bits1 = MinBits(RotationMaxBits)
            // bits2 = MinBits(MaxPositionChangePerSnapshot * MaxFrameDelta)
            // bits = max(bits1, bits2)

            std::array<Perodic_renomalisation<12, 16>, 2> bits_for_value;
            std::array<Binary_tree<Simple, 12>, 3> value;

            bool Reduce_vector_using_magnitude;
        };

        // Change my mind, for now, only have one global model.
        // Investigate multiple models later.
        Vector_model quant_delta;
        Vector_model quant;

        std::array<Simple, 2> position_changed;

        Vector_model position;

        // 1 bit == previous interactive
        // 1 bit == quant changed
        // 1 bit == position changed
        std::array<Simple, 8> interactive;
    };
}

// //////////////////////////////////////////////////////

struct IntVec3
{
    int x;
    int y;
    int z;
};

// //////////////////////////////////////////////////////

auto Largest_next_magnitude(
        unsigned magnitude,
        unsigned zig_zag_axis) -> unsigned
{
    // max magnitude for zigzag before we overflow
    // == sqrt(1 << 32) / 2. == 32768.
    assert(magnitude < 32768);

    // positive zig zag numbers are just number * 2.
    unsigned next = magnitude * 2;
    next *= next;

    unsigned axis = zig_zag_axis;
    axis *= axis;

    assert(next >= axis);

    next -= axis;

    // C++11 casts the int into a double for sqrt, would it be faster
    // if I manually cast to a float instead?
    auto root = sqrt(next);

    // +1 for rounding.
    return 1 + static_cast<unsigned>(root) / 2;
}

void Vector3Encode(
        IntVec3 vec,
        unsigned max_magnitude,
        Everything_model::Vector_model& model,
        Encoder& range,
        Binary_encoder& binary)
{
    // +1 for the sign bit.
    unsigned max_bits_required = 1 + MinBits(max_magnitude);

    // //////////////////////////////////////////

    assert(abs(vec.x) <= static_cast<int>(max_magnitude));
    assert(abs(vec.y) <= static_cast<int>(max_magnitude));
    assert(abs(vec.z) <= static_cast<int>(max_magnitude));

    // Sort from largest to smallest
    auto zx = ZigZag(vec.x);
    auto zy = ZigZag(vec.y);
    auto zz = ZigZag(vec.z);

    // default order, x,y,z.
    {
        using std::swap;

        auto top = 0;
        auto next = 0;

        if  (
                (zy > zx) &&
                (zy >= zz)
            )
        {
            swap(zx, zy);
            top = 1;
        }
        else
        {
            if  (
                    (zz > zx) &&
                    (zz >= zy)
                )
            {
                swap(zx, zz);
                swap(zy, zz);
                top = 2;
            }
        }

        assert(zx >= zy);
        assert(zx >= zz);

        if (zz > zy)
        {
            swap(zy, zz);
            next = 1;
        }

        assert(zy >= zz);

        model.largest_index.Encode(range, top);
        model.next_largest_index[top].Encode(binary, next);
    }

    // //////////////////////////////////////////

    auto Code = [&](unsigned zig_zag, unsigned model_index) -> bool
    {
        unsigned bits = MinBits(zig_zag);

        model.bits_for_value[model_index].Encode(range, bits);

        assert(max_bits_required >= bits);

        if (!bits)
        {
            // everything after is zero, we're done.
            return false;
        }

        // Don't need to send the top bit, as we know it's set since
        // otherwise bits would be smaller.
        if (bits > 1)
        {
            auto t = zig_zag & ((1 << (bits - 1)) - 1);

            model.value[model_index].Encode(binary, t, bits - 1);
        }

        if (model.Reduce_vector_using_magnitude)
        {
            max_magnitude = Largest_next_magnitude(max_magnitude, zig_zag);
        }

        max_bits_required = 1 + MinBits(max_magnitude);
        max_bits_required = std::min(max_bits_required, bits);

        return true;
    };

    // //////////////////////////////////////////

    if (!Code(zx, 0))
    {
        return;
    }

    if (!Code(zy, 1))
    {
        return;
    }

    // //////////////////////////////////////////

    assert(max_bits_required >= MinBits(zz));

    model.value[2].Encode(binary, zz, max_bits_required);
}

IntVec3 Vector3Decode(
        unsigned max_magnitude,
        Everything_model::Vector_model& model,
        Decoder& range,
        Binary_decoder& binary)
{
    IntVec3 result = {0,0,0};

    // +1 for the sign bit.
    unsigned max_bits_required = 1 + MinBits(max_magnitude);

    // //////////////////////////////////////////

    // Read the Order
    auto top = model.largest_index.Decode(range);
    auto next = model.next_largest_index[top].Decode(binary);

    auto ReturnSorted = [&top, &next](IntVec3 vec) -> IntVec3
    {
        using std::swap;

        if (next)
        {
            swap(vec.y, vec.z);
        }

        if (top)
        {
            if (top == 1)
            {
                swap(vec.x, vec.y);
            }
            else
            {
                swap(vec.x, vec.y);
                swap(vec.y, vec.z);
            }
        }

        return vec;
    };

    // //////////////////////////////////////////

    auto Code = [&](int& target, unsigned model_index) -> bool
    {
        auto bits = model.bits_for_value[model_index].Decode(range);

        if (!bits)
        {
            return false;
        }

        auto zig_zag = bits > 1 ?
            model.value[model_index].Decode(binary, bits - 1) :
            0;

        // Don't need to send the top bit, as we know it's set since
        // otherwise bitsForzx would be smaller.
        zig_zag |= 1 << (bits - 1);
        target = ZigZag(zig_zag);

        if (model.Reduce_vector_using_magnitude)
        {
            max_magnitude = Largest_next_magnitude(max_magnitude, zig_zag);
        }

        max_bits_required = 1 + MinBits(max_magnitude);
        max_bits_required = std::min(max_bits_required, bits);

        return true;
    };

    // //////////////////////////////////////////

    if (!Code(result.x, 0))
    {
        return ReturnSorted(result);
    }

    if (!Code(result.y, 1))
    {
        return ReturnSorted(result);
    }

    // //////////////////////////////////////////

    auto zz = model.value[2].Decode(binary, max_bits_required);
    result.z = ZigZag(zz);

    return ReturnSorted(result);
}



auto Encode_frames(
    const Frame& base,
    const Frame& target,
    unsigned frameDelta) -> Range_types::Bytes
{
    // for now use defaults for everything.
    Range_types::Bytes              data;

    {
        Everything_model                model;
        Range_coders::Encoder           range(data);
        Range_coders::Binary_encoder    binary(range);

        model.quant.    Reduce_vector_using_magnitude = false;
        model.position. Reduce_vector_using_magnitude = true;

        unsigned max_position_delta =
            1 + static_cast<unsigned>(
                frameDelta * MaxPositionChangePerSnapshot);

        auto size = base.size();
        for (unsigned i = 0; i < size; ++i)
        {
            auto quant_index_changed =
                base[i].orientation_largest != target[i].orientation_largest;

            auto quant_changed =
                (quant_index_changed) ||
                (base[i].orientation_a != target[i].orientation_a) ||
                (base[i].orientation_b != target[i].orientation_b) ||
                (base[i].orientation_c != target[i].orientation_c);

            model.quant_changed.Encode(binary, quant_changed);

            if (quant_changed)
            {
                model.largest_index_quant_changed.Encode(
                    binary,
                    quant_index_changed);

                if (quant_index_changed)
                {
                    // RAM: TODO: correlation between previous and current?
                    // RAM: TODO: Binary tree instead?
                    model.largest_index_quant.Encode(
                        range, target[i].orientation_largest);

                    IntVec3 not_delta
                    {
                        target[i].orientation_a,
                        target[i].orientation_b,
                        target[i].orientation_c
                    };

                    Vector3Encode(
                        not_delta,
                        (1u << RotationMaxBits) - 1,
                        model.quant,
                        range,
                        binary);
                }

                if (!quant_index_changed)
                {
                    IntVec3 delta
                    {
                        target[i].orientation_a - base[i].orientation_a,
                        target[i].orientation_b - base[i].orientation_b,
                        target[i].orientation_c - base[i].orientation_c
                    };

                    Vector3Encode(
                        delta,
                        (1u << RotationMaxBits) - 1,
                        model.quant_delta,
                        range,
                        binary);
                }
            }

            auto pos_changed =
                (base[i].position_x != target[i].position_x) ||
                (base[i].position_y != target[i].position_y) ||
                (base[i].position_z != target[i].position_z);

            model.position_changed[quant_changed].Encode(binary, pos_changed);

            if (pos_changed)
            {
                IntVec3 delta
                {
                    target[i].position_x - base[i].position_x,
                    target[i].position_y - base[i].position_y,
                    target[i].position_z - base[i].position_z
                };

                Vector3Encode(
                    delta,
                    max_position_delta,
                    model.position,
                    range,
                    binary);
            }

            unsigned interactive_index = (base[i].interacting == 1) << 2;
            interactive_index |= quant_changed << 1;
            interactive_index |= pos_changed;

            model.interactive[interactive_index].Encode(
                binary,
                target[i].interacting);
        }
    }

    return data;
}

auto Decode_frames(
    const Frame& base,
    const Range_types::Bytes& data,
    unsigned frameDelta) -> Frame
{
    Frame target;
    Everything_model                model;
    Range_coders::Decoder           range(data);
    Range_coders::Binary_decoder    binary(range);

    model.quant.    Reduce_vector_using_magnitude = false;
    model.position. Reduce_vector_using_magnitude = true;

    unsigned max_position_delta =
        1 + static_cast<unsigned>(
            frameDelta * MaxPositionChangePerSnapshot);

    auto size = base.size();
    for (unsigned i = 0; i < size; ++i)
    {
        auto quant_changed =
            model.quant_changed.Decode(binary);

        if (quant_changed)
        {
            auto quant_index_changed =
                model.largest_index_quant_changed.Decode(binary);

            if (quant_index_changed)
            {
                target[i].orientation_largest =
                    model.largest_index_quant.Decode(range);

                auto not_delta = Vector3Decode(
                    (1u << RotationMaxBits) - 1,
                    model.quant,
                    range,
                    binary);

                target[i].orientation_a = not_delta.x;
                target[i].orientation_b = not_delta.y;
                target[i].orientation_c = not_delta.z;
            }

            if (!quant_index_changed)
            {
                IntVec3 delta = Vector3Decode(
                    (1u << RotationMaxBits) - 1,
                    model.quant_delta,
                    range,
                    binary);

                target[i].orientation_largest =
                        base[i].orientation_largest;
                target[i].orientation_a = base[i].orientation_a + delta.x;
                target[i].orientation_b = base[i].orientation_b + delta.y;
                target[i].orientation_c = base[i].orientation_c + delta.z;
            }
        }

        if (!quant_changed)
        {
            target[i].orientation_largest =
                    base[i].orientation_largest;
            target[i].orientation_a = base[i].orientation_a;
            target[i].orientation_b = base[i].orientation_b;
            target[i].orientation_c = base[i].orientation_c;
        }

        auto pos_changed =
            model.position_changed[quant_changed].Decode(binary);

        if (pos_changed)
        {
            auto delta = Vector3Decode(
                max_position_delta,
                model.position,
                range,
                binary);

            target[i].position_x = base[i].position_x + delta.x;
            target[i].position_y = base[i].position_y + delta.y;
            target[i].position_z = base[i].position_z + delta.z;
        }

        if (!pos_changed)
        {
            target[i].position_x = base[i].position_x;
            target[i].position_y = base[i].position_y;
            target[i].position_z = base[i].position_z;
        }

        unsigned interactive_index = (base[i].interacting == 1) << 2;
        interactive_index |= quant_changed << 1;
        interactive_index |= pos_changed;

        target[i].interacting =
            model.interactive[interactive_index].Decode(binary);
    }

    return target;
}

auto Model_tests()
{
    {
        int const max = static_cast<int>((MaxPositionChangePerSnapshot) * 6 + 1);

        IntVec3 data =
        {
            -434,
            -90,
            0,
        };

        Range_coders::Bytes bytes;


        {
            Everything_model::Vector_model  model;
            Range_coders::Encoder           e_range(bytes);
            Range_coders::Binary_encoder    e_binary(e_range);

            model.Reduce_vector_using_magnitude = true;

            for (auto l = 0; l < 2; ++l)
            {
                Vector3Encode(
                    data,
                    max,
                    model,
                    e_range,
                    e_binary);
            }
        }

        {
            Everything_model::Vector_model  model;
            Range_coders::Decoder           d_range(bytes);
            Range_coders::Binary_decoder    d_binary(d_range);

            model.Reduce_vector_using_magnitude = true;

            for (auto l = 0; l < 2; ++l)
            {
                auto decoded = Vector3Decode(
                    max,
                    model,
                    d_range,
                    d_binary);

                assert(data.x == decoded.x);
                assert(data.y == decoded.y);
                assert(data.z == decoded.z);
            }
        }
    }
}

// //////////////////////////////////////////////////////

using RangeBits = std::array<uint8_t, 3>;

unsigned constexpr MaxRange(const RangeBits& ranges)
{
    return
        (1u << ranges[0]) +
        (1u << ranges[1]) +
        (1u << ranges[2]);
}

unsigned GaffersRangeEncode(
        RangeBits ranges,
        unsigned value,
        BitStream& target)
{
    unsigned bitsUsed = 0;

    struct MinCount
    {
        unsigned min;
        unsigned count;
    };

    using MinsCounts = std::array<MinCount, 3>;

    MinsCounts minCounts =
    {
        MinCount{0,                                     1u << ranges[0]},
        MinCount{1u << ranges[0],                       1u << ranges[1]},
        MinCount{(1u << ranges[0]) + (1u << ranges[1]), 1u << ranges[2]},
    };

    if (value < minCounts[0].count)
    {
        target.Write(1, 1);
        target.Write(value, ranges[0]);
        return bitsUsed + 1 + ranges[0];
    }

    target.Write(0, 1);
    ++bitsUsed;

    if (value < (minCounts[1].count + minCounts[1].min))
    {
        target.Write(1, 1);
        target.Write(value - minCounts[1].min, ranges[1]);
        return bitsUsed + 1 + ranges[1];
    }

    target.Write(0, 1);
    ++bitsUsed;

    assert (value < MaxRange(ranges));

    target.Write(value - minCounts[2].min, ranges[2]);

    return bitsUsed + 1 + ranges[2];
}

unsigned GaffersRangeDecode(
        RangeBits ranges,
        BitStream& source)
{
    auto first = source.Read(1);

    if (first)
    {
        return source.Read(ranges[0]);
    }

    auto second = source.Read(1);

    if (second)
    {
        auto result = source.Read(ranges[1]);
        return result + (1 << ranges[0]);
    }

    auto result = source.Read(ranges[2]);
    return result + (1 << ranges[0]) + (1 << ranges[1]);
}

void GaffersRangeTest()
{
    auto ranges = RangeBits
    {
        5,
        6,
        7
    };

    auto maxRange = MaxRange(ranges);
    for (unsigned i = 0; i < maxRange; ++i)
    {
        BitStream stream;

        GaffersRangeEncode(ranges, i, stream);
        stream.Reset();
        auto result = GaffersRangeDecode(ranges, stream);

        assert(i == result);
    }
}

// //////////////////////////////////////////////////////

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

    float next = static_cast<float>(maxMagnitude * maxMagnitude);
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

    next = static_cast<float>(maxMagnitude * maxMagnitude);
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

    float next = static_cast<float>(maxMagnitude * maxMagnitude);
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

    next = static_cast<float>(maxMagnitude * maxMagnitude);
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

    float next = static_cast<float>(maxMagnitude * maxMagnitude);
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

    next = static_cast<float>(maxMagnitude * maxMagnitude);
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

    float next = static_cast<float>(maxMagnitude * maxMagnitude);
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

    next = static_cast<float>(maxMagnitude * maxMagnitude);
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

// //////////////////////////////////////////////////////

using Quat = std::array<float, 4>;

struct Gaffer
{
    unsigned largest_index;
    int a;
    int b;
    int c;
};

float constexpr Magnitude_squared(const Quat& quat)
{
    return
        quat[0] * quat[0] +
        quat[1] * quat[1] +
        quat[2] * quat[2] +
        quat[3] * quat[3];
}

float q_to_g = 256.0 * 1.414213562373095048801688724209698078569671875;
float g_to_q = 1.0 / (256.0 * 1.414213562373095048801688724209698078569671875);

// stupid float precision and quantising
const float q_max_s = 256.49 * 256.49 * 2.0;
const float q_max = std::truncf(sqrt(q_max_s)) + 0.49;

// Again: did some research, -16, -16, -256 seems to be the worse.
// try using that as max sum?
const float gaffer_one_squared =
        (256 * 256) +
        (256 * 256) +
        (16 * 16) +
        (16 * 16);

// Ah thats because gaffer_one_squared == (sqrt(2) * 256.4995127)^2
// gaffer_one = sqrt(2) * 256.4995127

// Question is: When do we use 256, and when do we use 256.4995127?

float g_256 = 256.4995127;
float q_to_g2 = g_256 * 1.414213562373095048801688724209698078569671875;
float g_to_q2 = 1.0 / (g_256 * 1.414213562373095048801688724209698078569671875);

Quat2 ConvertGaffer2(const Gaffer& gaffer)
{
    Quat2 result
    {
        static_cast<float>(gaffer.a - 256),
        static_cast<float>(gaffer.b - 256),
        static_cast<float>(gaffer.c - 256),
        0.0,
    };

    auto largest_squared = Quat2
    {
        result[0] * result[0],
        result[1] * result[1],
        result[2] * result[2],
        0,
    };

    auto magnitude_squared =
            largest_squared[0] + largest_squared[1] + largest_squared[2];

    auto largest_value_squared = gaffer_one_squared - magnitude_squared;

    assert(largest_value_squared >= 0);

    // Deal with quantising errors.
    auto next_largest = largest_squared[0];
    if (next_largest < largest_squared[1])
    {
        next_largest = largest_squared[1];
    }
    if (next_largest < largest_squared[2])
    {
        next_largest = largest_squared[2];
    }

    if (next_largest > largest_value_squared)
    {
        // Do I need to add an offset to make sure it is still
        // the largest as opposed to equal?
        // RAM: TODO: YES!
        // RAM: TODO: figure out by what number ill use.
        largest_value_squared = next_largest + 20;
    }

    auto largest = sqrt(largest_value_squared);

    if (gaffer.largest_index == 0)
    {
        result[3] = result[2];
        result[2] = result[1];
        result[1] = result[0];
    }

    if (gaffer.largest_index == 1)
    {
        result[3] = result[2];
        result[2] = result[1];
    }

    if (gaffer.largest_index == 2)
    {
        result[3] = result[2];
    }

    result[gaffer.largest_index] = largest;
    result = Mul(result, g_to_q2);

    {
        auto mag = Magnitude_squared(result);
        auto quantised = 256 * (mag - 1.0f);
        assert(std::abs(quantised) < 0.5f);
    }

    result = Normalise(result);
    return result;
}

Gaffer ConvertGaffer2(const Quat2& quat)
{
    const auto size = quat.q.size();
    unsigned largest_index = 0;
    float largest = 0;

    auto squared = Quat2
    {
        quat[0] * quat[0],
        quat[1] * quat[1],
        quat[2] * quat[2],
        quat[3] * quat[3],
    };

    for (unsigned i = 0; i < size; ++i)
    {
        if (squared[i] >= largest)
        {
            if (squared[i] == largest)
            {
                if (quat[i] < 0)
                {
                    continue;
                }
            }

            largest = squared[i];
            largest_index = i;
        }
    }


    // If the largest is -ve then we need to negate the quat so it's positive.
    auto fixed_quat = quat;
    if (quat[largest_index] < 0)
    {
        fixed_quat = Mul(quat, -1.0f);
    }

    auto write_index = 0;
    Quat gaffer;

    for (unsigned i = 0; i < size; ++i)
    {
        if (i != largest_index)
        {
            gaffer[write_index++] = fixed_quat[i];
        }
    }

    // Screw it, if the quat would produce 256 or -257
    // then re adjust
    auto multiple = 1.0f;
    const auto max_positive = (255.0f / 256.0f) / sqrt(2.0f);
    const auto max_negative = -1.0f / sqrt(2.0f);

    if (gaffer[0] > max_positive)
    {
        multiple = max_positive / gaffer[0];
    }
    if (gaffer[1] > max_positive)
    {
        multiple = max_positive / gaffer[1];
    }
    if (gaffer[2] > max_positive)
    {
        multiple = max_positive / gaffer[2];
    }


    if (gaffer[0] < max_negative)
    {
        multiple = max_negative / gaffer[0];
    }
    if (gaffer[1] < max_negative)
    {
        multiple = max_negative / gaffer[1];
    }
    if (gaffer[2] < max_negative)
    {
        multiple = max_negative / gaffer[2];
    }

    auto result = Gaffer
    {
        largest_index,
        256 + static_cast<int>(round(gaffer[0] * q_to_g2 * multiple)),
        256 + static_cast<int>(round(gaffer[1] * q_to_g2 * multiple)),
        256 + static_cast<int>(round(gaffer[2] * q_to_g2 * multiple)),
    };

//    auto result = Gaffer
//    {
//        largest_index,
//        256 + static_cast<int>(gaffer[0] * q_to_g * multiple),
//        256 + static_cast<int>(gaffer[1] * q_to_g * multiple),
//        256 + static_cast<int>(gaffer[2] * q_to_g * multiple),
//    };

//    auto result = Gaffer
//    {
//        largest_index,
//        256 + static_cast<int>(round(gaffer[0] * q_to_g2)),
//        256 + static_cast<int>(round(gaffer[1] * q_to_g2)),
//        256 + static_cast<int>(round(gaffer[2] * q_to_g2)),
//    };

//    auto result = Gaffer
//    {
//        largest_index,
//        256 + static_cast<int>(gaffer[0] * q_to_g2),
//        256 + static_cast<int>(gaffer[1] * q_to_g2),
//        256 + static_cast<int>(gaffer[2] * q_to_g2),
//    };

    assert(result.a >= 0);
    assert(result.b >= 0);
    assert(result.c >= 0);

    assert(result.a < 512);
    assert(result.b < 512);
    assert(result.c < 512);

    return result;
}

Quat ConvertGaffer(const Gaffer& gaffer)
{
    // RAM: TODO: calculate missing item before
    // converting to 0-1.
    Quat result
    {
        (gaffer.a - 256) * g_to_q,
        (gaffer.b - 256) * g_to_q,
        (gaffer.c - 256) * g_to_q,
        0.0,
    };

    // RAM: TODO: This still is broken.
    auto largest = sqrt(1.0f - Magnitude_squared(result));

    assert(largest >= 0);

    if (gaffer.largest_index == 0)
    {
        result[3] = result[2];
        result[2] = result[1];
        result[1] = result[0];
    }

    if (gaffer.largest_index == 1)
    {
        result[3] = result[2];
        result[2] = result[1];
    }

    if (gaffer.largest_index == 2)
    {
        result[3] = result[2];
    }

    result[gaffer.largest_index] = largest;

    {
        auto mag = Magnitude_squared(result);
        auto quantised = 256 * (mag - 1.0f);
        assert(std::abs(quantised) < 0.5f);
    }

    return result;
}

Gaffer ConvertGaffer(const Quat& quat)
{
    const auto size = quat.size();
    unsigned largest_index = 0;
    float largest = 0;

    Quat squared
    {
        quat[0] * quat[0],
        quat[1] * quat[1],
        quat[2] * quat[2],
        quat[3] * quat[3],
    };

    for (unsigned i = 0; i < size; ++i)
    {
        if (squared[i] > largest)
        {
            largest = squared[i];
            largest_index = i;
        }
    }

    auto write_index = 0;
    Quat gaffer;

    for (unsigned i = 0; i < size; ++i)
    {
        if (i != largest_index)
        {
            gaffer[write_index++] = quat[i];
        }
    }

    return
    {
        largest_index,
        256 + static_cast<int>(round(gaffer[0] * q_to_g)),
        256 + static_cast<int>(round(gaffer[1] * q_to_g)),
        256 + static_cast<int>(round(gaffer[2] * q_to_g)),
    };
}

// //////////////////////////////////////////////////////

// RAM: Need to actuall figure out these values
//static const float      MAX_ANGULAR_VELOCITY_PER_FRAME = 10.0f;
static const unsigned   ROTOR_BITS = 12;
static const unsigned   ROTOR_MULTIPLE = 1 << ROTOR_BITS;
static const float      ROTOR_MULTIPLE_INV = 1.0 / ROTOR_MULTIPLE;

IntVec3 Chris_Doran(const Quat2& base, const Quat2& target)
{
    assert(!Equal(base, target));

    auto r = R(base, target);
    auto rotor = to_rotor(r);

    // RAM: Do i need to round or truncate?
    auto result = IntVec3
    {
        static_cast<int>(rotor[0] * ROTOR_MULTIPLE),
        static_cast<int>(rotor[1] * ROTOR_MULTIPLE),
        static_cast<int>(rotor[2] * ROTOR_MULTIPLE),
    };

    // RAM: TODO: find out max bits needed to encode the rotor!
//    assert (result.x < (int) ROTOR_MULTIPLE);
//    assert (result.y < (int) ROTOR_MULTIPLE);
//    assert (result.z < (int) ROTOR_MULTIPLE);

    return result;
}

Quat2 Chris_Doran(const Quat2& base, const IntVec3& encoded)
{
    Rotor rotor =
    {
        encoded.x * ROTOR_MULTIPLE_INV,
        encoded.y * ROTOR_MULTIPLE_INV,
        encoded.z * ROTOR_MULTIPLE_INV,
    };

    auto r = to_quat(rotor);
    auto result = Mul(r, base);

    return result;
}

IntVec3 Rotorify(const Gaffer& base, const Gaffer& target)
{
    auto base_quat = ConvertGaffer2(base);
    auto target_quat = ConvertGaffer2(target);

    return Chris_Doran(base_quat, target_quat);
}

Gaffer Rotorify(const Gaffer& base, const IntVec3& encoded)
{
    auto base_quat = ConvertGaffer2(base);
    auto target_quat = Chris_Doran(base_quat, encoded);

    return ConvertGaffer2(target_quat);
}

// //////////////////////////////////////////////////////

void Gaffer_tests()
{
    {
        auto gaffer = Gaffer
        {
            3,
            256,
            256,
            255,
        };

        auto quat = ConvertGaffer2(gaffer);
        auto result = ConvertGaffer2(quat);

        assert(result.largest_index == gaffer.largest_index);
        assert(result.a == gaffer.a);
        assert(result.b == gaffer.b);
        assert(result.c == gaffer.c);
    }
    {
        auto gaffer = Gaffer
        {
            3,
            240,
            0,
            240,
        };

        auto quat = ConvertGaffer2(gaffer);
        auto result = ConvertGaffer2(quat);

        assert(result.largest_index == gaffer.largest_index);
        assert(result.a == gaffer.a);
        assert(result.b == gaffer.b);
        assert(result.c == gaffer.c);
    }
    {
        auto b = Gaffer
        {
            0,
            54,
            326,
            214,
        };

        auto t = Gaffer
        {
            0,
            101,
            314,
            28
        };

        auto rotor = Rotorify(b, t);
        auto result = Rotorify(b, rotor);

        assert(result.largest_index == t.largest_index);
        assert(result.a == t.a);
        assert(result.b == t.b);
        assert(result.c == t.c);
    }

    for (int i = 0; i < 512; ++i)
    {
        for (int j = 0; j < 512; j += 11)
        {
            for (int k = 0; k < 512; k += 37)
            {
                // Ahah! I spot an oppertunity to compress
                // Not all permitations of the gaffer are valid.
                Quat g
                {
                    (i - 256) * g_to_q,
                    (j - 256) * g_to_q,
                    (k - 256) * g_to_q,
                    0.0,
                };

                auto largest = sqrt(1.0f - Magnitude_squared(g));

                bool valid = (largest >= 0);
                valid &= largest >= std::abs(g[0]);
                valid &= largest >= std::abs(g[1]);
                valid &= largest >= std::abs(g[2]);

                if (valid)
                {
                    for (unsigned l = 0; l < 4; ++l)
                    {
                        auto gaffer = Gaffer
                        {
                            l,
                            i,
                            j,
                            k
                        };

                        {

                            auto quat = ConvertGaffer(gaffer);
                            auto result = ConvertGaffer(quat);

                            assert(result.largest_index == gaffer.largest_index);
                            assert(result.a == gaffer.a);
                            assert(result.b == gaffer.b);
                            assert(result.c == gaffer.c);
                        }

                        {
                            auto base = Gaffer
                            {
                                l,
                                256,
                                256,
                                256,
                            };

                            auto rotor = Rotorify(base, gaffer);
                            auto result = Rotorify(base, rotor);

                            assert(result.largest_index == gaffer.largest_index);
                            assert(result.a == gaffer.a);
                            assert(result.b == gaffer.b);
                            assert(result.c == gaffer.c);
                        }
                    }
                }
            }
        }
    }
}

int Max_gaffer_value(int second_largest)
{
    auto second = second_largest;

    // next value is the largest if the lagest value is the smallest.
    // The smallest fir lagest can be is the same size as the second
    // largest.
    // And since the magnitude sums to one, we can get a max bound
    // for the third value.
    auto s = 2 * second * second;

    float mag = sqrt(q_max_s - s);
    auto quantised = mag;

    if (quantised < 0.5f)
    {
        return 0;
    }

    return static_cast<int>(quantised) + 1;
}

int Max_gaffer_value(int second_largest, int third_largest)
{
    auto second = second_largest;
    auto third = third_largest;
    auto s = (2 * second * second) + (third * third);

    float mag = sqrt(q_max_s - s);
    auto quantised = mag;

    if (quantised < 0.5f)
    {
        return 0;
    }

    return static_cast<int>(quantised) + 1;
}

void Max_gaffer_tests()
{
    // RAM: TODO: Some tests maybe?
//    {
//        auto should_be_zero = Max_gaffer_value(256);
//        assert(should_be_zero == 0);
//    }

//    for (int i = 0; i < 256; ++i)
//    {
//        auto result = std::min(i, Max_gaffer_value(i));
//        printf("%d: -> %d\n", i, result);
//    }

//    for (int i = 0; i < 256; ++i)
//    {
//        for (int j = 0; j < i; ++j)
//        {
//            auto result = std::min(j, Max_gaffer_value(i, j));
//            printf("%d, %d: -> %d\n", i, j, result);
//        }
//    }
}

// //////////////////////////////////////////////////////

enum class Use_magnitude_as
{
    Constant,
    Vector_Magnitude,
    Gaffer_Encode,
};

unsigned Sorted_no_bit_count_Encode(
        IntVec3 vec,
        unsigned maxMagnitude,
        Use_magnitude_as use_magnitude_as,
        BitStream& target)
{
    unsigned bitsUsed = 0;

    // //////////////////////////////////////////

    assert(abs(vec.x) <= static_cast<int>(maxMagnitude));
    assert(abs(vec.y) <= static_cast<int>(maxMagnitude));
    assert(abs(vec.z) <= static_cast<int>(maxMagnitude));

    // Sort from largest to smallest
    auto zx = ZigZag(vec.x);
    auto zy = ZigZag(vec.y);
    auto zz = ZigZag(vec.z);

    // default order, x,y,z.
    auto top = 0;
    auto next = 0;
    auto odd = -1;

    {
        if ((zx != zy) && (zy == zz))
        {
            odd = 0;
        }
        if ((zx != zy) && (zx == zz))
        {
            odd = 1;
        }
        if ((zx == zy) && (zy != zz))
        {
            odd = 2;
        }
    }

    {
        using std::swap;

        if  (
                (zy > zx) &&
                (zy >= zz)
            )
        {
            swap(zx, zy);
            top = 1;
        }
        else
        {
            if  (
                    (zz > zx) &&
                    (zz >= zy)
                )
            {
                swap(zx, zz);
                swap(zy, zz);
                top = 2;
            }
        }

        assert(zx >= zy);
        assert(zx >= zz);

        if (zz > zy)
        {
            swap(zy, zz);
            next = 1;
        }

        assert(zy >= zz);
    }

    // //////////////////////////////////////////

    auto zig_zag_max = ZigZag(-static_cast<int>(maxMagnitude));
    auto index_count = 0;

    // Ick
    auto x = ZigZag(zx);

    // //////////////////////////////////////////

    auto Code = [&](unsigned zig_zag)
    {
        assert(zig_zag <= zig_zag_max);

        bitsUsed += TruncateEncode(zig_zag, zig_zag_max, target);

        if (use_magnitude_as == Use_magnitude_as::Vector_Magnitude)
        {
            maxMagnitude = Largest_next_magnitude(maxMagnitude, zig_zag);
        }

        if (use_magnitude_as == Use_magnitude_as::Gaffer_Encode)
        {
            switch (index_count++)
            {
                case 0:
                {
                    maxMagnitude = Max_gaffer_value(x);
                    break;
                }
                case 1:
                {
                    maxMagnitude = Max_gaffer_value(x,ZigZag(zig_zag));
                    break;
                }

                default:
                case 2:
                {
                    maxMagnitude = 0;
                    break;
                }
            }
        }

        zig_zag_max = std::min(
            zig_zag,
            ZigZag(-static_cast<int>(maxMagnitude)));
    };

    // //////////////////////////////////////////

    Code(zx);

    if ((zx) && (zig_zag_max))
    {
        Code(zy);

        if ((zy) && (zig_zag_max))
        {
            Code(zz);
        }
    }

    // //////////////////////////////////////////

    bool all_same = (zx == zy) && (zy == zz);

    if (!all_same)
    {
        if (odd < 0)
        {
            bitsUsed += TruncateEncode(top, 3, target);
            target.Write(next, 1);
            ++bitsUsed;
        }
        else
        {
            bitsUsed += TruncateEncode(odd, 3, target);
        }
    }

    return bitsUsed;
}

IntVec3 Sorted_no_bit_count_Decode(
        unsigned maxMagnitude,
        Use_magnitude_as use_magnitude_as,
        BitStream& source)
{
    IntVec3 result = {0,0,0};

    // //////////////////////////////////////////

    // Read the Order
    auto largest = 0;
    auto nextLargest = 0;

    auto ReturnSorted = [&largest, &nextLargest](IntVec3 vec) -> IntVec3
    {
        using std::swap;

        if (nextLargest)
        {
            swap(vec.y, vec.z);
        }

        if (largest)
        {
            if (largest == 1)
            {
                swap(vec.x, vec.y);
            }
            else
            {
                swap(vec.x, vec.y);
                swap(vec.y, vec.z);
            }
        }

        return vec;
    };

    // //////////////////////////////////////////

    auto zig_zag_max = ZigZag(-static_cast<int>(maxMagnitude));
    auto index_count = 0;

    // //////////////////////////////////////////

    auto Code = [&](int& target)
    {
        auto zig_zag = TruncateDecode(zig_zag_max, source);

        target = ZigZag(zig_zag);

        if (use_magnitude_as == Use_magnitude_as::Vector_Magnitude)
        {
            maxMagnitude = Largest_next_magnitude(maxMagnitude, zig_zag);
        }

        if (use_magnitude_as == Use_magnitude_as::Gaffer_Encode)
        {
            switch (index_count++)
            {
                case 0:
                {
                    maxMagnitude = Max_gaffer_value(target);
                    break;
                }
                case 1:
                {
                    maxMagnitude = Max_gaffer_value(result.x,target);
                    break;
                }

                default:
                case 2:
                {
                    maxMagnitude = 0;
                    break;
                }
            }
        }

        zig_zag_max = std::min(
            zig_zag,
            ZigZag(-static_cast<int>(maxMagnitude)));
    };

    // //////////////////////////////////////////

    Code(result.x);

    if ((result.x) && (zig_zag_max))
    {
        Code(result.y);

        if ((result.y) && (zig_zag_max))
        {
            Code(result.z);
        }
    }

    // //////////////////////////////////////////

    if ((result.x != result.y) || (result.y != result.z))
    {
        if ((result.x != result.y) && (result.y != result.z))
        {
            largest = TruncateDecode(3, source);
            nextLargest = source.Read(1);
        }
        else
        {
            bool all_same = (result.x == result.y) && (result.y == result.z);

            if (!all_same)
            {
                auto odd_one = TruncateDecode(3, source);

                if (result.x != result.y)
                {
                    largest = odd_one;
                }
                else
                {
                    switch (odd_one)
                    {
                        default:
                        case 0:
                        {
                            largest = 1;
                            nextLargest = 1;
                            break;
                        }
                        case 1:
                        {
                            largest = 0;
                            nextLargest = 1;
                            break;
                        }
                        case 2:
                        {
                            break;
                        }
                    }
                }
            }
        }
    }

    return ReturnSorted(result);
}

// //////////////////////////////////////////////////////

unsigned BitVector3SortedEncode(
        IntVec3 vec,
        unsigned maxMagnitude,
        Use_magnitude_as use_magnitude_as,
        BitStream& target)
{
    unsigned bitsUsed = 0;

    // +1 for the sign bit.
    unsigned maxBitsRequired = 1 + MinBits(maxMagnitude);

    // //////////////////////////////////////////

    assert(abs(vec.x) <= static_cast<int>(maxMagnitude));
    assert(abs(vec.y) <= static_cast<int>(maxMagnitude));
    assert(abs(vec.z) <= static_cast<int>(maxMagnitude));

    // Sort from largest to smallest
    auto zx = ZigZag(vec.x);
    auto zy = ZigZag(vec.y);
    auto zz = ZigZag(vec.z);

    // default order, x,y,z.
    {
        using std::swap;

        auto top = 0;
        auto next = 0;

        if  (
                (zy > zx) &&
                (zy >= zz)
            )
        {
            swap(zx, zy);
            top = 1;
        }
        else
        {
            if  (
                    (zz > zx) &&
                    (zz >= zy)
                )
            {
                swap(zx, zz);
                swap(zy, zz);
                top = 2;
            }
        }

        assert(zx >= zy);
        assert(zx >= zz);

        if (zz > zy)
        {
            swap(zy, zz);
            next = 1;
        }

        assert(zy >= zz);

        bitsUsed += TruncateEncode(top, 3, target);
        target.Write(next, 1);
        ++bitsUsed;
    }

    // //////////////////////////////////////////

    auto Code = [&](unsigned zig_zag) -> bool
    {
        unsigned bits = MinBits(zig_zag);

        assert(maxBitsRequired >= bits);

        bitsUsed += TruncateEncode(bits, maxBitsRequired, target);

        if (!bits)
        {
            // everything after is zero, we're done.
            return false;
        }

        // Don't need to send the top bit, as we know it's set since
        // otherwise bits would be smaller.
        {
            auto t = zig_zag & ((1 << (bits - 1)) - 1);
            target.Write(t, bits - 1);
            bitsUsed += bits - 1;
        }

        if (use_magnitude_as == Use_magnitude_as::Vector_Magnitude)
        {
            maxMagnitude = Largest_next_magnitude(maxMagnitude, zig_zag);
        }

        maxBitsRequired = 1 + MinBits(maxMagnitude);
        maxBitsRequired = std::min(maxBitsRequired, bits);

        return true;
    };

    // //////////////////////////////////////////

    if (!Code(zx))
    {
        return bitsUsed;
    }

    if (!Code(zy))
    {
        return bitsUsed;
    }

    // //////////////////////////////////////////

    assert(maxBitsRequired >= MinBits(zz));

    target.Write(zz, maxBitsRequired);
    bitsUsed += maxBitsRequired;

    return bitsUsed;
}

IntVec3 BitVector3SortedDecode(
        unsigned maxMagnitude,
        Use_magnitude_as use_magnitude_as,
        BitStream& source)
{
    IntVec3 result = {0,0,0};

    // +1 for the sign bit.
    unsigned maxBitsRequired = 1 + MinBits(maxMagnitude);

    // //////////////////////////////////////////

    // Read the Order
    auto largest = TruncateDecode(3, source);
    auto nextLargest = source.Read(1);

    auto ReturnSorted = [&largest, &nextLargest](IntVec3 vec) -> IntVec3
    {
        using std::swap;

        if (nextLargest)
        {
            swap(vec.y, vec.z);
        }

        if (largest)
        {
            if (largest == 1)
            {
                swap(vec.x, vec.y);
            }
            else
            {
                swap(vec.x, vec.y);
                swap(vec.y, vec.z);
            }
        }

        return vec;
    };

    // //////////////////////////////////////////

    auto Code = [&](int& target) -> bool
    {
        auto bits = TruncateDecode(maxBitsRequired, source);

        if (!bits)
        {
            return false;
        }

        auto zig_zag = source.Read(bits - 1);

        // Don't need to send the top bit, as we know it's set since
        // otherwise bitsForzx would be smaller.
        zig_zag |= 1 << (bits - 1);
        target = ZigZag(zig_zag);

        if (use_magnitude_as == Use_magnitude_as::Vector_Magnitude)
        {
            maxMagnitude = Largest_next_magnitude(maxMagnitude, zig_zag);
        }

        maxBitsRequired = 1 + MinBits(maxMagnitude);
        maxBitsRequired = std::min(maxBitsRequired, bits);

        return true;
    };

    // //////////////////////////////////////////

    if (!Code(result.x))
    {
        return ReturnSorted(result);
    }

    if (!Code(result.y))
    {
        return ReturnSorted(result);
    }

    // //////////////////////////////////////////

    auto zz = source.Read(maxBitsRequired);
    result.z = ZigZag(zz);

    return ReturnSorted(result);
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

    float next = static_cast<float>(maxMagnitude * maxMagnitude);
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

    next = static_cast<float>(maxMagnitude * maxMagnitude);
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

    float next = static_cast<float>(maxMagnitude * maxMagnitude);
    next -= result.x * result.x;
    maxMagnitude = static_cast<unsigned>(sqrt(next) + 1);

    // //////////////////////////////////////////

    maxBitsRequired = 1 + MinBits(maxMagnitude);

    auto shiftY = source.Read(prefixBits);

    auto zy = source.Read(maxBitsRequired >> shiftY);
    result.y = ZigZag(zy);

    next = static_cast<float>(maxMagnitude * maxMagnitude);
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

// //////////////////////////////////////////////////////

unsigned BitVector3UnrelatedEncode(
        IntVec3 vec,
        unsigned maxMagnitude,
        BitStream& target)
{
    unsigned bitsUsed = 0;
    unsigned maxBitsPerComponent = 1 + MinBits(maxMagnitude);

    auto zx = ZigZag(vec.x);
    assert(zx < (1u << maxBitsPerComponent));

    unsigned bitsForzx = MinBits(zx);
    bitsUsed += TruncateEncode(bitsForzx, maxBitsPerComponent, target);

    // Don't need to send the top bit, as we know it's set since
    // otherwise bitsForzx would be smaller.
    if (bitsForzx)
    {
        zx &= (1 << (bitsForzx - 1)) - 1;
        target.Write(zx, bitsForzx - 1);
        bitsUsed += bitsForzx - 1;
    }

    // //////////////////////////////////////////

    auto zy = ZigZag(vec.y);
    assert(zy < (1u << maxBitsPerComponent));

    unsigned bitsForzy = MinBits(zy);
    bitsUsed += TruncateEncode(bitsForzy, maxBitsPerComponent, target);

    if (bitsForzy)
    {
        zy &= (1 << (bitsForzy - 1)) - 1;
        target.Write(zy, bitsForzy - 1);
        bitsUsed += bitsForzy - 1;
    }

    // //////////////////////////////////////////

    auto zz = ZigZag(vec.z);
    assert(zz < (1u << maxBitsPerComponent));

    unsigned bitsForzz = MinBits(zz);
    bitsUsed += TruncateEncode(bitsForzz, maxBitsPerComponent, target);

    if (bitsForzz)
    {
        zz &= (1 << (bitsForzz - 1)) - 1;
        target.Write(zz, bitsForzz - 1);
        bitsUsed += bitsForzz - 1;
    }

    return bitsUsed;
}

IntVec3 BitVector3UnrelatedDecode(
        unsigned maxMagnitude,
        BitStream& source)
{
    IntVec3 result = {0,0,0};

    unsigned maxBitsPerComponent = 1 + MinBits(maxMagnitude);

    auto bitsX = TruncateDecode(maxBitsPerComponent, source);

    if (bitsX)
    {
        auto zx = source.Read(bitsX - 1);
        // Don't need to send the top bit, as we know it's set since
        // otherwise bitsForzx would be smaller.
        zx |= 1 << (bitsX - 1);
        result.x = ZigZag(zx);
    }

    // //////////////////////////////////////////

    auto bitsY = TruncateDecode(maxBitsPerComponent, source);

    if (bitsY)
    {
        auto zy = source.Read(bitsY - 1);
        zy |= 1 << (bitsY - 1);
        result.y = ZigZag(zy);
    }

    // //////////////////////////////////////////

    auto bitsZ = TruncateDecode(maxBitsPerComponent, source);

    if (bitsZ)
    {
        auto zz = source.Read(bitsZ - 1);
        zz |= 1 << (bitsZ - 1);
        result.z = ZigZag(zz);
    }

    return result;
}

// //////////////////////////////////////////////////////

unsigned BitVector3BitCountEncode(
        IntVec3 vec,
        unsigned maxMagnitude,
        BitStream& target)
{
    unsigned bitsUsed = 0;

    // +1 for the sign bit.
    unsigned maxBitsPerComponent = 1 + MinBits(maxMagnitude);

    auto zx = ZigZag(vec.x);
    auto zy = ZigZag(vec.y);
    auto zz = ZigZag(vec.z);

    auto minBits = std::max(MinBits(zx), MinBits(zy));
    minBits = std::max(minBits, MinBits(zz));
    auto prefixBitCount = maxBitsPerComponent - minBits;
    maxBitsPerComponent = minBits;

    while (prefixBitCount--)
    {
        target.Write(0,1);
        bitsUsed++;
    }

    // we know  the prefix count in advance, so
    // no need to terminate if we reach max.
    if (minBits)
    {
        target.Write(1,1);
        bitsUsed++;
    }

    if (!maxBitsPerComponent)
    {
        return bitsUsed;
    }

    // //////////////////////////////////////////

    assert(zx < (1u << maxBitsPerComponent));

    target.Write(zx, maxBitsPerComponent);
    bitsUsed += maxBitsPerComponent;

    // //////////////////////////////////////////

    assert(zy < (1u << maxBitsPerComponent));

    target.Write(zy, maxBitsPerComponent);
    bitsUsed += maxBitsPerComponent;

    // //////////////////////////////////////////

    assert(zz < (1u << maxBitsPerComponent));

    target.Write(zz, maxBitsPerComponent);
    bitsUsed += maxBitsPerComponent;

    return bitsUsed;
}

IntVec3 BitVector3BitCountDecode(
        unsigned maxMagnitude,
        BitStream& source)
{
    IntVec3 result = {0,0,0};

    unsigned maxBitsPerComponent = 1 + MinBits(maxMagnitude);

    unsigned prefixCount = 0;
    while (!source.Read(1))
    {
        prefixCount++;

        if (prefixCount == maxBitsPerComponent)
        {
            break;
        }
    }

    maxBitsPerComponent -= prefixCount;

    if (!maxBitsPerComponent)
    {
        return result;
    }

    // //////////////////////////////////////////

    auto zx = source.Read(maxBitsPerComponent);
    result.x = ZigZag(zx);

    // //////////////////////////////////////////

    auto zy = source.Read(maxBitsPerComponent);
    result.y = ZigZag(zy);

    // //////////////////////////////////////////

    auto zz = source.Read(maxBitsPerComponent);
    result.z = ZigZag(zz);

    return result;
}

// //////////////////////////////////////////////////////

unsigned BitVector3BitCountZigZagEncode(
        IntVec3 toEncode,
        IntVec3 base,
        unsigned maxMagnitude,
        BitStream& target)
{
    unsigned bitsUsed = 0;
    unsigned maxBitsPerComponent = MinBits(maxMagnitude);

    auto zx = ZigZagEncode(toEncode.x, base.x, maxBitsPerComponent);
    auto zy = ZigZagEncode(toEncode.y, base.y, maxBitsPerComponent);
    auto zz = ZigZagEncode(toEncode.z, base.z, maxBitsPerComponent);

    auto minBits = std::max(MinBits(zx), MinBits(zy));
    minBits = std::max(minBits, MinBits(zz));
    auto prefixBitCount = maxBitsPerComponent - minBits;
    maxBitsPerComponent = minBits;

    while (prefixBitCount--)
    {
        target.Write(0,1);
        bitsUsed++;
    }
    target.Write(1,1);
    bitsUsed++;

    if (!maxBitsPerComponent)
    {
        return bitsUsed;
    }

    // //////////////////////////////////////////

    assert(zx < (1u << maxBitsPerComponent));

    target.Write(zx, maxBitsPerComponent);
    bitsUsed += maxBitsPerComponent;

    // //////////////////////////////////////////

    assert(zy < (1u << maxBitsPerComponent));

    target.Write(zy, maxBitsPerComponent);
    bitsUsed += maxBitsPerComponent;

    // //////////////////////////////////////////

    assert(zz < (1u << maxBitsPerComponent));

    target.Write(zz, maxBitsPerComponent);
    bitsUsed += maxBitsPerComponent;

    return bitsUsed;
}

IntVec3 BitVector3BitCountZigZagDecode(
        IntVec3 base,
        unsigned maxMagnitude,
        BitStream& source)
{
    IntVec3 result = {0,0,0};

    unsigned maxBitsPerComponent = MinBits(maxMagnitude);

    unsigned prefixCount = 0;
    while (!source.Read(1))
    {
        prefixCount++;
    }

    maxBitsPerComponent -= prefixCount;

    if (!maxBitsPerComponent)
    {
        result.x = ZigZagDecode(0, base.x, maxBitsPerComponent + prefixCount);
        result.y = ZigZagDecode(0, base.y, maxBitsPerComponent + prefixCount);
        result.z = ZigZagDecode(0, base.z, maxBitsPerComponent + prefixCount);
        return result;
    }

    // //////////////////////////////////////////

    auto zx = source.Read(maxBitsPerComponent);
    result.x = ZigZagDecode(zx, base.x, maxBitsPerComponent + prefixCount);

    // //////////////////////////////////////////

    auto zy = source.Read(maxBitsPerComponent);
    result.y = ZigZagDecode(zy, base.y, maxBitsPerComponent + prefixCount);

    // //////////////////////////////////////////

    auto zz = source.Read(maxBitsPerComponent);
    result.z = ZigZagDecode(zz, base.z, maxBitsPerComponent + prefixCount);

    return result;
}

// //////////////////////////////////////////////////////

unsigned GafferMinThresholdBits = 3;

unsigned BitVector3ModifiedGafferEncode(
        IntVec3 toEncode,
        IntVec3 base,
        unsigned maxMagnitude,
        RangeBits ranges,
        BitStream& target)
{
    unsigned max        = MaxRange(ranges);
    unsigned maxBits    = MinBits(maxMagnitude);
    unsigned bitsUsed   = 0;

    unsigned zx = ZigZagEncode(toEncode.x, base.x, maxBits);
    unsigned zy = ZigZagEncode(toEncode.y, base.y, maxBits);
    unsigned zz = ZigZagEncode(toEncode.z, base.z, maxBits);

    auto minBits = std::max(MinBits(zx), MinBits(zy));
    minBits = std::max(minBits, MinBits(zz));

    if (minBits < GafferMinThresholdBits)
    {
        target.Write(1 ,1);
        bitsUsed++;

        auto prefixBitCount = GafferMinThresholdBits - minBits - 1;
        auto maxBitsPerComponent = minBits;

        while (prefixBitCount--)
        {
            target.Write(0,1);
            bitsUsed++;
        }

        // we know  the prefix count in advance, so
        // no need to terminate if we reach max.
        if (minBits)
        {
            target.Write(1,1);
            bitsUsed++;
        }

        if (!maxBitsPerComponent)
        {
            return bitsUsed;
        }

        // //////////////////////////////////////////

        assert(zx < (1u << maxBitsPerComponent));

        target.Write(zx, maxBitsPerComponent);
        bitsUsed += maxBitsPerComponent;

        // //////////////////////////////////////////

        assert(zy < (1u << maxBitsPerComponent));

        target.Write(zy, maxBitsPerComponent);
        bitsUsed += maxBitsPerComponent;

        // //////////////////////////////////////////

        assert(zz < (1u << maxBitsPerComponent));

        target.Write(zz, maxBitsPerComponent);
        bitsUsed += maxBitsPerComponent;

        return bitsUsed;
    }

    auto WriteMax = [&]() -> unsigned
    {
        target.Write(0, 1);
        ++bitsUsed;

        target.Write(1, 1);
        ++bitsUsed;

        target.Write(zx, maxBits);
        target.Write(zy, maxBits);
        target.Write(zz, maxBits);

        bitsUsed += maxBits * 3;

        return bitsUsed;
    };

    if ((zx >= max) || (zy >= max) || (zz >= max))
    {
        return WriteMax();
    }

    BitStream testEncode;

    GaffersRangeEncode(ranges, zx, testEncode);
    GaffersRangeEncode(ranges, zy, testEncode);
    GaffersRangeEncode(ranges, zz, testEncode);

    auto size = testEncode.Bits() ;

    if (size < 3 * maxBits)
    {
        target.Write(0, 1);
        ++bitsUsed;

        target.Write(0, 1);
        ++bitsUsed;

        target.Write(testEncode);

        return size + bitsUsed;
    }

    return WriteMax();
}


IntVec3 BitVector3ModifiedGafferDecode(
        IntVec3 base,
        unsigned maxMagnitude,
        RangeBits ranges,
        BitStream& source)
{
    unsigned maxBits                = MinBits(maxMagnitude);
    unsigned maxBitsPerComponent    = GafferMinThresholdBits - 1;

    auto UnZigZag = [&base, &maxBits](unsigned x, unsigned y, unsigned z) -> IntVec3
    {
        IntVec3 vec;

        vec.x = ZigZagDecode(x, base.x, maxBits);
        vec.y = ZigZagDecode(y, base.y, maxBits);
        vec.z = ZigZagDecode(z, base.z, maxBits);

        return vec;
    };

    auto doMin = source.Read(1);

    if (doMin)
    {
        unsigned prefixCount = 0;
        while (!source.Read(1))
        {
            prefixCount++;

            if (prefixCount == (GafferMinThresholdBits - 1))
            {
                break;
            }
        }

        maxBitsPerComponent -= prefixCount;

        if (!maxBitsPerComponent)
        {
            return UnZigZag(0, 0, 0);
        }

        auto zx = source.Read(maxBitsPerComponent);
        auto zy = source.Read(maxBitsPerComponent);
        auto zz = source.Read(maxBitsPerComponent);

        return UnZigZag(zx, zy, zz);
    }

    auto doMax = source.Read(1);

    if (doMax)
    {
        auto zx = source.Read(maxBits);
        auto zy = source.Read(maxBits);
        auto zz = source.Read(maxBits);

        return UnZigZag(zx, zy, zz);
    }

    auto zx = GaffersRangeDecode(ranges, source);
    auto zy = GaffersRangeDecode(ranges, source);
    auto zz = GaffersRangeDecode(ranges, source);

    return UnZigZag(zx, zy, zz);
}

// //////////////////////////////////////////////////////

unsigned GafferEncode(
        IntVec3 deltaData,
        RangeBits ranges,
        BitStream& target)
{
    unsigned bitsUsed   = 0;

    unsigned zx = ZigZag(deltaData.x);
    unsigned zy = ZigZag(deltaData.y);
    unsigned zz = ZigZag(deltaData.z);

    auto minBits = std::max(MinBits(zx), MinBits(zy));
    minBits = std::max(minBits, MinBits(zz));

    if (minBits < GafferMinThresholdBits)
    {
        target.Write(1, 1);
        ++bitsUsed;

        target.Write(zx, GafferMinThresholdBits - 1);
        target.Write(zy, GafferMinThresholdBits - 1);
        target.Write(zz, GafferMinThresholdBits - 1);
        bitsUsed += 3 * (GafferMinThresholdBits - 1);

        return bitsUsed;
    }

    BitStream testEncode;

    GaffersRangeEncode(ranges, zx, testEncode);
    GaffersRangeEncode(ranges, zy, testEncode);
    GaffersRangeEncode(ranges, zz, testEncode);

    auto size = testEncode.Bits();

    target.Write(0, 1);
    ++bitsUsed;

    target.Write(testEncode);

    return size + bitsUsed;
}


IntVec3 GafferDecode(
        RangeBits ranges,
        BitStream& source)
{
    auto UnZigZag = [](unsigned x, unsigned y, unsigned z) -> IntVec3
    {
        IntVec3 vec;

        vec.x = ZigZag(x);
        vec.y = ZigZag(y);
        vec.z = ZigZag(z);

        return vec;
    };

    auto doMin = source.Read(1);

    if (doMin)
    {
        auto bitsToRead = GafferMinThresholdBits - 1;
        auto zx = source.Read(bitsToRead);
        auto zy = source.Read(bitsToRead);
        auto zz = source.Read(bitsToRead);

        return UnZigZag(zx, zy, zz);
    }

    auto zx = GaffersRangeDecode(ranges, source);
    auto zy = GaffersRangeDecode(ranges, source);
    auto zz = GaffersRangeDecode(ranges, source);

    return UnZigZag(zx, zy, zz);
}

// //////////////////////////////////////////////////////

void BitVector3Tests()
{
    using namespace std::placeholders;

    struct Test
    {
        std::function<unsigned(IntVec3, unsigned, BitStream&)> encode;
        std::function<IntVec3(unsigned, BitStream&)> decode;
        std::function<bool(int, int, int, int)> valid;

    };

    auto Valid_gaffer = [](int i, int j, int k, int) -> bool
    {
        IntVec3 qs
        {
            i*i,
            j*j,
            k*k,
        };

        if (qs.x < qs.z)
        {
            qs = {qs.z, qs.y, qs.x};
        }
        if (qs.y < qs.z)
        {
            qs = {qs.x, qs.z, qs.y};
        }

        auto second = Max_gaffer_value(sqrt(qs.x));
        auto third = Max_gaffer_value(sqrt(qs.x),sqrt(qs.y));

        if (sqrt(qs.y) > second)
        {
            return false;
        }
        if (sqrt(qs.z) > third)
        {
            return false;
        }

        auto largest = qs.x;
        if (qs.y > largest)
        {
            largest = qs.y;
        }
        if (qs.z > largest)
        {
            largest = qs.z;
        }

        auto mag = sqrt(qs.x + qs.y + qs.z + largest);
        auto max = q_max;

        return (mag < max);
    };

    auto Valid_position = [](int i, int j, int k, int max) -> bool
    {
        auto mag = sqrt(i*i + j*j + k*k);

        return (mag <= max);
    };

    std::vector<Test> tests;

    tests.push_back(
    {
        std::bind(Sorted_no_bit_count_Encode, _1, _2, Use_magnitude_as::Gaffer_Encode, _3),
        std::bind(Sorted_no_bit_count_Decode, _1, Use_magnitude_as::Gaffer_Encode, _2),
        Valid_gaffer,
    });

    tests.push_back(
    {
        std::bind(Sorted_no_bit_count_Encode, _1, _2, Use_magnitude_as::Vector_Magnitude, _3),
        std::bind(Sorted_no_bit_count_Decode, _1, Use_magnitude_as::Vector_Magnitude, _2),
        Valid_position,
    });

    tests.push_back(
    {
        std::bind(Sorted_no_bit_count_Encode, _1, _2, Use_magnitude_as::Constant, _3),
        std::bind(Sorted_no_bit_count_Decode, _1, Use_magnitude_as::Constant, _2),
        Valid_position,
    });

    tests.push_back(
    {
        std::bind(BitVector3SortedEncode, _1, _2, Use_magnitude_as::Vector_Magnitude, _3),
        std::bind(BitVector3SortedDecode, _1, Use_magnitude_as::Vector_Magnitude, _2),
        Valid_position,
    });

    tests.push_back(
    {
        std::bind(BitVector3SortedEncode, _1, _2, Use_magnitude_as::Constant, _3),
        std::bind(BitVector3SortedDecode, _1, Use_magnitude_as::Constant, _2),
        Valid_position,
    });

    tests.push_back(
    {
        BitVector3UnrelatedEncode,
        BitVector3UnrelatedDecode,
        Valid_position,
    });

    tests.push_back(
    {
        BitVector3BitCountEncode,
        BitVector3BitCountDecode,
        Valid_position,
    });

    tests.push_back(
    {
        BitVector3TruncatedEncode,
        BitVector3TruncatedDecode,
        Valid_position,
    });

    tests.push_back(
    {
        BitVector3Encode2BitExpPrefix,
        BitVector3Decode2BitExpPrefix,
        Valid_position,
    });

    tests.push_back(
    {
        BitVector3Encode,
        BitVector3Decode,
        Valid_position,
    });

    int const max = static_cast<int>((MaxPositionChangePerSnapshot) * 6 + 1);
    for (const auto& test : tests)
    {
        auto values =
        {
            IntVec3{   -256,   5,     -6},
            IntVec3{   -1,    -32,   -255},
            IntVec3{   -155,   58,   -228},
            IntVec3{    68,   -32,    68},
            IntVec3{    68,    68,   -32},
            IntVec3{   -1,    -1586,  0},
            IntVec3{   -1634, -106,   17},
            IntVec3{   -1634, -106,   0},
            IntVec3{    0,     0,     1},
            IntVec3{    0,     0,     0},
        };

        for (const auto& data : values)
        {
            BitStream encoded;

            if (!test.valid(data.x, data.y, data.z, max))
            {
                continue;
            }

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
                    if (!test.valid(i, j, k, max))
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

    // //////////////////////////////////////////////////

    {
        auto zzMax = 16;
        for (auto i = 0; i < zzMax; ++i)
        {
            for (auto j = 0; j < zzMax; ++j)
            {
                for (auto k = 0; k < zzMax; ++k)
                {
                    for (auto b = 0; b < zzMax; ++b)
                    {
                        IntVec3 base =
                        {
                            b,
                            (zzMax - b) - 1,
                            b / 2,
                        };

                        IntVec3 target =
                        {
                            i,
                            j,
                            k,
                        };

                        BitStream encoded;

                        auto maxMag = (1u << 4) - 1;

                        BitVector3BitCountZigZagEncode(
                            target,
                            base,
                            maxMag,
                            encoded);

                        encoded.Reset();

                        auto decoded = BitVector3BitCountZigZagDecode(
                            base,
                            maxMag,
                            encoded);

                        assert(i == decoded.x);
                        assert(j == decoded.y);
                        assert(k == decoded.z);
                    }
                }
            }
        }
    }

    // //////////////////////////////////////////////////

    {
        auto testValues = {0, 12, 66, 156, 256, 289, 511};
        auto ranges = RangeBits{5,6,7};

        for (const auto& a : testValues)
        {
            for (const auto& b : testValues)
            {
                for (const auto& c : testValues)
                {
                    IntVec3 base{a,b,c};
                    IntVec3 target{c,a,b};

                    {
                        BitStream encoded;

                        BitVector3ModifiedGafferEncode(
                            target,
                            base,
                            511,
                            ranges,
                            encoded);

                        encoded.Reset();

                        auto decoded = BitVector3ModifiedGafferDecode(
                            base,
                            511,
                            ranges,
                            encoded);

                        assert(decoded.x == target.x);
                        assert(decoded.y == target.y);
                        assert(decoded.z == target.z);
                    }

                    {
                        BitStream encoded;

                        IntVec3 deltaData
                        {
                            target.x - base.x,
                            target.y - base.y,
                            target.z - base.z,
                        };

                        auto zx = ZigZag(deltaData.x);
                        auto zy = ZigZag(deltaData.y);
                        auto zz = ZigZag(deltaData.z);

                        auto max = MaxRange(ranges);

                        if ((zx < max) && (zy < max) && (zz < max))
                        {
                            GafferEncode(
                                deltaData,
                                ranges,
                                encoded);

                            encoded.Reset();

                            auto decoded = GafferDecode(
                                ranges,
                                encoded);

                            decoded.x += base.x;
                            decoded.y += base.y;
                            decoded.z += base.z;

                            assert(decoded.x == target.x);
                            assert(decoded.y == target.y);
                            assert(decoded.z == target.z);
                        }
                    }
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

// stats.ones	241612	9.59%
// stats.zeros	2276683	90.41%

const unsigned Simple_probability_one =
        static_cast<unsigned>(
            round((1.0 - 0.9041) * Range_types::TOTAL_RANGE));

BitStream RangeEncodeSimpleEncode(BitStream data)
{
    auto size = data.Bits();

    if (size < 2)
    {
        return data;
    }

    Bytes result_buffer;

    {
        Range_coders::Encoder range(result_buffer);
        Range_coders::Binary_encoder coder(range);

        data.Reset();
        for (unsigned i = 0; i < size; ++i)
        {
            auto value = data.Read(1);

            coder.Encode(value, Simple_probability_one);
        }
    }

    BitStream done = {result_buffer};
    done.SetOffset(8 * result_buffer.size());

    return done;
}

BitStream RangeEncodeSimpleDecode(BitStream& data, unsigned targetBits = 0)
{
    auto raw_data = data.Data();

    Range_coders::Decoder range(raw_data);
    Range_coders::Binary_decoder coder(range);

    BitStream result;

    while (result.Bits() < targetBits)
    {
        auto bit = coder.Decode(Simple_probability_one);

        result.Write(bit, 1);
    }

    auto bitsUsed = 8 * coder.FlushAndGetBytesRead();

    data.SetOffset(bitsUsed);

    return result;
}

// //////////////////////////////////////////////////////

// stats.zero_zero	2175025	95.65%
// stats.zero_one	98923	4.35%
// stats.one_one	140024	57.97%
// stats.one_zero	101528	42.03%

const unsigned Smarter_zeros_probability_one =
        static_cast<unsigned>(
            round((1.0 - 0.9565) * Range_types::TOTAL_RANGE));

const unsigned Smarter_ones_probability_one =
        static_cast<unsigned>(
            round((1.0 - 0.4203) * Range_types::TOTAL_RANGE));

BitStream RangeEncodeSmarterEncode(BitStream data)
{
    auto size = data.Bits();

    if (size < 2)
    {
        return data;
    }

    Range_types::Bytes result_buffer;

    {
        Range_coders::Encoder range(result_buffer);
        Range_coders::Binary_encoder coder(range);

        data.Reset();
        unsigned last = 1;
        for (unsigned i = 0; i < size; ++i)
        {
            auto value = data.Read(1);

            if (last == 0)
            {
                coder.Encode(value, Smarter_zeros_probability_one);
            }
            else
            {
                coder.Encode(value, Smarter_ones_probability_one);
            }

            last = value;
        }
    }

    BitStream done = {result_buffer};
    done.SetOffset(8 * result_buffer.size());

    return done;
}

BitStream RangeEncodeSmarterDecode(BitStream& data, unsigned targetBits = 0)
{
    auto raw_data = data.Data();

    Range_coders::Decoder range(raw_data);
    Range_coders::Binary_decoder coder(range);

    BitStream result;

    unsigned last = 1;
    while (result.Bits() < targetBits)
    {
        auto probability_one =
            (last == 0) ?
            Smarter_zeros_probability_one :
            Smarter_ones_probability_one;

        auto bit = coder.Decode(probability_one);

        result.Write(bit, 1);

        last = bit;
    }

    auto bitsUsed = 8 * coder.FlushAndGetBytesRead();

    data.SetOffset(bitsUsed);

    return result;
}

// //////////////////////////////////////////////////////

static const unsigned Simple_inertia_bits = 2;

BitStream RangeEncodeSimpleAdaptiveEncode(BitStream data)
{
    auto size = data.Bits();

    if (size < 2)
    {
        return data;
    }

    Range_types::Bytes result_buffer;

    {
        Range_coders::Encoder range(result_buffer);
        Range_coders::Binary_encoder coder(range);
        Range_models::Binary
            model(Simple_inertia_bits, Simple_probability_one);

        data.Reset();
        for (unsigned i = 0; i < size; ++i)
        {
            auto value = data.Read(1);

            model.Encode(coder, value);
        }
    }

    BitStream done = {result_buffer};
    done.SetOffset(8 * result_buffer.size());

    return done;
}

BitStream RangeEncodeSimpleAdaptiveDecode(BitStream& data, unsigned targetBits = 0)
{
    auto raw_data = data.Data();

    Range_coders::Decoder range(raw_data);
    Range_coders::Binary_decoder coder(range);
    Range_models::Binary
        model(Simple_inertia_bits, Simple_probability_one);

    BitStream result;

    while (result.Bits() < targetBits)
    {
        auto bit = model.Decode(coder);

        result.Write(bit, 1);
    }

    auto bitsUsed = 8 * coder.FlushAndGetBytesRead();

    data.SetOffset(bitsUsed);

    return result;
}

//// //////////////////////////////////////////////////////


static const unsigned Smarter_inertia_bits_zero = 5;
static const unsigned Smarter_inertia_bits_one = 3;

BitStream RangeEncodeSmarterAdaptiveEncode(BitStream data)
{
    auto size = data.Bits();

    if (size < 2)
    {
        return data;
    }

    Range_types::Bytes result_buffer;

    {
        Range_coders::Encoder range(result_buffer);
        Range_coders::Binary_encoder coder(range);

        // Phoar, the models are getting complex now...
        Range_models::Binary_history<
            Range_models::Binary,
            Range_models::Binary> model(
                1,
                {Smarter_inertia_bits_zero, Smarter_zeros_probability_one},
                {Smarter_inertia_bits_one, Smarter_ones_probability_one});

        data.Reset();
        for (unsigned i = 0; i < size; ++i)
        {
            auto value = data.Read(1);

            model.Encode(coder, value);
        }
    }

    BitStream done = {result_buffer};
    done.SetOffset(8 * result_buffer.size());

    return done;
}

BitStream RangeEncodeSmarterAdaptiveDecode(BitStream& data, unsigned targetBits = 0)
{
    auto raw_data = data.Data();

    Range_coders::Decoder range(raw_data);
    Range_coders::Binary_decoder coder(range);
    Range_models::Binary_history<
        Range_models::Binary,
        Range_models::Binary> model(
            1,
            {Smarter_inertia_bits_zero, Smarter_zeros_probability_one},
            {Smarter_inertia_bits_one, Smarter_ones_probability_one});

    BitStream result;

    while (result.Bits() < targetBits)
    {
        auto bit = model.Decode(coder);

        result.Write(bit, 1);
    }

    auto bitsUsed = 8 * coder.FlushAndGetBytesRead();

    data.SetOffset(bitsUsed);

    return result;
}

// //////////////////////////////////////////////////////

void RunLengthTests()
{
    std::vector<std::function<void(BitStream)>> tests;

    tests.push_back([](BitStream data)
    {
        auto bits = data.Bits();
        auto encoded = RangeEncodeSimpleAdaptiveEncode(data);
        auto bits_encoded = encoded.Bits();
        auto decoded = RangeEncodeSimpleAdaptiveDecode(encoded, bits);
        auto bits_consumed = encoded.Bits();

        assert(bits_encoded == bits_consumed);
        assert(data == decoded);
    });

    tests.push_back([](BitStream data)
    {
        auto bits = data.Bits();
        auto encoded = RangeEncodeSmarterAdaptiveEncode(data);
        auto bits_encoded = encoded.Bits();
        auto decoded = RangeEncodeSmarterAdaptiveDecode(encoded, bits);
        auto bits_consumed = encoded.Bits();

        assert(bits_encoded == bits_consumed);
        assert(data == decoded);
    });

    tests.push_back([](BitStream data)
    {
        auto bits = data.Bits();
        auto encoded = RangeEncodeSmarterEncode(data);
        auto bits_encoded = encoded.Bits();
        auto decoded = RangeEncodeSmarterDecode(encoded, bits);
        auto bits_consumed = encoded.Bits();

        assert(bits_encoded == bits_consumed);
        assert(data == decoded);
    });

    tests.push_back([](BitStream data)
    {
        auto bits = data.Bits();
        auto encoded = RangeEncodeSimpleEncode(data);
        auto bits_encoded = encoded.Bits();
        auto decoded = RangeEncodeSimpleDecode(encoded, bits);
        auto bits_consumed = encoded.Bits();

        assert(bits_encoded == bits_consumed);
        assert(data == decoded);
    });

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

    {
        auto four_in_the_middle = ByteVector(113, 0);
        four_in_the_middle[27] = 4;
        testData.push_back(BitStream(four_in_the_middle, 901));
    }

    {
        auto one_zeros = ByteVector(113, 0);
        one_zeros[0] = 1;
        testData.push_back(BitStream(one_zeros, 901));
    }

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
    QuatPacker              quatPacker;
};

// //////////////////////////////////////////////////////

std::vector<uint8_t> EncodeStats(
    const Frame& base,
    const Frame& target,
    unsigned frameDelta,
    Stats& stats,
    Config config)
{
    const auto count = base.size();

    assert(count > 0);
    assert(count == target.size());

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

    int last = -1;

    for (size_t i = 0; i < count; ++i)
    {
        if (base[i] == target[i])
        {
            changed.Write(0, 1);
            ++runLength0;

            ++(stats.changed1DistanceRunHistogram[runLength1]);
            runLength1 = 0;

            {
                ++stats.zeros;

                if (last >= 0)
                {
                    if (last == 0)
                    {
                        ++stats.zero_zero;
                    }
                    else
                    {
                        ++stats.one_zero;
                    }
                }

                last = 0;
            }
        }
        else
        {
            {
                ++stats.ones;

                if (last >= 0)
                {
                    if (last == 1)
                    {
                        ++stats.one_one;
                    }
                    else
                    {
                        ++stats.zero_one;
                    }
                }

                last = 1;
            }


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

            if (!quatSame)
            {
                auto bg = Gaffer
                {
                    static_cast<unsigned>(base[i].orientation_largest),
                    base[i].orientation_a,
                    base[i].orientation_b,
                    base[i].orientation_c,
                };
                auto tg = Gaffer
                {
                    static_cast<unsigned>(target[i].orientation_largest),
                    target[i].orientation_a,
                    target[i].orientation_b,
                    target[i].orientation_c,
                };
                auto gb1 = ConvertGaffer2(bg);
                auto gt1 = ConvertGaffer2(tg);
                auto rb = ConvertGaffer2(gb1);
                auto rt = ConvertGaffer2(gt1);

                assert(bg.a == rb.a);
                assert(bg.b == rb.b);
                assert(bg.c == rb.c);
                assert(bg.largest_index == rb.largest_index);

                assert(tg.a == rt.a);
                assert(tg.b == rt.b);
                assert(tg.c == rt.c);
                assert(tg.largest_index == rt.largest_index);

                // RAM: now for the fun stuff.
                {
                    {
                        static float max = 0;
                        static float min = 10000000;

                        auto base_quat = ConvertGaffer2(bg);
                        auto target_quat = ConvertGaffer2(rt);

                        auto r = R(base_quat, target_quat);
                        auto rotor = to_rotor(r);
                        auto mag_squared = Magnitude_squared(rotor);
                        auto mag = sqrt(mag_squared);

                        // Hmm, seem we get the stupid mags due
                        // to rotating between itself (from q to -q roughly).
                        auto target_quat_neg = Mul(target_quat, -1.0f);
                        auto r2 = R(base_quat, target_quat_neg);
                        auto rotor2 = to_rotor(r2);
                        auto mag_squared2 = Magnitude_squared(rotor2);
                        auto mag2 = sqrt(mag_squared2);

                        bool other = false;
                        auto to_encode = rotor;

                        if (mag > mag2)
                        {
                            to_encode = rotor2;
                            mag = mag2;
                            other = true;
                        }

                        if (mag > max)
                        {
                            max = mag;
                            if (other)
                            {
                                printf("Max: %f ***\n", mag);
                            }
                            else
                            {
                                printf("Max: %f\n", mag);
                            }
                        }
                        if (mag < min)
                        {
                            min = mag;
                            printf("Min: %f\n", mag);
                        }

                        // Ok, lets get stats assuming our multiple is
                        // 2^14 == 8192.
//                        {
//                            static const float ROTOR_MULTIPLY = 10000.0;
//                            IntVec3 coded =
//                            {
//                                static_cast<int>(round(to_encode[0] * ROTOR_MULTIPLY)),
//                                static_cast<int>(round(to_encode[1] * ROTOR_MULTIPLY)),
//                                static_cast<int>(round(to_encode[2] * ROTOR_MULTIPLY)),
//                            };

//                            Rotor decoded =
//                            {
//                                coded.x / ROTOR_MULTIPLY,
//                                coded.y / ROTOR_MULTIPLY,
//                                coded.z / ROTOR_MULTIPLY,
//                            };

//                            auto decoded_quat = to_quat(decoded);
//                            auto result_target_quat = Mul(decoded_quat, base_quat);
//                            auto result = ConvertGaffer2(result_target_quat);

//                            assert(result.largest_index == tg.largest_index);
//                            assert(result.a == tg.a);
//                            assert(result.b == tg.b);
//                            assert(result.c == tg.c);

//                            bool wtf = false;
//                            if (MinBits(ZigZag(coded.x)) > 12)
//                            {
//                                wtf = true;
//                            }
//                            if (MinBits(ZigZag(coded.y)) > 12)
//                            {
//                                wtf = true;
//                            }
//                            if (MinBits(ZigZag(coded.z)) > 12)
//                            {
//                                wtf = true;
//                            }

//                            // Hint, it doesn't do well for us, having to
//                            // encode over 11 bits at times.
//                            stats.rotor_bits.Update(MinBits(ZigZag(coded.x)));
//                            stats.rotor_bits.Update(MinBits(ZigZag(coded.y)));
//                            stats.rotor_bits.Update(MinBits(ZigZag(coded.z)));

//                            if (wtf)
//                            {
//                                wtf = false;
//                            }
//                        }

                        // ok, to figure out min bits needed find the smallest
                        // value to be encoded by following the multiplication.
                        // NUP, fails just like th pervious method.
//                        auto min_r = 100000.0;
//                        for (auto e : r)
//                        {
//                            if (std::abs(e) > 0.00000001)
//                            {
//                                if (std::abs(e) < min_r)
//                                {
//                                    min_r = std::abs(e);
//                                }
//                            }
//                        }


                        // Lets try another approach.
//                        auto mults =
//                        {
//                            128.0f,
//                            256.0f,
//                            384.0f,
//                            512.0f,
//                            640.0f,
//                            768.0f,
//                            1024.0f,
//                            1152.0f,
//                            1180.0f,
//                            4352.0f,
//                        };

//                        auto mults =
//                        {
//                            256.0f,
//                            384.0f,
//                            512.0f,
//                            640.0f,
//                            768.0f,
//                            1180.0f,
//                            4352.0f,
//                        };

                        // RAM: Another suggestion. Get the min value
                        // of the target quat to encode, maybe that dicates
                        // bits?

                        // ok, lets see how much bits we need.
                        static float MAX_MULTIPLY = 1.0f;
                        bool not_enough = true;

                        assert(
                            (to_encode[0] != 0) ||
                            (to_encode[1] != 0) ||
                            (to_encode[2] != 0));

                        // Well, that didn't work. Always guessed min of 10 bits
                        // and even guessed 19 bits at one point.
//                        // Maybe I only need to multiply so that the smallest
//                        // value truncates to 1.
//                        float min_e = 1000.0f;
//                        for (auto e : to_encode)
//                        {
//                            if (std::abs(e) > 0.000001)
//                            {
//                                if (std::abs(e) < min_e)
//                                {
//                                    min_e = std::abs(e);
//                                }
//                            }
//                        }

//                        //float ROTOR_MULTIPLY = 1.0 / min_e;
//                        float ROTOR_MULTIPLY = 1.49 / min_e;
//                        //ROTOR_MULTIPLY = std::min(ROTOR_MULTIPLY, 512.0f);
//                        while (ROTOR_MULTIPLY < 256.0f)
//                        {
//                            //ROTOR_MULTIPLY = 256.0f;
//                            ROTOR_MULTIPLY *= 10.0f;
//                        }
//                        //ROTOR_MULTIPLY = std::min(ROTOR_MULTIPLY, 256.0f);

                        //float ROTOR_MULTIPLY = 1.49f / min_r;
                        float ROTOR_MULTIPLY = 128.0;
                        while (not_enough)
                        {
                            IntVec3 coded =
                            {
                                static_cast<int>(round(to_encode[0] * ROTOR_MULTIPLY)),
                                static_cast<int>(round(to_encode[1] * ROTOR_MULTIPLY)),
                                static_cast<int>(round(to_encode[2] * ROTOR_MULTIPLY)),
                            };

                            Rotor decoded =
                            {
                                coded.x / ROTOR_MULTIPLY,
                                coded.y / ROTOR_MULTIPLY,
                                coded.z / ROTOR_MULTIPLY,
                            };

                            auto decoded_quat = to_quat(decoded);
                            auto result_target_quat = Mul(decoded_quat, base_quat);
                            auto result = ConvertGaffer2(result_target_quat);

                            if  (
                                    (result.largest_index == tg.largest_index) &&
                                    (result.a == tg.a) &&
                                    (result.b == tg.b) &&
                                    (result.c == tg.c)
                                )
                            {
                                // stats
                                if (ROTOR_MULTIPLY <= 128)
                                {
                                    ++stats.rotor_bits_8;
                                }
                                else if (ROTOR_MULTIPLY <= 256)
                                {
                                    ++stats.rotor_bits_9;
                                }
                                else if (ROTOR_MULTIPLY <= 512)
                                {
                                    ++stats.rotor_bits_10;
                                }
                                else if (ROTOR_MULTIPLY <= 1024)
                                {
                                    ++stats.rotor_bits_11;
                                }
                                else if (ROTOR_MULTIPLY <= 2048)
                                {
                                    ++stats.rotor_bits_12;
                                }
                                else if (ROTOR_MULTIPLY <= 4096)
                                {
                                    ++stats.rotor_bits_13;
                                }
                                else if (ROTOR_MULTIPLY <= 8192)
                                {
                                    ++stats.rotor_bits_14;
                                }
                                else if (ROTOR_MULTIPLY <= 16384)
                                {
                                    ++stats.rotor_bits_15;
                                }
                                else
                                {
                                    ++stats.rotor_bits_wtf;
                                }

                                if (ROTOR_MULTIPLY > MAX_MULTIPLY)
                                {
                                    if (ROTOR_MULTIPLY < 2048)
                                    {
                                        printf("Multiple: %f\n", ROTOR_MULTIPLY);
                                        MAX_MULTIPLY = ROTOR_MULTIPLY;
                                    }
                                    else
                                    {
                                        printf("Mult ***: %f\n", ROTOR_MULTIPLY);
                                    }

                                }
                                not_enough = false;

                                stats.rotor_bits.Update(MinBits(ZigZag(coded.x)));
                                stats.rotor_bits.Update(MinBits(ZigZag(coded.y)));
                                stats.rotor_bits.Update(MinBits(ZigZag(coded.z)));
                            }
                            else
                            {
                                // RAM: rotoro multiply so that the smallest
                                // non-zero value rounds to 1 (eg > 0.55).
                                // treat anything less than 1 / 1 << 16 as zero
                                // Then if that fails, just multiply by 10/8/2?
                                //ROTOR_MULTIPLY *= 1.6;
                                ROTOR_MULTIPLY += 128.0f;
                                //ROTOR_MULTIPLY += 64.0f;
                                //ROTOR_MULTIPLY *= 2;
                            }
                        }
                    }


                    auto rotor = Rotorify(bg, tg);
                    auto result = Rotorify(bg, rotor);

                    float squared =
                        (float(rotor.x) * float(rotor.x)) +
                        (float(rotor.y) * float(rotor.y)) +
                        (float(rotor.z) * float(rotor.z));

                    assert(squared > 0);

                    auto mag = sqrt(squared);

                    assert(mag > 0);

                    stats.rotor.Update(static_cast<unsigned>(round(mag)));

                    assert(result.largest_index == tg.largest_index);
                    assert(result.a == tg.a);
                    assert(result.b == tg.b);
                    assert(result.c == tg.c);
                }



//                if (
//                        (!base[i].orientation_a) ||
//                        (!base[i].orientation_b) ||
//                        (!base[i].orientation_c) ||
//                        (!target[i].orientation_a) ||
//                        (!target[i].orientation_b) ||
//                        (!target[i].orientation_c) ||
//                        (511 == base[i].orientation_a) ||
//                        (511 == base[i].orientation_b) ||
//                        (511 == base[i].orientation_c) ||
//                        (511 == target[i].orientation_a) ||
//                        (511 == target[i].orientation_b) ||
//                        (511 == target[i].orientation_c)
//                        )
//                {
//                    printf("%d,%d,%d,%d,%d,%d\n",
//                           base[i].orientation_a - 256,
//                           base[i].orientation_b - 256,
//                           base[i].orientation_c - 256,
//                           target[i].orientation_a - 256,
//                           target[i].orientation_b - 256,
//                           target[i].orientation_c - 256);

                    // -16, -16, -256 is the max one.

//                }
//                auto b = Gaffer
//                {
//                    static_cast<unsigned>(base[i].orientation_largest),
//                    base[i].orientation_a,
//                    base[i].orientation_b,
//                    base[i].orientation_c,
//                };

//                auto t = Gaffer
//                {
//                    static_cast<unsigned>(target[i].orientation_largest),
//                    target[i].orientation_a,
//                    target[i].orientation_b,
//                    target[i].orientation_c,
//                };

//                auto rotor = Rotorify(b, t);
//                auto result = Rotorify(b, rotor);

//                assert(result.largest_index == static_cast<unsigned>(target[i].orientation_largest));
//                assert(result.a == target[i].orientation_a);
//                assert(result.b == target[i].orientation_b);
//                assert(result.c == target[i].orientation_c);
            }

            const auto vec3BitsUncompressed = (RotationMaxBits * 3);
            stats.quatDeltaUnpackedBitCount +=
                (vec3BitsUncompressed + RotationIndexMaxBits);

            ++stats.quatChanged;
            if (target[i].orientation_largest != base[i].orientation_largest)
            {
                ++stats.largestDifferent;
            }

            auto da = abs(target[i].orientation_a - base[i].orientation_a);
            auto db = abs(target[i].orientation_b - base[i].orientation_b);
            auto dc = abs(target[i].orientation_c - base[i].orientation_c);
            auto dSum = da + db + dc;

            stats.deltaA.Update(da);
            stats.deltaB.Update(db);
            stats.deltaC.Update(dc);

            // yea, ignore the case when only the largest 3 index chagnes.
            if (dSum)
            {
                stats.deltaABC.Update(dSum);
            }

            auto WriteFullQuant = [&]()
            {
                deltas.Write(target[i].orientation_largest, RotationIndexMaxBits);
                deltas.Write(target[i].orientation_a, RotationMaxBits);
                deltas.Write(target[i].orientation_b, RotationMaxBits);
                deltas.Write(target[i].orientation_c, RotationMaxBits);
            };

            switch (config.quatPacker)
            {
                default:
                case QuatPacker::None:
                {
                    if (quatSame)
                    {
                        deltas.Write(0, 1);
                    }
                    else
                    {
                        deltas.Write(1, 1);
                        WriteFullQuant();
                    }

                    break;
                }

                case QuatPacker::BitVector3Unrelated:
                case QuatPacker::BitVector3BitCount:
                case QuatPacker::BitVector3Sorted:
                case QuatPacker::Sorted_no_bit_count:
                {
                    unsigned bitsWritten = 2;
                    deltas.Write(target[i].orientation_largest, RotationIndexMaxBits);

                    auto vec = IntVec3
                    {
                        target[i].orientation_a - base[i].orientation_a,
                        target[i].orientation_b - base[i].orientation_b,
                        target[i].orientation_c - base[i].orientation_c,
                    };

                    BitStream encoded;
                    unsigned codedBits = 0;

                    if (config.quatPacker != QuatPacker::BitVector3BitCount)
                    {
                        switch (config.quatPacker)
                        {
                            default:
                            case QuatPacker::BitVector3Unrelated:
                            {
                                codedBits = BitVector3UnrelatedEncode(
                                    vec,
                                    (1 << RotationMaxBits) - 1,
                                    encoded);

                                break;
                            }

                            case QuatPacker::BitVector3Sorted:
                            {
                                codedBits = BitVector3SortedEncode(
                                    vec,
                                    (1 << RotationMaxBits) - 1,
                                    Use_magnitude_as::Constant,
                                    encoded);

                                break;
                            }

                            case QuatPacker::Sorted_no_bit_count:
                            {
                                codedBits = Sorted_no_bit_count_Encode(
                                    vec,
                                    (1 << RotationMaxBits) - 1,
                                    Use_magnitude_as::Constant,
                                    encoded);

                                break;
                            }
                        }

                        if (codedBits < vec3BitsUncompressed)
                        {
                            deltas.Write(0, 1);
                            deltas.Write(encoded);

                            bitsWritten += 1 + codedBits;
                        }
                        else
                        {
                            deltas.Write(1, 1);

                            if (config.quatPacker == QuatPacker::Sorted_no_bit_count)
                            {
                                IntVec3 full
                                {
                                    target[i].orientation_a - 256,
                                    target[i].orientation_b - 256,
                                    target[i].orientation_c - 256,
                                };

                                BitStream encoded2;

                                auto codedBits = Sorted_no_bit_count_Encode(
                                    full,
                                    (1 << RotationMaxBits) - 1,
                                    Use_magnitude_as::Gaffer_Encode,
                                    encoded2);

                                deltas.Write(encoded2);
                                bitsWritten += 1 + codedBits;
                            }
                            else
                            {
                                deltas.Write(target[i].orientation_a, RotationMaxBits);
                                deltas.Write(target[i].orientation_b, RotationMaxBits);
                                deltas.Write(target[i].orientation_c, RotationMaxBits);

                                bitsWritten += 1 + vec3BitsUncompressed;
                            }
                        }
                    }
                    else
                    {
                        auto a = ZigZag(vec.x);
                        auto b = ZigZag(vec.y);
                        auto c = ZigZag(vec.z);

                        auto ba = MinBits(a);
                        auto bb = MinBits(b);
                        auto bc = MinBits(c);

                        auto minBits = std::max(ba, bb);
                        minBits = std::max(minBits, bc);

                        // Even though I would save 3 bits for 8, it'll
                        // cost me and extra bit for smaller minBits
                        // prefixes, which taking into consideration the
                        // probabilites, ends up using more bits in total.
                        if (minBits >= QuanternionUncompressedBitThreshold)
                        {
                            deltas.Write(1, 1);
                            deltas.Write(target[i].orientation_a, RotationMaxBits);
                            deltas.Write(target[i].orientation_b, RotationMaxBits);
                            deltas.Write(target[i].orientation_c, RotationMaxBits);

                            bitsWritten += 1 + vec3BitsUncompressed;

                            {
                                static unsigned countD = 0;

                                if (countD < 100)
                                {
                                    countD++;

                                    printf("-------------------\n");

                                    printf(
                                        "%d - %d = %d\n",
                                        target[i].orientation_a,
                                        base[i].orientation_a,
                                        target[i].orientation_a - base[i].orientation_a);
                                    printf(
                                        "%d - %d = %d\n",
                                        target[i].orientation_b,
                                        base[i].orientation_b,
                                        target[i].orientation_b - base[i].orientation_b);
                                    printf(
                                        "%d - %d = %d\n",
                                        target[i].orientation_c,
                                        base[i].orientation_c,
                                        target[i].orientation_c - base[i].orientation_c);
                                }
                            }
                        }
                        else
                        {
                            deltas.Write(0, 1);

                            // -2 == -1 for the largest value, then
                            // -1 to account that we already count the sign
                            // bit.
                            codedBits = BitVector3BitCountEncode(
                                vec,
                                (1 << (QuanternionUncompressedBitThreshold - 2)) - 1,
                                encoded);

                            deltas.Write(encoded);
                            bitsWritten += 1 + codedBits;
                        }
                    }

                    // *sniff* worst case == vec3BitsUncompressed + 1.
                    //assert(bitsWritten < vec3BitsUncompressed);
                    stats.quatDeltaPackedBitCount.Update(bitsWritten);
                    break;
                }

                case QuatPacker::BitVector3ModifiedZigZag:
                {
                    unsigned bitsWritten = 2;
                    deltas.Write(target[i].orientation_largest, RotationIndexMaxBits);

                    BitStream encoded;
                    unsigned codedBits = 0;

                    auto a = ZigZagEncode(
                                target[i].orientation_a,
                                base[i].orientation_a,
                                RotationMaxBits);
                    auto b = ZigZagEncode(
                                target[i].orientation_b,
                                base[i].orientation_b,
                                RotationMaxBits);
                    auto c = ZigZagEncode(
                                target[i].orientation_c,
                                base[i].orientation_c,
                                RotationMaxBits);

                    auto ba = MinBits(a);
                    auto bb = MinBits(b);
                    auto bc = MinBits(c);

                    auto minBits = std::max(ba, bb);
                    minBits = std::max(minBits, bc);

                    // Even though I would save 3 bits for 8, it'll
                    // cost me and extra bit for smaller minBits
                    // prefixes, which taking into consideration the
                    // probabilites, ends up using more bits in total.
                    if (minBits >= QuanternionUncompressedBitThreshold)
                    {
                        deltas.Write(1, 1);
                        deltas.Write(target[i].orientation_a, RotationMaxBits);
                        deltas.Write(target[i].orientation_b, RotationMaxBits);
                        deltas.Write(target[i].orientation_c, RotationMaxBits);

                        bitsWritten += 1 + vec3BitsUncompressed;
                    }
                    else
                    {
                        deltas.Write(0, 1);

                        IntVec3 toEncode =
                        {
                            target[i].orientation_a,
                            target[i].orientation_b,
                            target[i].orientation_c,
                        };

                        IntVec3 baseVec =
                        {
                            base[i].orientation_a,
                            base[i].orientation_b,
                            base[i].orientation_c,
                        };

                        // -2 == -1 for the largest value, then
                        // -1 to account that we already count the sign
                        // bit.
                        codedBits = BitVector3BitCountZigZagEncode(
                            toEncode,
                            baseVec,
                            (1 << RotationMaxBits) - 1,
                            encoded);

                        deltas.Write(encoded);
                        bitsWritten += 1 + codedBits;
                    }

                    // *sniff* worst case == vec3BitsUncompressed + 1.
                    //assert(bitsWritten < vec3BitsUncompressed);
                    stats.quatDeltaPackedBitCount.Update(bitsWritten);
                    break;
                }

                case QuatPacker::BitVector3ModifiedGaffer:
                {
                    unsigned bitsWritten = 2;
                    deltas.Write(
                        target[i].orientation_largest,
                        RotationIndexMaxBits);

                    IntVec3 toEncode =
                    {
                        target[i].orientation_a,
                        target[i].orientation_b,
                        target[i].orientation_c,
                    };

                    IntVec3 baseVec =
                    {
                        base[i].orientation_a,
                        base[i].orientation_b,
                        base[i].orientation_c,
                    };

                    BitStream encoded;
                    auto codedBits = BitVector3ModifiedGafferEncode(
                            toEncode,
                            baseVec,
                            (1 << RotationMaxBits) - 1,
                            RangeBits {5, 6, 7},
                            encoded);

                    deltas.Write(encoded);
                    bitsWritten += 1 + codedBits;

                    stats.quatDeltaPackedBitCount.Update(bitsWritten);

                    break;
                }

                case QuatPacker::Gaffer:
                {
                    auto WriteFull = [&]()
                    {
                        deltas.Write(1, 1);

                        WriteFullQuant();

                        stats.quatDeltaPackedBitCount.Update(
                                    1 +
                                    (RotationMaxBits * 3) +
                                    RotationIndexMaxBits);
                    };

                    if  (
                            target[i].orientation_largest !=
                            base[i].orientation_largest
                        )
                    {
                        WriteFull();
                        break;
                    }

                    IntVec3 deltaData
                    {
                        target[i].orientation_a - base[i].orientation_a,
                        target[i].orientation_b - base[i].orientation_b,
                        target[i].orientation_c - base[i].orientation_c,
                    };

                    auto zx = ZigZag(deltaData.x);
                    auto zy = ZigZag(deltaData.y);
                    auto zz = ZigZag(deltaData.z);

                    auto ranges = RangeBits{4, 5, 7};

                    auto max = MaxRange(ranges);

                    if ((zx >= max) || (zy >= max) || (zz >= max))
                    {
                        WriteFull();
                        break;
                    }

                    // using gaffer encode.
                    deltas.Write(0, 1);

                    auto codedBits = GafferEncode(
                            deltaData,
                            ranges,
                            deltas);

                    stats.quatDeltaPackedBitCount.Update(1 + codedBits);
                    break;
                }
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
                    case PosVector3Packer::BitVector3Sorted:
                    case PosVector3Packer::Sorted_no_bit_count:
                    {
                        // GLEE, actually writing a delta this time!
                        auto vec = IntVec3
                        {
                            target[i].position_x - base[i].position_x,
                            target[i].position_y - base[i].position_y,
                            target[i].position_z - base[i].position_z,
                        };

                        unsigned maxMagnitude = 1 + static_cast<unsigned>(MaxPositionChangePerSnapshot * frameDelta);
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

                            case PosVector3Packer::BitVector3Sorted:
                            {
                                bitsWritten = BitVector3SortedEncode(
                                    vec,
                                    maxMagnitude,
                                    Use_magnitude_as::Vector_Magnitude,
                                    deltas);

                                break;
                            }

                            case PosVector3Packer::Sorted_no_bit_count:
                            {
                                bitsWritten = Sorted_no_bit_count_Encode(
                                    vec,
                                    maxMagnitude,
                                    Use_magnitude_as::Vector_Magnitude,
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
    auto range_simple = RangeEncodeSimpleEncode(changed);
    auto range_smarter = RangeEncodeSmarterEncode(changed);
    auto range_simple_adaptive = RangeEncodeSimpleAdaptiveEncode(changed);
    auto range_smarter_adaptive = RangeEncodeSmarterAdaptiveEncode(changed);

    assert((rle.size() * 8) < Cubes);
    assert((bitpack.size() * 8) < Cubes);
    assert((bitexprle.Bits()) < Cubes);
    assert((expbitpack.Bits()) < Cubes);
    assert((bitbitpackfull.Bits()) < Cubes);
    assert((range_simple.Bits()) < Cubes);
    assert((range_smarter.Bits()) < Cubes);
    assert((range_simple_adaptive.Bits()) < Cubes);
    assert((range_smarter_adaptive.Bits()) < Cubes);

    stats.rle.Update(rle.size() * 8);
    stats.bitpack.Update(bitpack.size() * 8);
    stats.bitexprle.Update(bitexprle.Bits());
    stats.bitbitpack.Update(expbitpack.Bits());
    stats.bitbitfullpack.Update(bitbitpackfull.Bits());
    stats.range_simple.Update(range_simple.Bits());
    stats.range_smarter.Update(range_smarter.Bits());
    stats.range_simple_adaptive.Update(range_simple_adaptive.Bits());
    stats.range_smarter_adaptive.Update(range_smarter_adaptive.Bits());

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

            case ChangedArrayEncoding::RangeSimple:
            {
                return RangeEncodeSimpleEncode(changed);
            }

            case ChangedArrayEncoding::RangeSmarter:
            {
                return RangeEncodeSmarterEncode(changed);
            }

            case ChangedArrayEncoding::RangeSimpleAdaptive:
            {
                return RangeEncodeSimpleAdaptiveEncode(changed);
            }

            case ChangedArrayEncoding::RangeSmarterAdaptive:
            {
                return RangeEncodeSmarterAdaptiveEncode(changed);
            }
        }
    }();

    assert(changedCompressed.Bits() <= 901);

    {
        // do some research to see if probabilites of
        // bits to bits are non-random.
        auto countBits = changedCompressed.Bits();

        changedCompressed.Reset();

        auto last = changedCompressed.Read(1);

        while (changedCompressed.Bits() < countBits)
        {
            if (last)
            {
                ++stats.bitStream1;

                last = changedCompressed.Read(1);
                if (last)
                {
                    ++stats.bitStream11;
                }
                else
                {
                    ++stats.bitStream10;
                }
            }
            else
            {
                ++stats.bitStream0;

                last = changedCompressed.Read(1);
                if (last)
                {
                    ++stats.bitStream01;
                }
                else
                {
                    ++stats.bitStream00;
                }
            }
        }
    }

    result.Write(changedCompressed);
    result.Write(deltas);

    result.TrimZerosFromBack();

    return result.Data();
}

// //////////////////////////////////////////////////////

std::vector<uint8_t> Encode(
    const Frame& base,
    const Frame& target,
    unsigned frameDelta,
    Config config =
    {
        ChangedArrayEncoding::None,
        PosVector3Packer::None,
        QuatPacker::None,
    })
{
    const auto count = base.size();

    assert(count > 0);
    assert(count == target.size());

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

            runLength1 = 0;
        }
        else
        {
            changed.Write(1, 1);
            ++runLength1;

            runLength0 = 0;

            assert(target[i].orientation_a >= 0);
            assert(target[i].orientation_b >= 0);
            assert(target[i].orientation_c >= 0);

            bool posSame =
                    (base[i].position_x == target[i].position_x) &&
                    (base[i].position_y == target[i].position_y) &&
                    (base[i].position_z == target[i].position_z);

            bool quatSame =
                    (base[i].orientation_largest == target[i].orientation_largest) &&
                    (base[i].orientation_a == target[i].orientation_a) &&
                    (base[i].orientation_b == target[i].orientation_b) &&
                    (base[i].orientation_c == target[i].orientation_c);

            const auto vec3BitsUncompressed =
                    (RotationMaxBits * 3);

            auto WriteFullQuant = [&]()
            {
                deltas.Write(target[i].orientation_largest, RotationIndexMaxBits);
                deltas.Write(target[i].orientation_a, RotationMaxBits);
                deltas.Write(target[i].orientation_b, RotationMaxBits);
                deltas.Write(target[i].orientation_c, RotationMaxBits);
            };

            switch (config.quatPacker)
            {
                default:
                case QuatPacker::None:
                {
                    if (quatSame)
                    {
                        deltas.Write(0, 1);
                    }
                    else
                    {
                        deltas.Write(1, 1);
                        WriteFullQuant();
                    }

                    break;
                }

                case QuatPacker::BitVector3Unrelated:
                case QuatPacker::BitVector3BitCount:
                case QuatPacker::BitVector3Sorted:
                case QuatPacker::Sorted_no_bit_count:
                {
                    unsigned bitsWritten = 2;
                    deltas.Write(target[i].orientation_largest, RotationIndexMaxBits);

                    auto vec = IntVec3
                    {
                        target[i].orientation_a - base[i].orientation_a,
                        target[i].orientation_b - base[i].orientation_b,
                        target[i].orientation_c - base[i].orientation_c,
                    };

                    BitStream encoded;
                    unsigned codedBits = 0;

                    if (config.quatPacker != QuatPacker::BitVector3BitCount)
                    {                        
                        switch (config.quatPacker)
                        {
                            default:
                            case QuatPacker::BitVector3Unrelated:
                            {
                                codedBits = BitVector3UnrelatedEncode(
                                    vec,
                                    (1 << RotationMaxBits) - 1,
                                    encoded);

                                break;
                            }

                            case QuatPacker::BitVector3Sorted:
                            {
                                codedBits = BitVector3SortedEncode(
                                    vec,
                                    (1 << RotationMaxBits) - 1,
                                    Use_magnitude_as::Constant,
                                    encoded);

                                break;
                            }

                            case QuatPacker::Sorted_no_bit_count:
                            {
                                codedBits = Sorted_no_bit_count_Encode(
                                    vec,
                                    (1 << RotationMaxBits) - 1,
                                    Use_magnitude_as::Constant,
                                    encoded);

                                break;
                            }
                        }

                        if (codedBits < vec3BitsUncompressed)
                        {
                            deltas.Write(0, 1);
                            deltas.Write(encoded);

                            bitsWritten += 1 + codedBits;
                        }
                        else
                        {
                            deltas.Write(1, 1);                            

                            if (config.quatPacker == QuatPacker::Sorted_no_bit_count)
                            {
                                IntVec3 full
                                {
                                    target[i].orientation_a - 256,
                                    target[i].orientation_b - 256,
                                    target[i].orientation_c - 256,
                                };

                                BitStream encoded2;

                                auto codedBits = Sorted_no_bit_count_Encode(
                                    full,
                                    (1 << RotationMaxBits) - 1,
                                    Use_magnitude_as::Gaffer_Encode,
                                    encoded2);

                                // Bother :-(
                                assert(codedBits <= vec3BitsUncompressed);

                                deltas.Write(encoded2);
                                bitsWritten += 1 + codedBits;
                            }
                            else
                            {
                                deltas.Write(target[i].orientation_a, RotationMaxBits);
                                deltas.Write(target[i].orientation_b, RotationMaxBits);
                                deltas.Write(target[i].orientation_c, RotationMaxBits);

                                bitsWritten += 1 + vec3BitsUncompressed;
                            }
                        }
                    }
                    else
                    {
                        auto a = ZigZag(vec.x);
                        auto b = ZigZag(vec.y);
                        auto c = ZigZag(vec.z);

                        auto ba = MinBits(a);
                        auto bb = MinBits(b);
                        auto bc = MinBits(c);

                        auto minBits = std::max(ba, bb);
                        minBits = std::max(minBits, bc);

                        // Even though I would save 3 bits for 8, it'll
                        // cost me and extra bit for smaller minBits
                        // prefixes, which taking into consideration the
                        // probabilites, ends up using more bits in total.
                        if (minBits >= QuanternionUncompressedBitThreshold)
                        {
                            deltas.Write(1, 1);
                            deltas.Write(target[i].orientation_a, RotationMaxBits);
                            deltas.Write(target[i].orientation_b, RotationMaxBits);
                            deltas.Write(target[i].orientation_c, RotationMaxBits);

                            bitsWritten += 1 + vec3BitsUncompressed;
                        }
                        else
                        {
                            deltas.Write(0, 1);

                            // -2 == -1 for the largest value, then
                            // -1 to account that we already count the sign
                            // bit.
                            codedBits = BitVector3BitCountEncode(
                                vec,
                                (1 << (QuanternionUncompressedBitThreshold - 2)) - 1,
                                encoded);

                            deltas.Write(encoded);
                            bitsWritten += 1 + codedBits;
                        }
                    }

                    // *sniff* worst case == vec3BitsUncompressed + 1.
                    //assert(bitsWritten < vec3BitsUncompressed);
                    break;
                }

                case QuatPacker::BitVector3ModifiedZigZag:
                {
                    unsigned bitsWritten = 2;
                    deltas.Write(target[i].orientation_largest, RotationIndexMaxBits);

                    BitStream encoded;
                    unsigned codedBits = 0;

                    auto a = ZigZagEncode(
                                target[i].orientation_a,
                                base[i].orientation_a,
                                RotationMaxBits);
                    auto b = ZigZagEncode(
                                target[i].orientation_b,
                                base[i].orientation_b,
                                RotationMaxBits);
                    auto c = ZigZagEncode(
                                target[i].orientation_c,
                                base[i].orientation_c,
                                RotationMaxBits);

                    auto ba = MinBits(a);
                    auto bb = MinBits(b);
                    auto bc = MinBits(c);

                    auto minBits = std::max(ba, bb);
                    minBits = std::max(minBits, bc);

                    // Even though I would save 3 bits for 8, it'll
                    // cost me and extra bit for smaller minBits
                    // prefixes, which taking into consideration the
                    // probabilites, ends up using more bits in total.
                    if (minBits >= QuanternionUncompressedBitThreshold)
                    {
                        deltas.Write(1, 1);
                        deltas.Write(target[i].orientation_a, RotationMaxBits);
                        deltas.Write(target[i].orientation_b, RotationMaxBits);
                        deltas.Write(target[i].orientation_c, RotationMaxBits);

                        bitsWritten += 1 + vec3BitsUncompressed;
                    }
                    else
                    {
                        deltas.Write(0, 1);

                        IntVec3 toEncode =
                        {
                            target[i].orientation_a,
                            target[i].orientation_b,
                            target[i].orientation_c,
                        };

                        IntVec3 baseVec =
                        {
                            base[i].orientation_a,
                            base[i].orientation_b,
                            base[i].orientation_c,
                        };

                        codedBits = BitVector3BitCountZigZagEncode(
                            toEncode,
                            baseVec,
                            (1 << RotationMaxBits) - 1,
                            encoded);

                        deltas.Write(encoded);
                        bitsWritten += 1 + codedBits;
                    }

                    break;
                }

                case QuatPacker::BitVector3ModifiedGaffer:
                {
                    unsigned bitsWritten = 2;
                    deltas.Write(
                        target[i].orientation_largest,
                        RotationIndexMaxBits);

                    IntVec3 toEncode =
                    {
                        target[i].orientation_a,
                        target[i].orientation_b,
                        target[i].orientation_c,
                    };

                    IntVec3 baseVec =
                    {
                        base[i].orientation_a,
                        base[i].orientation_b,
                        base[i].orientation_c,
                    };

                    BitStream encoded;
                    auto codedBits = BitVector3ModifiedGafferEncode(
                            toEncode,
                            baseVec,
                            (1 << RotationMaxBits) - 1,
                            RangeBits {5, 6, 7},
                            encoded);

                    deltas.Write(encoded);
                    bitsWritten += 1 + codedBits;

                    break;
                }

                case QuatPacker::Gaffer:
                {
                    auto WriteFull = [&]()
                    {
                        deltas.Write(1, 1);

                        WriteFullQuant();
                    };

                    if  (
                            target[i].orientation_largest !=
                            base[i].orientation_largest
                        )
                    {
                        WriteFull();
                        break;
                    }

                    IntVec3 deltaData
                    {
                        target[i].orientation_a - base[i].orientation_a,
                        target[i].orientation_b - base[i].orientation_b,
                        target[i].orientation_c - base[i].orientation_c,
                    };

                    auto zx = ZigZag(deltaData.x);
                    auto zy = ZigZag(deltaData.y);
                    auto zz = ZigZag(deltaData.z);

                    auto ranges = RangeBits{4, 5, 7};

                    auto max = MaxRange(ranges);

                    if ((zx >= max) || (zy >= max) || (zz >= max))
                    {
                        WriteFull();
                        break;
                    }

                    // using gaffer encode.
                    deltas.Write(0, 1);

                    GafferEncode(
                            deltaData,
                            ranges,
                            deltas);
                    break;
                }
            }

            if (posSame)
            {
                deltas.Write(0, 1);
            }
            else
            {
                deltas.Write(1, 1);

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
                    case PosVector3Packer::BitVector3Sorted:
                    case PosVector3Packer::Sorted_no_bit_count:
                    {
                        // GLEE, actually writing a delta this time!
                        auto vec = IntVec3
                        {
                            target[i].position_x - base[i].position_x,
                            target[i].position_y - base[i].position_y,
                            target[i].position_z - base[i].position_z,
                        };

                        unsigned maxMagnitude = 1 + static_cast<unsigned>(MaxPositionChangePerSnapshot * frameDelta);
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

                            case PosVector3Packer::BitVector3Sorted:
                            {
                                bitsWritten = BitVector3SortedEncode(
                                    vec,
                                    maxMagnitude,
                                    Use_magnitude_as::Vector_Magnitude,
                                    deltas);

                                break;
                            }

                            case PosVector3Packer::Sorted_no_bit_count:
                            {
                                bitsWritten = Sorted_no_bit_count_Encode(
                                    vec,
                                    maxMagnitude,
                                    Use_magnitude_as::Vector_Magnitude,
                                    deltas);

                                break;
                            }
                        }

                        assert(bitsWritten < ((MaxBitsXY * 2) + MaxBitsZ));
                        break;
                    }
                }
            }

            deltas.Write(target[i].interacting, 1);
        }
    }

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

            case ChangedArrayEncoding::RangeSimple:
            {
                return RangeEncodeSimpleEncode(changed);
            }

            case ChangedArrayEncoding::RangeSmarter:
            {
                return RangeEncodeSmarterEncode(changed);
            }

            case ChangedArrayEncoding::RangeSimpleAdaptive:
            {
                return RangeEncodeSimpleAdaptiveEncode(changed);
            }

            case ChangedArrayEncoding::RangeSmarterAdaptive:
            {
                return RangeEncodeSmarterAdaptiveEncode(changed);
            }
        }
    }();

    assert(changedCompressed.Bits() <= Cubes);

    result.Write(changedCompressed);
    result.Write(deltas);

    result.TrimZerosFromBack();

    return result.Data();
}

Frame Decode(
    const Frame& base,
    std::vector<uint8_t>& buffer,
    unsigned frameDelta,
    Config config =
    {
        ChangedArrayEncoding::None,
        PosVector3Packer::None,
        QuatPacker::None,
    })
{
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

            case ChangedArrayEncoding::RangeSimple:
            {
                return RangeEncodeSimpleDecode(bits, Cubes);
            }

            case ChangedArrayEncoding::RangeSmarter:
            {
                return RangeEncodeSmarterDecode(bits, Cubes);
            }

            case ChangedArrayEncoding::RangeSimpleAdaptive:
            {
                return RangeEncodeSimpleAdaptiveDecode(bits, Cubes);
            }

            case ChangedArrayEncoding::RangeSmarterAdaptive:
            {
                return RangeEncodeSmarterAdaptiveDecode(bits, Cubes);
            }
        }
    }();

    changed.Reset();

    for (size_t i = 0; i < Cubes; ++i)
    {
        if (changed.Read(1))
        {
            switch (config.quatPacker)
            {
                default:
                case QuatPacker::None:
                {
                    auto changed = bits.Read(1);

                    if (changed)
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

                    break;
                }

                case QuatPacker::BitVector3Unrelated:
                case QuatPacker::BitVector3BitCount:
                case QuatPacker::BitVector3Sorted:
                case QuatPacker::Sorted_no_bit_count:
                {
                    result[i].orientation_largest =
                            bits.Read(RotationIndexMaxBits);

                    auto notCompressed = bits.Read(1);

                    if (!notCompressed)
                    {
                        IntVec3 vec;

                        switch (config.quatPacker)
                        {
                            default:
                            case QuatPacker::BitVector3BitCount:
                            {
                                vec = BitVector3BitCountDecode(
                                    (1 << (QuanternionUncompressedBitThreshold - 2)) - 1,
                                    bits);

                                break;
                            }

                            case QuatPacker::BitVector3Unrelated:
                            {
                                vec = BitVector3UnrelatedDecode(
                                    (1 << RotationMaxBits) - 1,
                                    bits);

                                break;
                            }

                            case QuatPacker::BitVector3Sorted:
                            {
                                vec = BitVector3SortedDecode(
                                    (1 << RotationMaxBits) - 1,
                                    Use_magnitude_as::Constant,
                                    bits);

                                break;
                            }

                            case QuatPacker::Sorted_no_bit_count:
                            {
                                vec = Sorted_no_bit_count_Decode(
                                    (1 << RotationMaxBits) - 1,
                                    Use_magnitude_as::Constant,
                                    bits);

                                break;
                            }
                        }

                        result[i].orientation_a = vec.x + base[i].orientation_a;
                        result[i].orientation_b = vec.y + base[i].orientation_b;
                        result[i].orientation_c = vec.z + base[i].orientation_c;
                    }
                    else
                    {
                        if (config.quatPacker == QuatPacker::Sorted_no_bit_count)
                        {
                            auto g = Sorted_no_bit_count_Decode(
                                (1 << RotationMaxBits) - 1,
                                Use_magnitude_as::Gaffer_Encode,
                                bits);

                            result[i].orientation_a = g.x + 256;
                            result[i].orientation_b = g.y + 256;
                            result[i].orientation_c = g.z + 256;
                        }
                        else
                        {

                            result[i].orientation_a = bits.Read(RotationMaxBits);
                            result[i].orientation_b = bits.Read(RotationMaxBits);
                            result[i].orientation_c = bits.Read(RotationMaxBits);
                        }
                    }

                    break;
                }

                case QuatPacker::BitVector3ModifiedZigZag:
                {
                    result[i].orientation_largest =
                            bits.Read(RotationIndexMaxBits);

                    auto notCompressed = bits.Read(1);

                    if (!notCompressed)
                    {
                        IntVec3 vec;

                        if (config.quatPacker == QuatPacker::BitVector3Unrelated)
                        {
                            vec = BitVector3UnrelatedDecode(
                                (1 << RotationMaxBits) - 1,
                                bits);
                        }
                        else
                        {
                            IntVec3 baseVec =
                            {
                                base[i].orientation_a,
                                base[i].orientation_b,
                                base[i].orientation_c,
                            };

                            vec = BitVector3BitCountZigZagDecode(
                                baseVec,
                                (1 << RotationMaxBits) - 1,
                                bits);
                        }

                        result[i].orientation_a = vec.x;
                        result[i].orientation_b = vec.y;
                        result[i].orientation_c = vec.z;
                    }
                    else
                    {
                        result[i].orientation_a = bits.Read(RotationMaxBits);
                        result[i].orientation_b = bits.Read(RotationMaxBits);
                        result[i].orientation_c = bits.Read(RotationMaxBits);
                    }

                    break;
                }

                case QuatPacker::BitVector3ModifiedGaffer:
                {
                    result[i].orientation_largest =
                        bits.Read(RotationIndexMaxBits);

                    IntVec3 baseVec =
                    {
                        base[i].orientation_a,
                        base[i].orientation_b,
                        base[i].orientation_c,
                    };

                    auto vec = BitVector3ModifiedGafferDecode(
                            baseVec,
                            (1 << RotationMaxBits) - 1,
                            RangeBits{5, 6, 7},
                            bits);

                    result[i].orientation_a = vec.x;
                    result[i].orientation_b = vec.y;
                    result[i].orientation_c = vec.z;

                    break;
                }

                case QuatPacker::Gaffer:
                {
                    auto full = bits.Read(1);

                    if (full)
                    {
                        result[i].orientation_largest   = bits.Read(RotationIndexMaxBits);
                        result[i].orientation_a         = bits.Read(RotationMaxBits);
                        result[i].orientation_b         = bits.Read(RotationMaxBits);
                        result[i].orientation_c         = bits.Read(RotationMaxBits);
                        break;
                    }

                    result[i].orientation_largest =
                        base[i].orientation_largest;

                    auto vec = GafferDecode({4, 5, 7}, bits);

                    result[i].orientation_a = base[i].orientation_a + vec.x;
                    result[i].orientation_b = base[i].orientation_b + vec.y;
                    result[i].orientation_c = base[i].orientation_c + vec.z;

                    break;
                }
            }

            auto posChanged = bits.Read(1);

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
                    case PosVector3Packer::BitVector3Sorted:
                    case PosVector3Packer::Sorted_no_bit_count:
                    {
                        unsigned maxMagnitude =
                            1 +
                            static_cast<unsigned>(
                                MaxPositionChangePerSnapshot *
                                frameDelta);

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

                            case PosVector3Packer::BitVector3Sorted:
                            {
                                vec = BitVector3SortedDecode(
                                    maxMagnitude,
                                    Use_magnitude_as::Vector_Magnitude,
                                    bits);

                                break;
                            }

                            case PosVector3Packer::Sorted_no_bit_count:
                            {
                                vec = Sorted_no_bit_count_Decode(
                                    maxMagnitude,
                                    Use_magnitude_as::Vector_Magnitude,
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

void Tests()
{
    Gaffer_tests();
    Quat_tests();
    TruncateTest();
    BitVector3Tests();
    Max_gaffer_tests();
    Model_tests();
    Range_tests();
    RunLengthTests();
    ZigZagTest();
    GaffersRangeTest();
    BitStreamTest();
}

// //////////////////////////////////////////////////////

#define PRINT_INT(x) printf("%-32s\t%d\n", #x, x);
#define PRINT_FLOAT(x) printf("%-32s\t%f\n", #x, x);

void CalculateStats(std::vector<Frame>& frames, const Config& config)
{
    // Lets do some research
    auto size = frames.size();
    Stats stats
    {
        {1024*1024,0,0,0},
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
        0,
        {static_cast<unsigned>(1 << RotationMaxBits),0,0,0},
        {static_cast<unsigned>(1 << RotationMaxBits),0,0,0},
        {static_cast<unsigned>(1 << RotationMaxBits),0,0,0},
        {static_cast<unsigned>(1 << RotationMaxBits),0,0,0},
        0,
        {static_cast<unsigned>(MaxPositionChangePerSnapshot * 20),0,0,0},
        0,
        {static_cast<unsigned>(MaxPositionChangePerSnapshot * 20),0,0,0},
        0,
        0,
        0,
        0,
        0,
        0,

        0,
        0,
        0,
        0,
        0,
        0,

        {1000000,0,0,0},
        {1000000,0,0,0},

        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    };

    // Lets actually do the stuff.
    unsigned bytes = 0;
    unsigned packetsCoded = 0;

    for (size_t i = PacketDelta; i < size; ++i)
    {
        auto buffer = EncodeStats(
            frames[i-PacketDelta],
            frames[i],
            PacketDelta,
            stats,
            config);

        bytes += buffer.size();
        stats.bytesPerPacket.Update(buffer.size());

        packetsCoded++;
    }

    printf("== Statistics ================================\n\n");

    float rle                       = Average(stats.rle);
    float bitpack                   = Average(stats.bitpack);
    float bitexprle                 = Average(stats.bitexprle);
    float expbitpack                = Average(stats.bitbitpack);
    float bitbitpackfull            = Average(stats.bitbitfullpack);
    float range_simple              = Average(stats.range_simple);
    float range_smarter             = Average(stats.range_smarter);
    float range_simple_adaptive     = Average(stats.range_simple_adaptive);
    float range_smarter_adaptive    = Average(stats.range_smarter_adaptive);

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

    PRINT_FLOAT(range_simple)
    PRINT_INT(stats.range_simple.min)
    PRINT_INT(stats.range_simple.max)

    PRINT_FLOAT(range_smarter)
    PRINT_INT(stats.range_smarter.min)
    PRINT_INT(stats.range_smarter.max)

    PRINT_FLOAT(range_simple_adaptive)
    PRINT_INT(stats.range_simple_adaptive.min)
    PRINT_INT(stats.range_simple_adaptive.max)

    PRINT_FLOAT(range_smarter_adaptive)
    PRINT_INT(stats.range_smarter_adaptive.min)
    PRINT_INT(stats.range_smarter_adaptive.max)

//    for (const auto h : stats.changed0DistanceRunHistogram)
//    {
//        printf("%d, ", h);
//    }
//    printf("\n");

//    for (const auto h : stats.changed1DistanceRunHistogram)
//    {
//        printf("%d, ", h);
//    }

    printf("\n");
    printf("\n");

    auto dxa =  Average(stats.deltaX);
    auto dya =  Average(stats.deltaY);
    auto dza =  Average(stats.deltaZ);
    auto dta =  Average(stats.deltaTotal);
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

    auto posPackedAvg = Average(stats.posDeltaPackedBitCount);
    PRINT_FLOAT(posPackedAvg)
    PRINT_INT(stats.posDeltaPackedBitCount.min)
    PRINT_INT(stats.posDeltaPackedBitCount.max)

    printf("\n");
    PRINT_INT(stats.quatChanged)
    PRINT_INT(stats.largestDifferent)

    auto daa = Average(stats.deltaA);
    auto dba = Average(stats.deltaB);
    auto dca = Average(stats.deltaC);
    auto dalla = Average(stats.deltaABC);
    PRINT_FLOAT(daa)
    PRINT_INT(stats.deltaA.min)
    PRINT_INT(stats.deltaA.max)
    PRINT_FLOAT(dba)
    PRINT_INT(stats.deltaB.min)
    PRINT_INT(stats.deltaB.max)
    PRINT_FLOAT(dca)
    PRINT_INT(stats.deltaC.min)
    PRINT_INT(stats.deltaC.max)
    PRINT_FLOAT(dalla)
    PRINT_INT(stats.deltaABC.min)
    PRINT_INT(stats.deltaABC.max)

    printf("\n");
    auto quatPackedAvg = Average(stats.quatDeltaPackedBitCount);
    PRINT_FLOAT(quatPackedAvg)
    PRINT_INT(stats.quatDeltaPackedBitCount.min)
    PRINT_INT(stats.quatDeltaPackedBitCount.max)

    printf("\n");
    auto bytesAvg = Average(stats.bytesPerPacket);
    PRINT_FLOAT(bytesAvg)
    PRINT_INT(stats.bytesPerPacket.min)
    PRINT_INT(stats.bytesPerPacket.max)

    printf("\n");

    PRINT_INT(stats.bitStream0)
    PRINT_INT(stats.bitStream1)
    PRINT_INT(stats.bitStream00)
    PRINT_INT(stats.bitStream01)
    PRINT_INT(stats.bitStream10)
    PRINT_INT(stats.bitStream11)

    printf("\n");

    PRINT_INT(stats.ones);
    PRINT_INT(stats.zeros);
    PRINT_INT(stats.zero_zero);
    PRINT_INT(stats.one_one);
    PRINT_INT(stats.zero_one);
    PRINT_INT(stats.one_zero);

    printf("\n");

    auto rotor_average = Average(stats.rotor);
    PRINT_FLOAT(rotor_average);
    PRINT_INT(stats.rotor.min);
    PRINT_INT(stats.rotor.max);

    printf("\n");

    PRINT_INT(stats.rotor_bits_8);
    PRINT_INT(stats.rotor_bits_9);
    PRINT_INT(stats.rotor_bits_10);
    PRINT_INT(stats.rotor_bits_11);
    PRINT_INT(stats.rotor_bits_12);
    PRINT_INT(stats.rotor_bits_13);
    PRINT_INT(stats.rotor_bits_14);
    PRINT_INT(stats.rotor_bits_15);
    PRINT_INT(stats.rotor_bits_wtf);

    printf("\n");

    auto rotor_bits_average = Average(stats.rotor_bits);
    PRINT_FLOAT(rotor_bits_average);
    PRINT_INT(stats.rotor_bits.min);
    PRINT_INT(stats.rotor_bits.max);

    printf("\n");


    printf("\n==============================================\n");

//    auto histoSize = stats.quantCommonHistogram.size();
//    for (unsigned i = 0; i < histoSize; ++i)
//    {
//        printf("%d, %d\n", stats.quantCommonHistogram[i], i);
//    }

//    for (unsigned i = 0; i < histoSize; ++i)
//    {
//        unsigned mask = (1 << quantHistogramBitsPerComponent) - 1;
//        unsigned shift =
//                RotationMaxBits -
//                quantHistogramBitsPerComponent;
//        auto value = stats.quantCommonHistogramTooBig[i];
//        unsigned a = (i & mask) << shift;
//        unsigned b = ((i >> quantHistogramBitsPerComponent) & mask) << shift;
//        unsigned c = ((i >> (2 * quantHistogramBitsPerComponent)) & mask) << shift;
//        printf(
//            "%d, %d, %d, %d, %d\n",
//            value,
//            i,
//            a,
//            b,
//            c);
//    }

}

// //////////////////////////////////////////////////////

void Compress(std::vector<Frame>& frames, const Config& config)
{
    auto packets = frames.size();
    unsigned bytes = 0;
    unsigned packetsCoded = 0;

    for (size_t i = PacketDelta; i < packets; ++i)
    {
        auto buffer = Encode(
            frames[i-PacketDelta],
            frames[i],
            PacketDelta,
            config);

        bytes += buffer.size();

        auto back = Decode(
            frames[i-PacketDelta],
            buffer,
            PacketDelta,
            config);

        assert(back == frames[i]);

        packetsCoded++;
    }

    float packetSizeAverge = ((float) bytes) / packetsCoded;
    float bytesPerSecondAverage = packetSizeAverge * 60.0f;
    float kbps = bytesPerSecondAverage * 8 / 1000.0f;

    printf("\n");
    printf("== Compression ===============================\n\n");

    PRINT_INT(packetsCoded)
    PRINT_INT(bytes)
    PRINT_FLOAT(bytesPerSecondAverage)
    PRINT_FLOAT(packetSizeAverge)
    PRINT_FLOAT(kbps)

    printf("\n==============================================\n");
}

void Range_compress(std::vector<Frame>& frames)
{
    auto packets = frames.size();
    unsigned bytes = 0;
    unsigned packetsCoded = 0;

    for (size_t i = PacketDelta; i < packets; ++i)
    {
        auto buffer = Encode_frames(
            frames[i-PacketDelta],
            frames[i],
            PacketDelta);

        bytes += buffer.size();

        auto back = Decode_frames(
            frames[i-PacketDelta],
            buffer,
            PacketDelta);

        assert(back == frames[i]);

        packetsCoded++;
    }

    float packetSizeAverge = ((float) bytes) / packetsCoded;
    float bytesPerSecondAverage = packetSizeAverge * 60.0f;
    float kbps = bytesPerSecondAverage * 8 / 1000.0f;

    printf("\n");
    printf("== Compression (model) =======================\n\n");

    PRINT_INT(packetsCoded)
    PRINT_INT(bytes)
    PRINT_FLOAT(bytesPerSecondAverage)
    PRINT_FLOAT(packetSizeAverge)
    PRINT_FLOAT(kbps)

    printf("\n==============================================\n");
}

template<class MODEL>
struct Test
{
    int id;
    MODEL model;
};

void Range_search(std::vector<Frame>& frames)
{
    // RAM: TODO: Make encode/decode virtual so I can group all models into
    // one array. Print out id of model so I can see which one is best.
    using namespace Range_models;

    // Generate our test itesm
    std::vector<Test<Binary>> binarys;
    for (auto i = 1 ; i < 8; ++i)
    {
        binarys.push_back({i,Binary(i)});
    }

    std::vector<Test<Binary_two_speed>> binary_two_speeds;
    for (auto i = 1 ; i < 8; ++i)
    {
        for (auto j = i ; j < 8; ++j)
        {
            binary_two_speeds.push_back({i*10 + j, Binary_two_speed(i,j)});
        }
    }

    std::vector<Test<Binary_history<Binary, Binary>>> histories;
    for (auto a : binarys)
    {
        for (auto b : binarys)
        {
            // first been 0 or 1 only matters to about 0.0001 bits.
            histories.push_back({100 + a.id * 10 + b.id, Binary_history<Binary, Binary>(0,a.model,b.model)});
        }
    }

    std::vector<Test<Range_models::Binary_history<Binary_two_speed, Binary_two_speed>>>
        histories_two_speed;

    for (auto a : binary_two_speeds)
    {
        for (auto b : binary_two_speeds)
        {
            histories_two_speed.push_back({10000 + 100 * a.id + b.id,
                Range_models::Binary_history<Binary_two_speed, Binary_two_speed>(0,a.model,b.model)});
        }
    }

    // lets go wtf for a moment.
    std::vector<Test<Binary_history<Binary_history<Binary,Binary>, Binary_history<Binary,Binary>>>>
        histories_two_history_binary;

    for (auto a : binarys)
    {
        for (auto b : binarys)
        {
            for (auto a2 : binarys)
            {
                for (auto b2 : binarys)
                {
                    // first been 0 or 1 only matters to about 0.0001 bits.
                    histories_two_history_binary.push_back(
                    {
                        (100 * (a.id * 10 + b.id)) + a2.id * 10 + b2.id,
                        Binary_history<Binary_history<Binary,Binary>, Binary_history<Binary,Binary>>(
                            0,
                            {0,a.model,b.model},
                            {1,a2.model,b2.model})
                    });
                }
            }
        }
    }

    auto Run = [&](auto test) -> MinMaxSum
    {
        auto packets = frames.size();
        MinMaxSum result = {10000,0,0,0};

        for (size_t p = PacketDelta; p < packets; ++p)
        {
            auto& base = frames[p-PacketDelta];
            auto& target = frames[p];
            Range_types::Bytes data;
            data.reserve(70);

            {
                Range_coders::Encoder           range(data);
                Range_coders::Binary_encoder    binary(range);

                auto size = base.size();
                for (unsigned i = 0; i < size; ++i)
                {
                    auto quant_index_changed =
                        base[i].orientation_largest != target[i].orientation_largest;

                    auto quant_changed =
                        (quant_index_changed) ||
                        (base[i].orientation_a != target[i].orientation_a) ||
                        (base[i].orientation_b != target[i].orientation_b) ||
                        (base[i].orientation_c != target[i].orientation_c);

                    test.Encode(binary, quant_changed);
                }
            }

            result.Update(data.size());
        }

        return result;
    };

    auto Loop = [&](const auto& group, const char* name)
    {
        printf("%s\n", name);

        float best = 100000;
        int best_id = -1;

        for (auto model : group)
        {
            auto r = Run(model.model);

            auto a = Average(r);

            if (a < best)
            {
                best = a;
                best_id = model.id;
            }

            printf("%d: %d - %d - %f\n", model.id, r.min, r.max, Average(r));
        }

        printf("\nBest: %d -> %f\n", best_id, best);
        printf("\n");
    };

    Loop(binarys, "Binarys");
    Loop(binary_two_speeds, "binary_two_speeds");
    Loop(histories, "histories");
    Loop(histories_two_speed, "histories_two_speed");
    Loop(histories_two_history_binary, "histories_two_history_binary");

    printf("\n==============================================\n");
}


int main(int, char**)
{
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

    // //////////////////////////////////////////////////////

    if (doTests)
    {
        Tests();
    }

    Config config
    {
        ChangedArrayEncoding::Exp,
        PosVector3Packer::Sorted_no_bit_count,
        QuatPacker::BitVector3Unrelated,
    };

    if (doStats)
    {
        CalculateStats(frames, config);
    }

    if (doCompression)
    {
        Compress(frames, config);
    }

    if (doRangeCompression)
    {
        Range_compress(frames);
    }

    if (doRangeSearch)
    {
        Range_search(frames);
    }

    return 0;
}
