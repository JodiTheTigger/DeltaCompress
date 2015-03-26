// Copyright 2015 Richard Maxwell, all rights reserved.

// http://gafferongames.com/2015/03/14/the-networked-physics-data-compression-challenge/

// g++ -std=c++14 DeltaCompress.cpp -Wall -Wextra -Werror -g -o DeltaCompress

// //////////////////////////////////////////////////////

#include <stdint.h>
#include <stdio.h>
#include <assert.h>

#include <array>
#include <vector>

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
struct Stats
{
    unsigned itemCount;

    unsigned notChangedCount;

    unsigned interactingTotal;
    unsigned interactingNotChanged;
    unsigned notInteractingNotChanged;

    std::array<unsigned, Cubes> changed0DistanceRunHistogram;
    std::array<unsigned, Cubes> changed1DistanceRunHistogram;
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
    return (n >> 1) ^ (-(n & 1));
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
        return m_bitOffset + 1;
    }

    void Reset()
    {
        m_bitOffset = 0;
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

// http://michael.dipperstein.com/rle/

BitStream BitPack8Bit(BitStream toPack, unsigned runThreshold)
{
    assert(runThreshold >= 2);
    assert(runThreshold < 131);

    auto bits = toPack.Bits();
    toPack.Reset();
    BitStream result;
    BitStream hold;

    auto read = toPack.Read(1);
    auto last = read;
    auto running = false;
    unsigned run = 1;
    unsigned unRun = 1;
    hold.Write(read, 1);

    auto dumpRun = [&result](int value, unsigned last)
    {
        result.Write(ZigZag(value), 8);
        result.Write(last, 1);
    };

    auto dumpUnRun = [&result, &hold](int value)
    {
        if (!value)
        {
            return;
        }

        result.Write(ZigZag(value), 8);
        hold.Reset();

        for (int i = 0; i < value; ++i)
        {
            result.Write(hold.Read(1), 1);
        }
    };

    while (toPack.Bits() < bits)
    {
        auto read = toPack.Read(1);
        hold.Write(read, 1);

        if (running)
        {
            bool dump = true;

            if (last == read)
            {
                ++run;

                if (run < 130)
                {
                    dump = false;
                }
            }

            if (dump)
            {
                dumpRun(-(run - 2), last);

                hold.Reset();

                if (toPack.Bits() == bits)
                {
                    break;
                }

                read = toPack.Read(1);
                hold.Write(read, 1);
                running = false;
                run = 1;
            }
        }
        else
        {
            if (last == read)
            {
                ++run;
                ++unRun;

                if (run > runThreshold)
                {
                    dumpUnRun(unRun - (runThreshold + 1));

                    unRun = 1;
                    running = true;
                }
            }
            else
            {
                ++unRun;
                run = 1;
            }
        }
    }

    if (running)
    {
        dumpRun(-(run - 2), last);
    }
    else
    {
        if (hold.Bits() > 0)
        dumpUnRun(unRun - (runThreshold + 1));
    }

    return result;
}


BitStream BitPack8BitUnpack(BitStream toUnPack, unsigned totalExpectedBits)
{
    unsigned bitCount = 0;
    BitStream result;

    while (bitCount < totalExpectedBits)
    {
        auto run = ZigZag(toUnPack.Read(8));

        if (run < 0)
        {
            auto repeat = toUnPack.Read(1);

            auto count = -run + 2;

            for (auto i = 0; i < count; ++i)
            {
                result.Write(repeat, 1);
            }

            bitCount += count;
        }
        else
        {
            ++run;
            result.Write(toUnPack.ReadArray(run));

            bitCount += run;
        }
    }

    return result;
}

// //////////////////////////////////////////////////////

void BitPack8BitTest()
{
    BitStream start;

    start.Write(0xFF, 8);

    auto result = BitPack8Bit(start, 3);

    result.Reset();

    auto a = result.Read(8);
    auto az = ZigZag(a);
    auto b = result.Read(1);

    assert(-az == (8 - 2));
    assert(b == 1);

    result.Reset();

    auto back = BitPack8BitUnpack(result, 8);

    assert(start == back);
}

// //////////////////////////////////////////////////////

std::vector<uint8_t> Encode(
    const Frame& base,
    const Frame& target,
    Stats& stats)
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

            // Meh, first time delta, just write the target.
            deltas.Write(target[i].orientation_largest, RotationIndexMaxBits);
            deltas.Write(target[i].orientation_a, RotationMaxBits);
            deltas.Write(target[i].orientation_b, RotationMaxBits);
            deltas.Write(target[i].orientation_c, RotationMaxBits);
            deltas.Write(ZigZag(target[i].position_x), MaxBitsXY);
            deltas.Write(ZigZag(target[i].position_y), MaxBitsXY);
            deltas.Write(target[i].position_z, MaxBitsZ);
            deltas.Write(target[i].interacting, 1);
        }
    }

    ++(stats.changed0DistanceRunHistogram[runLength0]);
    ++(stats.changed1DistanceRunHistogram[runLength1]);

    result.Write(changed);
    result.Write(deltas);

    result.TrimZerosFromBack();

    return result.Data();
}

Frame Decode(const Frame& base, std::vector<uint8_t>& buffer)
{
    // A.
    if (buffer.empty())
    {
        return base;
    }

    BitStream bits(buffer);
    Frame result;

    auto changed = bits.ReadArray(Cubes);

    for (size_t i = 0; i < Cubes; ++i)
    {
        if (changed.Read(1))
        {
            result[i].orientation_largest   = bits.Read(RotationIndexMaxBits);
            result[i].orientation_a         = bits.Read(RotationMaxBits);
            result[i].orientation_b         = bits.Read(RotationMaxBits);
            result[i].orientation_c         = bits.Read(RotationMaxBits);
            result[i].position_x            = ZigZag(bits.Read(MaxBitsXY));
            result[i].position_y            = ZigZag(bits.Read(MaxBitsXY));
            result[i].position_z            = bits.Read(MaxBitsZ);
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

ByteVector RunLengthDecode(const ByteVector& data)
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

    return result;
}

void RunLengthTest()
{
    {
        auto data = ByteVector
        {
            0,1,3,3,4,4,4,4,4,5,6,6,6,5,4,3,3,3,3,4,
        };

        auto encoded = RunLengthEncode(data);
        auto decoded = RunLengthDecode(encoded);

        assert(data == decoded);
    }
    {
        auto data = ByteVector
        {
            0,1,3,3,4,4,4,4,4,5,6,6,6,5,4,3,3,3,3,
        };

        auto encoded = RunLengthEncode(data);
        auto decoded = RunLengthDecode(encoded);

        assert(data == decoded);
    }
    {
        auto data = ByteVector(300, 5);

        auto encoded = RunLengthEncode(data);
        auto decoded = RunLengthDecode(encoded);

        assert(data == decoded);
    }
}

// //////////////////////////////////////////////////////

// http://michael.dipperstein.com/rle/

ByteVector BitPackEncode(const ByteVector& data)
{
    auto size = data.size();

    if (size < 1)
    {
        return {};
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

            if ((startIndex + i) == size)
            {
                break;
            }

            ++i;
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

                if (data[index] != data[index - 1])
                {
                    break;
                }

                if ((startIndex + run) == size)
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
    }

    return result;
}

ByteVector BitPackDecode(const ByteVector& data)
{
    auto size = data.size();

    if (size < 2)
    {
        return {};
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
    }

    return result;
}

void BitPackTest()
{
    {
        auto data = ByteVector
        {
            0,1,3,3,4,4,4,4,4,5,6,6,6,5,4,3,3,3,3,4,
        };

        auto encoded = BitPackEncode(data);
        auto decoded = BitPackDecode(encoded);

        assert(data == decoded);
    }
    {
        auto data = ByteVector
        {
            0,1,3,3,4,4,4,4,4,5,6,6,6,5,4,3,3,3,3,
        };

        auto encoded = BitPackEncode(data);
        auto decoded = BitPackDecode(encoded);

        assert(data == decoded);
    }
    {
        auto data = ByteVector(300, 5);

        auto encoded = BitPackEncode(data);
        auto decoded = BitPackDecode(encoded);

        assert(data == decoded);
    }
}

// //////////////////////////////////////////////////////

static const std::array<unsigned, 8> exponentialBitLevelRunLengthLookup
{
    0, 2, 7, 27, 57, 90, 496, 898,
    //0, 2, 4, 12, 30, 62, 126, 254,
};

BitStream ExponentialBitLevelRunLengthEncode(BitStream data)
{
    auto size = data.Bits();

    if (size < 2)
    {
        return data;
    }

    unsigned count = 1;
    data.Reset();
    auto previous = data.Read(1);
    BitStream result;

    while (count < size)
    {
        auto current = data.Read(1);
        count++;

        unsigned run = 1;

        while (current == previous)
        {
            if (count == size)
            {
                break;
            }

            run++;
            count++;
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
            }

        }

        //result.Write(current, 1);
        previous = current;
    }

    return result;
}

BitStream ExponentialBitLevelRunLengthDecode(BitStream data)
{
    auto size = data.Bits();

    if (size < 2)
    {
        return data;
    }

    unsigned count = 1;
    data.Reset();
    auto previous = data.Read(1);
    BitStream result;

    while (count < size)
    {
        auto current = data.Read(1);
        count++;

        unsigned run = 1;

        while (current == previous)
        {
            if (count == size)
            {
                break;
            }

            run++;
            count++;
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

            if (count < size)
            {
                current = data.Read(1);
                count++;
            }
        }
        else
        {
            result.Write(previous, 1);
        }

        previous = current;
    }

    return result;
}

void ExponentialBitLevelRunLengthEncodingTest()
{
    {
        auto data = BitStream(ByteVector
        {
            0,1,3,3,0,0,0,0,0,5,6,6,6,5,4,3,3,3,3,4,
        }, 20 * 8);

        auto encoded = ExponentialBitLevelRunLengthEncode(data);
        auto decoded = ExponentialBitLevelRunLengthDecode(encoded);

        assert(data == decoded);
    }
    {
        auto data = BitStream(ByteVector
        {
            0,1,3,3,0,0,0,0,0,5,6,6,6,5,4,3,3,3,3,
        }, 19 * 8);

        auto encoded = ExponentialBitLevelRunLengthEncode(data);
        auto decoded = ExponentialBitLevelRunLengthDecode(encoded);

        assert(data == decoded);
    }
    {
        auto data = BitStream(ByteVector(100, 0), 100 * 8);

        auto encoded = ExponentialBitLevelRunLengthEncode(data);
        auto decoded = ExponentialBitLevelRunLengthDecode(encoded);

        assert(data == decoded);
    }
}

// //////////////////////////////////////////////////////

int main(int, char**)
{
    ExponentialBitLevelRunLengthEncodingTest();
    RunLengthTest();
    BitPackTest();
    BitPack8BitTest();
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

    for (size_t i = FirstBase; i < size; ++i)
    {
        auto buffer = Encode(frames[i-FirstBase], frames[i], stats);

        bytes += buffer.size();

        auto back = Decode(frames[i-FirstBase], buffer);

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

    return 0;
}
