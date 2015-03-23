// Copyright 2015 Richard Maxwell, all rights reserved.

// http://gafferongames.com/2015/03/14/the-networked-physics-data-compression-challenge/

// g++ -std=c++14 DeltaCompress.cpp -o DeltaCompress

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

static const size_t Cubes = 901;
static const size_t CubeBits = 10;

static_assert(Cubes < ((1 << CubeBits) - 1), "CubeBits is too small.");

// Using contraints from
// http://gafferongames.com/networked-physics/snapshot-compression/

// With this technique I�ve found that minimum sufficient precision for
// my simulation is 9 bits per-smallest component. This gives a result
// of 2 + 9 + 9 + 9 = 29 bits per-orientation (originally 128!).

static const unsigned RotationMaxBits = 9;
static const unsigned RotationIndexMaxBits = 2;

// I found a maximum speed of 32 meters per-second is a nice power of two
// and doesn�t negatively affect the player experience in the cube simulation.
static const unsigned MaxSpeedMetersPerSecond = 32;

// Now we have only position to compress. We�ll use the same trick that we used
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

    for (auto i = 0; i < size; ++i)
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
        auto bitOffset = m_bitOffset;
        auto byteOffset = bitOffset / 8u;
        auto byteOffsetAfter = (bitOffset + bitsToWrite) / 8u;
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

    std::vector<uint8_t> ReadArray(unsigned bitsToRead)
    {
        std::vector<uint8_t> result;

        while (bitsToRead > 7)
        {
            result.push_back(Read(8));
            bitsToRead -= 8;
        }

        if (bitsToRead)
        {
            result.push_back(Read(bitsToRead));
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

    BitStream bitsOutToSplit(bitsIn.Data());
    BitStream bitsSplit(bitsOutToSplit.ReadArray(6));

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

std::vector<uint8_t> Encode(const Frame& base, const Frame& target)
{
    const auto count = base.size();

    assert(count > 0);
    assert(count == target.size());

    // ////////////////////////////////
    // A. If nothing changed, don't send anything at all.

    bool same = true;
    auto firstChanged = 0;
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

    for (auto i = 0; i < count; ++i)
    {
        if (base[i] == target[i])
        {
            changed.Write(0, 1);
        }
        else
        {
            changed.Write(1, 1);

            // Meh, first time delta, just write the target.
            deltas.Write(target[i].orientation_largest, RotationIndexMaxBits);
            deltas.Write(ZigZag(target[i].orientation_a), RotationMaxBits);
            deltas.Write(ZigZag(target[i].orientation_b), RotationMaxBits);
            deltas.Write(ZigZag(target[i].orientation_c), RotationMaxBits);
            deltas.Write(ZigZag(target[i].position_x), MaxBitsXY);
            deltas.Write(ZigZag(target[i].position_y), MaxBitsXY);
            deltas.Write(target[i].position_z, MaxBitsZ);
            deltas.Write(target[i].interacting, 1);
        }
    }

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

    auto changed = BitStream(bits.ReadArray(Cubes));

    for (auto i = 0; i < Cubes; ++i)
    {
        if (changed.Read(1))
        {
            result[i].orientation_largest   = bits.Read(RotationIndexMaxBits);
            result[i].orientation_a         = ZigZag(bits.Read(RotationMaxBits));
            result[i].orientation_b         = ZigZag(bits.Read(RotationMaxBits));
            result[i].orientation_c         = ZigZag(bits.Read(RotationMaxBits));
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
}

// //////////////////////////////////////////////////////

int main(int, char**)
{
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
    Stats stats = {0};
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
        auto buffer = Encode(frames[i-FirstBase], frames[i]);

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

    return 0;
}
