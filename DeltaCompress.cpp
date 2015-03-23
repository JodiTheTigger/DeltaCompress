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
}

// //////////////////////////////////////////////////////

unsigned CrappyByteArrayWrite(
        std::vector<uint8_t>& target,
        unsigned bitOffset,
        unsigned value,
        unsigned bitsToWrite)
{
    auto byteOffset = bitOffset / 8u;
    auto byteOffsetAfter = (bitOffset + bitsToWrite) / 8u;
    auto index = bitOffset - (byteOffset * 8);
    auto mask = ((1 << index) - 1);

    while (target.size() <= byteOffsetAfter)
    {
        target.push_back(0);
    }

    while (bitsToWrite)
    {
        target[byteOffset] |= (value & 0xFF) << index;

        if ((index + bitsToWrite) > 8)
        {
            target[byteOffset + 1] |= (value >> (8 - index)) & mask;
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

    return bitOffset;
}

unsigned CrappyByteArrayRead(
        std::vector<uint8_t>& source,
        unsigned bitOffset,
        unsigned& value,
        unsigned bitsToRead)
{
    auto byteOffset = bitOffset / 8u;
    auto index = bitOffset - (byteOffset * 8);
    auto size = source.size();
    auto mask = ((1 << bitsToRead) - 1);
    auto shift = 0;

    value = 0;

    while (bitsToRead)
    {
        if (byteOffset < size)
        {
            value |= (source[byteOffset] >> index) << shift;
        }

        if (((index + bitsToRead) > 8) && ((byteOffset + 1) < size))
        {
            value |= ((source[byteOffset + 1] << (8 - index) & 0xFF)) << shift;
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

    return bitOffset;
}

void CrappyByteArrayTest()
{
    std::vector<uint8_t> array;

    unsigned bitOffset = 0;

    bitOffset = CrappyByteArrayWrite(array, bitOffset, 46, 6);
    bitOffset = CrappyByteArrayWrite(array, bitOffset, 666, 10);
    bitOffset = CrappyByteArrayWrite(array, bitOffset, 169, 8);

    unsigned a = 0;
    unsigned b = 0;
    unsigned c = 0;

    bitOffset = 0;

    bitOffset = CrappyByteArrayRead(array, bitOffset, a, 6);
    bitOffset = CrappyByteArrayRead(array, bitOffset, b, 10);
    bitOffset = CrappyByteArrayRead(array, bitOffset, c, 8);

    assert(a == 46);
    assert(b == 666);
    assert(c == 169);
}

// //////////////////////////////////////////////////////
/// \brief Encode
/// \param base
/// \param target
/// \param stats
/// \return
///
std::vector<uint8_t> Encode(const Frame& base, const Frame& target, Stats& stats)
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
}

Frame Decode(const Frame& base, std::vector<uint8_t>& buffer)
{
    // A.
    if (buffer.empty())
    {
        return base;
    }
}

// //////////////////////////////////////////////////////

int main(int, char**)
{
    CrappyByteArrayTest();
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

    return 0;
}
