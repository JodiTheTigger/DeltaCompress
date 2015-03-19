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

// Still making up my mind on the best way to assert.
#define ASSERT_TRUE(x) assert(x)

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

static const size_t Cubes = 901;

typedef std::array<DeltaData, Cubes> Frame;


// //////////////////////////////////////////////////////

int main(int, char**)
{
    auto frames = []() -> std::vector<Frame>
    {
        const auto fileHandle = fopen("delta_data.bin", "rb");
        if (!fileHandle)
        {
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
        printf("ERROR: Cannot open 'delta_data.bin' for reading.\n");
        return 1;
    }

    return 0;
}
