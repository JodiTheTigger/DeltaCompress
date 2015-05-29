// Copyright 2015 Richard Maxwell, all rights reserved.

// http://gafferongames.com/2015/03/14/the-networked-physics-data-compression-challenge/

// g++ -std=c++14 DeltaCompress.cpp -Wall -Wextra -Werror -g -o DeltaCompress

// //////////////////////////////////////////////////////

#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

#include "Range_coding.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <vector>
#include <cmath>
#include <type_traits>

// Note sure I need these.
//#include <stdio.h>
//#include <cstdlib>
//#include <functional>

// //////////////////////////////////////////////////////

bool do_tests       = true;
bool do_compression = true;

// //////////////////////////////////////////////////////

#ifdef WIN32
#define msvc_constexpr
#else
#define msvc_constexpr constexpr
#endif

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

inline constexpr bool operator==(const DeltaData& lhs, const DeltaData& rhs)
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

inline constexpr bool operator!=(const DeltaData& lhs, const DeltaData& rhs)
{
    return !operator==(lhs,rhs);
}

// //////////////////////////////////////////////////////

inline constexpr bool quat_equal(const DeltaData& lhs, const DeltaData& rhs)
{
    return
    (
            (lhs.orientation_a          == rhs.orientation_a)
        &&  (lhs.orientation_b          == rhs.orientation_b)
        &&  (lhs.orientation_c          == rhs.orientation_c)
        &&  (lhs.orientation_largest    == rhs.orientation_largest)
    );
}

inline constexpr bool pos_equal(const DeltaData& lhs, const DeltaData& rhs)
{
    return
    (
            (lhs.position_x             == rhs.position_x)
        &&  (lhs.position_y             == rhs.position_y)
        &&  (lhs.position_z             == rhs.position_z)
    );
}

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
//static const unsigned RotationIndexMaxBits = 2;

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

inline bool operator!=(const Frame& lhs, const Frame& rhs)
{
    return !operator==(lhs,rhs);
}

// //////////////////////////////////////////////////////

using Vec3i = std::array<int, 3>;
using Vec4f = std::array<float, 4>;
using Vec3f = std::array<float, 4>;

struct Gaffer
{
    unsigned orientation_largest;
    Vec3i vec;
};

struct Maxwell
{
    unsigned multiplier_index;
    Vec3i vec;
};

struct Quat
{
    Vec4f vec;

    float constexpr operator[](unsigned index) const
    {
        return vec[index];
    }

    float& operator[](unsigned index)
    {
        return vec[index];
    }
};

float* begin(Quat& q)
{
    return &(q.vec[0]);
}

float* end(Quat& q)
{
    return &(q.vec[4]);
}

struct Rotor
{
    Vec3f vec;

    float constexpr operator[](unsigned index) const
    {
        return vec[index];
    }

    float& operator[](unsigned index)
    {
        return vec[index];
    }
};

// //////////////////////////////////////////////////////

static_assert
(
    std::is_pod<Gaffer>::value,
    "Not a POD."
);

static_assert
(
    std::is_pod<Maxwell>::value,
    "Not a POD."
);

static_assert
(
    std::is_pod<Quat>::value,
    "Not a POD."
);

static_assert
(
    std::is_pod<Rotor>::value,
    "Not a POD."
);

// //////////////////////////////////////////////////////

inline constexpr bool operator==(const Gaffer& lhs, const Gaffer& rhs)
{
    return
    (
        (lhs.orientation_largest == rhs.orientation_largest) &&
        (lhs.vec == rhs.vec)
    );
}

inline constexpr bool operator!=(const Gaffer& lhs, const Gaffer& rhs)
{
    return !operator==(lhs,rhs);
}

// //////////////////////////////////////////////////////

auto constexpr mul(const Quat& lhs, float rhs) -> Quat
{
    return
    {
        lhs[0]*rhs,
        lhs[1]*rhs,
        lhs[2]*rhs,
        lhs[3]*rhs
    };
}

auto constexpr mul(const Quat& lhs, const Quat& rhs) -> Quat
{
    return
    {
        (lhs[0]*rhs[0] - lhs[1]*rhs[1] - lhs[2]*rhs[2] - lhs[3]*rhs[3]),
        (lhs[0]*rhs[1] + lhs[1]*rhs[0] + lhs[2]*rhs[3] - lhs[3]*rhs[2]),
        (lhs[0]*rhs[2] - lhs[1]*rhs[3] + lhs[2]*rhs[0] + lhs[3]*rhs[1]),
        (lhs[0]*rhs[3] + lhs[1]*rhs[2] - lhs[2]*rhs[1] + lhs[3]*rhs[0])
    };
}


auto constexpr magnitude_squared(const Quat& quat) -> float
{
    return
        quat[0] * quat[0] +
        quat[1] * quat[1] +
        quat[2] * quat[2] +
        quat[3] * quat[3];
}

auto msvc_constexpr normalise(const Quat& q) -> Quat
{
    return mul(q, 1.0f / std::sqrt(magnitude_squared(q)));
}

auto constexpr conjugate(const Quat& q) -> Quat
{
    return
    {
        q[0],
        -q[1],
        -q[2],
        -q[3]
    };
}

auto constexpr magnitude_squared(const Rotor& r) -> float
{
    return
        r[0] * r[0] +
        r[1] * r[1] +
        r[2] * r[2];
}

auto constexpr mul(const Rotor& lhs, float rhs) -> Rotor
{
    return
    {
        lhs[0]*rhs,
        lhs[1]*rhs,
        lhs[2]*rhs,
    };
}

// //////////////////////////////////////////////////////

Rotor constexpr to_rotor(const Quat& q)
{
    return
    {
        q[1] / (1.0f + q[0]),
        q[2] / (1.0f + q[0]),
        q[3] / (1.0f + q[0]),
    };
}

Quat constexpr to_quat(const Rotor& r)
{
    return
    {
        (1.0f - magnitude_squared(r)) / (1 + magnitude_squared(r)),
        (r[0] * 2.0f) / (1 + magnitude_squared(r)),
        (r[1] * 2.0f) / (1 + magnitude_squared(r)),
        (r[2] * 2.0f) / (1 + magnitude_squared(r)),
    };
}

// //////////////////////////////////////////////////////
// Converting between quantised quats and back again is really pissy.
// //////////////////////////////////////////////////////

static const float G_256    = 256.4995127f;
static const float Q_TO_G   = G_256 * M_SQRT2;
static const float G_TO_Q   = 1.0f / (G_256 * M_SQRT2);

// Did some research, -16, -16, -256 seems to be the worse.
// Use that as max sum.
static const float GAFFER_ONE_SQUARED =
        (256 * 256) +
        (256 * 256) +
        (16 * 16) +
        (16 * 16);

Quat to_quat(const Gaffer& gaffer)
{
    Quat result
    {
        static_cast<float>(gaffer.vec[0] - 256),
        static_cast<float>(gaffer.vec[1] - 256),
        static_cast<float>(gaffer.vec[2] - 256),
        0.0,
    };

    auto largest_squared = Quat
    {
        result[0] * result[0],
        result[1] * result[1],
        result[2] * result[2],
        0,
    };

    auto largest_value_squared =
        GAFFER_ONE_SQUARED -
        (
            largest_squared[0] +
            largest_squared[1] +
            largest_squared[2]
        );

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
        // I need to add an offset to make sure it is still
        // the largest as opposed to equal.
        // CBFd to figure out what number to actually use.
        largest_value_squared = next_largest + 20;
    }

    auto largest = sqrt(largest_value_squared);

    if (gaffer.orientation_largest == 0)
    {
        result[3] = result[2];
        result[2] = result[1];
        result[1] = result[0];
    }

    if (gaffer.orientation_largest == 1)
    {
        result[3] = result[2];
        result[2] = result[1];
    }

    if (gaffer.orientation_largest == 2)
    {
        result[3] = result[2];
    }

    result[gaffer.orientation_largest] = largest;
    result = mul(result, G_TO_Q);

    {
        auto mag = magnitude_squared(result);
        auto quantised = 256 * (mag - 1.0f);
        assert(std::abs(quantised) < 0.5f);
    }

    result = normalise(result);
    return result;
}

Gaffer to_gaffer(const Quat& quat)
{
    const auto size = quat.vec.size();
    unsigned largest_index = 0;
    float largest = 0;

    auto squared = Quat
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
        fixed_quat = mul(quat, -1.0f);
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
    // then re-adjust
    auto multiple = 1.0f;
    static const auto MAX_POSITIVE = (255.0f / 256.0f) * M_SQRT1_2;
    static const auto MAX_NEGATIVE = -1.0f * M_SQRT1_2;

    if (gaffer[0] > MAX_POSITIVE)
    {
        multiple = MAX_POSITIVE / gaffer[0];
    }
    if (gaffer[1] > MAX_POSITIVE)
    {
        multiple = MAX_POSITIVE / gaffer[1];
    }
    if (gaffer[2] > MAX_POSITIVE)
    {
        multiple = MAX_POSITIVE / gaffer[2];
    }


    if (gaffer[0] < MAX_NEGATIVE)
    {
        multiple = MAX_NEGATIVE / gaffer[0];
    }
    if (gaffer[1] < MAX_NEGATIVE)
    {
        multiple = MAX_NEGATIVE / gaffer[1];
    }
    if (gaffer[2] < MAX_NEGATIVE)
    {
        multiple = MAX_NEGATIVE / gaffer[2];
    }

    auto result = Gaffer
    {
        largest_index,
        {
            256 + static_cast<int>(round(gaffer[0] * Q_TO_G * multiple)),
            256 + static_cast<int>(round(gaffer[1] * Q_TO_G * multiple)),
            256 + static_cast<int>(round(gaffer[2] * Q_TO_G * multiple)),
        }
    };

    assert(result.vec[0] >= 0);
    assert(result.vec[1] >= 0);
    assert(result.vec[2] >= 0);

    assert(result.vec[0] < 512);
    assert(result.vec[1] < 512);
    assert(result.vec[2] < 512);

    return result;
}

Gaffer to_gaffer(const DeltaData& delta)
{
    return Gaffer
    {
        static_cast<unsigned>(delta.orientation_largest),
        {
            delta.orientation_a,
            delta.orientation_b,
            delta.orientation_c,
        }
    };
}

// //////////////////////////////////////////////////////

using Multipliers = std::array<float, 8>;

// For each rotor I bascially tried all multipliers starting from 1
// until the coverted back value matched the original Gaffer value.
// To convert the search space to just eight items I took the quartiles
// and some other values and played around Until I got an max value encoded
// being 12 bits, and the average below 5 bits.
// Pretty sure this would only work well when the frame delta is 6.
static const Multipliers DEFAULT_MUTLIPLES =
{
    130.0f,
    222.0f,
    329.0f,
    471.0f,
    512.0f,
    969.0f,
    4247.0f,
    65535.0f,
};

//static const Multipliers DEFAULT_MUTLIPLES =
//{
//    471.0f,
//    512.0f,
//    969.0f,
//    65535.0f,
//};

//static const float MULTIPLIER_MAX_VALUE_ADJUST = 0.6f;

// The for love of god I couldn't derive the formula for generating the
// multiplier. So I have to just search for it instead.
auto to_maxwell(
    const Rotor& to_encode,
    const Quat& base_quat,
    const Gaffer& target_gaffer,
    const Multipliers& multipliers) -> Maxwell
{
    unsigned index = 0;
    Rotor coded;

    for (auto rotor_multiply : multipliers)
    {
        coded = Rotor
        {
            std::round(to_encode[0] * rotor_multiply),
            std::round(to_encode[1] * rotor_multiply),
            std::round(to_encode[2] * rotor_multiply)
        };

        auto decoded            = mul(coded, 1.0f / rotor_multiply);
        auto decoded_quat       = to_quat(decoded);
        auto result_target_quat = mul(decoded_quat, base_quat);
        auto result             = to_gaffer(result_target_quat);

        if  (result == target_gaffer)
        {
            return
            {
                index,
                {
                    static_cast<int>(coded[0]),
                    static_cast<int>(coded[1]),
                    static_cast<int>(coded[2]),
                }
            };
        }

        ++index;
    }

    // :-(
    assert(false);

    return
    {
        index - 1,
        {
            static_cast<int>(coded[0]),
            static_cast<int>(coded[1]),
            static_cast<int>(coded[2]),
        }
    };
}

auto to_maxwell
(
    const Gaffer& base,
    const Gaffer& target,
    const Multipliers& multipliers = DEFAULT_MUTLIPLES
)
-> Maxwell
{
    auto base_quat      = to_quat(base);
    auto target_quat    = to_quat(target);

    // Hmm, seem we get stupid magnitudess due
    // to rotating near itself (from q to -q roughly).
    auto target_quat_neg = mul(target_quat, -1.0f);

    // http://www.geomerics.com/blogs/quaternions-rotations-and-compression/
    auto r                  = mul(target_quat, conjugate(base_quat));
    auto r_neg              = mul(target_quat_neg, conjugate(base_quat));
    auto rotor              = to_rotor(r);
    auto rotor_neg          = to_rotor(r_neg);
    auto mag_squared        = magnitude_squared(rotor);
    auto mag_squared_neg    = magnitude_squared(rotor_neg);
    auto mag                = sqrt(mag_squared);
    auto mag_neg            = sqrt(mag_squared_neg);

    auto result_rotor =
        (mag < mag_neg) ?
            rotor :
            rotor_neg;

    return to_maxwell
    (
        result_rotor,
        base_quat,
        target,
        multipliers
    );
}

auto to_gaffer
(
    const Gaffer& base,
    const Maxwell& coded,
    const Multipliers& multipliers = DEFAULT_MUTLIPLES
)
-> Gaffer
{
    assert(coded.multiplier_index < multipliers.size());

    auto base_quat = to_quat(base);
    auto inv_multiplier = 1.0f / multipliers[coded.multiplier_index];
    auto rotor = Rotor
    {
        coded.vec[0] * inv_multiplier,
        coded.vec[1] * inv_multiplier,
        coded.vec[2] * inv_multiplier,
    };

    auto r = to_quat(rotor);
    auto target_quat = mul(r, base_quat);

    return to_gaffer(target_quat);
}

// //////////////////////////////////////////////////////

inline constexpr uint32_t Zig_zag(int32_t n)
{
    return (n << 1) ^ (n >> 31);
}

inline constexpr int32_t Zig_zag(uint32_t n)
{
    return (n >> 1) ^ (-(static_cast<int>(n) & 1));
}

void zig_zag_test()
{
    assert(42 == Zig_zag(Zig_zag(42)));
    assert(-42 == Zig_zag(Zig_zag(-42)));
    assert(0 == Zig_zag(Zig_zag(0)));
    assert(-12345 == Zig_zag(Zig_zag(-12345)));
    assert(30654 == Zig_zag(Zig_zag(30654)));
    assert(-31654 == Zig_zag(Zig_zag(-31654)));
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

using namespace Range_models;

namespace Sorted_position
{
    struct Model
    {
        static const unsigned position_0_update = 16;
        static const unsigned position_1_update = 8;
        static const unsigned position_2_update = 5;
        // Quat changed:
        // Model adaption should be based on previous + 30 previous.
        // bit 0 = any of previous 5 set to 1
        // bit 1 = any of last 33 to 28 set to 1
        std::array<Binary, 4> quat_changed;

        // Position changed:
        // Model on if quat changed or not.
        std::array<Binary, 2> position_changed;

        // Rotor:
        // Encode signs seperatly (3 bit range).
        // ?: different models based on previous signs
        // ?: code sorted or not
        // ?: If bits > 256, only code top 8 bits
        // ?: Can we use max rotational velocity
        Periodic_update rotor_multiplier_lookup;;
        Periodic_update rotor_signs;

        // Um, this will be sloooow.
        Periodic_update_with_kernel rotor_magnitudes_low;
        Periodic_update rotor_magnitudes_high;

        // Position:
        // ?: Corrlation between max axis and our max axis?
        // Encode signs seperatly
        // Use max velocity
        // ?: different models based on previous signs
        // ?: code sorted or not
        // ?: If bits > 256, only code top 8 bits
        Periodic_update position_signs;

        // Since this is sorted, the max values are basically
        // [0]: 1 + MaxPositionChangePerSnapshot * PacketDelta
        // [1]: 1 + [0] / 2
        // [2]: 1 + [0] / 3
        Periodic_update_with_kernel position_0;
        Periodic_update_with_kernel position_1;
        Periodic_update_with_kernel position_2;

        // If all three items are differnt use largest_index and next_largest
        // index. Otherwise two items match so use different_index only (rare).
        Periodic_update  largest_index;
        Periodic_update  different_index;
        std::array<Binary, 3>               next_largest_index;

        // Interactive:
        // one bit
        // model on:
        // Previous was interactive
        // Anything has changed
        std::array<Binary, 2> interactive;
    };

    auto encode
    (
        const Frame& base,
        const Frame& target,
        unsigned frameDelta
    )
    -> Range_types::Bytes
    {
        auto                size = base.size();
        Range_types::Bytes  data;

        const unsigned max_position_0 =
            1 +
            static_cast<unsigned>
            (
                MaxPositionChangePerSnapshot *
                frameDelta
            );

        const unsigned max_position_1 = 1 + (max_position_0 >> 1);
        const unsigned max_position_2 = 1 + (max_position_0 / 3);

        Model model =
        {
            {},
            {},
            {8, 16},
            {8, 16},
            {256, 16*3},
            {32, 16*3},


            {8, 16},
            {max_position_0, Model::position_0_update},
            {max_position_1, Model::position_1_update},
            {max_position_2, Model::position_2_update},
            {3, 8},
            {3, 1},

            {},
            {},
        };

        {
            Range_coders::Encoder           range(data);
            Range_coders::Binary_encoder    binary(range);

            for (unsigned i = 0; i < size; ++i)
            {
                // //////////////////////////////////////////////////////

                unsigned quat_lookup = 0;

                auto any_5_changed =
                    [&base, &target](unsigned i, unsigned history)
                    -> bool
                {
                    if (i < history)
                    {
                        return false;
                    }

                    auto max = (i - history) + 5;
                    for (unsigned j = max - 5; j < max; ++j)
                    {
                        if (!quat_equal(base[j], target[j]))
                        {
                            return true;
                            break;
                        }
                    }

                    return false;
                };

                if (any_5_changed(i, 33))
                {
                    quat_lookup = 1;
                }

                if (any_5_changed(i, 5))
                {
                    quat_lookup |= 2;
                }

                auto quat_changed = !quat_equal(base[i], target[i]);

                model.quat_changed[quat_lookup].Encode(binary, quat_changed);

                auto get_signs = [](const Vec3i& v) -> unsigned
                {
                    unsigned result = 0;

                    if (v[0] < 0)
                    {
                        result |= 1;
                    }

                    if (v[1] < 0)
                    {
                        result |= 2;
                    }

                    if (v[2] < 0)
                    {
                        result |= 4;
                    }

                    return result;
                };

                auto strip_signs = [](const Vec3i& v) -> Vec3i
                {
                    return
                    {
                        std::abs(v[0]),
                        std::abs(v[1]),
                        std::abs(v[2]),
                    };
                };

                if (quat_changed)
                {
                    auto b = to_gaffer(base[i]);
                    auto t = to_gaffer(target[i]);
                    auto m = to_maxwell(b, t);
                    auto signs = get_signs(m.vec);
                    auto vec = strip_signs(m.vec);

                    model.rotor_multiplier_lookup.Encode
                    (
                        range,
                        m.multiplier_index
                    );

                    model.rotor_signs.Encode
                    (
                        range,
                        signs
                    );

                    for (auto v: vec)
                    {
                        model.rotor_magnitudes_high.Encode(range, v >> 8);
                        model.rotor_magnitudes_low.Encode(range, v & 0xFF);
                    }
                }

                // //////////////////////////////////////////////////////

                auto pos_changed = !pos_equal(base[i], target[i]);

                unsigned pos_lookup = quat_changed ? 1 : 0;

                model.position_changed[pos_lookup].Encode(binary, pos_changed);

                if (pos_changed)
                {
                    Vec3i delta
                    {
                        target[i].position_x - base[i].position_x,
                        target[i].position_y - base[i].position_y,
                        target[i].position_z - base[i].position_z
                    };

                    auto signs = get_signs(delta);
                    auto vec = strip_signs(delta);

                    model.position_signs.Encode(range, signs);

                    // Sort
                    int odd = -1;
                    {
                        if ((vec[0] != vec[1]) && (vec[1] == vec[2]))
                        {
                            odd = 0u;
                        }
                        if ((vec[0] != vec[1]) && (vec[0] == vec[2]))
                        {
                            odd = 1u;
                        }
                        if ((vec[0] == vec[1]) && (vec[1] != vec[2]))
                        {
                            odd = 2u;
                        }
                    }

                    unsigned top = 0;
                    unsigned next = 0;
                    {
                        using std::swap;

                        if  (
                                (vec[1] > vec[0]) &&
                                (vec[1] >= vec[2])
                            )
                        {
                            swap(vec[0], vec[1]);
                            top = 1;
                        }
                        else
                        {
                            if  (
                                    (vec[2] > vec[0]) &&
                                    (vec[2] >= vec[1])
                                )
                            {
                                swap(vec[0], vec[2]);
                                swap(vec[1], vec[2]);
                                top = 2;
                            }
                        }

                        assert(vec[0] >= vec[1]);
                        assert(vec[0] >= vec[2]);

                        if (vec[2] > vec[1])
                        {
                            swap(vec[1], vec[2]);
                            next = 1;
                        }

                        assert(vec[1] >= vec[2]);
                    }

                    // RAM: TODO: Use max magnitude to determine
                    // truncation values. Also, truncate
                    assert(vec[0] <= static_cast<int>(max_position_0));
                    assert(vec[1] <= static_cast<int>(max_position_1));
                    assert(vec[2] <= static_cast<int>(max_position_2));

                    model.position_0.Encode(range, vec[0]);

                    if (vec[0] > 0)
                    {
                        const auto x = static_cast<unsigned>(vec[0]);
                        const auto trunc_x = std::min(x, max_position_1 - 1);

                        model.position_1.Encode(range, vec[1], trunc_x);

                        if (vec[1] > 0)
                        {
                            const auto y = static_cast<unsigned>(vec[1]);
                            const auto trunc_y = std::min(y, max_position_2 - 1);

                            model.position_2.Encode(range, vec[2], trunc_y);
                        }
                    }

                    bool all_same = (vec[0] == vec[1]) && (vec[1] == vec[2]);

                    if (!all_same)
                    {
                        if (odd < 0)
                        {
                            model.largest_index.Encode(range, top);
                            model.next_largest_index[top].Encode(binary, next);
                        }
                        else
                        {
                            model.different_index.Encode(range, odd);
                        }
                    }
                }

                // //////////////////////////////////////////////////////

                unsigned interactive_lookup = 0;

                if ((pos_changed) || (quat_changed))
                {
                    interactive_lookup = 1;
                }

                model.interactive[interactive_lookup].Encode
                (
                    binary,
                    target[i].interacting
                );
            }
        }

        return data;
    }

    auto decode
    (
        const Frame& base,
        const Range_types::Bytes& data,
        unsigned frameDelta
    )
    -> Frame
    {
        auto    size = base.size();
        Frame   target;

        const unsigned max_position_0 =
            1 +
            static_cast<unsigned>
            (
                MaxPositionChangePerSnapshot *
                frameDelta
            );

        const unsigned max_position_1 = 1 + (max_position_0 >> 1);
        const unsigned max_position_2 = 1 + (max_position_0 / 3);

        Model model =
        {
            {},
            {},
            {8, 16},
            {8, 16},
            {256, 16*3},
            {32, 16*3},


            {8, 16},
            {max_position_0, Model::position_0_update},
            {max_position_1, Model::position_1_update},
            {max_position_2, Model::position_2_update},
            {3, 8},
            {3, 1},

            {},
            {},
        };

        {
            Range_coders::Decoder           range(data);
            Range_coders::Binary_decoder    binary(range);

            for (unsigned i = 0; i < size; ++i)
            {
                // //////////////////////////////////////////////////////

                unsigned quat_lookup = 0;

                auto any_5_changed =
                    [&base, &target](unsigned i, unsigned history)
                    -> bool
                {
                    if (i < history)
                    {
                        return false;
                    }

                    auto max = (i - history) + 5;
                    for (unsigned j = max - 5; j < max; ++j)
                    {
                        if (!quat_equal(base[j], target[j]))
                        {
                            return true;
                            break;
                        }
                    }

                    return false;
                };

                auto add_signs = [](unsigned signs, const Vec3i& v) -> Vec3i
                {
                    return
                    {
                        (signs & 1) ? -v[0] : v[0],
                        (signs & 2) ? -v[1] : v[1],
                        (signs & 4) ? -v[2] : v[2],
                    };
                };

                if (any_5_changed(i, 33))
                {
                    quat_lookup = 1;
                }

                if (any_5_changed(i, 5))
                {
                    quat_lookup |= 2;
                }

                auto quat_changed =
                    model.quat_changed[quat_lookup].Decode(binary);

                if (quat_changed)
                {
                    auto index = model.rotor_multiplier_lookup.Decode(range);
                    auto signs = model.rotor_signs.Decode(range);
                    auto v = Vec3i
                    {
                        static_cast<int>
                        (
                            model.rotor_magnitudes_high.Decode(range) << 8 |
                            model.rotor_magnitudes_low.Decode(range)
                        ),
                        static_cast<int>
                        (
                            model.rotor_magnitudes_high.Decode(range) << 8 |
                            model.rotor_magnitudes_low.Decode(range)
                        ),
                        static_cast<int>
                        (
                            model.rotor_magnitudes_high.Decode(range) << 8 |
                            model.rotor_magnitudes_low.Decode(range)
                        ),
                    };

                    auto m = Maxwell
                    {
                        index,
                        add_signs(signs, v)
                    };

                    auto b = to_gaffer(base[i]);
                    auto t = to_gaffer(b, m);

                    target[i].orientation_largest = t.orientation_largest;
                    target[i].orientation_a       = t.vec[0];
                    target[i].orientation_b       = t.vec[1];
                    target[i].orientation_c       = t.vec[2];
                }
                else
                {
                    target[i].orientation_largest = base[i].orientation_largest;
                    target[i].orientation_a       = base[i].orientation_a;
                    target[i].orientation_b       = base[i].orientation_b;
                    target[i].orientation_c       = base[i].orientation_c;
                }

                // //////////////////////////////////////////////////////

                auto pos_changed =
                    model.position_changed[quat_changed].Decode(binary);

                if (pos_changed)
                {
                    auto signs = model.position_signs.Decode(range);

                    auto x = model.position_0.Decode(range);

                    const auto trunc_x = std::min(x, max_position_1 - 1);

                    auto y = x ? model.position_1.Decode(range, trunc_x) : 0;

                    const auto trunc_y = std::min(y, max_position_2 - 1);

                    auto z = y ? model.position_2.Decode(range, trunc_y) : 0;

                    auto v = Vec3i
                    {
                        static_cast<int>(x),
                        static_cast<int>(y),
                        static_cast<int>(z),
                    };

                    // Read the Order
                    auto largest = 0;
                    auto nextLargest = 0;

                    if ((v[0] != v[1]) || (v[1] != v[2]))
                    {
                        if ((v[0] != v[1]) && (v[1] != v[2]))
                        {
                            largest =
                                model.largest_index.Decode(range);

                            nextLargest =
                                model.next_largest_index[largest].Decode
                                (
                                    binary
                                );
                        }
                        else
                        {
                            bool all_same = (v[0] == v[1]) && (v[1] == v[2]);

                            if (!all_same)
                            {
                                auto odd_one =
                                    model.different_index.Decode(range);

                                if (v[0] != v[1])
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

                    {
                        using std::swap;

                        if (nextLargest)
                        {
                            swap(v[1], v[2]);
                        }

                        if (largest)
                        {
                            if (largest == 1)
                            {
                                swap(v[0], v[1]);
                            }
                            else
                            {
                                swap(v[0], v[1]);
                                swap(v[1], v[2]);
                            }
                        }
                    };

                    auto vec = add_signs(signs, v);

                    target[i].position_x = vec[0] + base[i].position_x;
                    target[i].position_y = vec[1] + base[i].position_y;
                    target[i].position_z = vec[2] + base[i].position_z;

                }
                else
                {
                    target[i].position_x = base[i].position_x;
                    target[i].position_y = base[i].position_y;
                    target[i].position_z = base[i].position_z;
                }

                // //////////////////////////////////////////////////////

                target[i].interacting =
                    model.interactive[(quat_changed | pos_changed)].Decode
                    (
                        binary
                    );
            }
        }

        return target;
    }
};

// //////////////////////////////////////////////////////

namespace Naieve_rotor
{
    struct Model
    {
        // Quat changed:
        // Model adaption should be based on previous + 30 previous.
        // bit 0 = any of previous 5 set to 1
        // bit 1 = any of last 33 to 28 set to 1
        std::array<Binary, 4> quat_changed;

        // Position changed:
        // Model on if quat changed or not.
        std::array<Binary, 2> position_changed;

        // Rotor:
        // Encode signs seperatly (3 bit range).
        // ?: different models based on previous signs
        // ?: code sorted or not
        // ?: If bits > 256, only code top 8 bits
        // ?: Can we use max rotational velocity
        Periodic_update rotor_multiplier_lookup  = {8, 16};
        Periodic_update rotor_signs              = {8, 16};

        // Um, this will be sloooow.
        Periodic_update rotor_magnitudes_low     = {256, 16*3};
        Periodic_update rotor_magnitudes_high    = {32, 16*3};

        // Position:
        // ?: Corrlation between max axis and our max axis?
        // Encode signs seperatly
        // Use max velocity
        // ?: different models based on previous signs
        // ?: code sorted or not
        // ?: If bits > 256, only code top 8 bits
        Periodic_update position_signs = {8, 16};

        // Um, this will be sloooow.
        Periodic_update position_magnitudes_low  = {256, 16*3};
        Periodic_update position_magnitudes_high = {256, 16*3};

        // Interactive:
        // one bit
        // model on:
        // Previous was interactive
        // Anything has changed
        std::array<Binary, 2> interactive;
    };

    auto encode
    (
        const Frame& base,
        const Frame& target,
        unsigned// frameDelta
    )
    -> Range_types::Bytes
    {
        auto                size = base.size();
        Model         model;
        Range_types::Bytes  data;

        {
            Range_coders::Encoder           range(data);
            Range_coders::Binary_encoder    binary(range);

            for (unsigned i = 0; i < size; ++i)
            {
                // //////////////////////////////////////////////////////

                unsigned quat_lookup = 0;

                auto any_5_changed =
                    [&base, &target](unsigned i, unsigned history)
                    -> bool
                {
                    if (i < history)
                    {
                        return false;
                    }

                    auto max = (i - history) + 5;
                    for (unsigned j = max - 5; j < max; ++j)
                    {
                        if (!quat_equal(base[j], target[j]))
                        {
                            return true;
                            break;
                        }
                    }

                    return false;
                };

                if (any_5_changed(i, 33))
                {
                    quat_lookup = 1;
                }

                if (any_5_changed(i, 5))
                {
                    quat_lookup |= 2;
                }

                auto quat_changed = !quat_equal(base[i], target[i]);

                model.quat_changed[quat_lookup].Encode(binary, quat_changed);

                auto get_signs = [](const Vec3i& v) -> unsigned
                {
                    unsigned result = 0;

                    if (v[0] < 0)
                    {
                        result |= 1;
                    }

                    if (v[1] < 0)
                    {
                        result |= 2;
                    }

                    if (v[2] < 0)
                    {
                        result |= 4;
                    }

                    return result;
                };

                auto strip_signs = [](const Vec3i& v) -> Vec3i
                {
                    return
                    {
                        std::abs(v[0]),
                        std::abs(v[1]),
                        std::abs(v[2]),
                    };
                };

                if (quat_changed)
                {
                    auto b = to_gaffer(base[i]);
                    auto t = to_gaffer(target[i]);
                    auto m = to_maxwell(b, t);
                    auto signs = get_signs(m.vec);
                    auto vec = strip_signs(m.vec);

                    model.rotor_multiplier_lookup.Encode
                    (
                        range,
                        m.multiplier_index
                    );

                    model.rotor_signs.Encode
                    (
                        range,
                        signs
                    );

                    for (auto v: vec)
                    {
                        model.rotor_magnitudes_high.Encode(range, v >> 8);
                        model.rotor_magnitudes_low.Encode(range, v & 0xFF);
                    }
                }

                // //////////////////////////////////////////////////////

                auto pos_changed = !pos_equal(base[i], target[i]);

                unsigned pos_lookup = quat_changed ? 1 : 0;

                model.position_changed[pos_lookup].Encode(binary, pos_changed);

                if (pos_changed)
                {
                    Vec3i delta
                    {
                        target[i].position_x - base[i].position_x,
                        target[i].position_y - base[i].position_y,
                        target[i].position_z - base[i].position_z
                    };

                    auto signs = get_signs(delta);
                    auto vec = strip_signs(delta);

                    model.position_signs.Encode(range, signs);

                    for (auto v: vec)
                    {
                        model.position_magnitudes_high.Encode(range, v >> 8);
                        model.position_magnitudes_low.Encode(range, v & 0xFF);
                    }
                }

                // //////////////////////////////////////////////////////

                unsigned interactive_lookup = 0;

                if ((pos_changed) || (quat_changed))
                {
                    interactive_lookup = 1;
                }

                model.interactive[interactive_lookup].Encode
                (
                    binary,
                    target[i].interacting
                );
            }
        }

        return data;
    }

    auto decode
    (
        const Frame& base,
        const Range_types::Bytes& data,
        unsigned// frameDelta
    )
    -> Frame
    {
        auto                size = base.size();
        Model         model;
        Frame               target;

        {
            Range_coders::Decoder           range(data);
            Range_coders::Binary_decoder    binary(range);

            for (unsigned i = 0; i < size; ++i)
            {
                // //////////////////////////////////////////////////////

                unsigned quat_lookup = 0;

                auto any_5_changed =
                    [&base, &target](unsigned i, unsigned history)
                    -> bool
                {
                    if (i < history)
                    {
                        return false;
                    }

                    auto max = (i - history) + 5;
                    for (unsigned j = max - 5; j < max; ++j)
                    {
                        if (!quat_equal(base[j], target[j]))
                        {
                            return true;
                            break;
                        }
                    }

                    return false;
                };

                auto add_signs = [](unsigned signs, const Vec3i& v) -> Vec3i
                {
                    return
                    {
                        (signs & 1) ? -v[0] : v[0],
                        (signs & 2) ? -v[1] : v[1],
                        (signs & 4) ? -v[2] : v[2],
                    };
                };

                if (any_5_changed(i, 33))
                {
                    quat_lookup = 1;
                }

                if (any_5_changed(i, 5))
                {
                    quat_lookup |= 2;
                }

                auto quat_changed =
                    model.quat_changed[quat_lookup].Decode(binary);

                if (quat_changed)
                {
                    auto index = model.rotor_multiplier_lookup.Decode(range);
                    auto signs = model.rotor_signs.Decode(range);
                    auto v = Vec3i
                    {
                        static_cast<int>
                        (
                            model.rotor_magnitudes_high.Decode(range) << 8 |
                            model.rotor_magnitudes_low.Decode(range)
                        ),
                        static_cast<int>
                        (
                            model.rotor_magnitudes_high.Decode(range) << 8 |
                            model.rotor_magnitudes_low.Decode(range)
                        ),
                        static_cast<int>
                        (
                            model.rotor_magnitudes_high.Decode(range) << 8 |
                            model.rotor_magnitudes_low.Decode(range)
                        ),
                    };

                    auto m = Maxwell
                    {
                        index,
                        add_signs(signs, v)
                    };

                    auto b = to_gaffer(base[i]);
                    auto t = to_gaffer(b, m);

                    target[i].orientation_largest = t.orientation_largest;
                    target[i].orientation_a       = t.vec[0];
                    target[i].orientation_b       = t.vec[1];
                    target[i].orientation_c       = t.vec[2];
                }
                else
                {
                    target[i].orientation_largest = base[i].orientation_largest;
                    target[i].orientation_a       = base[i].orientation_a;
                    target[i].orientation_b       = base[i].orientation_b;
                    target[i].orientation_c       = base[i].orientation_c;
                }

                // //////////////////////////////////////////////////////

                auto pos_changed =
                    model.position_changed[quat_changed].Decode(binary);

                if (pos_changed)
                {
                    auto signs = model.position_signs.Decode(range);
                    auto v = Vec3i
                    {
                        static_cast<int>
                        (
                            model.position_magnitudes_high.Decode(range) << 8 |
                            model.position_magnitudes_low.Decode(range)
                        ),
                        static_cast<int>
                        (
                            model.position_magnitudes_high.Decode(range) << 8 |
                            model.position_magnitudes_low.Decode(range)
                        ),
                        static_cast<int>
                        (
                            model.position_magnitudes_high.Decode(range) << 8 |
                            model.position_magnitudes_low.Decode(range)
                        ),
                    };

                    auto vec = add_signs(signs, v);

                    target[i].position_x = vec[0] + base[i].position_x;
                    target[i].position_y = vec[1] + base[i].position_y;
                    target[i].position_z = vec[2] + base[i].position_z;

                }
                else
                {
                    target[i].position_x = base[i].position_x;
                    target[i].position_y = base[i].position_y;
                    target[i].position_z = base[i].position_z;
                }

                // //////////////////////////////////////////////////////

                target[i].interacting =
                    model.interactive[(quat_changed | pos_changed)].Decode
                    (
                        binary
                    );
            }
        }

        return target;
    }
};

// //////////////////////////////////////////////////////

namespace Naieve_gaffer
{
    using namespace Range_models;

    struct Everything_model
    {
        // Note: I == investigate if dependency exists.
        // code quat changed as adaptive history binary
        // {
        //   assume no history on largest quat changed, simple binary adaptive
        //   3 models dependent on largest changed (I)
        //   {
        //     dependent on largest changed to find largest vector index (I)
        //     dependent on largest, to get next largest. (I)
        //     3 models per index
        //     {
        //       model on bit count for first item (per max magnitude value)
        //       tree model on bits left, but with force 0 to dictate how many
        //       bits are actually used. To avoid coding them at all.
        //       don't need to code last bits needed as we know already.
        //     }
        //   }
        //
        // code position changed dependent on quat changed. Hmm, rygorous does
        // 8 models, 2x dependent on quat changed, then 2x per quat different
        // component. weird. (I)
        // Code position just like quat, but with max magnitude specifying bits
        // as well! yay! force zero.
        // interacting based on current interacting + quat/position changed.
        using Simple = Binary;

        std::array<Simple, 2> quat_changed =
        {{
            {4, 31},
            {},
        }};

        // Test range vs tree
        // test if needed multiple models depending on previous model.
        Simple largest_index_quant_changed;
        Periodic_update largest_index_quant = {4, 8};

        struct Vector_model
        {
            Periodic_update largest_index = {3, 8};
            std::array<Simple, 3> next_largest_index;

            // Note: even though I hard code 12, it shouldn't be hard coded.
            // bits1 = MinBits(RotationMaxBits)
            // bits2 = MinBits(MaxPositionChangePerSnapshot * MaxFrameDelta)
            // bits = max(bits1, bits2)

            std::array<Periodic_update, 2> bits_for_value =
            {{
                Periodic_update{12, 16},
                Periodic_update{12, 16},
            }};

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

    // //////////////////////////////////////////////////////

    auto Largest_next_magnitude(
            unsigned magnitude,
            unsigned zig_zag_axis) -> unsigned
    {
        // max magnitude for Zig_zag before we overflow
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
            Vec3i vec,
            unsigned max_magnitude,
            Everything_model::Vector_model& model,
            Encoder& range,
            Binary_encoder& binary)
    {
        // +1 for the sign bit.
        unsigned max_bits_required = 1 + MinBits(max_magnitude);

        // //////////////////////////////////////////

        assert(abs(vec[0]) <= static_cast<int>(max_magnitude));
        assert(abs(vec[1]) <= static_cast<int>(max_magnitude));
        assert(abs(vec[2]) <= static_cast<int>(max_magnitude));

        // Sort from largest to smallest
        auto zx = Zig_zag(vec[0]);
        auto zy = Zig_zag(vec[1]);
        auto zz = Zig_zag(vec[2]);

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

    Vec3i Vector3Decode(
            unsigned max_magnitude,
            Everything_model::Vector_model& model,
            Decoder& range,
            Binary_decoder& binary)
    {
        Vec3i result = {0,0,0};

        // +1 for the sign bit.
        unsigned max_bits_required = 1 + MinBits(max_magnitude);

        // //////////////////////////////////////////

        // Read the Order
        auto top = model.largest_index.Decode(range);
        auto next = model.next_largest_index[top].Decode(binary);

        auto ReturnSorted = [&top, &next](Vec3i vec) -> Vec3i
        {
            using std::swap;

            if (next)
            {
                swap(vec[1], vec[2]);
            }

            if (top)
            {
                if (top == 1)
                {
                    swap(vec[0], vec[1]);
                }
                else
                {
                    swap(vec[0], vec[1]);
                    swap(vec[1], vec[2]);
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
            target = Zig_zag(zig_zag);

            if (model.Reduce_vector_using_magnitude)
            {
                max_magnitude = Largest_next_magnitude(max_magnitude, zig_zag);
            }

            max_bits_required = 1 + MinBits(max_magnitude);
            max_bits_required = std::min(max_bits_required, bits);

            return true;
        };

        // //////////////////////////////////////////

        if (!Code(result[0], 0))
        {
            return ReturnSorted(result);
        }

        if (!Code(result[1], 1))
        {
            return ReturnSorted(result);
        }

        // //////////////////////////////////////////

        auto zz = model.value[2].Decode(binary, max_bits_required);
        result[2] = Zig_zag(zz);

        return ReturnSorted(result);
    }

    auto encode(
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

                auto last_quat_changed =
                    (i != 0) ?
                        (base[i - 1].orientation_largest != target[i - 1].orientation_largest) ||
                        (base[i - 1].orientation_a != target[i - 1].orientation_a) ||
                        (base[i - 1].orientation_b != target[i - 1].orientation_b) ||
                        (base[i - 1].orientation_c != target[i - 1].orientation_c) ?
                            1 :
                            0
                        :
                        0;

                model.quat_changed[last_quat_changed].Encode
                (
                    binary,
                    quant_changed
                );

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

                        Vec3i not_delta
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
                        Vec3i delta
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
                    Vec3i delta
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

    auto decode(
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
            auto last_quat_changed =
                (i != 0) ?
                    (base[i - 1].orientation_largest != target[i - 1].orientation_largest) ||
                    (base[i - 1].orientation_a != target[i - 1].orientation_a) ||
                    (base[i - 1].orientation_b != target[i - 1].orientation_b) ||
                    (base[i - 1].orientation_c != target[i - 1].orientation_c) ?
                        1 :
                        0
                    :
                    0;

            auto quant_changed =
                model.quat_changed[last_quat_changed].Decode(binary);

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

                    target[i].orientation_a = not_delta[0];
                    target[i].orientation_b = not_delta[1];
                    target[i].orientation_c = not_delta[2];
                }

                if (!quant_index_changed)
                {
                    Vec3i delta = Vector3Decode(
                        (1u << RotationMaxBits) - 1,
                        model.quant_delta,
                        range,
                        binary);

                    target[i].orientation_largest =
                            base[i].orientation_largest;
                    target[i].orientation_a = base[i].orientation_a + delta[0];
                    target[i].orientation_b = base[i].orientation_b + delta[1];
                    target[i].orientation_c = base[i].orientation_c + delta[2];
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

                target[i].position_x = base[i].position_x + delta[0];
                target[i].position_y = base[i].position_y + delta[1];
                target[i].position_z = base[i].position_z + delta[2];
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
}

// //////////////////////////////////////////////////////

#define PRINT_INT(x) printf("%-32s\t%d\n", #x, x);
#define PRINT_FLOAT(x) printf("%-32s\t%f\n", #x, x);

void range_compress(std::vector<Frame>& frames)
{
    auto packets = frames.size();

    auto test = [&](auto encoder, auto decoder)
    {
        unsigned bytes = 0;
        unsigned packetsCoded = 0;
        unsigned min = 10000000;
        unsigned max = 0;

        for (size_t i = PacketDelta; i < packets; ++i)
        {
            auto buffer = encoder(
                frames[i-PacketDelta],
                frames[i],
                PacketDelta);

            const unsigned size = buffer.size();
            bytes += size;
            min = std::min(min, size);
            max = std::max(max, size);

            auto back = decoder(
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

        PRINT_INT(packetsCoded);
        PRINT_INT(bytes);
        PRINT_FLOAT(bytesPerSecondAverage);
        PRINT_FLOAT(packetSizeAverge);
        PRINT_INT(min);
        PRINT_INT(max);
        PRINT_FLOAT(kbps);

        printf("\n==============================================\n");
    };

    printf("Sorted_position");
    test
    (
        Sorted_position::encode,
        Sorted_position::decode
    );

    printf("Naieve_rotor");
    test
    (
        Naieve_rotor::encode,
        Naieve_rotor::decode
    );

    printf("Naieve_gaffer");
    test
    (
        Naieve_gaffer::encode,
        Naieve_gaffer::decode
    );
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

        // I really hate how I cannot make an array
        // without initilising it first.
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

    if (do_tests)
    {
        range_tests();
        zig_zag_test();
    }

    if (do_compression)
    {
        // To beat:
        // total packed size 1258584
        // Min packet size 2
        // Max packet size 910
        // 444.73 bytes/frame
        // 213.47 kbps

        range_compress(frames);
    }

    return 0;
}
