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

bool do_tests       = false;
bool do_compression = true;
bool do_decompress  = true;

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
//static const unsigned MaxSpeedMetersPerSecond = 32;

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

static const unsigned   MaxBitsAll        =
    std::max(MaxBitsZ, std::max(MaxBitsXY, RotationMaxBits));

//static const unsigned   MaxSnapshotsPerSecond = 60;

// Currently gets stuck when the delta is larger than 26 packets.
static const unsigned   PacketDelta           = 6;

// This one is important
//static const float      MaxPositionChangePerSnapshot =
//        MaxSpeedMetersPerSecond * ValuesPerMeter /
//        (float) MaxSnapshotsPerSecond;

static_assert(
    ((1 << MaxBitsZ) - 1) == MaxZInclusize,
    "MaxBitZ doesn't match MaxZInclusize");

static_assert(
    ((1 << MaxBitsXY) - 1) == XYRange,
    "MaxBitsXY doesn't match XYRange");

// //////////////////////////////////////////////////////

typedef std::array<DeltaData, Cubes> Frame;

// //////////////////////////////////////////////////////

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

// //////////////////////////////////////////////////////

auto constexpr sub(const Vec3i& lhs, const Vec3i& rhs) -> Vec3i
{
    return
    {
        lhs[0] - rhs[0],
        lhs[1] - rhs[1],
        lhs[2] - rhs[2]
    };
};

auto constexpr add(const Vec3i& lhs, const Vec3i& rhs) -> Vec3i
{
    return
    {
        lhs[0] + rhs[0],
        lhs[1] + rhs[1],
        lhs[2] + rhs[2]
    };
};

auto constexpr div(const Vec3i& lhs, int denominator) -> Vec3i
{
    return
    {
        lhs[0] / denominator,
        lhs[1] / denominator,
        lhs[2] / denominator
    };
};

auto constexpr mul(const Vec3i& lhs, int multiple) -> Vec3i
{
    return
    {
        lhs[0] * multiple,
        lhs[1] * multiple,
        lhs[2] * multiple
    };
};

// //////////////////////////////////////////////////////

auto constexpr sub(const Vec3f& lhs, const Vec3f& rhs) -> Vec3f
{
    return
    {
        lhs[0] - rhs[0],
        lhs[1] - rhs[1],
        lhs[2] - rhs[2]
    };
};

auto constexpr add(const Vec3f& lhs, const Vec3f& rhs) -> Vec3f
{
    return
    {
        lhs[0] + rhs[0],
        lhs[1] + rhs[1],
        lhs[2] + rhs[2]
    };
};

auto constexpr div(const Vec3f& lhs, float denominator) -> Vec3f
{
    return
    {
        lhs[0] / denominator,
        lhs[1] / denominator,
        lhs[2] / denominator
    };
};

auto constexpr mul(const Vec3f& lhs, float multiple) -> Vec3f
{
    return
    {
        lhs[0] * multiple,
        lhs[1] * multiple,
        lhs[2] * multiple
    };
};

// //////////////////////////////////////////////////////

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

struct Position_and_quat
{
    Vec3i position;
    Quat quat;
};

// //////////////////////////////////////////////////////

struct Predictors
{
    Vec3f linear_velocity_per_frame;
    Vec3f linear_acceleration_per_frame;
    Vec3f angular_velocity_per_frame;
    Vec3f angular_acceleration_per_frame;
};

typedef std::array<Predictors, Cubes> Frame_predicitons;

// //////////////////////////////////////////////////////

auto constexpr position(const DeltaData& lhs) -> Vec3i
{
    return
    {
        lhs.position_x,
        lhs.position_y,
        lhs.position_z,
    };
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

unsigned MinBits(unsigned value)
{
    unsigned result = 0;

    while ((1u << result) <= value)
    {
        ++result;
    }

    return result;
}


static const constexpr unsigned bitmask_lookup[] =
{
    0x00000000,
    0x00000001,
    0x00000003,
    0x00000007,
    0x0000000F,

    0x0000001F,
    0x0000003F,
    0x0000007F,
    0x000000FF,

    0x000001FF,
    0x000003FF,
    0x000007FF,
    0x00000FFF,

    0x00001FFF,
    0x00003FFF,
    0x00007FFF,
    0x0000FFFF,

    0x0001FFFF,
    0x0003FFFF,
    0x0007FFFF,
    0x000FFFFF,

    0x001FFFFF,
    0x003FFFFF,
    0x007FFFFF,
    0x00FFFFFF,

    0x01FFFFFF,
    0x03FFFFFF,
    0x07FFFFFF,
    0x0FFFFFFF,

    0x1FFFFFFF,
    0x3FFFFFFF,
    0x7FFFFFFF,
    0xFFFFFFFF,
};

inline constexpr auto bit_mask(unsigned bits) -> unsigned
{
    return bitmask_lookup[bits];
}

// //////////////////////////////////////////////////////

using namespace Range_models;
struct Error_distance
{
    unsigned distance_squared;
    unsigned error;
};

std::vector<Error_distance> g_errors;

namespace Naive_error
{
    // Found emperically.
    // Changing restitution for cube 0 makes no difference.
    static const constexpr int LOWEST_POINT         = 38;
    static const constexpr int LOWEST_POINT_CUBE_0  = 367;
    static const constexpr float RESTITUTION        = 0.869;
    static const constexpr float DRAG               = 0.997;


    struct Model
    {
        // Ok, lets just get coding first before simplification
        Binary_two_speed has_error                      = {1, 7};
        Binary_two_speed has_quat_largest               = {1, 7};
        std::array<Binary_two_speed, 4> interactive     = {};

        // If I get error, send both pos and quat errors.
        Exp_update      quat_largest                    = {4, 8};
        Exp_update      error_signs                     = {8, 7};

        // Worst case == 10 bits of error.
        Periodic_update error_low_5_bits                = {32, 2};
        Exp_update error_high_5_bits                    = {32, 2};
    };

    auto predict
    (
        const Predictors& v_and_a,
        const Position_and_quat& base,
        int zero_height,
        unsigned frame_delta
    )
    -> Position_and_quat
    {
        // p = p0 + v0t + at^2/2
        // p = p0 + t(v0 + at/2)
        // p = p0 + tv
        // v = v0 + at/2
        auto at_2 =
            mul(v_and_a.linear_acceleration_per_frame, frame_delta / 2.0f);

        auto v = add(v_and_a.linear_velocity_per_frame, at_2);

        // RAM: TODO: angular!
        if ((std::abs(v[0]) > 0.0001f) || (std::abs(v[1]) > 0.0001f))
        {
            if (std::abs(v[2]) < 0.001f)
            {
                v = mul(v, DRAG);
            }
        }

        auto pos_delta = mul(v, frame_delta);

        auto pos = Vec3i
        {
            static_cast<int>(std::round(base.position[0] + pos_delta[0])),
            static_cast<int>(std::round(base.position[1] + pos_delta[1])),
            static_cast<int>(std::round(base.position[2] + pos_delta[2]))
        };

        // reflect z about lowest point.
        if (pos[2] < zero_height)
        {
            pos[2] = zero_height + RESTITUTION * (zero_height - pos[2]);
        }

        // Yay, this seems to work!
        auto wt_2 =
            mul(v_and_a.angular_acceleration_per_frame, frame_delta / 2.0f);

        auto w = add(v_and_a.angular_velocity_per_frame, wt_2);

        auto w_delta = mul(w, frame_delta);

        // Ok, need to convert to quat to do actual multiplications
        auto r = to_quat({w_delta});

        auto rotation = mul(r, base.quat);

        return
        {
            pos,
            rotation,
        };
    }

    auto update_prediciton
    (
        const Predictors& v_and_a,
        const Position_and_quat& base,
        const Position_and_quat& target,
        unsigned frame_delta
    )
    -> Predictors
    {
        float frame_delta_f = static_cast<float>(frame_delta);
        auto pos_delta = sub(target.position, base.position);
        auto v = Vec3f
        {
            pos_delta[0] / frame_delta_f,
            pos_delta[1] / frame_delta_f,
            pos_delta[2] / frame_delta_f,
        };

        auto v_delta = sub(v, v_and_a.linear_velocity_per_frame);
        auto a = div(v_delta, frame_delta);

        auto angle_delta = Rotor{0.0f,0.0f,0.0f};
        if
        (
            (base.quat[0] != target.quat[0])
            || (base.quat[1] != target.quat[1])
            || (base.quat[2] != target.quat[2])
            || (base.quat[3] != target.quat[3])
        )
        {
            // Hmm, seem we get stupid magnitudess due
            // to rotating near itself (from q to -q roughly).
            auto target_quat_neg = mul(target.quat, -1.0f);

            // http://www.geomerics.com/blogs/quaternions-rotations-and-compression/
            auto r                  = mul(target.quat, conjugate(base.quat));
            auto r_neg              = mul(target_quat_neg, conjugate(base.quat));
            auto rotor              = to_rotor(r);
            auto rotor_neg          = to_rotor(r_neg);
            auto mag_squared        = magnitude_squared(rotor);
            auto mag_squared_neg    = magnitude_squared(rotor_neg);

            angle_delta =
                (mag_squared < mag_squared_neg) ?
                    rotor :
                    rotor_neg;
        }

        auto w = div(angle_delta.vec, frame_delta);
        auto w_delta = sub(w, v_and_a.angular_velocity_per_frame);
        auto wa = div(w_delta, frame_delta);

        return
        {
            v,
            a,
            w,
            wa
        };
    }

    auto swap_orientation_largest
    (
        const Gaffer& quat,
        unsigned new_largest
    )
    -> Gaffer
    {
        // Due to errors, we cannot reliablly predict what the next
        // largest value is going to be. So just encode it with the
        // prediciton.
        auto max = 256;
        for (auto i = 0 ; i < 3; ++i)
        {
            if
            (
                std::abs(quat.vec[i] - 256)
                >
                std::abs(max - 256)
            )
            {
                max = quat.vec[i];
            }
        }

        int full_quat[4];

        {
            auto j = 0;
            for (unsigned i = 0; i < 3; ++i)
            {
                if (i == quat.orientation_largest)
                {
                    j++;
                }

                full_quat[j++] = quat.vec[i];
            }
        }

        full_quat[quat.orientation_largest] = max;

        auto result = Gaffer
        {
            new_largest,
            {0, 0, 0}
        };

        {
            auto j = 0;
            for (unsigned i = 0; i < 3; ++i)
            {
                if (i == new_largest)
                {
                    j++;
                }

                result.vec[i] = full_quat[j++];
            }
        }

        return result;
    }

    auto encode
    (
        const Frame& base,
        const Frame& target,
        Frame_predicitons& predicitons,
        unsigned frame_delta
    )
    -> Range_types::Bytes
    {
        auto                size = base.size();
        Range_types::Bytes  data;

        // //////////////////////////////////////////////////////

        Model model;

        {
            Range_coders::Encoder           range(data);
            Range_coders::Binary_encoder    binary(range);;

            // //////////////////////////////////////////////////////

            for (unsigned i = 0; i < size; ++i)
            {
                // Ok, get the predicted values, calculate the error
                // and if the error is non-zero, then we have some encoding
                // to do.
                auto g_b = to_gaffer(base[i]);
                auto g_t = to_gaffer(target[i]);
                auto q_b = to_quat(g_b);
                auto q_t = to_quat(g_t);

                // Right, for fun, lets see how good the predictor is.
                auto b = Position_and_quat
                {
                    position(base[i]),
                    q_b
                };

                auto t = Position_and_quat
                {
                    position(target[i]),
                    q_t
                };

                auto zero_height = i ? LOWEST_POINT : LOWEST_POINT_CUBE_0;

                auto calculated = predict
                (
                    predicitons[i],
                    b,
                    zero_height,
                    frame_delta
                );

                auto calculated_quat = to_gaffer(calculated.quat);

                // Update the predicitons for next time.
                predicitons[i] = update_prediciton
                (
                    predicitons[i],
                    b,
                    t,
                    frame_delta
                );

                auto error_pos = sub(position(target[i]), calculated.position);
                auto error_quat_largest =
                    calculated_quat.orientation_largest !=
                        static_cast<unsigned>(target[i].orientation_largest);

                if (error_quat_largest)
                {
                    calculated_quat =
                        swap_orientation_largest
                        (
                            calculated_quat,
                            static_cast<unsigned>
                            (
                                target[i].orientation_largest
                            )
                        );
                }

                auto error_quat = Vec3i
                {
                    target[i].orientation_a - calculated_quat.vec[0],
                    target[i].orientation_b - calculated_quat.vec[1],
                    target[i].orientation_c - calculated_quat.vec[2],
                };

                auto has_error_pos =
                    error_pos[0]
                    || error_pos[1]
                    || error_pos[2];

                auto has_error_quat =
                    error_quat[0]
                        || error_quat[1]
                        || error_quat[2];

                auto has_error =
                    has_error_pos
                    || has_error_quat
                    || error_quat_largest;

                // RAM: Store all the erros.
                if (i)
                {
                    if (has_error_pos)
                    {
                        // RAM: TODO: distance from line drawn by cube 0 please.
                        // you've done it before.
                        auto d = sub(position(base[i]), position(base[0]));
                        g_errors.push_back
                        ({
                            static_cast<unsigned>
                            (
                                d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
                            ),
                            static_cast<unsigned>
                            (
                                std::max
                                (
                                    std::abs(error_pos[0]),
                                    std::max
                                    (
                                        std::abs(error_pos[1]),
                                        std::abs(error_pos[2])
                                    )
                                )
                            )
                        });
                    }
                }

                // Encode!
                model.has_error.Encode(binary, has_error);

                if (has_error)
                {
                    model.has_quat_largest.Encode(binary, error_quat_largest);

                    if (error_quat_largest)
                    {
                        model.quat_largest.Encode
                        (
                            range,
                            target[i].orientation_largest
                        );
                    }

                    auto get_signs = [](const Vec3i& v) -> unsigned
                    {
                        unsigned result = 0;

                        if (v[0] < 0) { result |= 1; }
                        if (v[1] < 0) { result |= 2; }
                        if (v[2] < 0) { result |= 4; }

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

                    auto signs_pos  = get_signs(error_pos);
                    auto signs_quat = get_signs(error_quat);
                    auto vec_pos    = strip_signs(error_pos);
                    auto vec_quat   = strip_signs(error_quat);

                    for (auto v: vec_pos)
                    {
                        assert(v < (1 << 10));

                        model.error_high_5_bits.Encode(range, v >> 5);
                        model.error_low_5_bits.Encode
                        (
                            range,
                            v & ((1 << 5) - 1)
                        );
                    }

                    for (auto v: vec_quat)
                    {
                        assert(v < (1 << 10));

                        model.error_high_5_bits.Encode(range, v >> 5);
                        model.error_low_5_bits.Encode
                        (
                            range,
                            v & ((1 << 5) - 1)
                        );
                    }

                    if (vec_pos[0] || vec_pos[1] || vec_pos[2])
                    {
                        model.error_signs.Encode
                        (
                            range,
                            signs_pos
                        );
                    }

                    if (vec_quat[0] || vec_quat[1] || vec_quat[2])
                    {
                        model.error_signs.Encode
                        (
                            range,
                            signs_quat
                        );
                    }
                }

                // //////////////////////////////////////////////////////

                auto quat_changed = !quat_equal(base[i], target[i]);
                auto pos_changed = !pos_equal(base[i], target[i]);

                // //////////////////////////////////////////////////////

                // Note: You CAN get no interaction even if the quat or pos
                // changes.
                unsigned interact_lookup = base[i].interacting;
                if (pos_changed | quat_changed)
                {
                    interact_lookup |= 2;
                }

                model.interactive[interact_lookup].Encode
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
        Frame_predicitons& predicitons,
        const Range_types::Bytes& data,
        unsigned frame_delta
    )
    -> Frame
    {
        auto    size = base.size();
        Frame   target;        

        Model model;

        {
            Range_coders::Decoder           range(data);
            Range_coders::Binary_decoder    binary(range);

            for (unsigned i = 0; i < size; ++i)
            {
                auto g_b = to_gaffer(base[i]);
                auto q_b = to_quat(g_b);

                // Right, for fun, lets see how good the predictor is.
                auto b = Position_and_quat
                {
                    position(base[i]),
                    q_b
                };

                auto zero_height = i ? LOWEST_POINT : LOWEST_POINT_CUBE_0;

                auto calculated = predict
                (
                    predicitons[i],
                    b,
                    zero_height,
                    frame_delta
                );

                auto calculated_quat = to_gaffer(calculated.quat);

                auto has_error = model.has_error.Decode(binary);

                if (has_error)
                {
                    auto has_error_quat_largest =
                        model.has_quat_largest.Decode(binary);

                    if (has_error_quat_largest)
                    {
                        auto largest = model.quat_largest.Decode(range);

                        calculated_quat =
                            swap_orientation_largest(calculated_quat, largest);
                    }

                    // Get the errors and add them to things.
                    auto decode_vec = [&model, &range]() -> Vec3i
                    {
                        Vec3i result;

                        for (auto& v: result)
                        {
                            auto p = model.error_high_5_bits.Decode(range) << 5;
                            p += model.error_low_5_bits.Decode(range);
                            v = p;
                        }

                        return result;
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

                    auto vec_pos        = decode_vec();
                    auto vec_quat       = decode_vec();
                    unsigned signs_pos  = 0;
                    unsigned signs_quat = 0;

                    if (vec_pos[0] || vec_pos[1] || vec_pos[2])
                    {
                        signs_pos = model.error_signs.Decode(range);
                    }

                    if (vec_quat[0] || vec_quat[1] || vec_quat[2])
                    {
                        signs_quat = model.error_signs.Decode(range);
                    }

                    auto error_pos  = add_signs(signs_pos, vec_pos);
                    auto error_quat = add_signs(signs_quat, vec_quat);

                    calculated.position = add(calculated.position, error_pos);
                    calculated_quat.vec = add(calculated_quat.vec, error_quat);
                }

                target[i].orientation_largest =
                    calculated_quat.orientation_largest;

                target[i].orientation_a = calculated_quat.vec[0];
                target[i].orientation_b = calculated_quat.vec[1];
                target[i].orientation_c = calculated_quat.vec[2];

                target[i].position_x = calculated.position[0];
                target[i].position_y = calculated.position[1];
                target[i].position_z = calculated.position[2];

                // //////////////////////////////////////////////////////

                auto g_t = to_gaffer(target[i]);
                auto q_t = to_quat(g_t);

                auto t = Position_and_quat
                {
                    position(target[i]),
                    q_t
                };

                // Update the predicitons for next time.
                predicitons[i] = update_prediciton
                (
                    predicitons[i],
                    b,
                    t,
                    frame_delta
                );

                // //////////////////////////////////////////////////////

                auto quat_changed = !quat_equal(base[i], target[i]);
                auto pos_changed = !pos_equal(base[i], target[i]);

                // //////////////////////////////////////////////////////

                unsigned interact_lookup = base[i].interacting;
                if (pos_changed | quat_changed)
                {
                    interact_lookup |= 2;
                }

                target[i].interacting =
                    model.interactive[interact_lookup].Decode
                    (
                        binary
                    );
            }
        }

        return target;
    }
}

#define PRINT_INT(x) printf("%-32s\t%d\n", #x, x);
#define PRINT_FLOAT(x) printf("%-32s\t%f\n", #x, x);

void range_compress(std::vector<Frame>& frames)
{
    auto packets = frames.size();

    auto test = [&](auto encoder, auto decoder, const auto& title)
    {
        unsigned bytes = 0;
        unsigned packetsCoded = 0;
        unsigned min = 10000000;
        unsigned max = 0;

        Frame_predicitons predict_server = {};
        Frame_predicitons predict_client = {};

        for (size_t i = PacketDelta; i < packets; ++i)
        {
            auto buffer = encoder(
                frames[i-PacketDelta],
                frames[i],
                predict_server,
                PacketDelta);

            const unsigned size = buffer.size();
            bytes += size;
            if (bytes)
            {
                min = std::min(min, size);
                max = std::max(max, size);
            }

            if (do_decompress)
            {
                auto back = decoder(
                    frames[i-PacketDelta],
                    predict_client,
                    buffer,
                    PacketDelta);

                assert(back == frames[i]);
            }

            packetsCoded++;
        }

        float packetSizeAverge = ((float) bytes) / packetsCoded;
        float bytesPerSecondAverage = packetSizeAverge * 60.0f;
        float kbps = bytesPerSecondAverage * 8 / 1000.0f;

        printf("\n");
        printf("== Compression (%s) =======================\n\n", title);

        PRINT_INT(packetsCoded);
        PRINT_INT(bytes);
        PRINT_FLOAT(bytesPerSecondAverage);
        PRINT_FLOAT(packetSizeAverge);
        PRINT_INT(min);
        PRINT_INT(max);
        PRINT_FLOAT(kbps);

        printf("\n==============================================\n");
    };

    test
    (
        Naive_error::encode,
        Naive_error::decode,
        "Naive_error"
    );

    std::sort
    (
        begin(g_errors),
        end(g_errors),
        [](const Error_distance& lhs, const Error_distance& rhs)
        {
            return lhs.distance_squared < rhs.distance_squared;
        }
    );

    for (unsigned i = 0; i < 40; ++i)
    {
        printf
        (
            "Error: %d\tDistance: %f\n",
            g_errors[i].error,
            std::sqrt(g_errors[i].distance_squared)
        );
    }
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
    }

    if (do_compression)
    {
        range_compress(frames);
    }

    return 0;
}
