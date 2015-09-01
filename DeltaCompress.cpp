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
bool do_decompress  = false;

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

auto constexpr to_f(const Vec3i& i) -> Vec3f
{
    return
    {
        static_cast<float>(i[0]),
        static_cast<float>(i[1]),
        static_cast<float>(i[2])
    };
}

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

auto constexpr dot(const Vec3i& lhs, const Vec3i rhs) -> int
{
    return (lhs[0] * rhs[0]) + (lhs[1] * rhs[1]) + (lhs[2] * rhs[2]);
}

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

auto constexpr dot(const Vec3f& lhs, const Vec3f& rhs) -> float
{
    return
        lhs[0] * rhs[0] +
        lhs[1] * rhs[1] +
        lhs[2] * rhs[2];
};

auto constexpr normalise(const Vec3f& lhs) -> Vec3f
{
    return mul(lhs, 1.0f / std::sqrt(dot(lhs, lhs)));
}

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
    static const constexpr unsigned W_INDEX = 0;

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

// //////////////////////////////////////////////////////

struct Dual_quat
{
    Quat real;
    Quat dual;
};

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

auto constexpr add(const Quat& lhs, const Quat& rhs) -> Quat
{
    return
    {
        lhs[0] + rhs[0],
        lhs[1] + rhs[1],
        lhs[2] + rhs[2],
        lhs[3] + rhs[3]
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

auto constexpr dot(const Quat& lhs, const Quat& rhs) -> float
{
    return
        lhs[0]*rhs[0] +
        lhs[1]*rhs[1] +
        lhs[2]*rhs[2] +
        lhs[3]*rhs[3];
}

auto constexpr magnitude_squared(const Quat& quat) -> float
{
    return
        quat[0] * quat[0] +
        quat[1] * quat[1] +
        quat[2] * quat[2] +
        quat[3] * quat[3];
}

// NOTE: Interesting theory about normalisation from David Hammen
// http://stackoverflow.com/a/12934750 however need to convert that one to
// to use floats.
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

// //////////////////////////////////////////////////////

// Thanks Ben Kenwright
// http://www.xbdev.net/misc_demos/demos/dual_quaternions_beyond/paper.pdf

auto constexpr add(const Dual_quat& lhs, const Dual_quat& rhs) -> Dual_quat
{
    return
    {
        add(lhs.real, rhs.real),
        add(lhs.dual, rhs.dual)
    };
}

auto constexpr mul(const Dual_quat& lhs, const Dual_quat& rhs) -> Dual_quat
{
    return
    {
        mul(lhs.real, rhs.real),
        add(mul(lhs.real, rhs.dual), mul(lhs.dual, rhs.real))
    };
}

auto constexpr mul(const Dual_quat& lhs, float rhs) -> Dual_quat
{
    return
    {
        mul(lhs.real, rhs),
        mul(lhs.dual, rhs)
    };
}

auto constexpr conjugate(const Dual_quat& lhs) -> Dual_quat
{
    return
    {
        conjugate(lhs.real),
        conjugate(lhs.dual)
    };
}

// Norm(D) = sqrt(D.conjugate(D)) = ND.real, ND.dual = (1, 0)
// Thats why we only divide by the real part of the dual quat.
// I still think it's shakey becuase you also need the other constraint of:
// ND.real * conjugate(ND.dual) == conjugate(ND.real) * ND.dual.
// and I haven't seen proof anywhere for that :-(
auto constexpr normalise(const Dual_quat& lhs) -> Dual_quat
{
    return mul
    (
        lhs,
        1.0f / std::sqrt(dot(lhs.real, conjugate(lhs).real))
    );
}

auto constexpr to_dual(const Quat& rotation, const Vec3f& position) -> Dual_quat
{
    // RMA: TODO: Validate order of multiplication.
    return
    {
        rotation,
        mul({0.0f, position[0], position[1], position[2]}, mul(rotation, 0.5f))
    };
}

auto constexpr to_dual(const Vec3f& position) -> Dual_quat
{
    return
    {
        {1, 0, 0, 0},
        {0, position[0] * 0.5f, position[1] * 0.5f, position[2] * 0.5f}
    };
}

auto to_position(const Dual_quat& dual) -> Vec3f
{
    auto as_quat = mul(mul(dual.dual, conjugate(dual.real)), 2.0f).vec;

    return
    {
        as_quat[1],
        as_quat[2],
        as_quat[3]
    };
}

// //////////////////////////////////////////////////////

struct Screw
{
    float theta;
    Vec3f direction;

    float pitch;
    Vec3f moment;
};

// Oh hey, our friend cayley gets mentioned before.
// https://en.wikipedia.org/wiki/Rotation_matrix#Skew_parameters_via_Cayley.27s_formula

//static const constexpr float S_EPISLON  = 0.000001f;
static const constexpr float S_EPISLON2 = 0.00000001f;

auto to_screw(const Dual_quat& d) -> Screw
{
    // RAM: TODO: REad
    // http://www.seas.upenn.edu/~ladislav/kavan08geometric/kavan08geometric.pdf
    // To figure out how to screw a pure translation.
    // Appended A.3
    // I quote
    // (in the degenerate case with θ(0) = 0 or θ(0) = 2π, s(0) represents
    // the direction of the translation vector)
    // Note that for us θ(0) == theta, and s(0) = direction.
    // I assume that s(0) is a unit vec (is it?). so, would distance
    // still be valid?
    auto wr         = d.real[0];
    auto theta      = 2.0f * std::acos(wr);
    auto vd         = Vec3f{d.dual[1], d.dual[2], d.dual[3]};

    if (theta != 0.0f)
    {
        auto vr         = Vec3f{d.real[1], d.real[2], d.real[3]};
        auto vr_mag2    = dot(vr, vr);

        auto inv_vr_mag = (vr_mag2 > S_EPISLON2)
            ? 1.0f / std::sqrt(vr_mag2)
            : 0.0f;

        auto distance   = -2.0f * d.dual[0] * inv_vr_mag;
        auto direction  = mul(vr, inv_vr_mag);

        auto moment = mul
        (
            add(vd, mul(direction, -0.5f * distance * wr)),
            inv_vr_mag
        );

        return
        {
            theta,
            direction,
            distance,
            moment
        };
    }
    else
    {
        // Ok, this is where I have no idea what I'm doing :-(
        // I make sure the pitch is multiplied by -2 since the
        // pitch value is 1/2 the actual distance?
        auto vd_mag2 = dot(vd, vd);
        return
        {
            0.0f,
            vd_mag2 > S_EPISLON2 ? normalise(vd) : Vec3f{0.0f, 0.0f, 0.0f},
            vd_mag2 > S_EPISLON2 ? 2.0f * std::sqrt(vd_mag2) : 0,
            Vec3f{0.0f, 0.0f, 0.0f}
        };
    }
}

auto to_dual_quat(const Screw& s) -> Dual_quat
{
    auto wr             = std::cos(s.theta * 0.5f);
    auto sin_half_theta = std::sin(s.theta * 0.5f);
    auto vr             = mul(s.direction, sin_half_theta);
    auto wd             = -0.5f * s.pitch * sin_half_theta;

    auto vd = add
    (
        mul(s.moment,    sin_half_theta),
        mul(s.direction, 0.5f * s.pitch * wr)
    );

    return
    {
        {wr, vr[0], vr[1], vr[2]},
        {wd, vd[0], vd[1], vd[2]},
    };
}

// //////////////////////////////////////////////////////

// Ok, right, I have no idea if the maths here is even correct, but I'm
// going to give it a shot anyway.

struct Dual_angle_axis
{
    Vec3f real;
    Vec3f dual;
};

auto constexpr to_dual_angle_axis(const Screw& s) -> Dual_angle_axis
{
    // RAM: NO! angles are dual numbers, as are the vectors. Need to do dual mul.
    // Also, don't even try until you can figure out how to get the angles
    // back out of a dual number (what is the magnitude of a dual vector?).

    // https://en.wikipedia.org/wiki/Screw_theory states:
    // The multiplication of a screw S=(S, V) by the dual scalar â=(a, b)
    // is computed componentwise to be,
    // âS = (a, b)(S, V) = (aS, aV + bS) = (S', V')
    //
    // So, just working backwards. to get a and b back:
    // a = |S'|
    // S = S'/|S'|
    // b = how?
    // V = V' - (b/a)S
    //
    // Aha! We can get b by knowing that |v| == 1. reduce to a quadratic and
    // then solve that for b. Time for pen and paper :-)

    // RAM: MAybe this is what im looking for ?
    // http://rain.aa.washington.edu/@api/deki/files/401/=Dual_q_landing_CDC12.pdf
    // Equation 29
    // omega' = ( delta R, Delta D + cross(R, D) )
    //
    // omega = R + D, where R == real, D = dual part, and both are visors
    // i.e. have a zero w component. that the imaginary part is the
    // angular velocity vector.
    // only problem is that is D just bV, or is it aV + bS ?
    //
    // Another deveriation (for quanterions only atm) is from:
    // http://www.mbnexplorer.com/users-guide/5-molecular-dynamics-simulations/51-equations-motion/514-quaternions-application-rigid
    // quation 5.32 and 5.33
    // q'  = 0.5 * mul(q, R)
    // q'' = 0.5 * add(mul(q', R), mul(q, R'))
    // But then, how do I get R' ?

    // RAM: Note that theta and distance are not half values, but whole values.
    if (s.theta != 0)
    {
        return
        {
            mul(s.direction, s.theta),
            mul(s.moment, s.pitch),
        };
    }
    else
    {
        // to_screw explains this stupidity.
        return
        {
            Vec3f{0.0f, 0.0f, 0.0f},
            mul(s.direction, s.pitch),
        };
    }
}

auto to_screw(const Dual_angle_axis& a) -> Screw
{
    auto theta = std::sqrt(dot(a.real, a.real));

    if (theta != 0)
    {
        auto direction = normalise(a.real);

        // RAM: This is the dodgy bit - is it even valid?

        auto distance = std::sqrt(dot(a.dual, a.dual));

        // Should I return 0,0,0 for zero moment? or an axis of some descrption?
        auto moment =
            (distance * distance > 0.00000001f)
            ? normalise(a.dual)
            : Vec3f{0.0f, 0.0f, 0.0f};

        return
        {
            theta,
            direction,
            distance,
            moment
        };
    }
    else
    {
        auto distance = std::sqrt(dot(a.dual, a.dual));

        if (distance != 0.0f)
        {
            auto direction = normalise(a.dual);

            return
            {
                theta,
                direction,
                distance,
                Vec3f{0.0f, 0.0f, 0.0f}
            };
        }
        else
        {
            return
            {
                0.0f,
                Vec3f{0.0f, 0.0f, 0.0f},
                0.0f,
                Vec3f{0.0f, 0.0f, 0.0f},
            };
        }
    }
}


// RAM: TODO: DUAL paper with maybe useful info
// http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3576712/

// //////////////////////////////////////////////////////

struct Angle_and_axis
{
    float angle;
    Vec3f axis;
};

auto to_angle_and_axis(const Quat& q) -> Angle_and_axis
{
    auto theta      = 2.0f * std::acos(q[0]);
    auto vec_raw    = Vec3f{q[1],q[2],q[3]};
    auto mag2       = dot(vec_raw, vec_raw);

    // RAM: TODO: Use a proper epislon.
    if (mag2 < 0.000001)
    {
        return
        {
            0.0f,
            {1.0f, 0.0f, 0.0f}
        };
    }

    auto vec = mul(vec_raw, 1.0f / std::sqrt(mag2));

    // shortest arc.
    if (theta > M_PI)
    {
        theta -= 2.0f * M_PI;
    }

    return
    {
        theta,
        vec
    };
}

auto constexpr to_quat(const Angle_and_axis& aaa) -> Quat
{
    return
    {
        std::cos(aaa.angle * 0.5f),
        aaa.axis[0] * std::sin(aaa.angle * 0.5f),
        aaa.axis[1] * std::sin(aaa.angle * 0.5f),
        aaa.axis[2] * std::sin(aaa.angle * 0.5f)
    };
}

// //////////////////////////////////////////////////////

// Where the unit axis is multiplied by the angle itself (in radians).
struct Angle_axis
{
    Vec3f vec;
};

auto to_angle_axis(const Quat& q) -> Angle_axis
{
    auto theta  = 2.0f * std::acos(q[0]);
    auto vec    = normalise(Vec3f{q[1],q[2],q[3]});

    // shortest arc.
    if (theta > M_PI)
    {
        theta -= 2.0f * M_PI;
    }

    return
    {
        mul(vec, theta)
    };
}

auto to_quat(const Angle_axis& aa) -> Quat
{
    auto theta = std::sqrt(dot(aa.vec, aa.vec));

    if (theta * theta < 0.00000001f)
    {
        return {1,0,0,0};
    }

    auto vec   = normalise(aa.vec);

    return
    {
        std::cos(theta * 0.5f),
        vec[0] * std::sin(theta * 0.5f),
        vec[1] * std::sin(theta * 0.5f),
        vec[2] * std::sin(theta * 0.5f)
    };
}

// //////////////////////////////////////////////////////

auto delta_t(const Quat& q, float delta_t) -> Quat
{
    if ((q[0] == 1.0f) && (q[1] == 0.0f)  && (q[2] == 0.0f)  && (q[3] == 0.0f))
    {
        return {1.0f, 0.0f, 0.0f, 0.0f};
    }

    auto axis_angle = to_angle_axis(q);

    axis_angle.vec = mul(axis_angle.vec, 1.0f / delta_t);

    return to_quat(axis_angle);
}

// //////////////////////////////////////////////////////


// RAM: Ok, now I know. I need to do verlet integration
//      to calculate the new position. Secondly, that requires
//      accleration. My naieve way of calculating acceleration
//      didn't work. So I had to guess a formula that would.
//      I suspect that my guessed formula is related to verlet
//      somehow.
//      https://en.wikipedia.org/wiki/Verlet_integration
//
//      for acceleration calculation:
//      x0   = x at t0
//      x-1  = x at t-1
//      x-2  = x at t-2
//      dt0  = t0  - t-1
//      dt-1 = t-1 - t-2
//
//      calculation:
//      vt0  = (x0 - x-1)   / dt0
//      vt-1 = (x-1 - x-2)  / dt-1
//      dv   = (vt0 - vt-1) / dt0
//      a    = dv / ((dt0 + dt-1)
//
//      RAM: Argh, I cannot reproduce proper results with the above
//           calculation. I had to rederive one.
//
//      a    = (2 * (vt0 - vt-1)) / (t0 - t-2)
//
//
//      To calculate the next position:
//      x1  = next posisiton
//      x0  = current position
//      x-1 = previous position
//      dt1 = t1 - t0 (the time delta you're predicting over)
//      dt0 = t0 - t-1
//      a
//
//      calculation (verlet):
//      dx = x0 - x-1
//      a  = x0
//      b  = (dx * dt1) / dt0
//      c  = 0.5 * a * (dt0 + dt1) * dt1
//      x1 = a + b + c
//
//
//     struct Predict_data
//     {
//         float     acceleration;
//         Position  x-1;
//         float     dt;
//     }
//
//     auto Update(x0, x-1, x-2, dt0, dt-1) -> Predict_data;
//     auto Predict(Predict_data, x0, predict_t) -> x1;

struct Prediciton_data
{
    Vec3f acceleration;
    Vec3f x0;
    Vec3f xm1;
    float dt0;
};

auto verlet_acceleration
(
    Vec3f x0,
    Vec3f xm1,
    Vec3f xm2,
    float dt0,
    float dtm1
)
-> Prediciton_data
{
    if (dt0 == 0)
    {
        return
        {
            {0.0f, 0.0f, 0.0f},
            x0,
            xm1,
            dt0
        };
    }

    auto dx0  = sub(x0, xm1);
    auto dx1  = sub(xm1, xm2);
    auto vt0  = div(dx0, dt0);
    auto vtm1 = (dtm1 != 0.0f) ? div(dx1, dtm1) : Vec3f{0.0f, 0.0f, 0.0f};
    auto dvt  = sub(vt0, vtm1);
    auto dv   = div(dvt, dt0);

    auto a = div(dv, dt0 + dtm1);

    return
    {
        a,
        x0,
        xm1,
        dt0
    };
}

auto verlet_prediction
(
    Prediciton_data prediction,
    float dt1
)
-> Vec3f
{
    if ((dt1 == 0) || (prediction.dt0 == 0))
    {
        return prediction.x0;
    }

    auto dx = sub(prediction.x0, prediction.xm1);
    auto a  = prediction.x0;
    auto b  = div(mul(dx, dt1), prediction.dt0);
    auto c  = mul(prediction.acceleration, 0.5 * (prediction.dt0 + dt1) * dt1);
    auto x1 = add(a, add(b, c));

    return x1;
}

// //////////////////////////////////////////////////////

struct Predictors
{
    Vec3f linear_velocity_per_frame;
    Vec3f linear_acceleration_per_frame;
    Vec3f angular_velocity_per_frame;
    Vec3f angular_acceleration_per_frame;

    Prediciton_data verlet;
};

typedef std::array<Predictors, Cubes> Frame_predicitons;

// //////////////////////////////////////////////////////

void dual_tests()
{
    static const float EPISLON = 0.00001f;
    static const float FRAME_DELTA = 6;

    static const constexpr unsigned TEST_COUNT = 8;
    struct Posisiton_angle_axis
    {
        Vec3f pos;
        Angle_and_axis aa;
    };

    std::vector<std::vector<Posisiton_angle_axis>> tests;

    {
        // testing to_dual from a position and quat to make sure
        // multiplication order is correct - it is.
        auto q = normalise(Quat{0.3,0.5,0.5,0.0});

        auto dq = to_dual(q, {1.0f, 10.0f, -5.0f});

        auto p = to_position(dq);
        auto q2 = dq.real;

        assert(std::abs(p[0] - 1.0f)  < EPISLON);
        assert(std::abs(p[1] - 10.0f) < EPISLON);
        assert(std::abs(p[2] + 5.0f)  < EPISLON);


        assert(std::abs(q2[0] - q[0]) < EPISLON);
        assert(std::abs(q2[1] - q[1]) < EPISLON);
        assert(std::abs(q2[2] - q[2]) < EPISLON);
        assert(std::abs(q2[3] - q[3]) < EPISLON);
    }

    // constant velocity
    {
        tests.push_back({});
        for (unsigned i = 0; i < TEST_COUNT; ++i)
        {
            auto t = Posisiton_angle_axis
            {
                {i * 10.0f, 0.0f, 0.0f},
                {0.0f, {0.0f, 0.0f, 1.0f}}
            };

            tests.back().push_back(t);
        }
    }

    // constant acceleration
    {
        tests.push_back({});
        for (unsigned i = 0; i < TEST_COUNT; ++i)
        {
            // d = ut + 1/2 at^2
            auto position =
                0.5f * 10.0f * (i * FRAME_DELTA) * (i * FRAME_DELTA);

            auto t = Posisiton_angle_axis
            {
                {position, 0.0f, 0.0f},
                {0.0f, {0.0f, 0.0f, 1.0f}}
            };

            tests.back().push_back(t);
        }
    }

    // constant angular velocity, no offset
    {
        tests.push_back({});
        for (unsigned i = 0; i < TEST_COUNT; ++i)
        {
            auto t = Posisiton_angle_axis
            {
                {0.0f, 0.0f, 0.0f},
                {
                    static_cast<float>(2.0f * M_PI * i / 100.0f),
                    {0.0f, 0.0f, 1.0f}
                }
            };

            tests.back().push_back(t);
        }
    }

    // constant angular velocity
    {
        tests.push_back({});
        for (unsigned i = 0; i < TEST_COUNT; ++i)
        {
            auto t = Posisiton_angle_axis
            {
                {10.0f, 10.0f, 10.0f},
                {
                    static_cast<float>(2.0f * M_PI * i / 100.0f),
                    {0.0f, 0.0f, 1.0f}
                }
            };

            tests.back().push_back(t);
        }
    }

//    // constant angular acceleration
//    {
//        tests.push_back({});
//        for (unsigned i = 0; i < TEST_COUNT; ++i)
//        {
//            static const constexpr float wa = 2.0f * M_PI / 1000.0f;

//            // d = ut + 1/2 at^2
//            auto angle = 0.5f * wa * (i * FRAME_DELTA)*(i * FRAME_DELTA);

//            auto t = Posisiton_angle_axis
//            {
//                {10.0f, 10.0f, 10.0f},
//                {
//                    angle,
//                    {0.0f, 0.0f, 1.0f}
//                }
//            };

//            tests.back().push_back(t);
//        }
//    }

//    // //////////////////////////////////////////////////////

//    // constant velocity + angular velocity
//    {
//        tests.push_back({});
//        for (unsigned i = 0; i < TEST_COUNT; ++i)
//        {
//            auto t = Posisiton_angle_axis
//            {
//                {i * 10.0f, 0.0f, 0.0f},
//                {
//                    static_cast<float>(2.0f * M_PI * i / 100.0f),
//                    {0.0f, 0.0f, 1.0f}
//                }
//            };

//            tests.back().push_back(t);
//        }
//    }

//    // constant acceleration + angular acceleration
//    {
//        tests.push_back({});
//        for (unsigned i = 0; i < TEST_COUNT; ++i)
//        {
//            auto position =
//                    0.5f * 10.0f * (i * FRAME_DELTA) * (i * FRAME_DELTA);

//            static const constexpr float wa = 2.0f * M_PI / 1000.0f;

//            // d = ut + 1/2 at^2
//            auto angle = 0.5f * wa * (i * FRAME_DELTA)*(i * FRAME_DELTA);

//            auto t = Posisiton_angle_axis
//            {
//                {position, 0.0f, 0.0f},
//                {
//                    angle,
//                    {0.0f, 0.0f, 1.0f}
//                }
//            };

//            tests.back().push_back(t);
//        }
//    }


    // //////////////////////////////////////////////////////

    auto compare_q = []
    (
        const Quat& lhs,
        const Quat& rhs,
        const float EPISLON
    )
    {
        std::vector<float> errors;
        errors.reserve(4);

        errors.push_back(std::abs(lhs[0] - rhs[0]));
        errors.push_back(std::abs(lhs[1] - rhs[1]));
        errors.push_back(std::abs(lhs[2] - rhs[2]));
        errors.push_back(std::abs(lhs[3] - rhs[3]));

        assert(errors[0] < EPISLON);
        assert(errors[1] < EPISLON);
        assert(errors[2] < EPISLON);
        assert(errors[3] < EPISLON);
    };

    auto compare_dq = []
    (
        const Dual_quat& lhs,
        const Dual_quat& rhs,
        const float EPISLON
    )
    {
        std::vector<float> errors;
        errors.reserve(8);

        errors.push_back(std::abs(lhs.real[0] - rhs.real[0]));
        errors.push_back(std::abs(lhs.real[1] - rhs.real[1]));
        errors.push_back(std::abs(lhs.real[2] - rhs.real[2]));
        errors.push_back(std::abs(lhs.real[3] - rhs.real[3]));

        errors.push_back(std::abs(lhs.dual[0] - rhs.dual[0]));
        errors.push_back(std::abs(lhs.dual[1] - rhs.dual[1]));
        errors.push_back(std::abs(lhs.dual[2] - rhs.dual[2]));
        errors.push_back(std::abs(lhs.dual[3] - rhs.dual[3]));

        assert(errors[0] < EPISLON);
        assert(errors[1] < EPISLON);
        assert(errors[2] < EPISLON);
        assert(errors[3] < EPISLON);

        assert(errors[4] < EPISLON);
        assert(errors[5] < EPISLON);
        assert(errors[6] < EPISLON);
        assert(errors[7] < EPISLON);
    };


    auto compare = []
    (
        const Vec3f& lhs,
        const Vec3f& rhs,
        const float EPISLON
    )
    {
        std::vector<float> errors;
        errors.reserve(3);

        errors.push_back(std::abs(lhs[0] - rhs[0]));
        errors.push_back(std::abs(lhs[1] - rhs[1]));
        errors.push_back(std::abs(lhs[2] - rhs[2]));

        assert(errors[0] < EPISLON);
        assert(errors[1] < EPISLON);
        assert(errors[2] < EPISLON);
    };

    Dual_quat previous =
    {
        {1.0f, 0.0f, 0.0f, 0.0f},
        {0.0f, 0.0f, 0.0f, 0.0f}
    };


    Dual_quat previous_delta =
    {
        {1.0f, 0.0f, 0.0f, 0.0f},
        {0.0f, 0.0f, 0.0f, 0.0f}
    };

    const auto test_count = tests.size();
    for (unsigned test_index = 0; test_index < test_count; ++test_index)
    {
        const auto& test = tests[test_index];

        const auto item_count = test.size();
        for (unsigned item_index = 0; item_index < item_count; ++item_index)
        {
            const auto& item = test[item_index];

            Dual_quat dq = to_dual
            (
                to_quat(item.aa),
                item.pos
            );

            if (item_index > 2)
            {
                const auto& item_1 = test[item_index - 1];
                const auto& item_2 = test[item_index - 2];
                const auto& item_3 = test[item_index - 3];

                {
                    // Just just multiplying.
                    auto to_dq = [](const auto& i) -> Dual_quat
                    {
                        return to_dual
                        (
                            to_quat(i.aa),
                            i.pos
                        );
                    };

                    // returns velocity * time (or just delta)
                    auto delta = [](const auto &b, const auto& t) -> Dual_quat
                    {
                        return mul(t, conjugate(b));
                    };

                    auto b = to_dq(item_3);
                    auto c = to_dq(item_2);
                    auto t = to_dq(item_1);
                    //auto t2 = to_dq(test[item_index]);

                    auto vt = delta(b, c);

                    auto g = mul(vt, b);

                    compare_dq(g, c, 0.0001f);

                    // Did some maths
                    // d = t/2(u + v)
                    // so at t = 2: d = u + v. t==2==item t.
                    auto v1 = delta(b, c);
                    auto v2 = delta(c, t);

                    auto g2 = mul(mul(v2, v1), b);

                    compare_dq(g2, t, 0.0001f);

                    // RAM: TODO: Next step is getting ln and exp for a dual
                    // quat and using that.
                }

                // RAM: TODO: Get acceleration not by angle-axis vectors
                // but by quats as you cannot add angleaxis vectors apparently?
                {
                    // Hack in acceleration of 0.1
                    auto acc = [](float t) -> float
                    {
                        return 0.1f * 0.5f * (t * t);
                    };

                    auto values = std::vector<Quat>
                    {
                        to_quat(Angle_and_axis{acc(0.0f), {0.0f, 1.0f, 0.0f}}),
                        to_quat(Angle_and_axis{acc(1.0f), {0.0f, 1.0f, 0.0f}}),
                        to_quat(Angle_and_axis{acc(2.0f), {0.0f, 1.0f, 0.0f}}),
                        to_quat(Angle_and_axis{acc(3.0f), {0.0f, 1.0f, 0.0f}}),
                        to_quat(Angle_and_axis{acc(4.0f), {0.0f, 1.0f, 0.0f}})
                    };

                    auto v1t = mul(values[1], conjugate(values[0]));
                    auto v2t = mul(values[2], conjugate(values[1]));

                    auto v1t_aa = to_angle_and_axis(v1t);
                    auto v2t_aa = to_angle_and_axis(v2t);

                    // lets do this the proper way.
                    const auto t = 10.0f;

                    auto v1_aa = v1t_aa;
                    v1_aa.angle /= t;

                    auto v2_aa = v2t_aa;
                    v2_aa.angle /= t;

                    auto v1 = to_quat(v1_aa);
                    auto v2 = to_quat(v2_aa);

                    auto at = mul(v2, conjugate(v1));

                    auto a_aa = to_angle_and_axis(at);
                    a_aa.angle /= t;

                    auto a_calc = to_quat(a_aa);

                    // Great. Precision issues means acc is identity :-(
                    auto to_at2_2 = [](const Quat&a, float t) -> Quat
                    {
                        auto aa = to_angle_and_axis(a);
                        aa.angle *= t * t;
                        return to_quat(aa);
                    };

                    // Doh, I assume acceleration is relative to the first
                    // point? (assume point 0 is identity, so skip mul).
                    // anyway - these values are more or less correct
                    // however even then they have precision issue of up to 0.01.
                    // Assuming that q1a*q2b = ab(q1*q2) then just need to divide
                    // by ab to get a.
                    auto c_p1 = to_quat({a_aa.angle * 0.5f * 10.0f * 10.0f, a_aa.axis});
                    auto c_p2 = to_quat({a_aa.angle * 0.5f * 20.0f * 20.0f, a_aa.axis});
                    auto c_p3 = to_quat({a_aa.angle * 0.5f * 30.0f * 30.0f, a_aa.axis});
                    auto c_p4 = to_quat({a_aa.angle * 0.5f * 40.0f * 40.0f, a_aa.axis});





                    // LOL, wut? Can you even do that?                    
                    // RAM: YES! And more accurate than above!
                    // RAM: Downside is that it assumes t is constant
                    // RAM: for v2t and v1t.
                    auto at2 = mul(v2t, conjugate(v1t));

                    // Now, lets get some things striaght.
                    auto t_v1 = mul(v1t, values[0]);
                    compare_q(t_v1, values[1], 0.00001f);

                    auto t_v2 = mul(v2t, values[1]);
                    compare_q(t_v2, values[2], 0.00001f);

                    // d = delta change!
                    // d = ut + 0.5at^2
                    // d = d0 + v1t - 0.5at^2
                    // assume t = 10?
                    // do divide or mult by time, we need axis_angle version.
                    auto at2_aa = to_angle_and_axis(at2);
                    auto a = at2_aa;
                    a.angle /= 10.0f * 10.0f;


                    // Yup, these are correct.
                    auto d_p1 = to_quat({a.angle * 0.5f * 10.0f * 10.0f, a.axis});
                    auto d_p2 = to_quat({a.angle * 0.5f * 20.0f * 20.0f, a.axis});
                    auto d_p3 = to_quat({a.angle * 0.5f * 30.0f * 30.0f, a.axis});
                    auto d_p4 = to_quat({a.angle * 0.5f * 40.0f * 40.0f, a.axis});


                    // RAM: But we're not working with only acc, how to
                    // incoperate ut as well?


                    // now, get 0.5at^2
                    auto at2_2_aa = a;
                    at2_2_aa.angle *= 10 * 10 * 0.5;
                    auto at2_2 = to_quat(at2_2_aa);

                    // d = ut + 0.5at^2
                    // p = p0 + d

                    // first point, u == 0
                    auto d1 = at2_2;

                    // this should be right!
                    auto p1 = mul(d1, values[0]);

                    // second point
                    auto d2 = at2_2;
                    auto p2 = mul(d2, values[1]);

                    auto d3 = mul(v1t, at2_2);
                    auto p3 = mul(d2, values[1]);


                    compare_q(p1, values[1], 0.00001f);
                    compare_q(p2, values[2], 0.00001f);
                    compare_q(p3, values[3], 0.00001f);
                }

//                {
//                    // Test reflection to see if it works.
//                    // RAM: No it doesn't :-( I don't know what to believe now.
//                    Dual_quat dq_0 = to_dual
//                    (
//                        to_quat(item_1.aa),
//                        item_1.pos
//                    );

//                    Dual_quat dq_1 = to_dual
//                    (
//                        to_quat(item_2.aa),
//                        item_2.pos
//                    );

//                    // Reflect d0(d-1)d0 should give us d1
//                    auto dq_calc = mul(mul(dq_0, dq_1), dq_0);

//                    compare_dq(dq, dq_calc, EPISLON);
//                }

                // First deduce the "previous velocity"
                const auto velocity = []
                (
                    const Posisiton_angle_axis& i_0,
                    const Posisiton_angle_axis& i_1
                )
                -> Dual_angle_axis
                {
                    Dual_quat dq_0 = to_dual
                    (
                        to_quat(i_0.aa),
                        i_0.pos
                    );

                    Dual_quat dq_1 = to_dual
                    (
                        to_quat(i_1.aa),
                        i_1.pos
                    );

                    auto delta = mul(dq_1, conjugate(dq_0));

                    auto screw_delta = to_screw(delta);
                    auto dual_aa = to_dual_angle_axis(screw_delta);

                    // Ok, does this work?
                    // I mean either this maths works, or is 100%
                    // shonky.

                    // RAM: Ok, I think the values of the angles are for
                    // HALF the angle/distance since it's reflected.
                    // I also fear for dual quats, it's actually only 1/4
                    // since four reflections are involved? Ugh.
                    auto velocity = Dual_angle_axis
                    {
                        mul(dual_aa.real, 1.0f / FRAME_DELTA),
                        mul(dual_aa.dual, 1.0f / FRAME_DELTA)
                    };

                    return velocity;
                };

                auto v_0 = velocity(item_3, item_2);
                auto v_1 = velocity(item_2, item_1);

                Dual_angle_axis acc =
                {
                    mul(sub(v_1.real, v_0.real), 1.0f / FRAME_DELTA),
                    mul(sub(v_1.dual, v_0.dual), 1.0f / FRAME_DELTA)
                };

                // Ok, now calculate R, and see if we get the correct
                // predicted value.

                // p = p0 + v0t + at^2/2
                // p = p0 + t(v0 + at/2)
                // p = p0 + tv
                // v = v0 + at/2

                auto at_2 = Dual_angle_axis
                {
                    mul(acc.real, FRAME_DELTA / 1.0f),
                    mul(acc.dual, FRAME_DELTA / 1.0f),
                };

                auto v = Dual_angle_axis
                {
                    add(at_2.real, v_1.real),
                    add(at_2.dual, v_1.dual),
                };

                auto delta_aa = Dual_angle_axis
                {
                    mul(v.real, FRAME_DELTA),
                    mul(v.dual, FRAME_DELTA)
                };

                // now, how to get back to screw?
                auto screw_delta_new = to_screw(delta_aa);
                auto delta_new = to_dual_quat(screw_delta_new);

                Dual_quat previous_d = to_dual
                (
                    to_quat(item_1.aa),
                    item_1.pos
                );

                auto dq_calc = mul(delta_new, previous_d);

                {
                    // RAM: WTF, when I assert, the norm has a value > 1.0
                    // Is it a rounding issue, or a maths issue?
                    // Seems to be a precision issue :-(
                    // comparing a,b,c to spreadsheet it comes under 1.0
                    // but the spreadsheet says abouve 1.0.
                    // question is: will it affect the angle and axis?
                    auto dq_calc_norm = normalise(dq_calc);
                    auto dq_calc_norm2 = normalise(dq_calc_norm);


                    //auto mag2 = dot(dq_calc.real, conjugate(dq_calc).real);
                    auto a = dq_calc.real[0] * dq_calc.real[0];
                    auto b = dq_calc.real[3] * dq_calc.real[3];
                    auto mag2c = a+b;
    //                auto mag2b =
    //                    a +
    //                    dq_calc.real[1] * dq_calc.real[1] +
    //                    dq_calc.real[2] * dq_calc.real[2] +
    //                    b;


                    auto s = std::sqrt(mag2c);
                    auto s_1 = 1.0f / s;
                    auto dq_calc_norm_m = mul(dq_calc, s_1);
                    mul(dq_calc_norm, 1.0f);
                    mul(dq_calc_norm2, 1.0f);
                    mul(dq_calc_norm_m, 1.0f);
                }



                // convert back to angle axis to help debugging.
                auto aa_calc = to_angle_and_axis(dq_calc.real);
                assert(std::abs(aa_calc.angle - item.aa.angle) < EPISLON);
                auto pos_calc = to_position(dq_calc);
                assert(std::abs(pos_calc[0] - item.pos[0]) < EPISLON);
                assert(std::abs(pos_calc[1] - item.pos[1]) < EPISLON);
                assert(std::abs(pos_calc[2] - item.pos[2]) < EPISLON);

                compare_dq(dq, dq_calc, EPISLON);

            }

            // First, make sure that our "delta" can be converted
            // back and forward.
            {
                previous_delta = mul(dq, conjugate(previous));

                auto dq_calc = mul(previous_delta, previous);

                compare_dq(dq, dq_calc, EPISLON);

                auto pos = to_position(dq_calc);

                compare(item.pos, pos, EPISLON);
            }

            // Now try and predict the future!
            {
                auto dq_calc = mul(previous_delta, previous);

                auto pos = to_position(dq_calc);
                auto rotation = normalise(dq_calc).real;

                compare(item.pos, pos, EPISLON);

                if (rotation[0] < 1.0f)
                {
                    auto aa = to_angle_and_axis(rotation);
                    compare(item.aa.axis, aa.axis, EPISLON);
                    assert(std::abs(item.aa.angle - aa.angle) < EPISLON);
                }
            }
        }
    }
}



// //////////////////////////////////////////////////////

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

// //////////////////////////////////////////////////////

using namespace Range_models;
struct Error_distance
{
    unsigned distance_floor;
    unsigned distance_floor_calc;
    unsigned distance_squared_cube_0;
    float    shortest_distance_between_movements;
    unsigned velocity_bits;
    float    velocity;
    float    angular_v;
    unsigned error;
};

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

        // This seems to do the trick.
        Unsigned_golomb_range error_bits                = {10, 5};
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
        Quat r;

        r = to_quat(Angle_axis{w_delta});

        auto rotation = mul(r, base.quat);

        // RAM: How is verlet doing?
        // RAM: It sufferes from initial conditions, but seems much
        // RAM: better after that :-)
        // RAM: Pity the intial condition error > encode bits, so I
        // RAM: need to improve that.
        {
            // RAM: My interface is broken if this works.
            auto t = v_and_a.verlet;
            t.x0 = to_f(base.position);

            auto verlet = verlet_prediction
            (
                t,
                frame_delta
            );

            // RAM: DO IT!
            pos =
            {
                static_cast<int>(std::round(verlet[0])),
                static_cast<int>(std::round(verlet[1])),
                static_cast<int>(std::round(verlet[2])),
            };
        }

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

        auto angle_delta = Vec3f{0.0f,0.0f,0.0f};
        if
        (
            (base.quat[0] != target.quat[0])
            || (base.quat[1] != target.quat[1])
            || (base.quat[2] != target.quat[2])
            || (base.quat[3] != target.quat[3])
        )
        {
            auto r = mul(target.quat, conjugate(base.quat));
            angle_delta = to_angle_axis(r).vec;
        }

        auto w = div(angle_delta, frame_delta);
        auto w_delta = sub(w, v_and_a.angular_velocity_per_frame);
        auto wa = div(w_delta, frame_delta);

        // Do verlet predicition too.
        auto data = verlet_acceleration
        (
            to_f(target.position),
            v_and_a.verlet.x0,
            v_and_a.verlet.xm1,
            frame_delta,
            v_and_a.verlet.dt0
        );

        return
        {
            v,
            a,
            w,
            wa,

            data
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

    //===================================================================

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

            // RAM: TODO: Change size back please!
            for (unsigned i = 0; i < 1; ++i)
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

                    auto encode_error_bits = [&model, &range](const Vec3i& vec)
                    {
                        for (auto v: vec)
                        {
                            assert(v < (1 << 10));

                            model.error_bits.Encode(range,v);
                        }
                    };

                    encode_error_bits(vec_pos);
                    encode_error_bits(vec_quat);

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
                            v = model.error_bits.Decode(range);
                        }

                        return result;
                    };

                    Vec3i vec_pos;
                    Vec3i vec_quat;

                    auto add_signs = [](unsigned signs, const Vec3i& v) -> Vec3i
                    {
                        return
                        {
                            (signs & 1) ? -v[0] : v[0],
                            (signs & 2) ? -v[1] : v[1],
                            (signs & 4) ? -v[2] : v[2],
                        };
                    };

                    vec_pos        = decode_vec();
                    vec_quat       = decode_vec();

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

    // RAM: NOTES: For position, 6 frames error doesn't go above 145
    //      but for quats it gets as high as 450 (close and near to cube 0)
    //      also distance between segments and distance between point and
    //      cube 0 are pretty much the same looking distrubtuion, so just
    //      using the point to segment distance.
    //
    //      For quat errors, worst errors seem to happen when the angular
    //      velocity is near zero. range seems to be up to 0.12 / frame
    //      this is component velocity, not vector magnitude.
    //
    //      Revisited error when close to ground, but for predicted
    //      Z as opposed to base Z. Have 3 peaks. No errors under 23 (well, duh)
    //      Put peaks then toughs at 150. A second peak around 180 to 334, then
    //      the worsrt peak (but not as dense) at 900, before dropping off at
    //      1600 and stabilising to errors under 20. probably related to
    //      cube 0 as well.
    //
    //      For distance to bottom I'll use base distance from floor and
    //      anything under 700 for a simple 2 context model.
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
        dual_tests();
        range_tests();
    }

    if (do_compression)
    {
        range_compress(frames);
    }

    return 0;
}
