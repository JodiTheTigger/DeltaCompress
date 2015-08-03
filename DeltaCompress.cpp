// Copyright 2015 Richard Maxwell, all rights reserved.

// http://gafferongames.com/2015/03/14/the-networked-physics-data-compression-challenge/

// g++ -std=c++14 DeltaCompress.cpp -Wall -Wextra -Werror -g -o DeltaCompress

// //////////////////////////////////////////////////////

#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

#include <vector>
#include <array>
#include <cstdint>
#include <algorithm>
#include <cassert>
#include <type_traits>
#include <chrono>

// //////////////////////////////////////////////////////

bool do_tests       = true;
bool do_compression = true;
bool do_decompress  = true;

// //////////////////////////////////////////////////////

#ifdef WIN32
#define msvc_constexpr
#else
#define msvc_constexpr constexpr
#endif

// //////////////////////////////////////////////////////

inline unsigned constexpr min_bits(unsigned value)
{
    unsigned result = 0;

    while ((1u << result) <= value)
    {
        ++result;
    }

    return result;
}

// //////////////////////////////////////////////////////

namespace Coders
{
    // //////////////////////////////////////////////
    // Types
    // //////////////////////////////////////////////

    using Bytes = std::vector<uint8_t>;

    // //////////////////////////////////////////////
    // Coders
    // //////////////////////////////////////////////

    template<unsigned PROBABILITY_RANGE_BITS = 15>
    struct Fpaq0p_32bits
    {
        static const unsigned PROBABILITY_BITS  = PROBABILITY_RANGE_BITS;
        static const unsigned PROBABILITY_RANGE = 1 << PROBABILITY_RANGE_BITS;

        // w >= 2f + 1 comes from:
        // http://sachingarg.com/compression/entropy_coding/64bit/
        static_assert
        (
            PROBABILITY_RANGE_BITS < 16,
            "w >= 2f + 1. w = 32, PROBABILITY_RANGE_BITS maximum value is 15."
        );

        static_assert
        (
            PROBABILITY_RANGE_BITS > 1,
            "Comon, what's the point if PROBABILITY_RANGE_BITS < 2?"
        );

        class Encoder
        {
        public:
            Encoder(Bytes& bytes)
                : m_bytes(&bytes)
            {}

            void encode(unsigned bit, uint16_t one_probability)
            {
                assert(one_probability < PROBABILITY_RANGE);

                // cast to 64 bits for more precision.
                const uint32_t middle =
                    m_low
                    +
                    (
                        (
                            static_cast<uint64_t>(m_high - m_low)
                            * one_probability
                        )
                        >> PROBABILITY_BITS
                    );

                if (bit)
                {
                    m_high = middle;
                }
                else
                {
                    // +1 has something to do with squashing carry bits and
                    // underflow. I'm fuzzy with the maths.
                    m_low = middle + 1;
                }

                // flush while the top byte is all set for both low and high
                // i.e. normalise.
                while ((m_low ^ m_high) < (1u << 24))
                {
                    m_bytes->push_back(m_low >> 24);
                    m_low <<= 8;

                    // Some people ignore the 0xff at the end, can't remember
                    // why it was there. *shrug*
                    m_high = (m_high << 8) | 0xff;
                }
            }

            ~Encoder()
            {
                // Just copy Fabian's flush code.

                // Find shortest encoding that still decodes to the
                // right symbols. This is assume the decoder implicitly
                // zero-pads out of bounds array accesses.
                uint32_t round_up = 0xffffffu;
                while (round_up)
                {
                    if ((m_low | round_up) != ~0u)
                    {
                        const uint32_t rounded = (m_low + round_up) & ~round_up;

                        if (rounded <= m_high)
                        {
                            // inside interval, we're good!
                            m_low = rounded;
                            break;
                        }
                    }

                    round_up >>= 8;
                }

                while (m_low)
                {
                    m_bytes->push_back(m_low >> 24);
                    m_low <<= 8;
                }
            }

        private:
            Bytes*      m_bytes = 0;
            uint32_t    m_low   = 0;
            uint32_t    m_high  = 0xffffffff;
        };

        class Decoder
        {
        public:
            Decoder(const Bytes& bytes)
                : m_bytes(&bytes)
                , m_byte_count(bytes.size())
            {
                m_value = (m_value << 8) + read();
                m_value = (m_value << 8) + read();
                m_value = (m_value << 8) + read();
                m_value = (m_value << 8) + read();
            }

            unsigned decode(uint16_t one_probability)
            {
                assert(one_probability < (1 << PROBABILITY_RANGE_BITS));

                const uint32_t middle =
                    m_low
                    +
                    (
                        (
                            static_cast<uint64_t>(m_high - m_low)
                            * one_probability
                        )
                        >> PROBABILITY_RANGE_BITS
                    );

                auto bit = 0u;

                if (m_value <= middle)
                {
                    bit = 1;
                    m_high = middle;
                }
                else
                {
                    m_low = middle + 1;
                }

                // Renormaise
                while ((m_low ^ m_high) < (1u << 24))
                {
                    m_low <<= 8;
                    m_high = (m_high << 8) | 0xff;

                    m_value <<= 8;
                    m_value += read();
                }

                return bit;
            }

            uint32_t constexpr bytes_read() const
            {
                return m_read_index;
            }

        private:
            const Bytes*    m_bytes         = 0;
            uint32_t        m_value         = 0;
            uint32_t        m_low           = 0;
            uint32_t        m_high          = 0xffffffff;
            uint32_t        m_read_index    = 0;
            uint32_t        m_byte_count    = 0;

            uint8_t read()
            {
                if (m_read_index < m_byte_count)
                {
                    return (*m_bytes)[m_read_index++];
                }

                return 0;
            }
        };
    };

    // //////////////////////////////////////////////
    // Models
    // //////////////////////////////////////////////

    namespace Models
    {
        // Ideas initially from:
        // https://github.com/rygorous/gaffer_net/blob/master/main.cpp

        template<typename CODER>
        class Dual_exponential
        {
        public:
            typedef CODER                   Coder;
            typedef typename CODER::Encoder Encoder;
            typedef typename CODER::Decoder Decoder;

            Dual_exponential(
                    unsigned inertia_1 = 4,
                    unsigned inertia_2 = 4,
                    unsigned initial_probability_1 = QUARTER_RANGE,
                    unsigned initial_probability_2 = QUARTER_RANGE)
                : m_inertia_1(inertia_1)
                , m_inertia_2(inertia_2)
                , m_probabilities_1(initial_probability_1)
                , m_probabilities_2(initial_probability_2)
            {
                assert(m_inertia_1);
                assert(m_inertia_2);
            }

            void encode(Encoder& coder, unsigned value)
            {
                coder.encode(value, m_probabilities_1 + m_probabilities_2);
                adapt(value);
            }

            unsigned decode(Decoder& coder)
            {
                auto result =
                    coder.decode(m_probabilities_1 + m_probabilities_2);

                adapt(result);

                return result;
            }

        private:
            static const unsigned HALF_RANGE    = CODER::PROBABILITY_RANGE / 2;
            static const unsigned QUARTER_RANGE = HALF_RANGE / 2;

            unsigned m_inertia_1;
            unsigned m_inertia_2;
            uint16_t m_probabilities_1;
            uint16_t m_probabilities_2;

            void adapt(unsigned value)
            {
                if (value)
                {
                    m_probabilities_1 +=
                        (HALF_RANGE - m_probabilities_1) >> m_inertia_1;

                    m_probabilities_2 +=
                        (HALF_RANGE - m_probabilities_2) >> m_inertia_2;
                }
                else
                {
                    m_probabilities_1 -= m_probabilities_1 >> m_inertia_1;
                    m_probabilities_2 -= m_probabilities_2 >> m_inertia_2;
                }
            }
        };

        // //////////////////////////////////////////////

        template<typename CODER>
        class Bitstream
        {
        public:
            typedef CODER                   Coder;
            typedef typename CODER::Encoder Encoder;
            typedef typename CODER::Decoder Decoder;

            static const unsigned HALF_RANGE =  CODER::PROBABILITY_RANGE / 2;

            static void encode(Encoder& coder, unsigned value)
            {
                coder.encode(value, HALF_RANGE);
            }

            static unsigned decode(Decoder& coder)
            {
                auto result = coder.decode(HALF_RANGE);
                return result;
            }
        };

        // //////////////////////////////////////////////

        template
        <
            class BINARY_MODEL,
            unsigned BITS,
            unsigned MAX_VALUE = (1u << BITS) - 1
        >
        class Tree
        {
        public:
            typedef BINARY_MODEL Binary_model ;

            Tree() = default;

            // Forward the binary model constructor arguments.
            template<typename... Args>
            Tree(Args&&... args)
            {
                for(auto& model : m_models)
                {
                    model = BINARY_MODEL(std::forward<Args>(args)...);
                }
            }

            void encode(typename BINARY_MODEL::Encoder& coder, unsigned value)
            {
                assert(value < MODEL_COUNT);

                // Model the MSB first, then work our way down.
                // Seems adds are better than << 1.
                unsigned rebuilt = 1;
                unsigned mask = MAX_VALUE;

                while (rebuilt < MODEL_COUNT)
                {
                    unsigned bit = ((value & TOP_BIT) != 0);
                    value += value;

                    if (mask & TOP_BIT)
                    {
                        m_models[rebuilt - 1].encode(coder, bit);

                        // At any point we get a zero, then that means
                        // we are not constrained by max value anymore.
                        // RAM: TODO: algorithm is clunky - fix.
                        if (!bit)
                        {
                            mask = MODEL_COUNT - 1;
                        }
                    }

                    rebuilt += rebuilt + bit;
                    mask += mask;
                }
            }

            unsigned decode(typename BINARY_MODEL::Decoder& coder)
            {
                unsigned rebuilt    = 1;
                unsigned count      = MODEL_COUNT;
                unsigned mask       = MAX_VALUE;

                while (rebuilt < count)
                {
                    auto new_bit = 0u;

                    if (mask & TOP_BIT)
                    {
                        new_bit = m_models[rebuilt - 1].decode(coder);

                        if (!new_bit)
                        {
                            mask = MODEL_COUNT - 1;
                        }
                    }

                    rebuilt += rebuilt + new_bit;
                    mask    += mask;
                }

                // Clear the top bit due to starting rebuilt with 1.
                return rebuilt - count;
            }

        private:
            static const unsigned MODEL_COUNT   = 1 << BITS;
            static const unsigned TOP_BIT       = MODEL_COUNT / 2;

            std::array<BINARY_MODEL, MODEL_COUNT - 1> m_models;
        };

        // //////////////////////////////////////////////

        template
        <
            typename BINARY_MODEL,
            unsigned BITCOUNT_FOR_MAX_VALUE,
            unsigned MAX_VALUE=((1 << BITCOUNT_FOR_MAX_VALUE) - 1)
        >
        class Unsigned_golomb
        {
        public:
            Unsigned_golomb() = default;

            // Forward the binary model constructor arguments.
            template<typename... Args>
            Unsigned_golomb(Args&&... args)
            {
                m_bits =
                    Tree
                    <
                        BINARY_MODEL,
                        BITCOUNT_FOR_MAX_VALUE,
                        MAX_VALUE
                    >
                    (
                        std::forward<Args>(args)...
                    );
            }

            void encode(typename BINARY_MODEL::Encoder& coder, unsigned value)
            {
                // Enocde how many bits we are sending
                unsigned value_bits = min_bits(value);

                assert(value_bits <= MAX_VALUE);

                {
                    unsigned min_bits_2 = min_bits(value_bits);

                    assert(min_bits_2 <= BITCOUNT_FOR_MAX_VALUE);
                }

                m_bits.encode(coder, value_bits);

                if (value_bits)
                {
                    // No need to send the top bit
                    // as we can assume it's set due
                    // to value_bits.
                    --value_bits;

                    // Send the bits, no probabilities.
                    while (value_bits--)
                    {
                        const auto mask = (1 << value_bits);

                        Stream::encode
                        (
                            coder,
                            value & mask
                        );
                    }
                }
            }

            unsigned decode(typename BINARY_MODEL::Decoder& coder)
            {
                auto min_bits = m_bits.decode(coder);

                unsigned result = 0;

                if (min_bits)
                {
                    // top bit is always set.
                    result |= 1;

                    --min_bits;

                    while (min_bits)
                    {
                        result <<= 1;

                        result |= Stream::decode(coder);

                        --min_bits;
                    }
                }

                return result;
            }

        private:
            typedef Bitstream<typename BINARY_MODEL::Coder> Stream;

            Tree<BINARY_MODEL, BITCOUNT_FOR_MAX_VALUE, MAX_VALUE> m_bits;
        };
    }
}

void run_tests()
{
    using namespace Coders;
    using namespace Models;

    const auto binary_test_items =
    {
        0,1,0,1,0,1,1,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0
    };

    {
        Bytes data;
        const auto p = (Fpaq0p_32bits<>::PROBABILITY_RANGE * 2) / 3;

        {
            Fpaq0p_32bits<>::Encoder encoder(data);

            for (const unsigned t : binary_test_items)
            {
                encoder.encode(t, p);
            }
        }

        {
            Fpaq0p_32bits<>::Decoder decoder(data);

            for (const unsigned t : binary_test_items)
            {
                auto value = decoder.decode(p);

                assert(value == t);
            }

            auto read = decoder.bytes_read();
            assert(read == data.size());
        }
    }

    // Binary Model Tests
    {
        auto binary_test = [&binary_test_items](auto model_in, auto model_out)
        {
            Bytes data;

            {
                Fpaq0p_32bits<>::Encoder test_encoder(data);

                for (const auto t : binary_test_items)
                {
                    model_in.encode(test_encoder, t);
                }
            }
            {
                Fpaq0p_32bits<>::Decoder test_decoder(data);

                for (unsigned t : binary_test_items)
                {
                    auto value = model_out.decode(test_decoder);
                    assert(value == t);
                }

                auto read = test_decoder.bytes_read();
                assert(read == data.size());
            }
        };

        binary_test
        (
            Bitstream<Fpaq0p_32bits<>>(),
            Bitstream<Fpaq0p_32bits<>>()
        );
        binary_test
        (
            Dual_exponential<Fpaq0p_32bits<>>(5,2),
            Dual_exponential<Fpaq0p_32bits<>>(5,2)
        );
        binary_test
        (
            Dual_exponential<Fpaq0p_32bits<>>(6,1),
            Dual_exponential<Fpaq0p_32bits<>>(6,1)
        );
        binary_test
        (
            Dual_exponential<Fpaq0p_32bits<>>(3,4),
            Dual_exponential<Fpaq0p_32bits<>>(3,4)
        );
    }

    // Binary tree tests.
    {
        using Tree_model = Tree<Dual_exponential<Fpaq0p_32bits<>>, 3, 5>;
        Bytes data;
        Bytes tests
        {
            0,1,2,3,4,5,6,7, 0,1,2,3,4,5,6,7,
            0,1,2,3,4,5,6,7, 0,1,2
        };

        {
            Tree_model::Binary_model::Encoder test_encoder(data);
            Tree_model tree;

            for (auto t : tests)
            {
                tree.encode(test_encoder, t % 6);
            }
        }
        {
            Tree_model::Binary_model::Decoder test_decoder(data);
            Tree_model tree;

            for (unsigned t : tests)
            {
                auto value = tree.decode(test_decoder);
                assert(value == (t % 6));
            }

            auto read = test_decoder.bytes_read();
            assert(read == data.size());
        }
    }
}

// //////////////////////////////////////////////////////

struct Delta_data
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

// //////////////////////////////////////////////////////

inline constexpr bool quat_equal(const Delta_data& lhs, const Delta_data& rhs)
{
    return
    (
            (lhs.orientation_a          == rhs.orientation_a)
        &&  (lhs.orientation_b          == rhs.orientation_b)
        &&  (lhs.orientation_c          == rhs.orientation_c)
        &&  (lhs.orientation_largest    == rhs.orientation_largest)
    );
}

inline constexpr bool pos_equal(const Delta_data& lhs, const Delta_data& rhs)
{
    return
    (
            (lhs.position_x             == rhs.position_x)
        &&  (lhs.position_y             == rhs.position_y)
        &&  (lhs.position_z             == rhs.position_z)
    );
}

// //////////////////////////////////////////////////////

inline constexpr bool operator==(const Delta_data& lhs, const Delta_data& rhs)
{
    return (quat_equal(lhs, rhs) && pos_equal(lhs,rhs));
}

inline constexpr bool operator!=(const Delta_data& lhs, const Delta_data& rhs)
{
    return !operator==(lhs,rhs);
}

// //////////////////////////////////////////////////////

using ByteVector = std::vector<uint8_t>;

// //////////////////////////////////////////////////////

static const size_t     CUBES           = 901;
static const unsigned   PACKET_DELTA    = 6;

typedef std::array<Delta_data, CUBES> Frame;

// //////////////////////////////////////////////////////

inline bool constexpr operator==(const Frame& lhs, const Frame& rhs)
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

// //////////////////////////////////////////////////////

struct Gaffer
{
    unsigned orientation_largest;
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

typedef std::array<Predictors, CUBES> Frame_predicitons;

// //////////////////////////////////////////////////////

auto constexpr position(const Delta_data& lhs) -> Vec3i
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

Gaffer to_gaffer(const Delta_data& delta)
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

using namespace Coders;
using namespace Coders::Models;

// Found emperically.
// Changing restitution for cube 0 makes no difference.
static const constexpr int LOWEST_POINT         = 38;
static const constexpr int LOWEST_POINT_CUBE_0  = 367;
static const constexpr float RESTITUTION        = 0.869;
static const constexpr float DRAG               = 0.997;
static const constexpr float ANGULAR_CORRECTION = 0.999;
static const constexpr float DISTANCE_CUBE_0_SQ = 1742 * 1742;

struct Model
{
    using Coder         = Fpaq0p_32bits<13>;
    using Base_model    = Dual_exponential<Coder>;

    // Ok, lets just get coding first before simplification
    Base_model has_error                      = {1, 5};
    Base_model has_quat_largest               = {1, 4};
    std::array<Base_model, 4> interactive     =
    {{
        {1, 1},
        {1, 4},
        {1, 2},
        {1, 1}
    }};

    // If I get error, send both pos and quat errors.
    Tree<Base_model, 2> quat_largest = {5, 7};
    Tree<Base_model, 3> error_signs  = {5, 6};

    // This seems to do the trick.
    Unsigned_golomb<Base_model, 4, 10> error_bits =
    {
        2, 5
    };

    Unsigned_golomb<Base_model, 4, 10> error_bits_near_cube =
    {
        2, 4
    };
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

    const auto drag =
        ((std::abs(v[0]) > 0.0001f) || (std::abs(v[1]) > 0.0001f))
            ? (std::abs(v[2]) < 0.001f)
                ? true
                : false
            : false;

    if (drag)
    {
        v = mul(v, DRAG);
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

    if (!drag)
    {
        w = mul(w, ANGULAR_CORRECTION);
    }
    else
    {
        w = mul(w, ANGULAR_CORRECTION * ANGULAR_CORRECTION);
    }

    auto w_delta = mul(w, frame_delta);

    // Ok, need to convert to quat to do actual multiplications
    Quat r = to_quat(Rotor{w_delta});

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
    auto a = div(v_delta, frame_delta / 2.0f);

    auto angle_delta = Vec3f{0.0f,0.0f,0.0f};
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
                rotor.vec :
                rotor_neg.vec;
    }

    auto w = div(angle_delta, frame_delta);
    auto w_delta = sub(w, v_and_a.angular_velocity_per_frame);
    auto wa = div(w_delta, frame_delta / 2.0f);

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

// //////////////////////////////////////////////////////

auto encode
(
    const Frame& base,
    const Frame& target,
    Frame_predicitons& predicitons,
    unsigned frame_delta
)
-> Coders::Bytes
{
    auto                size = base.size();
    Coders::Bytes  data;

    // //////////////////////////////////////////////////////

    Model model;

    {
        Model::Coder::Encoder binary(data);

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
            model.has_error.encode(binary, has_error);

            if (has_error)
            {
                model.has_quat_largest.encode(binary, error_quat_largest);

                if (error_quat_largest)
                {
                    model.quat_largest.encode
                    (
                        binary,
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

                auto encode_error_bits = [&model, &binary](const Vec3i& vec)
                {
                    for (auto v: vec)
                    {
                        assert(v < (1 << 10));

                        model.error_bits.encode(binary,v);
                    }
                };
                auto encode_error_bits_near_cube =
                    [&model, &binary](const Vec3i& vec)
                {
                    for (auto v: vec)
                    {
                        assert(v < (1 << 10));

                        model.error_bits_near_cube.encode(binary,v);
                    }
                };

                auto previous_cube_0_distance =
                    sub(position(base[i]), position(base[0]));

                auto previous_cube_0_distance2 =
                    dot(previous_cube_0_distance, previous_cube_0_distance);

                if (previous_cube_0_distance2 < DISTANCE_CUBE_0_SQ)
                {
                    encode_error_bits_near_cube(vec_pos);
                    encode_error_bits_near_cube(vec_quat);
                }
                else
                {
                    encode_error_bits(vec_pos);
                    encode_error_bits(vec_quat);
                }

                if (vec_pos[0] || vec_pos[1] || vec_pos[2])
                {
                    model.error_signs.encode
                    (
                        binary,
                        signs_pos
                    );
                }

                if (vec_quat[0] || vec_quat[1] || vec_quat[2])
                {
                    model.error_signs.encode
                    (
                        binary,
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

            model.interactive[interact_lookup].encode
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
    const Coders::Bytes& data,
    unsigned frame_delta
)
-> Frame
{
    auto    size = base.size();
    Frame   target;

    Model model;

    {
        Model::Coder::Decoder binary(data);

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

            auto has_error = model.has_error.decode(binary);

            if (has_error)
            {
                auto has_error_quat_largest =
                    model.has_quat_largest.decode(binary);

                if (has_error_quat_largest)
                {
                    auto largest = model.quat_largest.decode(binary);

                    calculated_quat =
                        swap_orientation_largest(calculated_quat, largest);
                }

                // Get the errors and add them to things.
                auto decode_vec = [&model, &binary]() -> Vec3i
                {
                    Vec3i result;

                    for (auto& v: result)
                    {
                        v = model.error_bits.decode(binary);
                    }

                    return result;
                };
                auto decode_vec_near_cube = [&model, &binary]() -> Vec3i
                {
                    Vec3i result;

                    for (auto& v: result)
                    {
                        v = model.error_bits_near_cube.decode(binary);
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

                auto previous_cube_0_distance =
                    sub(position(base[i]), position(base[0]));

                auto previous_cube_0_distance2 =
                    dot(previous_cube_0_distance, previous_cube_0_distance);


                if (previous_cube_0_distance2 < DISTANCE_CUBE_0_SQ)
                {
                    vec_pos        = decode_vec_near_cube();
                    vec_quat       = decode_vec_near_cube();
                }
                else
                {
                    vec_pos        = decode_vec();
                    vec_quat       = decode_vec();
                }

                unsigned signs_pos  = 0;
                unsigned signs_quat = 0;

                if (vec_pos[0] || vec_pos[1] || vec_pos[2])
                {
                    signs_pos = model.error_signs.decode(binary);
                }

                if (vec_quat[0] || vec_quat[1] || vec_quat[2])
                {
                    signs_quat = model.error_signs.decode(binary);
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
                model.interactive[interact_lookup].decode
                (
                    binary
                );
        }
    }

    return target;
}

#define PRINT_INT(x) printf("%-32s\t%d\n", #x, x);
#define PRINT_LONG(x) printf("%-32s\t%ld\n", #x, x);
#define PRINT_FLOAT(x) printf("%-32s\t%f\n", #x, x);

void compress(std::vector<Frame>& frames)
{
    using namespace std::chrono;

    using Clock = high_resolution_clock;

    auto packets            = frames.size();
    unsigned bytes          = 0;
    unsigned packetsCoded   = 0;
    unsigned min            = 10000000;
    unsigned max            = 0;

    auto min_speed_encode = microseconds::max();
    auto max_speed_encode = microseconds::min();

    auto min_speed_decode = microseconds::max();
    auto max_speed_decode = microseconds::min();

    Frame_predicitons predict_server = {};
    Frame_predicitons predict_client = {};

    for (size_t i = PACKET_DELTA; i < packets; ++i)
    {
        auto encode_start = Clock::now();

        auto buffer = encode(
            frames[i-PACKET_DELTA],
            frames[i],
            predict_server,
            PACKET_DELTA);

        auto encode_time
            = duration_cast<microseconds>(Clock::now() - encode_start);

        min_speed_encode = std::min(min_speed_encode, encode_time);
        max_speed_encode = std::max(max_speed_encode, encode_time);

        const unsigned size = buffer.size();
        bytes += size;
        if (bytes)
        {
            min = std::min(min, size);
            max = std::max(max, size);
        }

        if (do_decompress)
        {
            auto decode_start = Clock::now();

            auto back = decode(
                frames[i-PACKET_DELTA],
                predict_client,
                buffer,
                PACKET_DELTA);

            assert(back == frames[i]);

            auto decode_time
                = duration_cast<microseconds>(Clock::now() - decode_start);

            min_speed_decode = std::min(min_speed_decode, decode_time);
            max_speed_decode = std::max(max_speed_decode, decode_time);
        }

        packetsCoded++;
    }

    float packetSizeAverge = ((float) bytes) / packetsCoded;
    float bytesPerSecondAverage = packetSizeAverge * 60.0f;
    float kbps = bytesPerSecondAverage * 8 / 1000.0f;

    printf("\n==============================================\n\n");

    PRINT_INT(packetsCoded);
    PRINT_INT(bytes);
    PRINT_FLOAT(bytesPerSecondAverage);
    PRINT_FLOAT(packetSizeAverge);

    PRINT_LONG(min_speed_encode.count());
    PRINT_LONG(max_speed_encode.count());
    PRINT_LONG(min_speed_decode.count());
    PRINT_LONG(max_speed_decode.count());

    PRINT_INT(min);
    PRINT_INT(max);
    PRINT_FLOAT(kbps);

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

        fclose(fileHandle);

        return frames;
    }();

    if (frames.empty())
    {
        return 1;
    }

    if (do_tests)
    {
        run_tests();
    }

    if (do_compression)
    {
        compress(frames);
    }

    return 0;
}
