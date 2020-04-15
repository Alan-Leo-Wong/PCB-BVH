#ifndef BVH_UTILITIES_HPP
#define BVH_UTILITIES_HPP

#include <cstring>
#include <cstdint>
#include <atomic>
#include <memory>
#include <algorithm>
#include <cmath>
#include <climits>

#include "bvh/bounding_box.hpp"

#ifdef _OPENMP
#include <omp.h>
#include <cassert>
#define bvh__get_num_threads()        omp_get_num_threads()
#define bvh__get_thread_num()         omp_get_thread_num()
#define bvh__assert_not_in_parallel() assert(omp_get_level() == 0);
#define bvh__assert_in_parallel()     assert(omp_get_level() > 0);
#else
#define bvh__get_num_threads()        1
#define bvh__get_thread_num()         0
#define bvh__assert_not_in_parallel() (void)0
#define bvh__assert_in_parallel()     (void)0
#endif

namespace bvh {

/// Safe function to reinterpret the bits of the given value as another type.
template <typename To, typename From>
To as(From from) {
    To to;
    std::memcpy(&to, &from, sizeof(from));
    return to;
}

/// Equivalent to copysign(x, x * y).
inline float product_sign(float x, float y) {
    return as<float>(as<uint32_t>(x) ^ (as<uint32_t>(y) & UINT32_C(0x80000000)));
}

/// Equivalent to copysign(x, x * y).
inline double product_sign(double x, double y) {
    return as<double>(as<uint64_t>(x) ^ (as<uint64_t>(y) & UINT64_C(0x8000000000000000)));
}

inline float multiply_add(float x, float y, float z) {
#ifdef FP_FAST_FMAF
    return std::fmaf(x, y, z);
#else
    return x * y + z;
#endif
}

inline double multiply_add(double x, double y, double z) {
#ifdef FP_FAST_FMA
    return std::fma(x, y, z);
#else
    return x * y + z;
#endif
}

template <typename Scalar>
void atomic_max(std::atomic<Scalar>& x, Scalar y) {
    auto z = x.load();
    while (z < y && !x.compare_exchange_weak(z, y)) ;
}

/// Shuffles primitives such that the primitive at index i is `primitives[indices[i]]`.
template <typename Primitive>
std::unique_ptr<Primitive[]> shuffle_primitives(const Primitive* primitives, const size_t* indices, size_t primitive_count) {
    auto primitives_copy = std::make_unique<Primitive[]>(primitive_count);
    for (size_t i = 0; i < primitive_count; ++i)
        primitives_copy[i] = primitives[indices[i]];
    return primitives_copy;
}

/// Computes the bounding box and the center of each primitive in given array.
template <typename Primitive, typename Scalar = typename Primitive::ScalarType>
std::pair<std::unique_ptr<BoundingBox<Scalar>[]>, std::unique_ptr<Vector3<Scalar>[]>>
compute_bounding_boxes_and_centers(const Primitive* primitives, size_t primitive_count) {
    auto bounding_boxes  = std::make_unique<BoundingBox<Scalar>[]>(primitive_count);
    auto centers         = std::make_unique<Vector3<Scalar>[]>(primitive_count);

    #pragma omp parallel for
    for (size_t i = 0; i < primitive_count; ++i) {
        bounding_boxes[i] = primitives[i].bounding_box();
        centers[i]        = primitives[i].center();
    }

    return std::make_pair(std::move(bounding_boxes), std::move(centers));
}

/// Computes the union of all the bounding boxes in the given array.
template <typename Scalar>
BoundingBox<Scalar> compute_bounding_boxes_union(const BoundingBox<Scalar>* bboxes, size_t count) {
    auto bbox = BoundingBox<Scalar>::empty();

    #pragma omp declare reduction \
        (bbox_extend:BoundingBox<Scalar>:omp_out.extend(omp_in)) \
        initializer(omp_priv = BoundingBox<Scalar>::empty())

    #pragma omp parallel for reduction(bbox_extend: bbox)
    for (size_t i = 0; i < count; ++i)
        bbox.extend(bboxes[i]);

    return bbox;
}

/// Templates that contains signed and unsigned integer types of the given number of bits.
template <size_t Bits>
struct SizedIntegerType {
    static_assert(Bits <= 8);
    using Signed   = int8_t;
    using Unsigned = uint8_t;
};

template <>
struct SizedIntegerType<64> {
    using Signed   = int64_t;
    using Unsigned = uint64_t;
};

template <>
struct SizedIntegerType<32> {
    using Signed   = int32_t;
    using Unsigned = uint32_t;
};

template <>
struct SizedIntegerType<16> {
    using Signed   = int16_t;
    using Unsigned = uint16_t;
};

/// Computes the (rounded-up) compile-time log in base-2 of an unsigned integer.
template <size_t P, size_t I = 0, bool Found = P == 0>
struct RoundUpLog2 {};

template <size_t P, size_t I>
struct RoundUpLog2<P, I, false> {
    static constexpr size_t value = RoundUpLog2<P, I + 1, (1 << (I + 1)) >= P>::value;
};

template <size_t P, size_t I>
struct RoundUpLog2<P, I, true> {
    static constexpr size_t value = I;
};

// Returns the number of bits that are equal to zero,
// starting from the most significant one.
template <typename T>
size_t count_leading_zeros(T value) {
    static constexpr size_t bit_count = sizeof(T) * CHAR_BIT;
    size_t a = 0;
    size_t b = bit_count;
    auto all = T(-1);
    for (size_t i = 0; i < RoundUpLog2<bit_count>::value; i++) {
        auto m = (a + b) / 2;
        auto mask = all << m;
        if (value & mask) a = m + 1;
        else              b = m;
    }
    return bit_count - b;
}

} // namespace bvh

#endif
