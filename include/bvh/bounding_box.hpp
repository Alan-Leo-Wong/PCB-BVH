#ifndef BVH_BOUNDING_BOX_HPP
#define BVH_BOUNDING_BOX_HPP

#include "bvh/vector.hpp"

namespace bvh {

/// A bounding box, represented with two extreme points.
template <typename Scalar>
struct BoundingBox {
    Vector3<Scalar> min, max;

    BoundingBox() = default;
    BoundingBox(const Vector3<Scalar>& v) : min(v), max(v) {}
    BoundingBox(const Vector3<Scalar>& min, const Vector3<Scalar>& max) : min(min), max(max) {}

    BoundingBox& extend(const BoundingBox& bbox) {
        min = bvh::min(min, bbox.min);
        max = bvh::max(max, bbox.max);
        return *this;
    }

    BoundingBox& extend(const Vector3<Scalar>& v) {
        min = bvh::min(min, v);
        max = bvh::max(max, v);
        return *this;
    }

    Vector3<Scalar> diagonal() const {
        return max - min;
    }

    Scalar half_area() const {
        Vector3<Scalar> d = diagonal();
        return (d[0] + d[1]) * d[2] + d[0] * d[1];
    }

    static BoundingBox full() {
        return BoundingBox(
            Vector3<Scalar>(-std::numeric_limits<Scalar>::max()),
            Vector3<Scalar>(std::numeric_limits<Scalar>::max()));
    }

    static BoundingBox empty() {
        return BoundingBox(
            Vector3<Scalar>(std::numeric_limits<Scalar>::max()),
            Vector3<Scalar>(-std::numeric_limits<Scalar>::max()));
    }
};

} // namespace bvh

#endif