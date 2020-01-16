/// Distributed under the MIT license:
///
/// Copyright (C) 2020 Arsène Pérard-Gayot
///
/// Permission is hereby granted, free of charge, to any person obtaining a
/// copy of this software and associated documentation files (the "Software"), to
/// deal in the Software without restriction, including without limitation the
/// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in
/// all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
/// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
/// SOFTWARE.

#ifndef BVH_HPP
#define BVH_HPP

#include <cstdint>
#include <cstddef>
#include <cassert>
#include <cstring>
#include <cmath>
#include <climits>

#include <algorithm>
#include <numeric>
#include <optional>
#include <memory>
#include <limits>
#include <array>
#include <stack>

namespace bvh {

#ifdef BVH_DOUBLE
using Scalar = double;
#else
using Scalar = float;
#endif

/// A simple 3D vector class.
struct Vec3 {
    Scalar x, y, z;

    Vec3() = default;
    Vec3(Scalar s) : x(s), y(s), z(s) {}
    Vec3(Scalar x, Scalar y, Scalar z) : x(x), y(y), z(z) {}

    Vec3 operator - () const {
        return Vec3(-x, -y, -z);
    }

    Vec3 inverse() const {
        return Vec3(Scalar(1.0) / x, Scalar(1.0) / y, Scalar(1.0) / z);
    }

    Scalar& operator [] (size_t i) {
        return (&x)[i];
    }

    Scalar operator [] (size_t i) const {
        return (&x)[i];
    }
};

inline Vec3 operator + (const Vec3& a, const Vec3& b) {
    return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline Vec3 operator - (const Vec3& a, const Vec3& b) {
    return Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline Vec3 operator * (const Vec3& a, const Vec3& b) {
    return Vec3(a.x * b.x, a.y * b.y, a.z * b.z);
}

inline Vec3 operator * (const Vec3& a, Scalar s) {
    return Vec3(a.x * s, a.y * s, a.z * s);
}

inline Vec3 operator * (Scalar s, const Vec3& b) {
    return b * s;
}

inline Scalar dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

/// A ray, defined by an origin and a direction, with minimum and maximum distances along the direction from the origin.
struct Ray {
    Vec3 origin;
    Vec3 direction;
    Scalar tmin;
    Scalar tmax;

    Ray() = default;
    Ray(const Vec3& origin, const Vec3& direction, Scalar tmin = Scalar(0), Scalar tmax = std::numeric_limits<Scalar>::max())
        : origin(origin), direction(direction), tmin(tmin), tmax(tmax)
    {}

    unsigned classify() const {
        return (direction.x > 0 ? 1 : 0) |
               (direction.y > 0 ? 2 : 0) |
               (direction.z > 0 ? 4 : 0);
    }
};

/// A bounding box, represented with two extreme points.
struct BBox {
    Vec3 min, max;

    BBox() = default;
    BBox(const Vec3& v) : min(v), max(v) {}
    BBox(const Vec3& min, const Vec3& max) : min(min), max(max) {}

    BBox& extend(const BBox& bbox) {
        min.x = std::min(min.x, bbox.min.x);
        min.y = std::min(min.y, bbox.min.y);
        min.z = std::min(min.z, bbox.min.z);
        max.x = std::max(max.x, bbox.max.x);
        max.y = std::max(max.y, bbox.max.y);
        max.z = std::max(max.z, bbox.max.z);
        return *this;
    }

    BBox& extend(const Vec3& v) {
        return extend(BBox(v));
    }

    Vec3 diagonal() const {
        return max - min;
    }

    float half_area() const {
        Vec3 d = diagonal();
        return (d.x + d.y) * d.z + d.x * d.y;
    }

    static BBox full() {
        return BBox(Vec3(-std::numeric_limits<Scalar>::max()), Vec3(std::numeric_limits<Scalar>::max()));
    }

    static BBox empty() {
        return BBox(Vec3(std::numeric_limits<Scalar>::max()), Vec3(-std::numeric_limits<Scalar>::max()));
    }
};

#pragma omp declare reduction \
    (bbox_extend:BBox:omp_out.extend(omp_in)) \
    initializer(omp_priv = BBox::empty())

template <typename T, size_t Size>
struct TraversalStack {
    T elements[Size];
    size_t size = 0;

    void push(const T& t) { elements[size++] = t; }
    T pop() { return elements[--size]; }
    bool empty() const { return size == 0; }
};

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

inline Vec3 multiply_add(Vec3 x, Vec3 y, Vec3 z) {
    return Vec3(
        multiply_add(x.x, y.x, z.x),
        multiply_add(x.y, y.y, z.y),
        multiply_add(x.z, y.z, z.z)
    );
}

/// This structure represents a BVH with a list of nodes and primitives indices.
/// This API is very low level but offers full control of the algorithm, including
/// the number of bins used during building or the maximum tree depth.
/// For a higher-level API, use the Accel structure.
template <size_t BinCount = 32, size_t MaxDepth = 64, bool SweepAllAxes = true, size_t ParallelThreshold = 1024>
struct BVH {
    // The size of this structure is 32 bytes in single precision and 64 bytes in double precision.
    struct Node {
#ifdef BVH_DOUBLE
        using Size = uint64_t;
#else
        using Size = uint32_t;
#endif
        BBox bbox;
        bool leaf : 1;
        Size primitive_count : sizeof(Size) * CHAR_BIT - 1;
        Size first_child_or_primitive;

        std::pair<Scalar, Scalar> intersect(const Vec3& inverse_origin, const Vec3& inverse_direction, Scalar tmin, Scalar tmax, unsigned ix, unsigned iy, unsigned iz) const {
            auto values = reinterpret_cast<const float*>(this);
            Scalar entry_x = multiply_add(values[    ix], inverse_direction.x, inverse_origin.x);
            Scalar entry_y = multiply_add(values[    iy], inverse_direction.y, inverse_origin.y);
            Scalar entry_z = multiply_add(values[    iz], inverse_direction.z, inverse_origin.z);
            Scalar exit_x  = multiply_add(values[3 - ix], inverse_direction.x, inverse_origin.x);
            Scalar exit_y  = multiply_add(values[5 - iy], inverse_direction.y, inverse_origin.y);
            Scalar exit_z  = multiply_add(values[7 - iz], inverse_direction.z, inverse_origin.z);
            return std::make_pair(
                std::max(std::max(entry_x, entry_y), std::max(entry_z, tmin)),
                std::min(std::min(exit_x,  exit_y),  std::min(exit_z,  tmax))
            );
        }
    };

    struct Bin {
        BBox   bbox;
        size_t primitive_count;
        Scalar right_cost;
    };

    static constexpr size_t bin_count = BinCount;
    static constexpr size_t max_depth = MaxDepth;

    static constexpr size_t parallel_threshold = ParallelThreshold;

    struct BuildTask {
        struct WorkItem {
            size_t node_index;
            size_t begin;
            size_t end;
            size_t depth;

            WorkItem() = default;
            WorkItem(size_t node_index, size_t begin, size_t end, size_t depth)
                : node_index(node_index), begin(begin), end(end), depth(depth)
            {}

            size_t work_size() const { return end - begin; }
        };

        BVH* bvh;

        const BBox* bboxes;
        const Vec3* centers;

        BuildTask(BVH* bvh, const BBox* bboxes, const Vec3* centers)
            : bvh(bvh), bboxes(bboxes), centers(centers)
        {}

        std::optional<std::pair<WorkItem, WorkItem>> build(const WorkItem& item) {
            Node& node = bvh->nodes[item.node_index];
            auto make_leaf = [] (Node& node, size_t begin, size_t end) {
                node.first_child_or_primitive = begin;
                node.primitive_count          = end - begin;
                node.leaf                     = true;
            };

            if (item.work_size() <= 1 || item.depth >= max_depth) {
                make_leaf(node, item.begin, item.end);
                return std::nullopt;
            }

            auto primitive_indices = bvh->primitive_indices.get();

            // Compute the bounding box of the centers of the primitives in this node
            BBox center_bbox = BBox::empty();
            for (size_t i = item.begin; i < item.end; ++i)
                center_bbox.extend(centers[primitive_indices[i]]);

            int first_axis = 0, last_axis = 3;
            if (!SweepAllAxes) {
                // Compute largest axis
                Vec3 extents = node.bbox.diagonal();
                int axis = 0;
                if (extents.x < extents.y)
                    axis = 1;
                if (extents[axis] < extents.z)
                    axis = 2;
                first_axis = axis;
                last_axis  = axis + 1;
            }

            int    best_axis = -1;
            size_t best_split = 0;
            Scalar best_cost  = std::numeric_limits<Scalar>::max();
            std::array<Bin, bin_count> bins_per_axis[SweepAllAxes ? 3 : 1];

            auto inverse = Vec3(bin_count) * center_bbox.diagonal().inverse();
            auto base    = -center_bbox.min * inverse;
            auto bin_index = [=] (const Vec3& center, int axis) {
                return std::min(size_t(center[axis] * inverse[axis] + base[axis]), size_t(bin_count - 1));
            };

            for (int axis = first_axis; axis < last_axis; ++axis) {
                auto& bins = bins_per_axis[SweepAllAxes ? axis : 0];

                // Setup bins
                for (auto& bin : bins) {
                    bin.bbox = BBox::empty();
                    bin.primitive_count = 0;
                }

                // Fill bins
                for (size_t i = item.begin; i < item.end; ++i) {
                    size_t primitive_index = primitive_indices[i];
                    Bin& bin = bins[bin_index(centers[primitive_index], axis)];
                    bin.primitive_count++;
                    bin.bbox.extend(bboxes[primitive_index]);
                }

                // Right sweep to compute partial SAH
                BBox   current_bbox  = BBox::empty();
                size_t current_count = 0;
                for (size_t i = bin_count - 1; i > 0; --i) {
                    current_bbox.extend(bins[i].bbox);
                    current_count += bins[i].primitive_count;
                    bins[i].right_cost = current_bbox.half_area() * current_count;
                }

                // Left sweep to compute full cost and find minimum
                current_bbox  = BBox::empty();
                current_count = 0;
                for (size_t i = 0; i < bin_count - 1; ++i) {
                    current_bbox.extend(bins[i].bbox);
                    current_count += bins[i].primitive_count;
                    float cost = current_bbox.half_area() * current_count + bins[i + 1].right_cost;
                    if (cost < best_cost) {
                        best_axis  = axis;
                        best_split = i + 1;
                        best_cost  = cost;
                    }
                }
            }

            size_t total_primitives = item.end - item.begin;
            float  half_total_area  = node.bbox.half_area();

            // Check that the split is useful
            if (best_axis != -1 && best_cost + bvh->traversal_cost * half_total_area < total_primitives * half_total_area) {
                // Split primitives according to split position
                size_t begin_right = std::partition(primitive_indices + item.begin, primitive_indices + item.end, [&] (size_t i) {
                    return bin_index(centers[i], best_axis) < best_split;
                }) - primitive_indices;

                // Check that the split does not leave one side empty
                if (begin_right > item.begin && begin_right < item.end) {
                    // Allocate two nodes
                    size_t left_index;
                    #pragma omp atomic capture
                    { left_index = bvh->node_count; bvh->node_count += 2; }
                    Node& left  = bvh->nodes[left_index + 0];
                    Node& right = bvh->nodes[left_index + 1];
                    node.first_child_or_primitive = left_index;
                    node.primitive_count          = 0;
                    node.leaf                     = false;
                    
                    // Compute the bounding boxes of each node
                    auto& bins = bins_per_axis[SweepAllAxes ? best_axis : 0];
                    left.bbox = BBox::empty();
                    for (size_t i = 0; i < best_split; ++i)
                        left.bbox.extend(bins[i].bbox);
                    right.bbox = BBox::empty();
                    for (size_t i = best_split; i < bin_count; ++i)
                        right.bbox.extend(bins[i].bbox);

                    // Return new work items
                    WorkItem first_item (left_index + 0, item.begin, begin_right, item.depth + 1);
                    WorkItem second_item(left_index + 1, begin_right, item.end,   item.depth + 1);
                    return std::make_optional(std::make_pair(first_item, second_item));
                }
            }

            make_leaf(node, item.begin, item.end);
            return std::nullopt;
        }

        template <typename... Args>
        void run(Args&&... args) {
            std::stack<WorkItem> stack;
            stack.emplace(std::forward<Args&&>(args)...);
            while (!stack.empty()) {
                auto work_item = stack.top();
                stack.pop();

                auto more_work = build(work_item);
                if (more_work) {
                    auto [first_item, second_item] = *more_work;
                    if (first_item.work_size() > second_item.work_size())
                        std::swap(first_item, second_item);

                    stack.push(second_item);
                    if (first_item.work_size() > parallel_threshold) {
                        BuildTask task(*this);
                        #pragma omp task firstprivate(task)
                        { task.run(first_item); }
                    } else {
                        stack.push(first_item);
                    }
                }
            }
        }
    };

    void build(const BBox* bboxes, const Vec3* centers, size_t primitive_count) {
        // Allocate buffers
        nodes.reset(new Node[2 * primitive_count + 1]);
        primitive_indices.reset(new size_t[primitive_count]);
        std::iota(primitive_indices.get(), primitive_indices.get() + primitive_count, 0);

        // Initialize root node
        BBox& root_bbox = nodes[0].bbox;
        root_bbox = BBox::empty();
        node_count = 1;

        #pragma omp parallel
        {
            #pragma omp for reduction(bbox_extend: root_bbox)
            for (size_t i = 0; i < primitive_count; ++i)
                root_bbox.extend(bboxes[i]);
            #pragma omp single
            { BuildTask(this, bboxes, centers).run(0, 0, primitive_count, 0); }
        }
    }

    /// Intersects the BVH with the given ray and intersector.
    template <bool AnyHit, typename Intersector>
    std::optional<std::pair<size_t, typename Intersector::Result>> intersect(Ray ray, const Intersector& intersector) const {
        auto best_hit       = std::optional<std::pair<size_t, typename Intersector::Result>>(std::nullopt);
        auto intersect_leaf = [&] (size_t first_primitive, size_t primitive_count) {
            for (size_t i = 0; i < primitive_count; ++i) {
                auto hit = intersector(first_primitive + i, ray);
                if (hit) {
                    best_hit = std::make_optional(std::make_pair(primitive_indices[first_primitive + i], *hit));
                    if (AnyHit)
                        return true;
                    ray.tmax = hit->distance;
                }
            }
            return false;
        };

        // If the root is a leaf, intersect it and return
        if (nodes[0].leaf) {
            intersect_leaf(nodes[0].first_child_or_primitive, nodes[0].primitive_count);
            return best_hit;
        }

        // Precompute the inverse direction to avoid divisions and refactor
        // the computation to allow the use of FMA instructions (when available).
        auto inverse_direction = ray.direction.inverse();
        auto inverse_origin    = -ray.origin * inverse_direction;

        static constexpr size_t stack_size = max_depth + 3;

        // Indices into the node bounding box values are precomputed based on the ray octant
        auto ix = ray.direction.x > Scalar(0) ? 0 : 3;
        auto iy = ray.direction.y > Scalar(0) ? 1 : 4;
        auto iz = ray.direction.z > Scalar(0) ? 2 : 5;

        // This traversal loop is eager, because it immediately processes leaves instead of pushing them on the stack.
        // This is generally beneficial for performance because intersections will likely be found which will
        // allow to cull more subtrees with the ray-box test of the traversal loop.
        TraversalStack<size_t, stack_size> stack;
        stack.push(0);
        while (!stack.empty()) {
            auto& node = nodes[stack.pop()];
            auto first_child = node.first_child_or_primitive;

            auto& left  = nodes[first_child + 0];
            auto& right = nodes[first_child + 1];
            auto distance_left  = left .intersect(inverse_origin, inverse_direction, ray.tmin, ray.tmax, ix, iy, iz);
            auto distance_right = right.intersect(inverse_origin, inverse_direction, ray.tmin, ray.tmax, ix, iy, iz);
            bool hit_left  = distance_left.first  <= distance_left.second;
            bool hit_right = distance_right.first <= distance_right.second;

            if (hit_left && left.leaf) {
                if (intersect_leaf(left.first_child_or_primitive, left.primitive_count) && AnyHit)
                    break;
                hit_left = false;
            }

            if (hit_right && right.leaf) {
                if (intersect_leaf(right.first_child_or_primitive, right.primitive_count) && AnyHit)
                    break;
                hit_right = false;
            }

            if (hit_left && hit_right) {
                int order = distance_left.first < distance_right.first ? 0 : 1;
                stack.push(first_child + (1 - order));
                stack.push(first_child + order);
            } else if (hit_left ^ hit_right) {
                stack.push(first_child + (hit_left ? 0 : 1));
            }
        }

        return best_hit;
    }

    std::unique_ptr<Node[]>   nodes;
    std::unique_ptr<size_t[]> primitive_indices;
    size_t                    node_count = 0;
    float                     traversal_cost = 1.0f;
};

template <typename To, typename From>
To as(From from) {
    To to;
    std::memcpy(&to, &from, sizeof(from));
    return to;
}

inline float prodsign(float x, float y) {
    return as<float>(as<uint32_t>(x) ^ (as<uint32_t>(y) & UINT32_C(0x80000000)));
}

inline double prodsign(double x, double y) {
    return as<double>(as<uint64_t>(x) ^ (as<uint64_t>(y) & UINT64_C(0x8000000000000000)));
}

/// Triangle primitive, compatible with the Accel structure.
struct Triangle {
    Vec3 p0, e1, e2, n;

    Triangle() = default;
    Triangle(const Vec3& p0, const Vec3& p1, const Vec3& p2)
        : p0(p0), e1(p0 - p1), e2(p2 - p0)
    {
        n = cross(e1, e2);
    }

    Vec3 p1() const { return p0 - e1; }
    Vec3 p2() const { return p0 + e2; }

    BBox bounding_box() const {
        BBox bbox(p0);
        bbox.extend(p1());
        bbox.extend(p2());
        return bbox;
    }

    Vec3 center() const {
        return (p0 + p1() + p2()) * (Scalar(1.0) / Scalar(3.0));
    }

    struct Intersection {
        Scalar distance; ///< Distance along the ray.
        Scalar u;        ///< First barycentric coordinate.
        Scalar v;        ///< Second barycentric coordinate.

        Intersection() = default;
        Intersection(Scalar distance, Scalar u, Scalar v)
            : distance(distance), u(u), v(v)
        {}
    };

    std::optional<Intersection> intersect(const Ray& ray) const {
        static constexpr Scalar tolerance = Scalar(1e-9);

        auto c = p0 - ray.origin;
        auto r = cross(ray.direction, c);
        auto det = dot(n, ray.direction);
        auto abs_det = std::fabs(det);

        auto u = prodsign(dot(r, e2), det);
        auto v = prodsign(dot(r, e1), det);
        auto w = abs_det - u - v;

        if (u >= -tolerance && v >= -tolerance && w >= -tolerance) {
            auto t = prodsign(dot(n, c), det);
            if (t >= abs_det * ray.tmin && abs_det * ray.tmax > t) {
                auto inv_det = Scalar(1.0) / abs_det;
                return std::make_optional(Intersection(t * inv_det, u * inv_det, v * inv_det));
            }
        }

        return std::nullopt;
    }
};

/// High-level API entry point to create BVHs from primitive arrays.
/// This API is simpler but it may be wasteful of memory in certain cases.
/// For those cases, use the lower-level BVH structure.
template <typename Primitive>
struct Accel {
    Accel() = default;
    Accel(Primitive* primitives, size_t primitive_count)
        : primitives(primitives), primitive_count(primitive_count)
    {}

    /// Builds the acceleration structure. 
    void build() {
        assert(primitive_count > 0);
        std::unique_ptr<BBox[]> bboxes(new BBox[primitive_count]);
        std::unique_ptr<Vec3[]> centers(new Vec3[primitive_count]);
        #pragma omp parallel for
        for (size_t i = 0; i < primitive_count; ++i) {
            bboxes[i]  = primitives[i].bounding_box();
            centers[i] = primitives[i].center();
        }
        bvh.build(bboxes.get(), centers.get(), primitive_count);

        // Remap primitives to avoid indirect access through primitive indices
        std::unique_ptr<Primitive[]> primitives_copy(new Primitive[primitive_count]);
        std::copy(primitives, primitives + primitive_count, primitives_copy.get());
        for (size_t i = 0; i < primitive_count; ++i)
            primitives[i] = primitives_copy[bvh.primitive_indices[i]];
    }

    using Intersection = typename Primitive::Intersection;

    struct ClosestIntersector {
        using Result = Intersection;

        ClosestIntersector() = default;
        ClosestIntersector(const Accel* accel)
            : primitives(accel->primitives)
        {}

        std::optional<Intersection> operator () (size_t index, const Ray& ray) const {
            return primitives[index].intersect(ray);
        }

        const Primitive* primitives = nullptr;
    };

    /// Returns the closest intersection with a given ray.
    /// The intersection type is a pair made of the primitive index and the intersection
    /// value returned by the primitive, which contains at least a distance.
    std::optional<std::pair<size_t, Intersection>> intersect_closest(const Ray& ray) const {
        return bvh.intersect<false>(ray, ClosestIntersector(this));
    }

    struct AnyIntersector {
        struct Result {
            Scalar distance;

            Result() = default;
            Result(Scalar distance)
                : distance(distance)
            {}
        };

        AnyIntersector() = default;
        AnyIntersector(const Accel* accel)
            : primitives(accel->primitives)
        {}

        std::optional<Result> operator () (size_t index, const Ray& ray) const {
            auto hit = primitives[index].intersect(ray);
            if (hit)
                return std::make_optional(Result(hit->distance));
            return std::nullopt;
        }

        const Primitive* primitives = nullptr;
    };

    /// Returns true if there is any intersection with the given ray, otherwise false.
    bool intersect_any(const Ray& ray) const {
        return bvh.intersect<true>(ray, AnyIntersector(this)).has_value();
    }

    Primitive* primitives = nullptr;
    size_t primitive_count = 0;
    BVH<> bvh;
};

}

#endif