#include <bvh/v2/bvh.h>
#include <bvh/v2/vec.h>
#include <bvh/v2/ray.h>
#include <bvh/v2/node.h>
#include <bvh/v2/default_builder.h>
#include <bvh/v2/thread_pool.h>
#include <bvh/v2/executor.h>
#include <bvh/v2/stack.h>
#include <bvh/v2/tri.h>

#include <iostream>

using Scalar = float;
using Vec3 = bvh::v2::Vec<Scalar, 3>;
using BBox = bvh::v2::BBox<Scalar, 3>;
using Tri = bvh::v2::Tri<Scalar, 3>;
using Node = bvh::v2::Node<Scalar, 3>;
using Bvh = bvh::v2::Bvh<Node>;
using Ray = bvh::v2::Ray<Scalar, 3>;

using PrecomputedTri = bvh::v2::PrecomputedTri<Scalar, 3>;

int main() {
    // This is the original data, which may come in some other data type/structure.
    std::vector<Tri> tris;
    tris.emplace_back(
            Vec3(0.0, 0.0, 0.0),
            Vec3(1.0, 0.0, 0.0),
            Vec3(0.0, 0.0, 1.0)
    );
    tris.emplace_back(
            Vec3(2.0, 2.0, 2.0),
            Vec3(3.0, 3.0, 3.0),
            Vec3(3.0, 3.0, 1.0)
    );

    bvh::v2::ThreadPool thread_pool;
    bvh::v2::ParallelExecutor executor(thread_pool);

    // Get triangle centers and bounding boxes (required for BVH builder)
    std::vector<BBox> bboxes(tris.size());
    std::vector<Vec3> centers(tris.size());
    executor.for_each(0, tris.size(), [&](size_t begin, size_t end) {
        for (size_t i = begin; i < end; ++i) {
            bboxes[i] = tris[i].get_bbox();
            centers[i] = tris[i].get_center();
        }
    });

    typename bvh::v2::DefaultBuilder<Node>::Config config;
    config.quality = bvh::v2::DefaultBuilder<Node>::Quality::High;
    auto bvh = bvh::v2::DefaultBuilder<Node>::build(thread_pool, bboxes, centers, config);

    // Permuting the primitive data allows to remove indirections during traversal, which makes it faster.
    static constexpr bool should_permute = true;

    // This precomputes some data to speed up traversal further.
    std::vector<PrecomputedTri> precomputed_tris(tris.size());
    executor.for_each(0, tris.size(), [&](size_t begin, size_t end) {
        for (size_t i = begin; i < end; ++i) {
            auto j = should_permute ? bvh.prim_ids[i] : i;
            precomputed_tris[i] = tris[j];
        }
    });

    auto ray = Ray{
            Vec3(0., 0., 0.), // Ray origin
            Vec3(0., 0., 1.), // Ray direction
            0.,               // Minimum intersection distance
            100.              // Maximum intersection distance
    };

    static constexpr size_t invalid_id = std::numeric_limits<size_t>::max();
    static constexpr size_t stack_size = 64;
    static constexpr bool use_robust_traversal = false;

    auto prim_id = invalid_id;
    Scalar u, v;

    // Traverse the BVH and get the u, v coordinates of the closest intersection.
    bvh::v2::SmallStack<Bvh::Index, stack_size> stack;
    bvh.intersect<false, use_robust_traversal>(ray, bvh.get_root().index, stack,
                                               [&](size_t begin, size_t end) {
                                                   for (size_t i = begin; i < end; ++i) {
                                                       size_t j = should_permute ? i : bvh.prim_ids[i];
                                                       if (auto hit = precomputed_tris[j].intersect(ray)) {
                                                           prim_id = i;
                                                           std::tie(u, v) = *hit;
                                                       }
                                                   }
                                                   return prim_id != invalid_id;
                                               });

    Vec3 q = Vec3(0.2, 0, 0.5);
    Scalar min_sqr_dis = std::numeric_limits<Scalar>::max();
    Vec3 closest_point;
    bvh.closest_point(q, bvh.get_root().index, stack,
                      [&](size_t begin, size_t end) {
                          for (size_t i = begin; i < end; ++i) {
                              size_t j = should_permute ? i : bvh.prim_ids[i];
                              auto res = precomputed_tris[j].point_sqr_dis(q);
                              if (min_sqr_dis > res->first) {
                                  prim_id = i;
                                  std::tie(min_sqr_dis, closest_point) = *res;
                              }
                          }
                          return prim_id != invalid_id;
                      });

    if (prim_id != invalid_id) {
        std::cout
                << "Intersection found\n"
                << "  primitive: " << prim_id << "\n"
                << "  distance: " << ray.tmax << "\n"
                << "  barycentric coords.: " << u << ", " << v << std::endl;

        std::cout << "===================\n";

        std::cout
                << "Closest distance\n"
                << "  primitive: " << prim_id << "\n"
                << "  square distance: " << min_sqr_dis << "\n"
                << "  closest_point: ("
                << closest_point[0] << ", " << closest_point[1] << ", " << closest_point[2]
                << ")" << std::endl;
        return 0;
    } else {
        std::cout << "No intersection found" << std::endl;
        return 1;
    }
}