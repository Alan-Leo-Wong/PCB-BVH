#ifndef BVH_PARALLEL_REINSERTION_OPTIMIZATION_HPP
#define BVH_PARALLEL_REINSERTION_OPTIMIZATION_HPP

#include <cassert>

#include "bvh/bvh.hpp"

namespace bvh {

template <typename Bvh>
class ParallelReinsertionOptimization {
public:
    ParallelReinsertionOptimization(Bvh* bvh)
        : bvh(bvh), parents(new size_t[bvh->node_count])
    {
        parents[0] = std::numeric_limits<size_t>::max();
        #pragma omp parallel for
        for (size_t i = 0; i < bvh->node_count; ++i) {
            if (bvh->nodes[i].is_leaf)
                continue;
            auto first_child = bvh->nodes[i].first_child_or_primitive;
            parents[first_child + 0] = i;
            parents[first_child + 1] = i;
        }
    }

private:
    Bvh* bvh;
    std::unique_ptr<size_t[]> parents;

    using Scalar    = typename Bvh::ScalarType;
    using Insertion = std::pair<size_t, Scalar>;

    Scalar cost() {
        Scalar cost(0);
        #pragma omp parallel for reduction(+: cost)
        for (size_t i = 0; i < bvh->node_count; ++i) {
            if (bvh->nodes[i].is_leaf)
                cost += bvh->nodes[i].bounding_box_proxy().half_area() * bvh->nodes[i].primitive_count;
            else
                cost += bvh->traversal_cost * bvh->nodes[i].bounding_box_proxy().half_area();
        }
        return cost;
    }

    void refit(size_t child) {
        auto bbox = bvh->nodes[child].bounding_box_proxy().to_bounding_box();
        while (child != 0) {
            auto parent = parents[child];
            bvh->nodes[parent].bounding_box_proxy() = bbox.extend(bvh->nodes[sibling(child)].bounding_box_proxy());
            child = parent;
        }
    }

    std::array<size_t, 6> conflicts(size_t in, size_t out) {
        auto parent_in = parents[in];
        return std::array<size_t, 6> {
            in,
            sibling(in),
            parent_in,
            parent_in == 0 ? in : parents[parent_in],
            out,
            out == 0 ? out : parents[out],
        };
    }

    void reinsert(size_t in, size_t out) {
        auto sibling_in   = sibling(in);
        auto parent_in    = parents[in];
        auto sibling_node = bvh->nodes[sibling_in];
        auto out_node     = bvh->nodes[out];

        // Re-insert it into the destination
        bvh->nodes[out].bounding_box_proxy().extend(bvh->nodes[in].bounding_box_proxy());
        bvh->nodes[out].first_child_or_primitive = std::min(in, sibling_in);
        bvh->nodes[out].is_leaf = false;
        bvh->nodes[sibling_in] = out_node;
        bvh->nodes[parent_in] = sibling_node;

        // Update parent-child indices
        if (!out_node.is_leaf) {
            parents[out_node.first_child_or_primitive + 0] = sibling_in;
            parents[out_node.first_child_or_primitive + 1] = sibling_in;
        }
        if (!sibling_node.is_leaf) {
            parents[sibling_node.first_child_or_primitive + 0] = parent_in;
            parents[sibling_node.first_child_or_primitive + 1] = parent_in;
        }
        parents[sibling_in] = out;
        parents[in] = out;
    }

    size_t sibling(size_t index) const {
        assert(index != 0);
        return index % 2 == 1 ? index + 1 : index - 1;
    }

    Insertion search(size_t in) {
        bool   down  = true;
        size_t pivot = parents[in];
        size_t out   = sibling(in);
        size_t out_best = out;

        auto bbox_in = bvh->nodes[in].bounding_box_proxy();
        auto bbox_parent = bvh->nodes[pivot].bounding_box_proxy();
        auto bbox_pivot = BoundingBox<Scalar>::empty();

        Scalar d = 0;
        Scalar d_best = 0;
        const Scalar d_bound = bbox_parent.half_area() - bbox_in.half_area();
        while (true) {
            auto bbox_out = bvh->nodes[out].bounding_box_proxy().to_bounding_box();
            auto bbox_merged = BoundingBox<Scalar>(bbox_in).extend(bbox_out);
            if (down) {
                auto d_direct = bbox_parent.half_area() - bbox_merged.half_area();
                if (d_best < d_direct + d) {
                    d_best = d_direct + d;
                    out_best = out;
                }
                d = d + bbox_out.half_area() - bbox_merged.half_area();
                if (bvh->nodes[out].is_leaf || d_bound + d <= d_best)
                    down = false;
                else
                    out = bvh->nodes[out].first_child_or_primitive;
            } else {
                d = d - bbox_out.half_area() + bbox_merged.half_area();
                if (pivot == parents[out]) {
                    bbox_pivot.extend(bbox_out);
                    out = pivot;
                    bbox_out = bvh->nodes[out].bounding_box_proxy();
                    if (out != parents[in]) {
                        bbox_merged = BoundingBox<Scalar>(bbox_in).extend(bbox_pivot);
                        auto d_direct = bbox_parent.half_area() - bbox_merged.half_area();
                        if (d_best < d_direct + d) {
                            d_best = d_direct + d;
                            out_best = out;
                        }
                        d = d + bbox_out.half_area() - bbox_pivot.half_area();
                    }
                    if (out == 0)
                        break;
                    out = sibling(pivot);
                    pivot = parents[out];
                    down = true;
                } else {
                    // If the node is the left sibling, go down
                    if (out % 2 == 1) {
                        down = true;
                        out = sibling(out);
                    } else {
                        out = parents[out];
                    }
                }
            }
        }

        if (in == out_best || sibling(in) == out_best || parents[in] == out_best)
            return Insertion { 0, 0 };
        return Insertion { out_best, d_best };
    }

public:
    void optimize(size_t u = 9, Scalar threshold = 0.1) {
        auto locks = std::make_unique<std::atomic<int64_t>[]>(bvh->node_count);
        auto outs  = std::make_unique<Insertion[]>(bvh->node_count);

        auto old_cost = cost();
        for (size_t iteration = 0; ; ++iteration) {
            size_t first_node = iteration % u + 1;

            #pragma omp parallel
            {
                // Clear the locks
                #pragma omp for nowait
                for (size_t i = 0; i < bvh->node_count; i++)
                    locks[i] = 0;

                // Search for insertion candidates
                #pragma omp for
                for (size_t i = first_node; i < bvh->node_count; i += u)
                    outs[i] = search(i);

                // Resolve topological conflicts with locking
                #pragma omp for
                for (size_t i = first_node; i < bvh->node_count; i += u) {
                    if (outs[i].second <= 0)
                        continue;
                    // Encode locks into 64bits using the highest 32 bits for the cost and
                    // the lowest 32 bits for the index of the node requesting the re-insertion
                    auto lock = (int64_t(as<int32_t>(float(outs[i].second))) << 32) | (int64_t(i) & INT64_C(0xFFFFFFFF));
                    // This takes advantage of the fact that IEEE-754 floats can be compared with regular integer comparisons
                    for (auto c : conflicts(i, outs[i].first))
                        atomic_max(locks[c], lock);
                }

                // Check the locks to disable conflicting re-insertions
                #pragma omp for
                for (size_t i = first_node; i < bvh->node_count; i += u) {
                    if (outs[i].second <= 0)
                        continue;
                    auto conflict_list = conflicts(i, outs[i].first);
                    // Make sure that this node owns all the locks for each and every conflicting node
                    bool is_conflict_free = std::all_of(conflict_list .begin(), conflict_list .end(), [&] (size_t j) {
                        return (locks[j] & INT64_C(0xFFFFFFFF)) == i;
                    });
                    if (!is_conflict_free)
                        outs[i] = Insertion { 0, 0 };
                }

                // Perform the reinsertions
                #pragma omp for
                for (size_t i = first_node; i < bvh->node_count; i += u) {
                    if (outs[i].second > 0)
                        reinsert(i, outs[i].first);
                }

                // Refit the nodes that have changed
                #pragma omp for
                for (size_t i = first_node; i < bvh->node_count; i += u) {
                    if (outs[i].second > 0) {
                        refit(i);
                        refit(outs[i].first);
                    }
                }
            }

            auto new_cost = cost();
            if (std::abs(new_cost - old_cost) <= threshold || iteration >= u) {
                if (u <= 1)
                    break;
                u = u - 1;
                iteration = 0;
            }
            old_cost = new_cost;
        }
    }
};

} // namespace bvh

#endif
