#ifndef BVH_V2_TOP_DOWN_SAH_BUILDER_H
#define BVH_V2_TOP_DOWN_SAH_BUILDER_H

#include "bvh/v2/bvh.h"
#include "bvh/v2/vec.h"
#include "bvh/v2/bbox.h"
#include "bvh/v2/split_heuristic.h"

#include <stack>
#include <span>
#include <algorithm>
#include <optional>
#include <numeric>
#include <cassert>

namespace bvh::v2 {

/// Base class for all SAH-based, top-down builders.
    template<typename Node>
    class TopDownSahBuilder {
    protected:
        using Scalar = typename Node::Scalar;
        using Vec = bvh::v2::Vec<Scalar, Node::dimension>;
        using BBox = bvh::v2::BBox<Scalar, Node::dimension>;

    public:
        struct Config {
            /// SAH heuristic parameters that control how primitives are partitioned.
            SplitHeuristic<Scalar> sah;

            /// Nodes containing less than this amount of primitives will not be split.
            /// This is mostly to speed up BVH construction, and using large values may lead to lower
            /// quality BVHs.
            size_t min_leaf_size = 1;

            /// Nodes that cannot be split based on the SAH and have a number of primitives larger than
            /// this will be split using a fallback strategy. This should not happen often, but may
            /// happen in worst-case scenarios or poorly designed scenes.
            size_t max_leaf_size = 8;
        };

    protected:
        struct WorkItem {
            size_t node_id;
            size_t begin; // prim begin
            size_t end;

            BVH_ALWAYS_INLINE size_t size() const { return end - begin; }
        };

        std::span<const BBox> bboxes_;
        std::span<const Vec> centers_;
        const Config &config_;

        BVH_ALWAYS_INLINE TopDownSahBuilder(
                std::span<const BBox> bboxes,
                std::span<const Vec> centers,
                const Config &config)
                : bboxes_(bboxes), centers_(centers), config_(config) {
            assert(bboxes.size() == centers.size());
            assert(config.min_leaf_size <= config.max_leaf_size);
        }

        virtual std::vector<size_t> &get_prim_ids() = 0;

        virtual std::optional<size_t> try_split(const BBox &bbox, size_t begin, size_t end) = 0;

        BVH_ALWAYS_INLINE const std::vector<size_t> &get_prim_ids() const {
            return const_cast<TopDownSahBuilder *>(this)->get_prim_ids();
        }

        /// 核心构建(自上向下)步骤
        Bvh<Node> build() {
            const auto prim_count = bboxes_.size();

            Bvh<Node> bvh;
            bvh.nodes.reserve((2 * prim_count) / config_.min_leaf_size);
            bvh.nodes.emplace_back();
            bvh.nodes.back().set_bbox(compute_bbox(0, prim_count));

            std::stack<WorkItem> stack;
            stack.push(WorkItem{0, 0, prim_count});
            while (!stack.empty()) {
                auto item = stack.top();
                stack.pop();

                auto &node = bvh.nodes[item.node_id];
//                if (item.node_id == 212)
                /*{
                    std::cout << std::boolalpha << (item.size() > config_.min_leaf_size) << std::endl;
                    std::cout << "bbox min: " << node.get_bbox().min[0] << ", "
                              << node.get_bbox().min[1] << std::endl;
                    std::cout << "bbox max: " << node.get_bbox().max[0] << ", "
                              << node.get_bbox().max[1] << std::endl;
                    std::cout << "item.size(): " << item.size() << std::endl;
                    std::cout << "item.begin: " << item.begin << std::endl;
                    std::cout << "item.end: " << item.end << std::endl;
                    std::cout << "==============\n";
                    system("pause");
                }*/
                if (item.size() > config_.min_leaf_size) {
                    /// 不同的bvh有不同的节点划分策略
                    auto split_pos1 = try_split(node.get_bbox(), item.begin, item.end);
                    /*if (split_pos1 == std::nullopt) {
                        std::cout << "my god" << std::endl;
                    }*/
                    if (auto split_pos = try_split(node.get_bbox(), item.begin, item.end)) {
                        auto first_child = bvh.nodes.size();
                        node.index = Node::Index::make_inner(first_child);

                        bvh.nodes.resize(first_child + 2);

//                        std::cout << "split_pos: " << *split_pos << std::endl;

                        auto first_bbox = compute_bbox(item.begin, *split_pos);
                        auto second_bbox = compute_bbox(*split_pos, item.end);
                        auto first_range = std::make_pair(item.begin, *split_pos);
                        auto second_range = std::make_pair(*split_pos, item.end);

                        /// 根据 sah 策略，始终让包裹 pri 数量最有可能多的 node 置为 left node
                        // For "any-hit" queries, the left child is chosen first, so we make sure that
                        // it is the child with the largest area, as it is more likely to contain an
                        // an occluder. See "SATO: Surface Area Traversal Order for Shadow Ray Tracing",
                        // by J. Nah and D. Manocha.
                        if (first_bbox.get_half_area() < second_bbox.get_half_area()) {
                            std::swap(first_bbox, second_bbox);
                            std::swap(first_range, second_range);
                        }

                        /*if (first_child == 212 || first_child == 211) {
                            std::cout << "yes first_child: " << first_child << std::endl;
                        }*/

                        auto first_item = WorkItem{first_child + 0, first_range.first, first_range.second};
                        auto second_item = WorkItem{first_child + 1, second_range.first, second_range.second};
                        bvh.nodes[first_child + 0].set_bbox(first_bbox);
                        bvh.nodes[first_child + 1].set_bbox(second_bbox);

                        // Process the largest child item first, in order to minimize the stack size.
                        if (first_item.size() < second_item.size())
                            std::swap(first_item, second_item);

                        /*{
                            std::cout << "first_child: " << first_child << std::endl;
                            std::cout << "first_item: " << first_item.begin << ", " << first_item.end << std::endl;
                            std::cout << "second_item: " << second_item.begin << ", " << second_item.end << std::endl;
                            std::cout << "==============\n";
                            system("pause");
                        }*/

                        stack.push(first_item);
                        stack.push(second_item);
                        continue;
                    }
                }

//                std::cout << "item.node_id: " << item.node_id << std::endl;
                node.index = Node::Index::make_leaf(item.begin, item.size());
                /*if (item.node_id == 212) {
                    std::cout << node.index.first_id() << std::endl;
                }*/
            }

            bvh.prim_ids = std::move(get_prim_ids());
            bvh.nodes.shrink_to_fit();
            return bvh;
        }

        BVH_ALWAYS_INLINE BBox compute_bbox(size_t begin, size_t end) const {
            const auto &prim_ids = get_prim_ids();
            auto bbox = BBox::make_empty();
            for (size_t i = begin; i < end; ++i) {
//                std::cout << "i: " << i << ", prim_ids[i]: " << prim_ids[i] << std::endl;
                bbox.extend(bboxes_[prim_ids[i]]);
            }
            return bbox;
        }
    };

} // namespace bvh::v2

#endif
