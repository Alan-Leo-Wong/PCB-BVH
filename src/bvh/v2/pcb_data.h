#ifndef PCB_BVH_PCB_DATA_H
#define PCB_BVH_PCB_DATA_H

#include <bvh/v2/bbox.h>
#include <bvh/v2/ray.h>
#include <bvh/v2/vec.h>


#include <array>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <utility>
#include <vector>


namespace bvh::v2 {

#ifndef M_PI
    constexpr double M_PI = 3.1415926535897932384626433832795;
#endif

    template<typename T, size_t N>
    struct PCBData {
        /// type
        using Point = Vec<T, N>;

        /// data
        Point p0, p1;
        bool is_arc = false;

        /// Constructors
        PCBData() = default;

        PCBData(const Point &_p0, const Point &_p1) : p0(_p0), p1(_p1) {}

        /**
         *
         * @return
         */
        BVH_ALWAYS_INLINE std::pair<Point, Point> get_ed() const {
            return std::pair<Point, Point>{p0, p1};
        }

        /**
         *
         * @return
         */
        virtual BVH_ALWAYS_INLINE BBox<T, N> get_bbox() const = 0;

        /**
         *
         * @return
         */
        virtual BVH_ALWAYS_INLINE Point get_bbox_center() const = 0;

        /**
         *
         * @return
         */
        virtual BVH_ALWAYS_INLINE std::pair<T, Point>
        get_closest_dis(const Point &q) const = 0;

        /**
         * Check for intersection with a given box, which is used for collision
         * detection.
         * @param box
         * @return
         */
        virtual BVH_ALWAYS_INLINE bool is_intersect(const BBox<T, N> &bbox) const = 0;
    };

    template<typename T, size_t N>
    struct PCBSeg : public PCBData<T, N> {
        using typename PCBData<T, N>::Point;

        /// Constructors
        PCBSeg() = default;

        BVH_ALWAYS_INLINE PCBSeg(const Point &_p0, const Point &_p1)
                : PCBData<T, N>{_p0, _p1} {}

        /// Override functions
        BVH_ALWAYS_INLINE BBox<T, N> get_bbox() const override {
            return BBox<T, N>(this->p0).extend(this->p1);
        }

        BVH_ALWAYS_INLINE Point get_bbox_center() const override {
            return (this->p0 + this->p1) * static_cast<T>(0.5);
        }

        BVH_ALWAYS_INLINE std::pair<T, Point>
        get_closest_dis(const Point &q) const override {
            Vec<T, 2> p01 = this->p1 - this->p0;
            Vec<T, 2> p0q = q - this->p0;

            T c1 = dot(p0q, p01);
            if (c1 <= 0) {
                T dis = (T) dot(p0q, p0q);
                return std::pair{dis, this->p0}; // point is closest to p0
            }

            T c2 = dot(p01, p01);
            if (c2 <= c1) {
                T dis = (T) dot(q - this->p1, q - this->p1);
                return std::pair{dis, this->p1}; // point is closest to p1
            }

            T b = c1 / c2;
            Vec<T, 2> p_c =
                    this->p0 + b * p01; // Projection of point onto the line segment
            T dis = (T) dot(q - p_c, q - p_c); /// squared distance
            return std::pair{dis, p_c};       // Distance to projection
        }

        BVH_ALWAYS_INLINE bool is_intersect(const BBox<T, N> &bbox) const override {
            /*std::cout << "p1 = (" << this->p0[0] << ", " << this->p0[1] <<
                      ")" << std::endl;
            std::cout << "p2 = (" << this->p1[0] << ", " <<
                      this->p1[1] << ")" << std::endl;*/

            T min_t = 0.0, max_t = 1.0;
            Point _p0 = this->p0;
            Point _p1 = this->p1;
            /*std::cout << "seg_x_delta: " << _p1[0] - _p0[0] << std::endl;
            std::cout << "seg_y_delta: " << _p1[1] - _p0[1] << std::endl;*/

            if (bbox.min[0] == bbox.max[0]) {
                if (_p1[0] < bbox.min[0] || _p0[0] > bbox.max[0])
                    return false;
            } else {
                T bbox_x_min_delta = (bbox.min[0] - _p0[0]);
                T bbox_x_max_delta = (bbox.max[0] - _p0[0]);
                /*std::cout << "bbox_x_min_delta: " << bbox_x_min_delta << std::endl;
                std::cout << "bbox_x_max_delta: " << bbox_x_max_delta << std::endl;*/

                T p_x_delta = _p1[0] - _p0[0];
//                 std::cout << "p_x_delta: " << p_x_delta << std::endl;
                if (p_x_delta != .0) {
                    min_t = bbox_x_min_delta / p_x_delta;
                    max_t = bbox_x_max_delta / p_x_delta;
                    if (min_t > max_t)
                        std::swap(min_t, max_t);
                    /*std::cout << "min_t: " << min_t << std::endl;
                    std::cout << "max_t: " << max_t << std::endl;*/
                    if (min_t > 1.0 || max_t < 0.0) return false;
                }
            }

            {
                T bbox_y_min_delta = (bbox.min[1] - _p0[1]);
                T bbox_y_max_delta = (bbox.max[1] - _p0[1]);

                /*std::cout << "bbox_y_min_delta: " << bbox_y_min_delta << std::endl;
                std::cout << "bbox_y_max_delta: " << bbox_y_max_delta << std::endl;*/

                T p_y_delta = _p1[1] - _p0[1];
                if (p_y_delta >= 0) {
                    if (min_t * p_y_delta > bbox_y_max_delta ||
                        max_t * p_y_delta < bbox_y_min_delta)
                        return false;
                } else {
                    if (min_t * p_y_delta < bbox_y_min_delta ||
                        max_t * p_y_delta > bbox_y_max_delta)
                        return false;
                }
            }

            return true;
        }
    };

    template<typename T, size_t N>
    struct ArcData {
        /// data
        Vec<T, N> center;
        T radius = (T) NAN;
        T theta_0 = (T) NAN, theta_1 = (T) NAN, delta_theta = (T) NAN;

        /// Constructors
        ArcData() = default;

        ArcData(const Vec<T, N> &_center, T _radius, T _theta_0, T _theta_1)
                : center(_center), radius(_radius), theta_0(_theta_0), theta_1(_theta_1) {
            delta_theta = theta_1 - theta_0;
        }
    };

    template<typename T, size_t N>
    struct PCBArc : public PCBData<T, N> {
        using typename PCBData<T, N>::Point;

        static constexpr double INF_THETA = 1e9;

        /// data
        ArcData<T, N> arc_data;
        static constexpr std::array<T, 11> sp_angles = {
                -2 * M_PI, -1.5 * M_PI, -M_PI, -0.5 * M_PI, 0, 0.5 * M_PI,
                M_PI, 1.5 * M_PI, 2 * M_PI, 2.5 * M_PI, 3.0 * M_PI}; // to optimize

        static constexpr std::array<T, 11> cos_sp_angles = {
                1, 0, -1, 0, 1, 0, -1, 0, 1, 0, -1}; // to optimize

        static constexpr std::array<T, 11> sin_sp_angles = {
                0, 1, 0, -1, 0, 1, 0, -1, 0, 1, 0}; // to optimize

        /// Constructors
        PCBArc() = default;

        BVH_ALWAYS_INLINE PCBArc(const Point &center, const Point &_p0,
                                 const Point &_p1)
                : PCBData<T, N>{_p0, _p1} {
            T radius = std::sqrt(dot(center - _p0, center - _p0));

            T theta_0 = std::atan2(_p0[1] - center[1], _p0[0] - center[0]);
            T theta_1 = std::atan2(_p1[1] - center[1], _p1[0] - center[0]);

            if (theta_1 <= theta_0)
                theta_1 += 2 * M_PI;

            arc_data = ArcData(center, radius, theta_0, theta_1);
        }

        /// Override functions
        BVH_ALWAYS_INLINE BBox<T, N> get_bbox() const override {
            auto bbox = BBox<T, N>(this->p0).extend(this->p1);

            for (int i = 0; i < 11; ++i) {
                double angle = sp_angles[i];
                if (angle > arc_data.theta_1)
                    break;
                if (arc_data.theta_0 <= angle && angle <= arc_data.theta_1) {
                    Point point = {
                            (T) (arc_data.center[0] + arc_data.radius * std::cos(angle)),
                            (T) (arc_data.center[1] + arc_data.radius * std::sin(angle))};
                    bbox.extend(point);
                }
            }

            return bbox;
        }

        BVH_ALWAYS_INLINE Point get_bbox_center() const override {
            return (this->p0 + this->p1) * static_cast<T>(0.5);
        }

        BVH_ALWAYS_INLINE std::pair<T, Point>
        get_closest_dis(const Point &q) const override {
            Vec<T, 2> c_q = q - arc_data.center;
            T q_angle = std::atan2(c_q[1], c_q[0]);

            /*std::cout << get_bbox().min[0] << ", " << get_bbox().min[1] << std::endl;
            std::cout << get_bbox().max[0] << ", " << get_bbox().max[1] << std::endl;*/
            /*std::cout << "q_angle: " << q_angle / M_PI * 180.0 << std::endl;
            std::cout << "theta_0: " << arc_data.theta_0 / M_PI * 180.0 << std::endl;
            std::cout << "theta_1: " << arc_data.theta_1 / M_PI * 180.0 << std::endl;
            std::cout << "=============\n";
            system("pause");*/

            bool is_in_range_1 =
                    (q_angle >= arc_data.theta_0 && q_angle <= arc_data.theta_1);
            bool is_in_range_2 = (q_angle + 2 * M_PI >= arc_data.theta_0 &&
                                  q_angle + 2 * M_PI <= arc_data.theta_1);
            if (is_in_range_1 || is_in_range_2) {
                Point closest;
                T dis_to_c = length(c_q);
                /// squared distance
                T dis_to_arc =
                        (dis_to_c - arc_data.radius) * (dis_to_c - arc_data.radius);
                /*std::cout << "q: " << q[0] << ", " << q[1] << std::endl;
                std::cout << "center: " << arc_data.center[0] << ", " <<
                arc_data.center[1] << std::endl; std::cout << "dis_to_c: " << dis_to_c <<
                std::endl; std::cout << "radius: " << arc_data.radius << std::endl;
                std::cout << "dis_to_arc: " << dis_to_arc << std::endl;
                system("pause");*/
                if (is_in_range_1)
                    closest = get_radium_pos(q_angle);
                else
                    closest = get_radium_pos(q_angle + 2 * M_PI);
                return std::pair{dis_to_arc, closest};
            } else {
                /// squared distance
                T dis_to_p0 = dot(q - this->p0, q - this->p0);
                T dis_to_p1 = dot(q - this->p1, q - this->p1);

                if (dis_to_p0 < dis_to_p1)
                    return std::pair{dis_to_p0, this->p0};
                else
                    return std::pair{dis_to_p1, this->p1};
            }
        }

        BVH_ALWAYS_INLINE bool check_cover(double bbox_x_min_theta,
                                           double bbox_x_max_theta,
                                           double bbox_y_min_delta,
                                           double bbox_y_max_delta) const {
            T valid_min_theta = std::max(bbox_x_min_theta, (double)arc_data.theta_0);
            T valid_max_theta = std::min(bbox_x_max_theta, (double)arc_data.theta_1);
            /*std::cout << "==========\nvalid_x_min_theta = " << valid_min_theta / M_PI
                                                               * 180.0 << std::endl;
            std::cout << "valid_x_max_theta = " << valid_max_theta
                                                   / M_PI * 180.0 << std::endl;*/
            if (valid_min_theta > bbox_x_max_theta || valid_max_theta < bbox_x_min_theta) return false; // 没有任何区间重叠

            T valid_y_min_delta, valid_y_max_delta;
            valid_y_min_delta = std::min(std::sin(valid_min_theta), std::sin(valid_max_theta));
            valid_y_max_delta = std::max(std::sin(valid_min_theta), std::sin(valid_max_theta));
            for (int i = 0; i < 11; ++i) {
                T angle = sp_angles[i];
                if (angle <= valid_min_theta) continue;
                if (angle >= valid_max_theta) break;

                valid_y_min_delta = std::min(valid_y_min_delta, sin_sp_angles[i]);
                valid_y_max_delta = std::max(valid_y_max_delta, sin_sp_angles[i]);
            }

            if (valid_y_min_delta * arc_data.radius > bbox_y_max_delta ||
                valid_y_max_delta * arc_data.radius < bbox_y_min_delta)
                return false;
            return true;
        }

//           BVH_ALWAYS_INLINE bool is_intersect(const BBox<T, N> &bbox) const override {
//               T bbox_x_min_delta = (bbox.min[0] - arc_data.center[0]);
//               T bbox_x_max_delta = (bbox.max[0] - arc_data.center[0]);
//
//               T bbox_y_min_delta = (bbox.min[1] - arc_data.center[1]);
//               T bbox_y_max_delta = (bbox.max[1] - arc_data.center[1]);
//               /*std::cout << "center: " << arc_data.center[0] << ", " <<
//                         arc_data.center[1] << std::endl;
//
//               std::cout << "p1 = (" << this->p0[0] << ", " << this->p0[1] << ")" <<
//                         std::endl;
//               std::cout << "p2 = (" << this->p1[0] << ", " << this->p1[1] <<
//                         ")" << std::endl;
//               std::cout << "bbox_x_min_delta = " << bbox_x_min_delta <<
//                         std::endl;
//               std::cout << "bbox_x_max_delta = " << bbox_x_max_delta <<
//                         std::endl;
//
//               std::cout << "arc_data.theta_0 = " << arc_data.theta_0 / M_PI * 180.0 <<
//                         std::endl;
//               std::cout << "arc_data.theta_1 = " << arc_data.theta_1 / M_PI *
//                                                     180.0 << std::endl;
//               std::cout << "R = " << arc_data.radius << std::endl;*/
//
//               T R = arc_data.radius;
//               if (bbox_x_min_delta > R || bbox_y_min_delta > R || bbox_x_max_delta < -R ||
//                   bbox_y_max_delta < -R)
//                   return false;
//
//               T bbox_x_min_theta[4], bbox_x_max_theta[4];
//               if (bbox_x_min_delta >= -R && bbox_x_min_delta <= R) {
//                   bbox_x_min_theta[0] = std::acos(bbox_x_min_delta / R); // [0, pi]
//                   bbox_x_min_theta[1] = -bbox_x_min_theta[0];            // [-pi, 0]
//                   bbox_x_min_theta[2] = 2 * M_PI - bbox_x_min_theta[0];  // [pi, 2 * pi]
//                   bbox_x_min_theta[3] = 2 * M_PI + bbox_x_min_theta[0];  // [2 * pi, 3 * pi]
//               } else {
//                   bbox_x_min_theta[0] = std::min((double)arc_data.theta_0, .0);
//                   bbox_x_min_theta[1] = std::min((double)arc_data.theta_0, -M_PI);
//                   bbox_x_min_theta[2] = std::min((double)arc_data.theta_0, M_PI);
//                   bbox_x_min_theta[3] = std::min((double)arc_data.theta_0, 2 * M_PI);
//               }
//               if (bbox_x_max_delta >= -R && bbox_x_max_delta <= R) {
//                   bbox_x_max_theta[0] = std::acos(bbox_x_max_delta / R); // [0, pi]
//                   bbox_x_max_theta[1] = -bbox_x_max_theta[0];            // [-pi, 0]
//                   bbox_x_max_theta[2] = 2 * M_PI - bbox_x_max_theta[0];  // [pi, 2 * pi]
//                   bbox_x_max_theta[3] = 2 * M_PI + bbox_x_max_theta[0];  // [2 * pi, 3 * pi]
//               } else {
//                   bbox_x_max_theta[0] = std::max((double)arc_data.theta_1, M_PI);
//                   bbox_x_max_theta[1] = std::max((double)arc_data.theta_1, .0);
//                   bbox_x_max_theta[2] = std::max((double)arc_data.theta_1, 2 * M_PI);
//                   bbox_x_max_theta[3] = std::max((double)arc_data.theta_1, 3 * M_PI);
//               }
//
//               for (int i = 0; i < 4; ++i) {
//                   if (bbox_x_min_theta[i] > bbox_x_max_theta[i])
//                       std::swap(bbox_x_min_theta[i], bbox_x_max_theta[i]);
//                   /*std::cout << "==========\nbbox_x_min_theta = " << bbox_x_min_theta[i] /
//                                                                     M_PI * 180.0 << std::endl;
//                   std::cout << "bbox_x_max_theta = " <<
//                             bbox_x_max_theta[i] / M_PI * 180.0 << std::endl;
//                   system("pause");*/
//
//                   if (check_cover(bbox_x_min_theta[i], bbox_x_max_theta[i],
//                                   bbox_y_min_delta, bbox_y_max_delta))
//                       return true;
//               }
//
//               return false;
//           }

        BVH_ALWAYS_INLINE bool is_onArc(T theta) const {
            if ((theta >= arc_data.theta_0 && theta <= arc_data.theta_1) ||
                (theta - 2 * M_PI >= arc_data.theta_0 && theta - 2 * M_PI <= arc_data.theta_1) ||
                (theta + 2 * M_PI >= arc_data.theta_0 && theta + 2 * M_PI <= arc_data.theta_1))
                return true;
            return false;
        }

        BVH_ALWAYS_INLINE bool is_intersect(const BBox<T, N> &bbox) const override {
            T bbox_x_min_delta = (bbox.min[0] - arc_data.center[0]);
            T bbox_x_max_delta = (bbox.max[0] - arc_data.center[0]);
            T bbox_y_min_delta = (bbox.min[1] - arc_data.center[1]);
            T bbox_y_max_delta = (bbox.max[1] - arc_data.center[1]);
            T R = arc_data.radius;
            // 若bbox和整圆的包围盒无交集则直接return false
            // if (bbox_x_min_delta > R || bbox_y_min_delta > R ||
            //     bbox_x_max_delta < -R || bbox_y_max_delta < -R)
            //     return false;
            // 若bbox和圆弧的包围盒无交集则直接return false
            /*const auto arc_bbox = get_bbox();
            if (bbox.max[0] < arc_bbox.min[0] || bbox.max[1] < arc_bbox.min[1] ||
                bbox.min[0] > arc_bbox.max[0] || bbox.min[1] > arc_bbox.max[1])
                return false;*/
            // 若bbox完全包住整圆的包围盒则直接return true
            if (bbox_x_min_delta < -R && bbox_y_min_delta < -R && bbox_x_max_delta > R && bbox_y_max_delta > R)
                return true;
            // 思路：先求bbox四条边所在直线与整个圆的交点，再判断交点是否在包围盒及圆弧上
            bool has_arc_extend_intersect = false;  // 圆弧上是否存在与延长线的交点，若这种交点存在，圆弧却没有与bbox相交，则return false
            T temp_theta; // 暂存交点的弧度角，方便作对称
            T temp_coord; // 若求的是直线x=C与圆交点，则该变量表示交点的y坐标，以判断交点是否在bbox上
            // if (bbox_x_min_delta >= -R && bbox_x_min_delta <= R) {
            if (bbox_x_min_delta >= -R) { // bbox的左侧边是否在圆包围盒以内
                temp_theta = std::acos(bbox_x_min_delta / R);
                temp_coord = arc_data.center[1] + arc_data.radius * std::sin(temp_theta);
                // 第一个交点，若该交点在bbox左侧边所在线段(而非延长线)上
                if (bbox.min[1] <= temp_coord && temp_coord <= bbox.max[1]) {
                    if (is_onArc(temp_theta)) return true;
                }
                else if (is_onArc(temp_theta))  // 交点在延长线上
                    has_arc_extend_intersect = true;
                // 上一个交点在圆上关于x轴的对称点的y坐标
                if (bbox.min[1] <= -temp_coord + 2 * arc_data.center[1] &&
                    -temp_coord + 2 * arc_data.center[1] <= bbox.max[1]) {
                    if (is_onArc(-temp_theta)) return true;
                }
                else if (!has_arc_extend_intersect) if (is_onArc(-temp_theta)) has_arc_extend_intersect = true;
            }
            // if (bbox_x_max_delta >= -R && bbox_x_max_delta <= R) {
            if (bbox_x_max_delta <= R) { // bbox的右侧边
                temp_theta = std::acos(bbox_x_max_delta / R);
                temp_coord = arc_data.center[1] + arc_data.radius * std::sin(temp_theta);
                if (bbox.min[1] <= temp_coord && temp_coord <= bbox.max[1]) {
                    if (is_onArc(temp_theta)) return true;
                }
                else if (!has_arc_extend_intersect) if (is_onArc(temp_theta)) has_arc_extend_intersect = true;
                if (bbox.min[1] <= -temp_coord + 2 * arc_data.center[1] &&
                    -temp_coord + 2 * arc_data.center[1] <= bbox.max[1]) {
                    if (is_onArc(-temp_theta)) return true;
                }
                else if (!has_arc_extend_intersect) if (is_onArc(-temp_theta)) has_arc_extend_intersect = true;
            }
            // if (bbox_y_min_delta >= -R && bbox_y_min_delta <= R) {
            if (bbox_y_min_delta >= -R) { // bbox的下侧边
                temp_theta = M_PI / 2.0 + std::acos(bbox_y_min_delta / R);
                temp_coord = arc_data.center[0] + arc_data.radius * std::cos(temp_theta);
                if (bbox.min[0] <= temp_coord && temp_coord <= bbox.max[0]) {
                    if (is_onArc(temp_theta)) return true;
                }
                else if (!has_arc_extend_intersect) if (is_onArc(temp_theta)) has_arc_extend_intersect = true;
                if (bbox.min[0] <= -temp_coord + 2 * arc_data.center[0] &&
                    -temp_coord + 2 * arc_data.center[0] <= bbox.max[0]) {
                    if (is_onArc(M_PI - temp_theta)) return true;
                }
                else if (!has_arc_extend_intersect) if (is_onArc(M_PI - temp_theta)) has_arc_extend_intersect = true;
            }
            // if (bbox_y_max_delta >= -R && bbox_y_max_delta <= R) {
            if (bbox_y_max_delta <= R) { // bbox的上侧边
                temp_theta = M_PI / 2.0 + std::acos(bbox_y_max_delta / R);
                temp_coord = arc_data.center[0] + arc_data.radius * std::cos(temp_theta);
                if (bbox.min[0] <= temp_coord && temp_coord <= bbox.max[0]) {
                    if (is_onArc(temp_theta)) return true;
                }
                else if (!has_arc_extend_intersect) if (is_onArc(temp_theta)) has_arc_extend_intersect = true;
                if (bbox.min[0] <= -temp_coord + 2 * arc_data.center[0] &&
                    -temp_coord + 2 * arc_data.center[0] <= bbox.max[0]) {
                    if (is_onArc(M_PI - temp_theta)) return true;
                }
                else if (!has_arc_extend_intersect) if (is_onArc(M_PI - temp_theta)) has_arc_extend_intersect = true;
            }
            if (has_arc_extend_intersect) return false;
            return true;
        }
        //
        // BVH_ALWAYS_INLINE bool is_intersect3(const BBox<T, N> &bbox) const {
        //     // return true;
        //     T bbox_x_min_delta = (bbox.min[0] - arc_data.center[0]);
        //     T bbox_x_max_delta = (bbox.max[0] - arc_data.center[0]);
        //
        //     T bbox_y_min_delta = (bbox.min[1] - arc_data.center[1]);
        //     T bbox_y_max_delta = (bbox.max[1] - arc_data.center[1]);
        //     std::cout << "center: " << arc_data.center[0] << ", " << arc_data.center[1]
        //               << std::endl;
        //
        //     std::cout << "p1 = (" << this->p0[0] << ", " << this->p0[1] << ")"
        //               << std::endl;
        //     std::cout << "p2 = (" << this->p1[0] << ", " << this->p1[1] << ")"
        //               << std::endl;
        //     std::cout << "bbox_x_min_delta = " << bbox_x_min_delta << std::endl;
        //     std::cout << "bbox_x_max_delta = " << bbox_x_max_delta << std::endl;
        //
        //     std::cout << "arc_data.theta_0 = " << arc_data.theta_0 / M_PI * 180.0
        //               << std::endl;
        //     std::cout << "arc_data.theta_1 = " << arc_data.theta_1 / M_PI * 180.0
        //               << std::endl;
        //     std::cout << "R = " << arc_data.radius << std::endl;
        //
        //     T R = arc_data.radius;
        //     if (bbox_x_min_delta > R || bbox_y_min_delta > R || bbox_x_max_delta < -R ||
        //         bbox_y_max_delta < -R)
        //         return false;
        //
        //     // std::cout << "yes\n";
        //     T bbox_x_min_theta, bbox_x_max_theta;
        //     T bbox_y_min_theta, bbox_y_max_theta;
        //     // acos 是单调递减函数
        //     if (bbox_x_min_delta >= -R && bbox_x_min_delta <= R) {
        //         bbox_x_min_theta = std::acos(bbox_x_min_delta / R);
        //     } else {
        //         std::cout << "yes\n";
        //         bbox_x_min_theta = std::max(M_PI, arc_data.theta_1);
        //         //                bbox_x_min_theta = M_PI;
        //     }
        //     if (bbox_x_max_delta >= -R && bbox_x_max_delta <= R) {
        //         bbox_x_max_theta = std::acos(bbox_x_max_delta / R);
        //     } else {
        //         bbox_x_max_theta = std::min(.0, arc_data.theta_0);
        //         //                bbox_x_max_theta = .0;
        //     }
        //     if (bbox_x_min_theta > bbox_x_max_theta)
        //         std::swap(bbox_x_min_theta, bbox_x_max_theta);
        //
        //     std::cout << "==========\nbefore:\nbbox_x_min_theta = "
        //               << bbox_x_min_theta / M_PI * 180.0 << std::endl;
        //     std::cout << "bbox_x_max_theta = " << bbox_x_max_theta / M_PI * 180.0
        //               << std::endl;
        //
        //     // TODO: 直接将acos的结果替换成所有可能的区间情况([max(-pi, theta_0),
        //     // min(theta_1, 3*pi)])，然后逐一判断
        //     bool double_check = false;
        //     if (arc_data.theta_0 < 0 &&
        //         (arc_data.theta_1 < 0 || bbox_y_max_delta < 0)) {
        //         double t_bbox_min_theta = bbox_x_min_theta;
        //         bbox_x_min_theta = -bbox_x_max_theta;
        //         bbox_x_max_theta = -t_bbox_min_theta;
        //     } else if (arc_data.theta_1 >= M_PI ||
        //                std::abs(arc_data.theta_0 - M_PI) < 1e-5) { // TODO: optimize
        //         if (arc_data.theta_1 >= M_PI && arc_data.theta_1 < 2 * M_PI) {
        //             bbox_x_min_theta =
        //                     -bbox_x_min_theta + 2 * M_PI; // 0 < arc_data.theta_0 < PI
        //             bbox_x_max_theta =
        //                     -bbox_x_max_theta + 2 * M_PI; // 0 < arc_data.theta_0 < PI
        //         } else if (arc_data.theta_1 > 2 * M_PI) {
        //             bbox_x_min_theta += 2 * M_PI; // 0 < arc_data.theta_0 < PI
        //             bbox_x_max_theta += 2 * M_PI; // 0 < arc_data.theta_0 < PI
        //         }
        //     } else if (arc_data.theta_0 < 0 && arc_data.theta_1 < M_PI &&
        //                bbox_y_min_delta < 0) {
        //         double_check = true;
        //     }
        //     if (bbox_x_min_theta > bbox_x_max_theta)
        //         std::swap(bbox_x_min_theta, bbox_x_max_theta);
        //     std::cout << "==========\nafter:\nbbox_x_min_theta = "
        //               << bbox_x_min_theta / M_PI * 180.0 << std::endl;
        //     std::cout << "bbox_x_max_theta = " << bbox_x_max_theta / M_PI * 180.0
        //               << std::endl;
        //
        //     /*bool cover = check_cover(bbox_x_min_theta, bbox_x_max_theta,
        //                              bbox_x_min_theta, bbox_x_max_theta,
        //                              bbox_y_min_delta, bbox_y_max_delta);
        //     std::cout << std::boolalpha << "cover: " << cover << std::endl;
        //     std::cout << std::boolalpha << "double_check: " << double_check <<
        //     std::endl; if (!cover) { if (double_check) return
        //     check_cover(-bbox_x_max_theta, -bbox_x_min_theta, bbox_x_min_theta,
        //     bbox_x_max_theta, bbox_y_min_delta, bbox_y_max_delta); else return false;
        //     }*/
        //
        //     return true;
        // }

        /// Specialized functions
        BVH_ALWAYS_INLINE Point get_radium_pos(double theta_radium) const {
            T x = T(arc_data.center[0] + arc_data.radius * std::cos(theta_radium));
            T y = T(arc_data.center[1] + arc_data.radius * std::sin(theta_radium));
            return Point(x, y);
        }

        BVH_ALWAYS_INLINE std::vector<Point> adaptive_sample() {
            std::vector<Point> sample_pts;

            int base_samples = 10;
            float growth_rate = 0.5f;
            float angle_scale = 0.5f;
            int nb_pts = base_samples *
                         std::exp(growth_rate * angle_scale * (arc_data.delta_theta));
            double theta_sep = arc_data.delta_theta / ((nb_pts + 1) * 1.0);

            sample_pts.emplace_back(this->p0);
            for (int i = 1; i <= nb_pts; ++i) {
                double theta_radium = (arc_data.theta_0 + i * theta_sep);
                Point cur = get_radium_pos(theta_radium);
                sample_pts.emplace_back(cur);
            }
            sample_pts.emplace_back(this->p1);

            return sample_pts;
        }

        BVH_ALWAYS_INLINE Vec<T, N> side_proj(const Vec<T, N> &side_dir) const {
            T min_arc_proj = std::numeric_limits<T>::max();
            T max_arc_proj = std::numeric_limits<T>::lowest();

            T p0_dot = dot(this->p0, side_dir);
            min_arc_proj = std::min(min_arc_proj, p0_dot);
            max_arc_proj = std::max(max_arc_proj, p0_dot);

            T p1_dot = dot(this->p1, side_dir);
            min_arc_proj = std::min(min_arc_proj, p1_dot);
            max_arc_proj = std::max(max_arc_proj, p1_dot);

            for (int i = 0; i < 9; ++i) {
                double angle = sp_angles[i];
                if (angle <= arc_data.theta_0)
                    continue;
                if (angle >= arc_data.theta_1)
                    break;

                Point point = {
                        (T) (arc_data.center[0] + arc_data.radius * std::cos(angle)),
                        (T) (arc_data.center[1] + arc_data.radius * std::sin(angle))};
                T p_dot = dot(point, side_dir);
                min_arc_proj = std::min(min_arc_proj, p_dot);
                max_arc_proj = std::max(max_arc_proj, p_dot);
            }

            return {min_arc_proj, max_arc_proj};
        }
    };

} // namespace bvh::v2

#endif // PCB_BVH_PCB_DATA_H
