#include <bvh/v2/bvh.h>
#include <bvh/v2/stack.h>
#include <bvh/v2/executor.h>
#include <bvh/v2/pcb_data.h>
#include <bvh/v2/thread_pool.h>
#include <bvh/v2/default_builder.h>

#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <unordered_map>

typedef uint64_t index_t;
enum class ERROR_CODE : int32_t {
    SUCCESS = 0,                      // 成功
    ERROR_FILE_NOT_FOUND = -1,         // 文件未找到
    ERROR_INVALID_PARAMETER = -2,      // 无效参数
    ERROR_OUT_OF_MEMORY = -3,          // 内存不足
    ERROR_PERMISSION_DENIED = -4,      // 权限拒绝
    ERROR_IO_FAILURE = -5,             // IO 操作失败
    ERROR_TIMEOUT = -6,                // 超时
    ERROR_UNSUPPORTED_OPERATION = -7,  // 不支持的操作
    ERROR_NETWORK_UNAVAILABLE = -8,    // 网络不可用
    ERROR_RESOURCE_BUSY = -9,          // 资源正忙
    ERROR_DATA_CORRUPTION = -10,       // 数据损坏
    ERROR_OVERFLOW = -11,              // 数据溢出
    ERROR_UNDERFLOW = -12,             // 数据下溢
    ERROR_NOT_IMPLEMENTED = -13,       // 未实现的功能
    ERROR_UNKNOWN = -100,              // 未知错误

    // 警告代码（正数）
    WARNING_DEPRECATED = 1,            // 警告：功能已过时
    WARNING_UNEXPECTED_BEHAVIOR = 2    // 警告：非预期行为
};

using
enum ERROR_CODE;

using Scalar = double;
using Vec2 = bvh::v2::Vec<Scalar, 2>;
using BBox = bvh::v2::BBox<Scalar, 2>;
using Node = bvh::v2::Node<Scalar, 2>;
using Bvh = bvh::v2::Bvh<Node>;
using Ray = bvh::v2::Ray<Scalar, 2>;
using PCBData = bvh::v2::PCBData<Scalar, 2>;
using PCBSeg = bvh::v2::PCBSeg<Scalar, 2>;
using PCBArc = bvh::v2::PCBArc<Scalar, 2>;

class PCBScene {
    using Point = PCBData::Point;

private:
    /// data
    std::unordered_map<index_t, Point> P_coord; //
    std::unordered_map<index_t, Point> C_coord; //

    std::vector<std::unique_ptr<PCBData>> pcb_data; //

private:
    /**
     *
     * @param line
     * @return
     */
    BVH_ALWAYS_INLINE ERROR_CODE
    read_pcb_points(const std::string &line) {
        index_t left_paren_pos = line.find('(');
        index_t right_paren_pos = line.find(')');
        index_t comma_pos = line.find(',');
        if (left_paren_pos == std::string::npos ||
            right_paren_pos == std::string::npos ||
            comma_pos == std::string::npos)
            return ERROR_CODE::ERROR_IO_FAILURE;

        double x = std::stod(line.substr(left_paren_pos + 1, comma_pos - left_paren_pos - 1));
        double y = std::stod(line.substr(comma_pos + 1, right_paren_pos - comma_pos - 1));

        index_t equal_pos = line.find('=');
        if (line[0] == 'P') {
            index_t p_index = std::stoull(line.substr(1, equal_pos - 1));
            P_coord[p_index] = Point(x, y);
        } else {
            index_t c_index = std::stoull(line.substr(1, equal_pos - 1));
            C_coord[c_index] = Point(x, y);
        }

        return ERROR_CODE::SUCCESS;
    }

    /**
     *
     * @param token
     * @return
     */
    BVH_ALWAYS_INLINE ERROR_CODE
    read_pcb_segs(const std::vector<std::wstring> &token) {
        if (token.size() != 2) {
            std::cerr << "token.size() != 2\n";
            return ERROR_CODE::ERROR_IO_FAILURE;
        }

        index_t pos_0 = std::stoull(token[0].substr(1));
        index_t pos_1 = std::stoull(token[1].substr(1));

        if (!P_coord.count(pos_0) || !P_coord.count(pos_1))
            return ERROR_CODE::ERROR_IO_FAILURE;

        Point p0 = P_coord.at(pos_0);
        Point p1 = P_coord.at(pos_1);

        pcb_data.push_back(std::make_unique<PCBSeg>(p0, p1));

        return ERROR_CODE::SUCCESS;
    }

    /**
     *
     * @param token
     * @return
     */
    BVH_ALWAYS_INLINE ERROR_CODE
    read_pcb_arcs(const std::vector<std::wstring> &token) {
        if (token.size() != 3)
            return ERROR_CODE::ERROR_IO_FAILURE;

        index_t pos_c = std::stoull(token[0].substr(1));
        index_t pos_0 = std::stoull(token[1].substr(1));
        index_t pos_1 = std::stoull(token[2].substr(1));

        if (!C_coord.count(pos_c) || !P_coord.count(pos_0) || !P_coord.count(pos_1))
            return ERROR_CODE::ERROR_IO_FAILURE;

        Point center = C_coord.at(pos_c);
        Point p0 = P_coord.at(pos_0);
        Point p1 = P_coord.at(pos_1);

        pcb_data.push_back(std::make_unique<PCBArc>(center, p0, p1));
        pcb_data.back()->is_arc = true;

        return ERROR_CODE::SUCCESS;
    }

public:
    const std::vector<std::unique_ptr<PCBData>> &get_data() const { return pcb_data; }

    /**
     *
     * @param in_file
     * @return
     */
    BVH_ALWAYS_INLINE ERROR_CODE
    read_data(const std::string &in_file) {
        std::ifstream in(in_file);
        if (!in.is_open()) return ERROR_CODE::ERROR_IO_FAILURE;

        std::string line;

        // preprocessing
        std::string _seg_name = "线段";
        std::string _arc_name = "圆弧";
        std::wstringstream wss1;
        wss1 << _seg_name.c_str();
        std::wstring seg_name = wss1.str();
        std::wstringstream wss2;
        wss2 << _arc_name.c_str();
        std::wstring arc_name = wss2.str();

        while (getline(in, line)) {
            if (line[0] == 'P' || line[0] == 'C') { /// points
                if (read_pcb_points(line) != ERROR_CODE::SUCCESS)
                    return ERROR_CODE::ERROR_IO_FAILURE;
            } else {
                std::wstringstream wss;
                wss << line.c_str();
                std::wstring c_line = wss.str();

                index_t equal_pos = c_line.find(L'=');
                index_t left_paren_pos = c_line.find(L'(');
                index_t right_paren_pos = c_line.find(L')');
                if (equal_pos == std::wstring::npos ||
                    left_paren_pos == std::wstring::npos ||
                    right_paren_pos == std::wstring::npos)
                    return ERROR_CODE::ERROR_IO_FAILURE;

                std::wstring type = c_line.substr(equal_pos + 2, left_paren_pos - equal_pos - 2);
                std::wstring data = c_line.substr(left_paren_pos + 1, right_paren_pos - left_paren_pos - 1);

                index_t comma_pos = 0;
                std::vector<std::wstring> token;
                while ((comma_pos = data.find(L',')) != std::wstring::npos) {
                    token.emplace_back(data.substr(0, comma_pos));
                    data.erase(0, comma_pos + 1);
                }
                token.emplace_back(data);

                if (type.compare(seg_name) == 0) { /// segments
                    if (read_pcb_segs(token) != ERROR_CODE::SUCCESS)
                        return ERROR_CODE::ERROR_IO_FAILURE;
                } else if (type.compare(arc_name) == 0) { /// arcs
                    if (read_pcb_arcs(token) != ERROR_CODE::SUCCESS)
                        return ERROR_CODE::ERROR_IO_FAILURE;
                } else {
                    return ERROR_CODE::ERROR_IO_FAILURE;
                }
            }
        }

        in.close();
        return ERROR_CODE::SUCCESS;
    }

    /**
     *
     * @param out_file
     * @return
     */
    BVH_ALWAYS_INLINE ERROR_CODE
    output_scene(const std::string &out_file) {
        std::ofstream out(out_file);
        if (!out) return ERROR_CODE::ERROR_IO_FAILURE;
        out << std::setprecision(15);

        index_t cnt = 1;
        for (const auto &pri: pcb_data) {
            if (pri->is_arc) {
                auto arc = dynamic_cast<PCBArc *>(pri.get());

                std::vector<Point> sample_points = arc->adaptive_sample();
                assert(sample_points.size() > 1);

                for (index_t i = 0; i < sample_points.size(); ++i) {
                    const Point &p = sample_points.at(i);
                    out << "v " << p[0] << " " << p[1] << " 0 " << " 0.71 0.49 0.86" << std::endl;
                    if (i >= 1) {
                        out << "l " << cnt - 1 << " " << cnt << std::endl;
                    }
                    ++cnt;
                }
            } else {
                auto [p0, p1] = pri->get_ed();
                out << "v " << p0[0] << " " << p0[1] << " 0 " << " 0.53 0.81 0.98" << std::endl;
                out << "v " << p1[0] << " " << p1[1] << " 0 " << " 0.53 0.81 0.98" << std::endl;

                out << "l " << cnt << " " << cnt + 1 << std::endl;

                cnt += 2;
            }
        }

        return ERROR_CODE::SUCCESS;
    }
};

int main(int argc, char **argv) {
    PCBScene pcb;

    if (argc >= 2) pcb.read_data(argv[1]);
    else pcb.read_data(R"(D:\VSProjects\PCB-BVH\test\pcb_data\data5.txt)");

    if (argc >= 3) pcb.read_data(argv[2]);
    else pcb.output_scene(R"(D:\VSProjects\PCB-BVH\test\pcb_data\scene5.obj)");

    const auto &pcb_data = pcb.get_data();
    size_t num_pris = pcb_data.size();

    /*PCBSeg pcbSeg(Vec2(0, 0), Vec2(1, 0));
    auto [dis, p_c] = pcbSeg.get_closest_dis(Vec2(0.5, 0.5));
    std::cout << "dis: " << dis << "\np_c: " << p_c[0] << ", " << p_c[1] << std::endl;*/

    /*PCBArc pcbArc(Vec2(0, 0), Vec2(1, 0), Vec2(0, 1));
    auto [dis, p_c] = pcbArc.get_closest_dis(Vec2(1, 1));
    std::cout << "dis: " << dis << "\np_c: " << p_c[0] << ", " << p_c[1] << std::endl;*/

    // construct bvh
    bvh::v2::ThreadPool thread_pool;
    bvh::v2::ParallelExecutor executor(thread_pool);

    std::cout << "num_pris: " << num_pris << std::endl;
    std::vector<BBox> bboxes(num_pris);
    std::vector<Vec2> centers(num_pris);
    executor.for_each(0, num_pris, [&](size_t begin, size_t end) {
        for (size_t i = begin; i < end; ++i) {
            bboxes[i] = pcb_data[i]->get_bbox();
            centers[i] = pcb_data[i]->get_bbox_center();
        }
    });

    typename bvh::v2::DefaultBuilder<Node>::Config config;
    config.quality = bvh::v2::DefaultBuilder<Node>::Quality::High;
    auto bvh = bvh::v2::DefaultBuilder<Node>::build(thread_pool, bboxes, centers, config);

    // Permuting the primitive data allows to remove indirections during traversal, which makes it faster.
    static constexpr bool should_permute = true;

    static constexpr size_t invalid_id = std::numeric_limits<size_t>::max();
    static constexpr size_t stack_size = 64;
    static constexpr bool use_robust_traversal = false;
    bvh::v2::SmallStack<Bvh::Index, stack_size> stack;

    auto prim_id = invalid_id;

    Vec2 q = Vec2(1.0, -2.0);
    Scalar min_dis = std::numeric_limits<Scalar>::max();
    Vec2 closest_point;
    bvh.closest_point(q, bvh.get_root().index, stack,
                      [&](size_t begin, size_t end) {
                          for (size_t i = begin; i < end; ++i) {
                              size_t j = should_permute ? i : bvh.prim_ids[i];
                              auto res = pcb_data[j]->get_closest_dis(q);
                              if (min_dis > res.first) {
                                  prim_id = i;
                                  std::tie(min_dis, closest_point) = res;
                              }
                          }
                          return prim_id != invalid_id;
                      });

    if (prim_id != invalid_id) {
        std::cout
                << "Closest distance test\n"
                << "  primitive: " << prim_id << "\n"
                << "  minimum distance: " << min_dis << "\n"
                << "  closest_point: ("
                << closest_point[0] << ", " << closest_point[1] << ")" << std::endl;
        return 0;
    } else {
        std::cout << "No closest point found" << std::endl;
        return 1;
    }

    return 0;
}