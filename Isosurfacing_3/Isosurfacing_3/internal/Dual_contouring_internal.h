#ifndef CGAL_DUAL_CONTOURING_3_INTERNAL_DUAL_CONTOURING_3_H
#define CGAL_DUAL_CONTOURING_3_INTERNAL_DUAL_CONTOURING_3_H

#include <array>
#include <eigen3/Eigen/SVD>
#include <map>

#include "Tables.h"

namespace CGAL {
namespace Isosurfacing {
namespace internal {

namespace Positioning {
template <bool use_bbox = false>
class QEM_SVD {
public:
    /// <summary>
    /// Compute vertex position for Dual Contouring
    /// </summary>
    /// <typeparam name="Domain_"></typeparam>
    /// <param name="domain"></param>
    /// <param name="iso_value"></param>
    /// <param name="i"></param>
    /// <param name="j"></param>
    /// <param name="k"></param>
    /// <returns> true, if there is a point in the cell</returns>
    template <class Domain_>
    bool position(const Domain_& domain, const typename Domain_::FT iso_value, const typename Domain_::Cell_handle& vh,
                  typename Domain_::Point_3& point) const {
        typedef typename Domain_::Point_3 Point_3;
        typedef typename Domain_::Geom_traits::Vector_3 Vector_3;
        typedef typename Domain_::FT FT;

        typename Domain_::Cell_vertices vertices = domain.cell_vertices(vh);

        namespace Tables = internal::Cube_table;

        std::array<typename Domain_::FT, Tables::N_VERTICES> s;
        std::transform(vertices.begin(), vertices.end(), s.begin(), [&](const auto& v) { return domain.value(v); });

        std::array<bool, Tables::N_VERTICES> b;
        std::transform(s.begin(), s.end(), b.begin(), [iso_value](const auto& e) { return e <= iso_value; });

        unsigned int cubeindex = 0;
        // set bit if corresponding corner is below iso
        for (int i = 0; i < Tables::N_VERTICES; ++i) {
            cubeindex |= b[i] << i;
        }

        if (cubeindex == 0 || cubeindex == 255) {
            return false;
        }

        std::array<Vector_3, Tables::N_VERTICES> pos;
        std::transform(vertices.begin(), vertices.end(), pos.begin(),
                       [&](const auto& v) { return domain.position(v) - CGAL::ORIGIN; });

        point = CGAL::ORIGIN + (pos[0] + 0.5 * (pos[7] - pos[0]));  // set point to voxel center

        std::array<Vector_3, Tables::N_VERTICES> normals;
        std::transform(vertices.begin(), vertices.end(), normals.begin(),
                       [&](const auto& v) { return domain.gradient(v); });

        // compute edge intersections
        std::vector<Point_3> edge_intersections;
        std::vector<Vector_3> edge_intersection_normals;

        for (int i = 0; i < Tables::N_EDGES; ++i) {
            const auto& v0 = Tables::edge_to_vertex[i][0];
            const auto& v1 = Tables::edge_to_vertex[i][1];

            if (b[v0] != b[v1]) {  // e0
                const FT u = (s[v0] - iso_value) / (s[v0] - s[v1]);
                const Point_3 p_lerp = CGAL::ORIGIN + ((1 - u) * pos[v0] + u * pos[v1]);
                edge_intersections.push_back(p_lerp);
                const Vector_3 n_lerp = (1 - u) * normals[v0] + u * normals[v1];
                edge_intersection_normals.push_back(n_lerp);
            }
        }

        // MC Polygon Center of Mass
        if (false) {
            Vector_3 com_vec(0, 0, 0);

            for (int i = 0; i < edge_intersections.size(); ++i) {
                com_vec += edge_intersections[i] - CGAL::ORIGIN;
            }

            Point_3 p = CGAL::ORIGIN + com_vec / edge_intersections.size();
            point = p;
        }

        // SVD QEM
        if (true) {
            Eigen::Matrix3d A;
            A.setZero();
            Eigen::Vector3d b;
            b.setZero();
            for (std::size_t i = 0; i < edge_intersections.size(); ++i) {
                Eigen::Vector3d n_k = {edge_intersection_normals[i].x(), edge_intersection_normals[i].y(),
                                       edge_intersection_normals[i].z()};
                Eigen::Vector3d p_k = {edge_intersections[i].x(), edge_intersections[i].y(), edge_intersections[i].z()};
                double d_k = n_k.transpose() * p_k;

                Eigen::Matrix3d A_k = n_k * n_k.transpose();
                Eigen::Vector3d b_k = d_k * n_k;
                A += A_k;
                b += b_k;
            }

            Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
            // set threshold as in Peter Lindstrom's paper, "Out-of-Core
            // Simplification of Large Polygonal Models"
            svd.setThreshold(1e-3);

            // Init x hat
            Eigen::Vector3d x_hat;
            x_hat << point.x(), point.y(), point.z();

            // Lindstrom formula for QEM new position for singular matrices
            Eigen::VectorXd v_svd = x_hat + svd.solve(b - A * x_hat);

            point = Point_3(v_svd[0], v_svd[1], v_svd[2]);
        }

        // bbox
        if constexpr (use_bbox) {
            CGAL::Bbox_3 bbox = (CGAL::ORIGIN + pos[0]).bbox() + (CGAL::ORIGIN + pos[7]).bbox();

            FT x = std::min<FT>(std::max<FT>(point.x(), bbox.xmin()), bbox.xmax());
            FT y = std::min<FT>(std::max<FT>(point.y(), bbox.ymin()), bbox.ymax());
            FT z = std::min<FT>(std::max<FT>(point.z(), bbox.zmin()), bbox.zmax());
            point = Point_3(x, y, z);
        }

        return true;
    }
};

class Voxel_center {
public:
    /// <summary>
    /// Compute vertex position for Dual Contouring
    /// </summary>
    /// <typeparam name="Domain_"></typeparam>
    /// <param name="domain"></param>
    /// <param name="iso_value"></param>
    /// <param name="i"></param>
    /// <param name="j"></param>
    /// <param name="k"></param>
    /// <returns> true, if there is a point in the cell</returns>
    template <class Domain_>
    bool position(const Domain_& domain, const typename Domain_::FT iso_value, const typename Domain_::Cell_handle& vh,
                  typename Domain_::Point_3& point) const {
        typedef typename Domain_::Point_3 Point_3;
        typedef typename Domain_::Vector_3 Vector_3;

        namespace Tables = internal::Cube_table;

        std::array<typename Domain_::FT, Tables::N_VERTICES> s = domain.voxel_values(vh);

        std::array<bool, Tables::N_VERTICES> b;
        std::transform(s.begin(), s.end(), b.begin(), [iso_value](const auto& e) { return e <= iso_value; });

        unsigned int cubeindex = 0;
        // set bit if corresponding corner is below iso
        for (int i = 0; i < Tables::N_VERTICES; ++i) {
            cubeindex |= b[i] << i;
        }

        if (cubeindex == 0 || cubeindex == 255) {
            return false;
        }

        std::array<Point_3, Tables::N_VERTICES> p = domain.voxel_vertex_positions(vh);
        std::array<Vector_3, Tables::N_VERTICES> pos;
        std::transform(p.begin(), p.end(), pos.begin(), [](const auto& e) { return e - CGAL::ORIGIN; });

        point = CGAL::ORIGIN + (pos[0] + 0.5 * (pos[7] - pos[0]));  // set point to voxel center

        return true;
    }
};

class MC_polygon_center {
public:
    /// <summary>
    /// Compute vertex position for Dual Contouring
    /// </summary>
    /// <typeparam name="Domain_"></typeparam>
    /// <param name="domain"></param>
    /// <param name="iso_value"></param>
    /// <param name="i"></param>
    /// <param name="j"></param>
    /// <param name="k"></param>
    /// <returns> true, if there is a point in the cell</returns>
    template <class Domain_>
    bool position(const Domain_& domain, const typename Domain_::FT iso_value, const typename Domain_::Cell_handle& vh,
                  typename Domain_::Point_3& point) const {
        typedef typename Domain_::Point_3 Point_3;
        typedef typename Domain_::Vector_3 Vector_3;
        typedef typename Domain_::FT FT;

        namespace Tables = internal::Cube_table;

        std::array<typename Domain_::FT, Tables::N_VERTICES> s = domain.voxel_values(vh);

        std::array<bool, Tables::N_VERTICES> b;
        std::transform(s.begin(), s.end(), b.begin(), [iso_value](const auto& e) { return e <= iso_value; });

        unsigned int cubeindex = 0;
        // set bit if corresponding corner is below iso
        for (int i = 0; i < Tables::N_VERTICES; ++i) {
            cubeindex |= b[i] << i;
        }

        if (cubeindex == 0 || cubeindex == 255) {
            return false;
        }

        std::array<Point_3, Tables::N_VERTICES> p = domain.voxel_vertex_positions(vh);
        std::array<Vector_3, Tables::N_VERTICES> pos;
        std::transform(p.begin(), p.end(), pos.begin(), [](const auto& e) { return e - CGAL::ORIGIN; });

        point = CGAL::ORIGIN + (pos[0] + 0.5 * (pos[7] - pos[0]));  // set point to voxel center

        std::array<Vector_3, Tables::N_VERTICES> normals = domain.gradient(vh);

        // compute edge intersections
        std::vector<Point_3> edge_intersections;
        std::vector<Vector_3> edge_intersection_normals;

        for (int i = 0; i < Tables::N_EDGES; ++i) {
            const auto& v0 = Tables::edge_to_vertex[i][0];
            const auto& v1 = Tables::edge_to_vertex[i][1];

            if (b[v0] != b[v1]) {  // e0
                const FT u = (s[v0] - iso_value) / (s[v0] - s[v1]);
                const Point_3 p_lerp = CGAL::ORIGIN + ((1 - u) * pos[v0] + u * pos[v1]);
                edge_intersections.push_back(p_lerp);
                const Vector_3 n_lerp = (1 - u) * normals[v0] + u * normals[v1];
                edge_intersection_normals.push_back(n_lerp);
            }
        }

        // MC Polygon Center of Mass
        Vector_3 com_vec(0, 0, 0);

        for (int i = 0; i < edge_intersections.size(); ++i) {
            com_vec += edge_intersections[i] - CGAL::ORIGIN;
        }

        point = CGAL::ORIGIN + com_vec / edge_intersections.size();

        return true;
    }
};
}  // namespace Positioning

template <class Domain_, class Positioning_>
class Dual_contouring_position_functor {
private:
    typedef Domain_ Domain;
    typedef Positioning_ Positioning;

    typedef typename Domain::FT FT;
    typedef typename Domain::Point_3 Point_3;
    typedef typename Domain::Cell_handle Cell_handle;

public:
    Dual_contouring_position_functor(const Domain& domain, FT iso_value, const Positioning& positioning)
        : domain(domain), iso_value(iso_value), positioning(positioning), points_counter(0) {}

    void operator()(const Cell_handle& v) {
        // compute dc-vertices
        Point_3 p;
        if (positioning.position(domain, iso_value, v, p)) {
            map_voxel_to_point[v] = p;
            map_voxel_to_point_id[v] = points_counter++;
        }
    }

    // private:
    const Domain& domain;
    FT iso_value;
    const Positioning& positioning;

    std::map<Cell_handle, std::size_t> map_voxel_to_point_id;
    std::map<Cell_handle, Point_3> map_voxel_to_point;
    std::size_t points_counter;
};

template <class Domain_>
class Dual_contouring_quads_functor {
private:
    typedef Domain_ Domain;

    typedef typename Domain_::FT FT;
    typedef typename Domain_::Edge_handle Edge_handle;
    typedef typename Domain_::Cell_handle Cell_handle;

public:
    Dual_contouring_quads_functor(const Domain& domain, FT iso_value) : domain(domain), iso_value(iso_value) {}

    void operator()(const Edge_handle& e) {
        // save all quads
        const auto& vertices = domain.edge_vertices(e);
        const FT s0 = domain.value(vertices[0]);
        const FT s1 = domain.value(vertices[1]);

        if (s0 <= iso_value && s1 > iso_value) {
            const auto voxels = domain.cells_incident_to_edge(e);
            quads[e] = {voxels[0], voxels[1], voxels[2], voxels[3]};
        } else if (s1 <= iso_value && s0 > iso_value) {
            const auto voxels = domain.cells_incident_to_edge(e);
            quads[e] = {voxels[0], voxels[3], voxels[2], voxels[1]};
        }
    }

    // private:
    std::map<Edge_handle, std::array<Cell_handle, 4>> quads;

    const Domain& domain;
    FT iso_value;
};

}  // namespace internal
}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_DUAL_CONTOURING_3_INTERNAL_DUAL_CONTOURING_3_H