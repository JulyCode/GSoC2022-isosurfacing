#ifndef CGAL_MARCHING_TETS_3_INTERNAL_H
#define CGAL_MARCHING_TETS_3_INTERNAL_H

#include <array>
#include <map>
#include <mutex>

#include "Tables.h"
#include "Marching_cubes_3_internal.h"

namespace CGAL {
namespace Isosurfacing {
namespace internal {

template <class Domain_, class PointRange, class PolygonRange>
class Marching_tets_functor {
private:
    typedef Domain_ Domain;
    typedef PointRange Point_range;
    typedef PolygonRange Polygon_range;

    typedef typename Domain::FT FT;
    typedef typename Domain::Point Point;
    typedef typename Domain::Vertex_handle Vertex_handle;
    typedef typename Domain::Edge_handle Edge_handle;
    typedef typename Domain::Cell_handle Cell_handle;

public:
    Marching_tets_functor(const Domain& domain, const FT iso_value, Point_range& points, Polygon_range& polygons)
        : domain(domain), iso_value(iso_value), points(points), polygons(polygons) {}


    void operator()(const Cell_handle& cell) {

        FT values[4];
        Point corners[4];
        const int i_case = get_cell_corners(domain, cell, iso_value, corners, values);

        const int all_bits_set = (1 << (4 + 1)) - 1;  // last 4 bits are 1
        if (i_case == 0 || i_case == all_bits_set) {
            return;
        }

        std::array<Point, 6> vertices;
        
        for (int i = 0; i < 4; i++) {
            if (i_case & (1 << i)) {

            }
        }

        vertices[e_id] = vertex_interpolation(corners[v0], corners[v1], values[v0], values[v1], iso_value);

        std::lock_guard<std::mutex> lock(mutex);
    }

private:
    const Domain& domain;
    FT iso_value;

    Point_range& points;
    Polygon_range& polygons;

    // compute a unique global index for vertices
    // use as key the unique edge number
    std::map<Edge_handle, std::size_t> vertex_map;

    std::mutex mutex;
};

}  // namespace internal
}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_MARCHING_TETS_3_INTERNAL_H
