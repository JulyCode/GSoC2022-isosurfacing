#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <math.h>

#include "../../Implicit_domain.h"
#include "../../Marching_cubes_3.h"

typedef CGAL::Simple_cartesian<float> Kernel;
typedef typename Kernel::Vector_3 Vector;
typedef typename Kernel::Point_3 Point;

typedef std::vector<Point> Point_range;
typedef std::vector<std::vector<std::size_t>> Polygon_range;

int main() {
    auto sphere_function = [](const Point& point) {
        return std::sqrt(point.x() * point.x() + point.y() * point.y() + point.z() * point.z());
    };

    CGAL::Isosurfacing::Implicit_domain<Kernel, decltype(sphere_function)> domain(
        sphere_function, {-1, -1, -1, 1, 1, 1}, Vector(0.02f, 0.02f, 0.02f));

    Point_range points;
    Polygon_range polygons;

    CGAL::Isosurfacing::make_triangle_mesh_using_marching_cubes(domain, 0.8f, points, polygons);

    CGAL::IO::write_OFF("result.off", points, polygons);
}
