#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <math.h>

#include <fstream>
#include <iostream>

#include "Cartesian_grid_3.h"
#include "Cartesian_grid_domain.h"
#include "Implicit_domain.h"
#include "Marching_cubes_3.h"
#include "Timer.h"

typedef CGAL::Simple_cartesian<float> Kernel;
typedef typename Kernel::Vector_3 Vector_3;
typedef typename Kernel::Point_3 Point_3;

typedef CGAL::Surface_mesh<Point_3> Mesh;
typedef CGAL::Cartesian_grid_3<Kernel> Grid;

typedef std::vector<Point_3> Point_range;
typedef std::vector<std::vector<std::size_t>> Polygon_range;

int main() {
    auto sphere_function = [](const Point_3& point) {
        return std::sqrt(point.x() * point.x() + point.y() * point.y() + point.z() * point.z());
    };
    CGAL::Implicit_domain<Kernel, decltype(sphere_function)> implicit_domain(sphere_function, {-1, -1, -1, 1, 1, 1},
                                                                             Vector_3(0.002f, 0.002f, 0.002f));


    // Grid grid(implicit_domain.size_x(), implicit_domain.size_y(), implicit_domain.size_z(), {-1, -1, -1, 1, 1, 1});

    // for (std::size_t x = 0; x < grid.xdim(); x++) {
    //     for (std::size_t y = 0; y < grid.ydim(); y++) {
    //         for (std::size_t z = 0; z < grid.zdim(); z++) {

    //             grid.value(x, y, z) = implicit_domain.value(x, y, z);
    //         }
    //     }
    // }

    const std::string fname = "../images/skull_2.9.inr";
    // Load image
    CGAL::Image_3 image;
    if (!image.read(fname)) {
        std::cerr << "Error: Cannot read file " << fname << std::endl;
        return EXIT_FAILURE;
    }
    Grid grid(image);

    CGAL::Cartesian_grid_domain<Kernel> grid_domain(grid);

    Point_range points;
    Polygon_range polygons;

    {
        ScopeTimer timer;
        CGAL::Isosurfacing::make_triangle_mesh_using_marching_cubes(implicit_domain, 0.8f, points, polygons);
    }

    // TODO: compare results with mesh_3

    Mesh mesh;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);

    CGAL::IO::write_OFF("result.off", mesh);
}
