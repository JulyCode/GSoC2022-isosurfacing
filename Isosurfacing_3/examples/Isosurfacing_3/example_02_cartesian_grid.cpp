#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/IO/OFF.h>

#include "../../Cartesian_grid_3.h"
#include "../../Cartesian_grid_domain.h"
#include "../../Marching_cubes_3.h"

typedef CGAL::Simple_cartesian<float> Kernel;
typedef typename Kernel::FT FT;
typedef typename Kernel::Point_3 Point;

typedef CGAL::Cartesian_grid_3<Kernel> Grid;

typedef std::vector<Point> Point_range;
typedef std::vector<std::vector<std::size_t>> Polygon_range;

int main() {

    Grid grid(100, 100, 100, {-1, -1, -1, 1, 1, 1});

    for (std::size_t x = 0; x < grid.xdim(); x++) {
        for (std::size_t y = 0; y < grid.ydim(); y++) {
            for (std::size_t z = 0; z < grid.zdim(); z++) {

                const FT pos_x = x * grid.voxel_x() + grid.offset_x();
                const FT pos_y = y * grid.voxel_y() + grid.offset_y();
                const FT pos_z = z * grid.voxel_z() + grid.offset_z();

                grid.value(x, y, z) = std::sqrt(pos_x * pos_x + pos_y * pos_y + pos_z * pos_z);
            }
        }
    }

    CGAL::Isosurfacing::Cartesian_grid_domain<Kernel> domain(grid);

    Point_range points;
    Polygon_range polygons;

    CGAL::Isosurfacing::make_triangle_mesh_using_marching_cubes(domain, 0.8f, points, polygons);

    CGAL::IO::write_OFF("result.off", points, polygons);
}
