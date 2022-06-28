
#ifndef CGAL_ISOSURFACE_3_H
#define CGAL_ISOSURFACE_3_H

#include <CGAL/Bbox_3.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include "Cartesian_grid_3.h"

namespace CGAL {

template <class InputFunction, class OutputMesh, template<class, class, class> class Algorithm_, class Traits>
class Isosurface_grid_3 {
public:
    typedef typename Traits::FT FT;

    typedef InputFunction Input_function;
    typedef OutputMesh Output_mesh;

    typedef std::vector<typename Traits::Point_3> Point_range;
    typedef std::vector<std::vector<std::size_t>> Polygon_range;
    
    typedef Algorithm_<Point_range, Polygon_range, Traits> Algorithm;

public:
    Isosurface_grid_3(const Input_function& function, const Bbox_3& domain, const FT resolution)
        : grid(domain.x_span() / resolution, domain.y_span() / resolution, domain.z_span() / resolution) {

        CGAL_static_assertion((std::is_floating_point<FT>::value));

        // TODO: set resolution

        for (std::size_t x = 0; x < grid.xdim(); x++) {
            const FT x_pos = x * resolution + domain.min(0);

            for (std::size_t y = 0; y < grid.ydim(); y++) {
                const FT y_pos = y * resolution + domain.min(1);

                for (std::size_t z = 0; z < grid.zdim(); z++) {
                    const FT z_pos = z * resolution + domain.min(2);

                    grid.value(x, y, z) = function(x_pos, y_pos, z_pos);
                }
            }
        }
    }

    void compute(const FT iso_value, Output_mesh& mesh) {
        Point_range points;
        Polygon_range polygons;

        Algorithm algorithm(grid, iso_value, points, polygons);
        algorithm.compute();

        Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);
    }

private:
    Cartesian_grid_3<FT> grid;
};


template <class Oracle, class OutputMesh, template<class, class, class, class> class Algorithm_, class Traits>
class Isosurface_implicit_3 {
public:
    typedef typename Traits::FT FT;

    typedef OutputMesh Output_mesh;

    typedef std::vector<typename Traits::Point_3> Point_range;
    typedef std::vector<std::vector<std::size_t>> Polygon_range;
    
    typedef Algorithm_<Oracle, Point_range, Polygon_range, Traits> Algorithm;

public:
    Isosurface_implicit_3(const Oracle& oracle) : oracle(oracle) {

        CGAL_static_assertion((std::is_floating_point<FT>::value));
    }

    Output_mesh compute(const FT iso_value) {
        Point_range points;
        Polygon_range polygons;

        Algorithm algorithm(oracle, iso_value, points, polygons);
        algorithm.compute();

        Output_mesh mesh;
        Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);
        return mesh;
    }

private:
    const Oracle& oracle;
};

} // namespace CGAL

#endif // CGAL_ISOSURFACE_3_H