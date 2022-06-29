#include <iostream>
#include <fstream>
#include <math.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/IO/OFF.h>

#include "Isosurface_3.h"
#include "Marching_cubes_3.h"

typedef CGAL::Simple_cartesian<float> Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;

template<class Traits>
class Sphere_3 {
public:
    typedef typename Traits::FT FT;
    typedef typename Traits::Vector_3 Vector_3;

public:
    Sphere_3(const FT resolution, const FT radius)
        : resolution(resolution), bbox(-radius, -radius, -radius, radius, radius, radius),
          size_x(bbox.x_span() / resolution), size_y(bbox.y_span() / resolution), size_z(bbox.z_span() / resolution) {}

    FT operator()(FT x, FT y, FT z) const {
        return std::sqrt(x * x + y * y + z * z);
    }

    std::size_t xdim() const { return size_x; }
    std::size_t ydim() const { return size_y; }
    std::size_t zdim() const { return size_z; }

    Vector_3 pos(std::size_t x, std::size_t y, std::size_t z) const {
        return Vector_3(x * resolution + bbox.xmin(), y * resolution + bbox.ymin(), z * resolution + bbox.zmin());
    }

    float value(std::size_t x, std::size_t y, std::size_t z) const {
        const Vector_3 p = pos(x, y, z);
        return (*this)(p.x(), p.y(), p.z());
    }

public:
    const float resolution;
    const CGAL::Bbox_3 bbox;

private:
    const std::size_t size_x;
    const std::size_t size_y;
    const std::size_t size_z;
};

typedef Sphere_3<Kernel> Sphere;
typedef CGAL::Cartesian_grid_3<Kernel> Grid;

int main() {
    Sphere sphere(0.002f, 1);


    Grid grid(sphere.xdim(), sphere.ydim(), sphere.zdim(), sphere.bbox);

    for (std::size_t x = 0; x < grid.xdim(); x++) {
        const float x_pos = x * sphere.resolution + sphere.bbox.min(0);

        for (std::size_t y = 0; y < grid.ydim(); y++) {
            const float y_pos = y * sphere.resolution + sphere.bbox.min(1);

            for (std::size_t z = 0; z < grid.zdim(); z++) {
                const float z_pos = z * sphere.resolution + sphere.bbox.min(2);

                grid.value(x, y, z) = sphere(x_pos, y_pos, z_pos);
            }
        }
    }

    //Mesh mesh;
    //CGAL::Isosurface_grid_3<Sphere, Mesh, CGAL::Marching_cubes_grid_3, Kernel> iso(sphere, sphere.bbox, sphere.resolution);
    //iso.compute(0.8f, mesh);

    CGAL::Isosurface_implicit_3<Sphere, Mesh, CGAL::Marching_cubes_implicit_3, Kernel> iso(sphere);
    Mesh mesh = iso.compute(0.8f);

    CGAL::IO::write_OFF("result.off", mesh);

    //write_off("result.off", result);
}
