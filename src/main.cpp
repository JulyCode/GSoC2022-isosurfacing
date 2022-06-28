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

int main() {
    struct Sphere {

        const float resolution = 0.002f;
        const CGAL::Bbox_3 bbox = CGAL::Bbox_3(-1, -1, -1, 1, 1, 1);

        float operator()(float x, float y, float z) const {
            return std::sqrt(x * x + y * y + z * z);
        }

        std::size_t xdim() const { return bbox.x_span() / resolution; }
        std::size_t ydim() const { return bbox.y_span() / resolution; }
        std::size_t zdim() const { return bbox.z_span() / resolution; }

        Kernel::Vector_3 pos(std::size_t x, std::size_t y, std::size_t z) const {
            return Kernel::Vector_3(x * resolution + bbox.xmin(), y * resolution + bbox.ymin(), z * resolution + bbox.zmin());
        }

        float value(std::size_t x, std::size_t y, std::size_t z) const {
            Kernel::Vector_3 p = pos(x, y, z);
            return (*this)(p.x(), p.y(), p.z());
        }
    };

    Sphere sphere;

    //Mesh mesh;
    //CGAL::Isosurface_grid_3<Sphere, Mesh, CGAL::Marching_cubes_grid_3, Kernel> iso(sphere, sphere.bbox, sphere.resolution);
    CGAL::Isosurface_implicit_3<Sphere, Mesh, CGAL::Marching_cubes_implicit_3, Kernel> iso(sphere);

    Mesh mesh = iso.compute(0.8f);

    CGAL::IO::write_OFF("result.off", mesh);

    //write_off("result.off", result);
}
