#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <math.h>

#include <functional>
#include <iostream>
#include <limits>
#include <unordered_map>

#include "../Cartesian_grid_3.h"
#include "../Cartesian_grid_domain.h"
#include "../Implicit_domain.h"
#include "../Marching_cubes_3.h"
#include "../Timer.h"

typedef CGAL::Simple_cartesian<double> Kernel;
typedef typename Kernel::FT FT;
typedef typename Kernel::Vector_3 Vector;
typedef typename Kernel::Point_3 Point;

typedef CGAL::Cartesian_grid_3<Kernel> Grid;

typedef std::vector<Point> Point_range;
typedef std::vector<std::vector<std::size_t>> Polygon_range;


template<typename Function, class Domain_>
int64_t run_func(Function func, const Domain_& domain, const typename Domain_::FT iso) {
    Point_range points;
    Polygon_range polygons;

    ScopeTimer timer;
    func(domain, iso, points, polygons);

    const int64_t ms = timer.stop();

    if (points.size() > std::numeric_limits<std::size_t>::max() - 2) {
        std::cout << "This should never print and only prevents optimizations" << std::endl;
    }
    return ms;
}

template<class Point_>
struct SphereFunction {
    FT operator()(const Point_& point) const {
        return std::sqrt(point.x() * point.x() + point.y() * point.y() + point.z() * point.z());
    }
};

template <class GeomTraits>
    CGAL::Isosurfacing::Implicit_domain < GeomTraits,
    SphereFunction<GeomTraits>> implicit_sphere(const std::size_t N) {

    const FT resolution = 2.0 / N;
    auto domain = CGAL::Isosurfacing::create_implicit_domain(SphereFunction < GeomTraits>(), {-1, -1, -1, 1, 1, 1},
                                                             Vector(resolution, resolution, resolution));
    return domain;
}

    template <class GeomTraits>
CGAL::Isosurfacing::Cartesian_grid_domain<GeomTraits> grid_sphere(const std::size_t N) {

    const FT resolution = 2.0 / N;

    Grid grid(N, N, N, {-1, -1, -1, 1, 1, 1});

#pragma omp parallel for
    for (std::size_t x = 0; x < grid.xdim(); x++) {
        const FT xp = x * resolution - 1.0;

        for (std::size_t y = 0; y < grid.ydim(); y++) {
            const FT yp = y * resolution - 1.0;

            for (std::size_t z = 0; z < grid.zdim(); z++) {
                const FT zp = z * resolution - 1.0;

                grid.value(x, y, z) = std::sqrt(xp * xp + yp * yp + zp * zp);
            }
        }
    }
    CGAL::Isosurfacing::Cartesian_grid_domain<Kernel> domain(grid);
    return domain;
    }

int main(int argc, char* argv[]) {
    std::unordered_map<std::string, std::function<int64_t(std::size_t)>> scenarios;
    scenarios["implicit_sphere"] = implicit_sphere;
    scenarios["grid_sphere"] = grid_sphere;

    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " <scenario> <N>" << std::endl;
        std::cout << "Available scenarios:" << std::endl;
        for (auto& s : scenarios) {
            std::cout << "  " << s.first << std::endl;
        }
        return 0;
    }

    const std::size_t N = std::stoull(argv[2]);

    auto it = scenarios.find(argv[1]);
    if (it == scenarios.end()) {
        std::cout << "Invalid scenario!" << std::endl;
        std::cout << "Available scenarios:" << std::endl;
        for (auto& s : scenarios) {
            std::cout << "  " << s.first << std::endl;
        }
        return 0;
    }

    auto& benchmark = it->second;
    int64_t ms = benchmark(N);

    std::cout << ms << std::endl;
}
