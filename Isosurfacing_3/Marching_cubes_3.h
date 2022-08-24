#ifndef CGAL_MARCHING_CUBES_3_H
#define CGAL_MARCHING_CUBES_3_H

#include <mutex>

#include "Isosurfacing_3/internal/Marching_cubes_3_internal.h"

namespace CGAL {
namespace Isosurfacing {

template <class Domain_, class PointRange, class PolygonRange>
void make_triangle_mesh_using_marching_cubes_old(const Domain_& domain, const typename Domain_::FT iso_value,
                                                 PointRange& points, PolygonRange& polygons) {

    std::mutex mutex;

    const std::size_t size_k = domain.size_z();
    const std::size_t size_j = domain.size_y();
    const std::size_t size_i = domain.size_x();

    const std::size_t blocking_size = 100;

    // TODO: look at polygon mesh processing for tbb (also linking)

    // compute a unique global index for vertices
    // use as key the unique edge number
    std::unordered_map<std::size_t, std::size_t> v_map;


    for (std::size_t bj = 0; bj < size_j - 1; bj += blocking_size) {
        //#pragma omp parallel for
        for (std::size_t k = 0; k < size_k - 1; k++) {

            const std::size_t j_start = bj;
            const std::size_t j_end = std::min(size_j - 1, bj + blocking_size);

            for (std::size_t j = j_start; j < j_end; j++) {
                for (std::size_t i = 0; i < size_i - 1; i++) {
                    // internal::marching_cubes_cell_old(i, j, k, domain, iso_value, points, polygons, mutex);
                    internal::marching_cubes_cell_RG(i, j, k, domain, iso_value, points, polygons, mutex, v_map);
                }
            }
        }
    }
}

template <typename Concurrency_tag = Sequential_tag, class Domain_, class PointRange, class PolygonRange>
void make_triangle_mesh_using_marching_cubes(const Domain_& domain, const typename Domain_::FT iso_value,
                                             PointRange& points, PolygonRange& polygons) {

    internal::Marching_cubes_functor<Domain_, PointRange, PolygonRange> functor(domain, iso_value, points, polygons);
    domain.iterate_cells(functor, Concurrency_tag());
}

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_MARCHING_CUBES_3_H
