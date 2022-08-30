#ifndef CGAL_MARCHING_TETS_3_H
#define CGAL_MARCHING_TETS_3_H

#include "Cell_type.h"
#include "Isosurfacing_3/internal/Marching_cubes_3_internal.h"

namespace CGAL {
namespace Isosurfacing {

template <typename Concurrency_tag = Sequential_tag, class Domain_, class PointRange, class PolygonRange>
void make_triangle_mesh_using_marching_tets(const Domain_& domain, const typename Domain_::FT iso_value,
                                            PointRange& points, PolygonRange& polygons) {

    // static_assert(Domain_::CELL_TYPE & TETRAHEDRAL_CELL);

    internal::Marching_tets_functor<Domain_, PointRange, PolygonRange> functor(domain, iso_value, points, polygons);
    domain.iterate_cells(functor, Concurrency_tag());
}

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_MARCHING_TETS_3_H
