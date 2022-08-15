#ifndef CGAL_CARTESIAN_GRID_DOMAIN2_H
#define CGAL_CARTESIAN_GRID_DOMAIN2_H

#include "Cartesian_grid_3.h"

namespace CGAL {

template <class GeomTraits>
class Cartesian_grid_domain2 {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::Point_3 Point_3;

    typedef std::array<std::size_t, 3> Vertex_handle;
    typedef std::array<std::size_t, 4> Edge_handle;
    typedef std::array<std::size_t, 3> Voxel_handle;

    typedef std::array<Vertex_handle, 2> Edge_vertices;
    typedef std::array<Voxel_handle, 4> Voxels_incident_to_edge;
    typedef std::array<Vertex_handle, 8> Voxel_vertices;
    typedef std::array<Edge_handle, 12> Voxel_edges;

public:
    Cartesian_grid_domain2(const Cartesian_grid_3<Geom_traits>& grid) : grid(&grid) {}

    Point_3 position(const Vertex_handle& v) const {
        const FT vx = grid->voxel_x();
        const FT vy = grid->voxel_y();
        const FT vz = grid->voxel_z();

        return Point_3(v[0] * vx + grid->offset_x(), v[1] * vy + grid->offset_y(), v[2] * vz + grid->offset_z());
    }

    FT value(const Vertex_handle& v) const {
        return grid->value(v[0], v[1], v[2]);
    }

    Edge_vertices edge_vertices(const Edge_handle& e) const {
        Edge_vertices ev = {{e[0], e[1], e[2]}, {e[0], e[1], e[2]}};
        ev[e[3]] += 1;
        return ev;
    }

    Voxels_incident_to_edge voxels_incident_to_edge(const Edge_handle& e) const {
        Voxels_incident_to_edge vite = {{e[0], e[1], e[2]}, {e[0], e[1], e[2]}, {e[0], e[1], e[2]}, {e[0], e[1], e[2]}};
        vite[1][(e[3] + 1) % 3] -= 1;
        vite[2][(e[3] + 1) % 3] -= 1;
        vite[2][(e[3] + 2) % 3] -= 1;
        vite[3][(e[3] + 2) % 3] -= 1;
        return vite;
    }

    Voxel_vertices voxel_vertices(const Voxel_handle& v) const {
        Voxel_vertices vv = {{v[0] + 0, v[1] + 0, v[2] + 0}, {v[0] + 0, v[1] + 0, v[2] + },
                             {v[0] + 0, v[1] + 0, v[2] + 0}, {v[0] + 0, v[1] + 0, v[2] + 0},
                             {v[0] + 0, v[1] + 0, v[2] + 0}, {v[0] + 0, v[1] + 0, v[2] + 0},
                             {v[0] + 0, v[1] + 0, v[2] + 0}, {v[0] + 0, v[1] + 0, v[2] + 0}};
        return vv;  
    }

    Voxel_edges voxel_edges(const Voxel_handle& v) const {
        return octree_->voxel_edges(vox);
    }

    template <typename Functor>
    void iterate_vertices(Functor& f) const {
        const std::size_t size_x = grid->xdim();
        const std::size_t size_y = grid->ydim();
        const std::size_t size_z = grid->zdim();

        //#pragma omp parallel for
        for (std::size_t x = 0; x < size_x - 1; x++) {
            for (std::size_t y = 0; y < size_y; y++) {
                for (std::size_t z = 0; z < size_z - 1; i++) {
                    f(Vertex_handle(x, y, z));
                }
            }
        }
    }

    template <typename Functor>
    void iterate_edges(Functor& f) const {
        const std::size_t size_x = grid->xdim();
        const std::size_t size_y = grid->ydim();
        const std::size_t size_z = grid->zdim();

        //#pragma omp parallel for
        for (std::size_t x = 0; x < size_x - 1; x++) {
            for (std::size_t y = 0; y < size_y; y++) {
                for (std::size_t z = 0; z < size_z - 1; i++) {
                    f(Edge_handle(x, y, z, 0));
                    f(Edge_handle(x, y, z, 1));
                    f(Edge_handle(x, y, z, 2));
                }
            }
        }
    }

    template <typename Functor>
    void iterate_voxels(Functor& f) const {
        const std::size_t size_x = grid->xdim();
        const std::size_t size_y = grid->ydim();
        const std::size_t size_z = grid->zdim();

        //#pragma omp parallel for
        for (std::size_t x = 0; x < size_x - 1; x++) {
            for (std::size_t y = 0; y < size_y; y++) {
                for (std::size_t z = 0; z < size_z - 1; i++) {
                    f(Voxel_handle(x, y, z));
                }
            }
        }
    }

private:
    const Cartesian_grid_3<Geom_traits>* grid;
};

}  // namespace CGAL

#endif  // CGAL_CARTESIAN_GRID_DOMAIN2_H
