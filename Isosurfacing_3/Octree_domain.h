#ifndef CGAL_OCTREE_GRID_DOMAIN_H
#define CGAL_OCTREE_GRID_DOMAIN_H

#include "Octree_wrapper.h"
#include "types.h"

#include <array>

namespace CGAL {

    template<typename GeomTraits>
    class Octree_domain {
      public:
        typedef GeomTraits Geom_traits;
        typedef typename Geom_traits::FT FT;
        typedef typename Geom_traits::Point_3 Point_3;
        typedef typename Geom_traits::Vector_3 Vector_3;

        typedef Octree_wrapper<Geom_traits> Octree;
        typedef typename Octree::Vertex_handle Vertex_handle;
        typedef typename Octree::Edge_handle Edge_handle;
        typedef typename Octree::Voxel_handle Voxel_handle;

        typedef std::array<Vertex_handle, 2> Edge_vertices;
        typedef std::array<Voxel_handle, 4> Voxels_incident_to_edge;
        typedef std::array<Vertex_handle, 8> Voxel_vertices;
        typedef std::array<Edge_handle, 12> Voxel_edges;

      public:
        Octree_domain( const Octree& octree ) : octree_( &octree ) {}

        Point_3 position( const Vertex_handle& v ) const { return octree_->point( v ); }

        Vector_3 gradient( const Vertex_handle& v ) const { return octree_->gradient( v ); }

        FT value( const Vertex_handle& v ) const { return octree_->vertex_value( v ); }

        Edge_vertices edge_vertices( const Edge_handle& e_id ) const { return octree_->edge_vertices( e_id ); }

        Voxels_incident_to_edge voxels_incident_to_edge( const Edge_handle& e_id ) const { return octree_->edge_voxels( e_id ); }

        Voxel_vertices voxel_vertices( const Voxel_handle& vox ) const { return octree_->voxel_vertices( vox ); }

        Voxel_edges voxel_edges( const Voxel_handle& vox ) const { return octree_->voxel_edges( vox ); }

        template<typename Functor>
        void iterate_vertices( Functor& f ) const {
            for( const Vertex_handle& v: octree_->leaf_vertices() ) {
                f( v );
            }
        }

        template<typename Functor>
        void iterate_edges( Functor& f ) const {
            for( const Edge_handle& e: octree_->leaf_edges() ) {
                f( e );
            }
        }

        template<typename Functor>
        void iterate_voxels( Functor& f ) const {
            for( const Voxel_handle& v: octree_->leaf_voxels() ) {
                f( v );
            }
        }

      private:
        const Octree* octree_;
    };

}    // namespace CGAL

#endif    // CGAL_OCTREE_GRID_DOMAIN_H
