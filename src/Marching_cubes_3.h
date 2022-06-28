#ifndef CGAL_MARCHING_CUBES_3_H
#define CGAL_MARCHING_CUBES_3_H

#include <array>
#include <vector>
#include <mutex>



#include "Cartesian_grid_3.h"

// temporary
#include "Tables.h"
#include "util.h"

namespace CGAL {

template <class PointRange, class PolygonRange, class Traits>
class Marching_cubes_grid_3 {
public:
    typedef typename Traits::FT FT;
    typedef typename Traits::Point_3 Point_3;
    typedef typename Traits::Vector_3 Vector_3;

    typedef Cartesian_grid_3<FT> Grid;

    typedef PointRange Point_range;
    typedef PolygonRange Polygon_range;
    typedef typename Polygon_range::value_type Polygon;

public:
    Marching_cubes_grid_3(const Grid& grid, const FT iso_value, Point_range& points, Polygon_range& polygons)
        : grid(grid), iso_value(iso_value), points(points), polygons(polygons) {}

    void compute() {
        ScopeTimer timer;

        #pragma omp parallel for
        for (std::size_t x =  0; x < grid.xdim(); x++) {
            for (std::size_t y =  0; y < grid.ydim(); y++) {
                for (std::size_t z =  0; z < grid.zdim(); z++) {
                    marche_cell(x, y, z);
                }
            }
        }
    }

private:
    void marche_cell(const std::size_t x, const std::size_t y, const std::size_t z) {
        if (x == 0 || y == 0 || z == 0) {
            return;
        }
        if (x >= grid.xdim() - 2 || y >= grid.ydim() - 2 || z >= grid.zdim() - 2) {
            return;
        }

        const std::array<std::size_t, 3> idx0 = {x + 0, y + 1, z + 0};
        const std::array<std::size_t, 3> idx1 = {x + 1, y + 1, z + 0};
        const std::array<std::size_t, 3> idx2 = {x + 1, y + 0, z + 0};
        const std::array<std::size_t, 3> idx3 = {x + 0, y + 0, z + 0};
        const std::array<std::size_t, 3> idx4 = {x + 0, y + 1, z + 1};
        const std::array<std::size_t, 3> idx5 = {x + 1, y + 1, z + 1};
        const std::array<std::size_t, 3> idx6 = {x + 1, y + 0, z + 1};
        const std::array<std::size_t, 3> idx7 = {x + 0, y + 0, z + 1};

        const FT dist0 = grid.value(idx0[0], idx0[1], idx0[2]);
        const FT dist1 = grid.value(idx1[0], idx1[1], idx1[2]);
        const FT dist2 = grid.value(idx2[0], idx2[1], idx2[2]);
        const FT dist3 = grid.value(idx3[0], idx3[1], idx3[2]);
        const FT dist4 = grid.value(idx4[0], idx4[1], idx4[2]);
        const FT dist5 = grid.value(idx5[0], idx5[1], idx5[2]);
        const FT dist6 = grid.value(idx6[0], idx6[1], idx6[2]);
        const FT dist7 = grid.value(idx7[0], idx7[1], idx7[2]);

        unsigned int cubeindex = 0;
        //TODO compute cubeindex & voxel positions
        // set bit if corresponding corner is below iso
        cubeindex |= (dist0 < iso_value) << 0;
        cubeindex |= (dist1 < iso_value) << 1;
        cubeindex |= (dist2 < iso_value) << 2;
        cubeindex |= (dist3 < iso_value) << 3;
        cubeindex |= (dist4 < iso_value) << 4;
        cubeindex |= (dist5 < iso_value) << 5;
        cubeindex |= (dist6 < iso_value) << 6;
        cubeindex |= (dist7 < iso_value) << 7;

        const Vector_3 spacing(grid.voxel_x(), grid.voxel_y(), grid.voxel_z());

        const Vector_3 pos0 = Vector_3(idx0[0] * spacing.x(), idx0[1] * spacing.y(), idx0[2] * spacing.z());
        const Vector_3 pos1 = Vector_3(idx1[0] * spacing.x(), idx1[1] * spacing.y(), idx1[2] * spacing.z());
        const Vector_3 pos2 = Vector_3(idx2[0] * spacing.x(), idx2[1] * spacing.y(), idx2[2] * spacing.z());
        const Vector_3 pos3 = Vector_3(idx3[0] * spacing.x(), idx3[1] * spacing.y(), idx3[2] * spacing.z());
        const Vector_3 pos4 = Vector_3(idx4[0] * spacing.x(), idx4[1] * spacing.y(), idx4[2] * spacing.z());
        const Vector_3 pos5 = Vector_3(idx5[0] * spacing.x(), idx5[1] * spacing.y(), idx5[2] * spacing.z());
        const Vector_3 pos6 = Vector_3(idx6[0] * spacing.x(), idx6[1] * spacing.y(), idx6[2] * spacing.z());
        const Vector_3 pos7 = Vector_3(idx7[0] * spacing.x(), idx7[1] * spacing.y(), idx7[2] * spacing.z());

        Point_3 vertlist[12];
        //TODO generate vertexlist
        // bitmap of edges that intersect the iso-plane(s)
        const int edges = edgeTable[cubeindex];

        if (edges == 0) {
            return;
        }

        // interpolation of vertices incident to the edge, according to diagram
        if (edges & (1 << 0))  { vertlist[0]  = vertex_interpolation(pos0, pos1, dist0, dist1); }
        if (edges & (1 << 1))  { vertlist[1]  = vertex_interpolation(pos1, pos2, dist1, dist2); }
        if (edges & (1 << 2))  { vertlist[2]  = vertex_interpolation(pos2, pos3, dist2, dist3); }
        if (edges & (1 << 3))  { vertlist[3]  = vertex_interpolation(pos3, pos0, dist3, dist0); }
        if (edges & (1 << 4))  { vertlist[4]  = vertex_interpolation(pos4, pos5, dist4, dist5); }
        if (edges & (1 << 5))  { vertlist[5]  = vertex_interpolation(pos5, pos6, dist5, dist6); }
        if (edges & (1 << 6))  { vertlist[6]  = vertex_interpolation(pos6, pos7, dist6, dist7); }
        if (edges & (1 << 7))  { vertlist[7]  = vertex_interpolation(pos7, pos4, dist7, dist4); }
        if (edges & (1 << 8))  { vertlist[8]  = vertex_interpolation(pos0, pos4, dist0, dist4); }
        if (edges & (1 << 9))  { vertlist[9]  = vertex_interpolation(pos1, pos5, dist1, dist5); }
        if (edges & (1 << 10)) { vertlist[10] = vertex_interpolation(pos2, pos6, dist2, dist6); }
        if (edges & (1 << 11)) { vertlist[11] = vertex_interpolation(pos3, pos7, dist3, dist7); }

        std::lock_guard<std::mutex> lock(mutex);

        for (int i = 0; triTable[cubeindex][i] != -1; i += 3) {

            points.push_back(vertlist[triTable[cubeindex][i + 0]]);
            points.push_back(vertlist[triTable[cubeindex][i + 1]]);
            points.push_back(vertlist[triTable[cubeindex][i + 2]]);

            Polygon t;
            // TODO: this is not allowed
            t.push_back(points.size() - 3);
            t.push_back(points.size() - 2);
            t.push_back(points.size() - 1);

            polygons.push_back(t);
        }
    }

    Point_3 vertex_interpolation(const Vector_3& r1, const Vector_3& r2, const FT d1, const FT d2) {
        FT mu = -1.f;

        // don't divide by 0
        if (abs(d2 - d1) < 0.000001) {
            mu = 0.5f;  // if both points have the same value, assume isolevel is in the middle
        } else {
            mu = (iso_value - d1) / (d2 - d1);
        }

        if (mu < 0.f || mu > 1.f) {
            printf("ERROR: isolevel is not between points\n");  // TODO
        }

        const Vector_3 res = r2 * mu + r1 * (1 - mu);
        return Point_3(res.x(), res.y(), res.z());
    }

private:
    const Grid& grid;
    const FT iso_value;

    Point_range& points;
    Polygon_range& polygons;

    std::mutex mutex;
};


template <class Oracle, class PointRange, class PolygonRange, class Traits>
class Marching_cubes_implicit_3 {
public:
    typedef typename Traits::FT FT;
    typedef typename Traits::Point_3 Point_3;
    typedef typename Traits::Vector_3 Vector_3;

    typedef PointRange Point_range;
    typedef PolygonRange Polygon_range;
    typedef typename Polygon_range::value_type Polygon;

public:
    Marching_cubes_implicit_3(const Oracle& oracle, const FT iso_value, Point_range& points, Polygon_range& polygons)
        : oracle(oracle), iso_value(iso_value), points(points), polygons(polygons) {}

    void compute() {
        ScopeTimer timer;

        #pragma omp parallel for
        for (std::size_t x =  0; x < oracle.xdim(); x++) {
            for (std::size_t y =  0; y < oracle.ydim(); y++) {
                for (std::size_t z =  0; z < oracle.zdim(); z++) {
                    marche_cell(x, y, z);
                }
            }
        }
    }

private:
    typedef std::array<std::size_t, 3> Idx_3;

private:
    void marche_cell(const std::size_t x, const std::size_t y, const std::size_t z) {
        if (x == 0 || y == 0 || z == 0) {
            return;
        }
        if (x >= oracle.xdim() - 2 || y >= oracle.ydim() - 2 || z >= oracle.zdim() - 2) {
            return;
        }

        const Idx_3 idx0 = {x + 0, y + 1, z + 0};
        const Idx_3 idx1 = {x + 1, y + 1, z + 0};
        const Idx_3 idx2 = {x + 1, y + 0, z + 0};
        const Idx_3 idx3 = {x + 0, y + 0, z + 0};
        const Idx_3 idx4 = {x + 0, y + 1, z + 1};
        const Idx_3 idx5 = {x + 1, y + 1, z + 1};
        const Idx_3 idx6 = {x + 1, y + 0, z + 1};
        const Idx_3 idx7 = {x + 0, y + 0, z + 1};

        const FT dist0 = oracle.value(idx0[0], idx0[1], idx0[2]);
        const FT dist1 = oracle.value(idx1[0], idx1[1], idx1[2]);
        const FT dist2 = oracle.value(idx2[0], idx2[1], idx2[2]);
        const FT dist3 = oracle.value(idx3[0], idx3[1], idx3[2]);
        const FT dist4 = oracle.value(idx4[0], idx4[1], idx4[2]);
        const FT dist5 = oracle.value(idx5[0], idx5[1], idx5[2]);
        const FT dist6 = oracle.value(idx6[0], idx6[1], idx6[2]);
        const FT dist7 = oracle.value(idx7[0], idx7[1], idx7[2]);

        unsigned int cubeindex = 0;
        //TODO compute cubeindex & voxel positions
        // set bit if corresponding corner is below iso
        cubeindex |= (dist0 < iso_value) << 0;
        cubeindex |= (dist1 < iso_value) << 1;
        cubeindex |= (dist2 < iso_value) << 2;
        cubeindex |= (dist3 < iso_value) << 3;
        cubeindex |= (dist4 < iso_value) << 4;
        cubeindex |= (dist5 < iso_value) << 5;
        cubeindex |= (dist6 < iso_value) << 6;
        cubeindex |= (dist7 < iso_value) << 7;

        Point_3 vertlist[12];
        //TODO generate vertexlist
        // bitmap of edges that intersect the iso-plane(s)
        const int edges = edgeTable[cubeindex];

        if (edges == 0) {
            return;
        }

        // interpolation of vertices incident to the edge, according to diagram
        if (edges & (1 << 0))  { vertlist[0]  = vertex_interpolation(idx0, idx1, dist0, dist1); }
        if (edges & (1 << 1))  { vertlist[1]  = vertex_interpolation(idx1, idx2, dist1, dist2); }
        if (edges & (1 << 2))  { vertlist[2]  = vertex_interpolation(idx2, idx3, dist2, dist3); }
        if (edges & (1 << 3))  { vertlist[3]  = vertex_interpolation(idx3, idx0, dist3, dist0); }
        if (edges & (1 << 4))  { vertlist[4]  = vertex_interpolation(idx4, idx5, dist4, dist5); }
        if (edges & (1 << 5))  { vertlist[5]  = vertex_interpolation(idx5, idx6, dist5, dist6); }
        if (edges & (1 << 6))  { vertlist[6]  = vertex_interpolation(idx6, idx7, dist6, dist7); }
        if (edges & (1 << 7))  { vertlist[7]  = vertex_interpolation(idx7, idx4, dist7, dist4); }
        if (edges & (1 << 8))  { vertlist[8]  = vertex_interpolation(idx0, idx4, dist0, dist4); }
        if (edges & (1 << 9))  { vertlist[9]  = vertex_interpolation(idx1, idx5, dist1, dist5); }
        if (edges & (1 << 10)) { vertlist[10] = vertex_interpolation(idx2, idx6, dist2, dist6); }
        if (edges & (1 << 11)) { vertlist[11] = vertex_interpolation(idx3, idx7, dist3, dist7); }

        std::lock_guard<std::mutex> lock(mutex);

        for (int i = 0; triTable[cubeindex][i] != -1; i += 3) {

            points.push_back(vertlist[triTable[cubeindex][i + 0]]);
            points.push_back(vertlist[triTable[cubeindex][i + 1]]);
            points.push_back(vertlist[triTable[cubeindex][i + 2]]);

            Polygon t;
            // TODO: this is not allowed
            t.push_back(points.size() - 3);
            t.push_back(points.size() - 2);
            t.push_back(points.size() - 1);

            polygons.push_back(t);
        }
    }

    Point_3 vertex_interpolation(const Idx_3& idx1, const Idx_3& idx2, const FT d1, const FT d2) {
        FT mu = -1.f;

        // don't divide by 0
        if (abs(d2 - d1) < 0.000001) {
            mu = 0.5f;  // if both points have the same value, assume isolevel is in the middle
        } else {
            mu = (iso_value - d1) / (d2 - d1);
        }

        if (mu < 0.f || mu > 1.f) {
            printf("ERROR: isolevel is not between points\n");  // TODO
        }

        const Vector_3 p1 = oracle.pos(idx1[0], idx1[1], idx1[2]);
        const Vector_3 p2 = oracle.pos(idx2[0], idx2[1], idx2[2]);

        const Vector_3 res = p1 * mu + p2 * (1 - mu);
        return Point_3(res.x(), res.y(), res.z());
    }

private:
    const Oracle& oracle;
    const FT iso_value;

    Point_range& points;
    Polygon_range& polygons;

    std::mutex mutex;
};

} // end namespace CGAL

#endif // CGAL_MARCHING_CUBES_3_H