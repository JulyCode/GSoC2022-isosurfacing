#ifndef CGAL_TMC_INTERNAL_TMC_H
#define CGAL_TMC_INTERNAL_TMC_H

#include <array>
#include <map>
#include <mutex>

#include "Tables.h"
#include "Tables_old.h"

namespace CGAL {
namespace Isosurfacing {
namespace internal {

void p_slice(const int i_index, const int j_index, const int k_index, const double i0, double* F,
                                 Point& p, Vector& n, const int i_case) {
    // there are 12 edges, assign to each vertex three edges, the global edge numbering
    // consist of 3*global_vertex_id + edge_offset.
    const unsigned long long gei_pattern_ = 670526590282893600ull;
    const uint axis_mask = 1;
    const uint offs_mask = 3;

    // code edge end vertices for each of the 12 edges
    const unsigned char l_edges_[12] = {16, 49, 50, 32, 84, 117, 118, 100, 64, 81, 115, 98};
    auto get_edge_vertex = [](const int e, unsigned int& v0, unsigned int& v1, const unsigned char l_edges_[12]) {
        v0 = (unsigned int)(l_edges_[e] & 0xF);
        v1 = (unsigned int)(l_edges_[e] >> 4) & 0xF;
    };

    // A hexahedron has twelve edges, save the intersection of the isosurface with the edge
    // save global edge and global vertex index of isosurface
    std::vector<Vertex> vertices(12);
    // save coordinates of intersection points in 3D
    std::vector<Point> ip(12);
    // save normals of intersection points from scalar data
    std::vector<Normal> in(12);
    // save loca coordinate along the edge of intersection point
    std::vector<double> ecoord{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    // collect vertices
    ushort flag{1};
    for (int eg = 0; eg < 12; eg++) {
        vertices[eg].g_edg = -1;
        vertices[eg].g_idx = -1;
        if (flag & e_pattern[i_case]) {
            // the edge global index is given by the vertex global index + the edge offset
            uint shift = 5 * eg;
            const int ix = i_index + (int)((gei_pattern_ >> shift) & 1);        // global_edge_id[eg][0];
            const int iy = j_index + (int)((gei_pattern_ >> (shift + 1)) & 1);  // global_edge_id[eg][1];
            const int iz = k_index + (int)((gei_pattern_ >> (shift + 2)) & 1);  // global_edge_id[eg][2];
            const int off_val = (int)((gei_pattern_ >> (shift + 3)) & 3);

            vertices[eg].g_edg = int(m_cell_shift_factor * m_ugrid.global_index(ix, iy, iz) + off_val);

            // generate vertex here, do not care at this point if vertex already exist
            // int* vert = l_edges[eg];
            // interpolation weight
            // uint v0 = (l_edges_[eg]&0xF);
            // uint v1 = (l_edges_[eg]>>4)&0xF;
            //            int v0 = vert[0];
            //            int v1 = vert[1];
            uint v0, v1;
            get_edge_vertex(eg, v0, v1, l_edges_);

            double l = (i0 - F[v0]) / (F[v1] - F[v0]);
            ecoord[eg] = l;
            // interpolate vertex
            ip[eg][0] = (1 - l) * p[v0][0] + l * p[v1][0];
            ip[eg][1] = (1 - l) * p[v0][1] + l * p[v1][1];
            ip[eg][2] = (1 - l) * p[v0][2] + l * p[v1][2];

            // interpolate normal
            in[eg][0] = (1 - l) * n[v0][0] + l * n[v1][0];
            in[eg][1] = (1 - l) * n[v0][1] + l * n[v1][1];
            in[eg][2] = (1 - l) * n[v0][2] + l * n[v1][2];

            const double n_size = std::sqrt(in[eg][0] * in[eg][0] + in[eg][1] * in[eg][1] + in[eg][2] * in[eg][2]);
            in[eg][0] /= n_size;
            in[eg][1] /= n_size;
            in[eg][2] /= n_size;

            // set vertex in map
            // set vertex index
            auto s_index = m_vertices.find(vertices[eg].g_edg);
            if (s_index == m_vertices.end()) {
                const int g_idx = (int)m_points.size();
                vertices[eg].g_idx = g_idx;
                m_vertices[vertices[eg].g_edg] = g_idx;
                m_points.push_back(ip[eg]);
                m_pnorms.push_back(in[eg]);
            } else {
                vertices[eg].g_idx = s_index->second;
            }
        }
        /*else {
        e_set[eg] = false;
        }*/
        // next edge
        flag <<= 1;
    }

    // compute oriented contours
    // A countour consists of segment at the faces connecting the intersection of the
    // iso-surface with the edges. For each edge we store the edge to which the segment
    // is outgoing and the edge from which the segment in comming. Therefore a contour
    // cab be reconstructed by connecting the edges in the direccion of the outgoing.
    // The contour is oriented in such a way, that the positive vertices are outside.
    // 1. build segments
    // 2. connect segments
    // build up segments
    // set segments map
    unsigned char segm_[12] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF};
    auto set_segm = [](const int e, const int pos, const int val, unsigned char segm_[12]) {
        if (pos == 0) {
            segm_[e] &= 0xF0;
            segm_[e] |= (unsigned char)val & 0xF;
        } else if (pos == 1) {
            segm_[e] &= 0xF;
            segm_[e] |= val << 4;
        }
    };
    auto get_segm = [](const int e, const int pos, unsigned char segm_[12]) {
        if (pos == 0)
            return (int)(segm_[e] & 0xF);
        else
            return (int)((segm_[e] >> 4) & 0xF);
    };
    auto is_segm_set = [](const int e, unsigned char segm_[12]) { return (segm_[e] != 0xFF); };
    auto unset_segm = [](const int e, unsigned char segm_[12]) { segm_[e] = 0xFF; };
    // In order to compute oriented segments, the hexahedron has to be flatten.
    // The insides of the faces of the hexahedron have to be all at the same
    // side of the flattend hexa. This requires changing the order of the
    // edges when reading from the faces
    // code edges at face
    // unsigned short face_e_[6] = { 12816, 30292, 33936, 46754, 34739, 38305 };
    std::array<unsigned short, 6> e_face_{{291, 18277, 18696, 10859, 33719, 38305}};
    // code vertices at face
    // unsigned short face_v_[6] = { 12816, 30292, 21520, 30258, 25632, 30001 };
    std::array<unsigned short, 6> v_face_{{12576, 25717, 5380, 29538, 8292, 30001}};

    // reading edge from face
    auto get_face_e = [e_face_](const int f, const int e) { return ((e_face_[f] >> (4 * e)) & 0xF); };
    auto get_face_v = [v_face_](const int f, const int e) { return ((v_face_[f] >> (4 * e)) & 0xF); };
    // compute oriented segments using the isoline scheme at the faces
    const unsigned int BIT_1 = 1;
    const unsigned int BIT_2 = 2;
    const unsigned int BIT_3 = 4;
    const unsigned int BIT_4 = 8;
    auto asymptotic_decider = [](const double f0, const double f1, const double f2, const double f3) {
        return (f0 * f3 - f1 * f2) / (f0 + f3 - f1 - f2);
    };
    std::vector<bool> f_flag(6, false);
    for (int f = 0; f < 6; f++) {
        // classify face
        unsigned int f_case{0};
        uint v0 = get_face_v(f, 0);
        uint v1 = get_face_v(f, 1);
        uint v2 = get_face_v(f, 2);
        uint v3 = get_face_v(f, 3);
        uint e0 = get_face_e(f, 0);
        uint e1 = get_face_e(f, 1);
        uint e2 = get_face_e(f, 2);
        uint e3 = get_face_e(f, 3);
        double f0 = F[v0];
        double f1 = F[v1];
        double f2 = F[v2];
        double f3 = F[v3];
        if (f0 >= i0) f_case |= BIT_1;
        if (f1 >= i0) f_case |= BIT_2;
        if (f2 >= i0) f_case |= BIT_3;
        if (f3 >= i0) f_case |= BIT_4;
        switch (f_case) {
            case 1:
                set_segm(e0, 0, e3, segm_);
                set_segm(e3, 1, e0, segm_);
                break;
            case 2:
                set_segm(e1, 0, e0, segm_);
                set_segm(e0, 1, e1, segm_);
                break;
            case 3:
                set_segm(e1, 0, e3, segm_);
                set_segm(e3, 1, e1, segm_);
                break;
            case 4:
                set_segm(e3, 0, e2, segm_);
                set_segm(e2, 1, e3, segm_);
                break;
            case 5:
                set_segm(e0, 0, e2, segm_);
                set_segm(e2, 1, e0, segm_);
                break;
            case 6: {
                const double val = asymptotic_decider(f0, f1, f2, f3);
                if (val > i0) {
                    set_segm(e3, 0, e0, segm_);
                    set_segm(e0, 1, e3, segm_);
                    set_segm(e1, 0, e2, segm_);
                    set_segm(e2, 1, e1, segm_);
                } else if (val < i0) {
                    set_segm(e1, 0, e0, segm_);
                    set_segm(e0, 1, e1, segm_);
                    set_segm(e3, 0, e2, segm_);
                    set_segm(e2, 1, e3, segm_);
                } else {
                    f_flag[f] = true;
                    // singular case val == i0, there are no asymptotes
                    // check if there is a reasonable triangulation of the face
                    unsigned short e_flag = e_flag = 0x218;
                    unsigned short bit_1 = 0x1;
                    unsigned short bit_2 = 0x2;
                    double ec0 = ecoord[e0];
                    double ec1 = ecoord[e1];
                    double ec2 = ecoord[e2];
                    double ec3 = ecoord[e3];
                    if ((e_flag >> (f * 2)) & bit_1) {
                        ec0 = 1 - ec0;
                        ec2 = 1 - ec2;
                    }
                    if ((e_flag >> (f * 2)) & bit_2) {
                        ec1 = 1 - ec1;
                        ec3 = 1 - ec3;
                    }
                    if (ec1 < ec3 && ec0 > ec2) {
                        set_segm(e1, 0, e0, segm_);
                        set_segm(e0, 1, e1, segm_);
                        set_segm(e3, 0, e2, segm_);
                        set_segm(e2, 1, e3, segm_);
                    } else if (ec1 > ec3 && ec0 < ec2) {
                        set_segm(e3, 0, e0, segm_);
                        set_segm(e0, 1, e3, segm_);
                        set_segm(e1, 0, e2, segm_);
                        set_segm(e2, 1, e1, segm_);
                    } else {
                        std::cerr << "ERROR: can't correctly triangulate cell's face\n";
                        return;
                    }
                }
            } break;
            case 7:
                set_segm(e1, 0, e2, segm_);
                set_segm(e2, 1, e1, segm_);
                break;
            case 8:
                set_segm(e2, 0, e1, segm_);
                set_segm(e1, 1, e2, segm_);
                break;
            case 9: {
                const double val = asymptotic_decider(f0, f1, f2, f3);
                if (val > i0) {
                    set_segm(e0, 0, e1, segm_);
                    set_segm(e1, 1, e0, segm_);
                    set_segm(e2, 0, e3, segm_);
                    set_segm(e3, 1, e2, segm_);
                } else if (val < i0) {
                    set_segm(e0, 0, e3, segm_);
                    set_segm(e3, 1, e0, segm_);
                    set_segm(e2, 0, e1, segm_);
                    set_segm(e1, 1, e2, segm_);
                } else {
                    f_flag[f] = true;
                    // singular case val == i0, there are no asymptotes
                    // check if there is a reasonable triangulation of the face
                    unsigned short e_flag = e_flag = 0x218;
                    unsigned short bit_1 = 0x1;
                    unsigned short bit_2 = 0x2;
                    double ec0 = ecoord[e0];
                    double ec1 = ecoord[e1];
                    double ec2 = ecoord[e2];
                    double ec3 = ecoord[e3];
                    if ((e_flag >> (f * 2)) & bit_1) {
                        ec0 = 1 - ec0;
                        ec2 = 1 - ec2;
                    }
                    if ((e_flag >> (f * 2)) & bit_2) {
                        ec1 = 1 - ec1;
                        ec3 = 1 - ec3;
                    }
                    if (ec1 < ec3 && ec0 > ec2) {
                        set_segm(e0, 0, e1, segm_);
                        set_segm(e1, 1, e0, segm_);
                        set_segm(e2, 0, e3, segm_);
                        set_segm(e3, 1, e2, segm_);
                    } else if (ec1 > ec3 && ec0 < ec2) {
                        set_segm(e0, 0, e3, segm_);
                        set_segm(e3, 1, e0, segm_);
                        set_segm(e2, 0, e1, segm_);
                        set_segm(e1, 1, e2, segm_);
                    } else {
                        std::cerr << "ERROR: can't correctly triangulate cell's face\n";
                        return;
                    }
                }
            } break;
            case 10:
                set_segm(e2, 0, e0, segm_);
                set_segm(e0, 1, e2, segm_);

                break;
            case 11:
                set_segm(e2, 0, e3, segm_);
                set_segm(e3, 1, e2, segm_);

                break;
            case 12:
                set_segm(e3, 0, e1, segm_);
                set_segm(e1, 1, e3, segm_);

                break;
            case 13:
                set_segm(e0, 0, e1, segm_);
                set_segm(e1, 1, e0, segm_);

                break;
            case 14:
                set_segm(e3, 0, e0, segm_);
                set_segm(e0, 1, e3, segm_);
                break;
            default:
                break;
        }
    }

    // connect oriented segments into oriented contours
    // closed contours are coded in 64 bit unsigned long long
    // 1) Each entry has 4 bits
    // 2) The first 4 entries are reserved for the size of the contours
    // 3) The next 12 entries are the indices of the edges constituting the contorus
    //    The indices are numbers from 0 to 12
    unsigned long long c_ = 0xFFFFFFFFFFFF0000;
    // in the 4 first bits store size of contours
    auto get_cnt_size = [](const int cnt, unsigned long long& c_) {
        return (size_t)((c_ & (0xF << 4 * cnt)) >> 4 * cnt);
    };
    auto set_cnt_size = [](const int cnt, const int size, unsigned long long& c_) {
        // unset contour size
        c_ &= ~(0xF << 4 * cnt);
        c_ |= (size << 4 * cnt);
    };
    // set corresponging edge
    auto set_c = [](const int cnt, const int pos, const int val, unsigned long long& c_) {
        const uint mask[4] = {0x0, 0xF, 0xFF, 0xFFF};
        const uint c_sz = c_ & mask[cnt];
        const uint e = 16 + 4 * ((c_sz & 0xF) + ((c_sz & 0xF0) >> 4) + ((c_sz & 0xF00) >> 8) + pos);
        c_ &= ~(((unsigned long long)0xF) << e);
        c_ |= (((unsigned long long)val) << e);
    };
    // read edge from contour
    auto get_c = [](const int cnt, const int pos, unsigned long long c_) {
        const uint mask[4] = {0x0, 0xF, 0xFF, 0xFFF};
        const uint c_sz = (uint)(c_ & mask[cnt]);
        const uint e = 16 + 4 * ((c_sz & 0xF) + ((c_sz & 0xF0) >> 4) + ((c_sz & 0xF00) >> 8) + pos);
        return (int)((c_ >> e) & 0xF);
    };


    // connect oriented contours
    uint cnt_{0};
    for (uint e = 0; e < 12; e++) {
        if (is_segm_set(e, segm_)) {
            uint eTo = get_segm(e, 0, segm_);
            uint eIn = get_segm(e, 1, segm_);
            uint eStart = e;
            uint pos = 0;
            set_c(cnt_, pos, eStart, c_);
            while (eTo != eStart) {
                pos = pos + 1;
                set_c(cnt_, pos, eTo, c_);
                eIn = eTo;
                eTo = get_segm(eIn, 0, segm_);
                unset_segm(eIn, segm_);
            }
            // set contour length
            set_cnt_size(cnt_, pos + 1, c_);
            // update number of contours
            cnt_ = cnt_ + 1;
        }
    }

    // compute intersection of opposite faces
    // It is enough to compute a pair of solutions for one face
    // The other solutions are obtained by evaluating the equations
    // for the common variable
    double ui[2]{};
    double vi[2]{};
    double wi[2]{};
    unsigned char q_sol{0};
    const double a = (F[0] - F[1]) * (-F[6] + F[7] + F[4] - F[5]) - (F[4] - F[5]) * (-F[2] + F[3] + F[0] - F[1]);
    const double b = (i0 - F[0]) * (-F[6] + F[7] + F[4] - F[5]) + (F[0] - F[1]) * (F[6] - F[4]) -
                     (i0 - F[4]) * (-F[2] + F[3] + F[0] - F[1]) - (F[4] - F[5]) * (F[2] - F[0]);
    const double c = (i0 - F[0]) * (F[6] - F[4]) - (i0 - F[4]) * (F[2] - F[0]);
    ;
    double d = b * b - 4 * a * c;
    if (d > 0) {
        d = std::sqrt(d);
        // compute u-coord of solutions
        ui[0] = (-b - d) / (2 * a);
        ui[1] = (-b + d) / (2 * a);
        // compute v-coord of solutions
        double g1 = F[0] * (1 - ui[0]) + F[1] * ui[0];
        double g2 = F[2] * (1 - ui[0]) + F[3] * ui[0];
        vi[0] = (i0 - g1) / (g2 - g1);
        if (std::isnan(vi[0]) || std::isinf(vi[0])) vi[0] = -1.f;
        g1 = F[0] * (1 - ui[1]) + F[1] * ui[1];
        g2 = F[2] * (1 - ui[1]) + F[3] * ui[1];
        vi[1] = (i0 - g1) / (g2 - g1);
        if (std::isnan(vi[1]) || std::isinf(vi[1])) vi[1] = -1.f;
        // compute w-coordinates of solutions
        g1 = F[0] * (1 - ui[0]) + F[1] * ui[0];
        g2 = F[4] * (1 - ui[0]) + F[5] * ui[0];
        wi[0] = (i0 - g1) / (g2 - g1);
        if (std::isnan(wi[0]) || std::isinf(wi[0])) wi[0] = -1.f;
        g1 = F[0] * (1 - ui[1]) + F[1] * ui[1];
        g2 = F[4] * (1 - ui[1]) + F[5] * ui[1];
        wi[1] = (i0 - g1) / (g2 - g1);
        if (std::isnan(wi[1]) || std::isinf(wi[1])) wi[1] = -1.f;
        // correct values for roots of quadratic equations
        // in case the asymptotic decider has failed
        if (f_flag[0] == true) {  // face 1, w = 0;
            if (wi[0] < wi[1])
                wi[0] = 0;
            else
                wi[1] = 0;
        }
        if (f_flag[1] == true) {  // face 2, w = 1
            if (wi[0] > wi[1])
                wi[1] = 1;
            else
                wi[1] = 1;
        }
        if (f_flag[2] == true) {  // face 3, v = 0
            if (vi[0] < vi[1])
                vi[0] = 0;
            else
                vi[1] = 0;
        }
        if (f_flag[3] == true) {  // face 4, v = 1
            if (vi[0] > vi[1])
                vi[0] = 1;
            else
                vi[1] = 1;
        }
        if (f_flag[4] == true) {  // face 5, u = 0
            if (ui[0] < ui[1])
                ui[0] = 0;
            else
                ui[1] = 0;
        }
        if (f_flag[5] == true) {  // face 6, u = 1
            if (ui[0] > ui[1])
                ui[0] = 1;
            else
                ui[1] = 1;
        }

        // check solution intervals
        if (0 < ui[0] && ui[0] < 1) {
            q_sol |= 1;
        }
        if (0 < ui[1] && ui[1] < 1) {
            q_sol |= 2;
        }
        if (0 < vi[0] && vi[0] < 1) {
            q_sol |= 4;
        }
        if (0 < vi[1] && vi[1] < 1) {
            q_sol |= 8;
        }
        if (0 < wi[0] && wi[0] < 1) {
            q_sol |= 16;
        }
        if (0 < wi[1] && wi[1] < 1) {
            q_sol |= 32;
        }
    }

    //
    // count the number of set bits
    auto numberOfSetBits = [](const unsigned char n) {
        // C or C++: use uint32_t
        uint b = (uint)n;
        b = b - ((b >> 1) & 0x55555555);
        b = (b & 0x33333333) + ((b >> 2) & 0x33333333);
        return (((b + (b >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
    };
    // compute the number of solutions to the quadratic equation for a given face
    auto nrQSolFace = [](const uint f, const unsigned char n) {
        uint nr{0};
        switch (f) {
            case 0:
                if ((n & 0x5) == 0x5) nr = nr + 1;
                if ((n & 0xA) == 0xA) nr = nr + 1;
                break;
            case 1:
                if ((n & 0x11) == 0x11) nr = nr + 1;
                if ((n & 0x22) == 0x22) nr = nr + 1;
                break;
            case 2:
                if ((n & 0x18) == 0x18) nr = nr + 1;
                if ((n & 0x24) == 0x24) nr = nr + 1;
                break;
        }
        return nr;
    };


    // triangulate contours
    // if all bits are set, then there are three pairs of nontrivial solutions
    // to the quadratic equations. In this case, there is a tunnel or a contour
    // with 12 vertices. If there are three contours, then there is a tunnel and
    // one of the contorus with only three vertices is not part of it.
    if (numberOfSetBits(q_sol) == 6) {
        m_ccases_tunnel += 1;
        // there are at most three contours
        // Possible cases:
        //  1) a single contour with 12 vertices
        //  2) two contours which build a tunnel
        //  3) three contours, one has only 3 vertices and does not belong to the tunnel

        // construct the six vertices of the inner hexagon
        double hvt[6][3];
        hvt[0][0] = ui[0];
        hvt[0][1] = vi[0];
        hvt[0][2] = wi[0];
        hvt[1][0] = ui[0];
        hvt[1][1] = vi[0];
        hvt[1][2] = wi[1];
        hvt[2][0] = ui[1];
        hvt[2][1] = vi[0];
        hvt[2][2] = wi[1];
        hvt[3][0] = ui[1];
        hvt[3][1] = vi[1];
        hvt[3][2] = wi[1];
        hvt[4][0] = ui[1];
        hvt[4][1] = vi[1];
        hvt[4][2] = wi[0];
        hvt[5][0] = ui[0];
        hvt[5][1] = vi[1];
        hvt[5][2] = wi[0];

        // construct vertices at intersections with the edges
        auto e_vert = [&ecoord](const int e, const int i) {
            const unsigned int l_coord[3]{1324855, 5299420, 16733440};
            unsigned char flag = (l_coord[i] >> (2 * e)) & 3;
            if (flag == 3)
                return ecoord[e];
            else
                return (double)(flag);
        };

        // if there are three contours, then there is a tunnel and one
        // of the contours is not part of it.
        unsigned char _not_tunnel = 0xF;
        if (cnt_ == 3) {
            // loop over the contorus
            // triangulate the contour which is not part of
            // the tunnel
            const double uc_min = (ui[0] < ui[1]) ? ui[0] : ui[1];
            const double uc_max = (ui[0] < ui[1]) ? ui[1] : ui[0];
            for (int t = 0; t < (int)cnt_; t++) {
                if (get_cnt_size(t, c_) == 3) {
                    double umin = 2;
                    double umax = -2;
                    uint e0 = get_c(t, 0, c_);
                    uint e1 = get_c(t, 1, c_);
                    uint e2 = get_c(t, 2, c_);
                    const double u_e0 = e_vert(e0, 0);
                    const double u_e1 = e_vert(e1, 0);
                    const double u_e2 = e_vert(e2, 0);
                    umin = (u_e0 < umin) ? u_e0 : umin;
                    umin = (u_e1 < umin) ? u_e1 : umin;
                    umin = (u_e2 < umin) ? u_e2 : umin;
                    umax = (u_e0 > umax) ? u_e0 : umax;
                    umax = (u_e1 > umax) ? u_e1 : umax;
                    umax = (u_e2 > umax) ? u_e1 : umax;
                    if (uc_min > umax || uc_max < umin) {
                        // this contour is not part of the tunnel
                        _not_tunnel = t;
                        Triangle tr;
                        tr.v[0] = vertices[e0].g_idx;
                        tr.v[1] = vertices[e1].g_idx;
                        tr.v[2] = vertices[e2].g_idx;
                        m_triangles.push_back(tr);
                        m_ccases_3++;
                    }
                }
            }
        }

        // compute vertices of inner hexagon, save new vertices in list and compute and keep
        // global vertice index to build triangle connectivity later on.
        uint tg_idx[6];
        Point p_ih[6];  // remember the six points of inner hexagon for debugging
        for (int i = 0; i < 6; i++) {
            Point hp;
            Normal hn;
            const double u = hvt[i][0];
            const double v = hvt[i][1];
            const double w = hvt[i][2];
            hp[0] =
                (1 - w) * ((1 - v) * (p[0][0] + u * (p[1][0] - p[0][0])) + v * (p[2][0] + u * (p[3][0] - p[2][0]))) +
                w * ((1 - v) * (p[4][0] + u * (p[5][0] - p[4][0])) + v * (p[6][0] + u * (p[7][0] - p[6][0])));
            hp[1] =
                (1 - w) * ((1 - v) * (p[0][1] + u * (p[1][1] - p[0][1])) + v * (p[2][1] + u * (p[3][1] - p[2][1]))) +
                w * ((1 - v) * (p[4][1] + u * (p[5][1] - p[4][1])) + v * (p[6][1] + u * (p[7][1] - p[6][1])));
            hp[2] =
                (1 - w) * ((1 - v) * (p[0][2] + u * (p[1][2] - p[0][2])) + v * (p[2][2] + u * (p[3][2] - p[2][2]))) +
                w * ((1 - v) * (p[4][2] + u * (p[5][2] - p[4][2])) + v * (p[6][2] + u * (p[7][2] - p[6][2])));
            hn[0] =
                (1 - w) * ((1 - v) * (n[0][0] + u * (n[1][0] - n[0][0])) + v * (n[2][0] + u * (n[3][0] - n[2][0]))) +
                w * ((1 - v) * (n[4][0] + u * (n[5][0] - n[4][0])) + v * (n[6][0] + u * (n[7][0] - n[6][0])));
            hn[1] =
                (1 - w) * ((1 - v) * (n[0][1] + u * (n[1][1] - n[0][1])) + v * (n[2][1] + u * (n[3][1] - n[2][1]))) +
                w * ((1 - v) * (n[4][1] + u * (n[5][1] - n[4][1])) + v * (n[6][1] + u * (n[7][1] - n[6][1])));
            hn[2] =
                (1 - w) * ((1 - v) * (n[0][2] + u * (n[1][2] - n[0][2])) + v * (n[2][2] + u * (n[3][2] - n[2][2]))) +
                w * ((1 - v) * (n[4][2] + u * (n[5][2] - n[4][2])) + v * (n[6][2] + u * (n[7][2] - n[6][2])));
            // normalize normal
            const double factor = std::sqrt(hn[0] * hn[0] + hn[1] * hn[1] + hn[2] * hn[2]);
            hn[0] = hn[0] / factor;
            hn[1] = hn[1] / factor;
            hn[2] = hn[2] / factor;
            tg_idx[i] = (uint)m_points.size();
            m_points.push_back(hp);
            m_pnorms.push_back(hn);
            p_ih[i] = hp;
        }

        // triangulate contours with inner hexagon
        std::vector<Triangle> lt_;  // remember triangles for debugging
        unsigned char tcon_[12];
        for (int i = 0; i < (int)cnt_; i++) {
            if (_not_tunnel != i) {  // contour belongs to tunnel
                const int cnt_sz = (int)get_cnt_size(i, c_);
                for (int r = 0; r < cnt_sz; r++) {
                    uint index = -1;
                    double dist = 1000.;
                    uint ci = get_c(i, r, c_);
                    const double u_edge = e_vert(ci, 0);
                    const double v_edge = e_vert(ci, 1);
                    const double w_edge = e_vert(ci, 2);
                    for (int s = 0; s < 6; s++) {
                        const double uval = u_edge - hvt[s][0];
                        const double vval = v_edge - hvt[s][1];
                        const double wval = w_edge - hvt[s][2];
                        double val = uval * uval + vval * vval + wval * wval;
                        if (dist > val) {
                            index = s;
                            dist = val;
                        }
                    }
                    tcon_[ci] = (unsigned char)index;
                }

                // correspondence between vertices found
                // create triangles
                // needs some functions
                auto distanceRingIntsModulo = [](const int d1, const int d2) {
                    const int r = (d1 - d2) < 0 ? d2 - d1 : d1 - d2;
                    return (r > 2 ? 6 - r : r);
                };
                auto midpointRingIntModulo = [](const int d1, const int d2) {
                    const int dmax = (d1 > d2) ? d1 : d2;
                    const int dmin = (d1 < d2) ? d1 : d2;
                    return ((dmax + 2) % 6 == dmin) ? (dmax + 1) % 6 : (dmax + dmin) / 2;
                };

                for (int r = 0; r < cnt_sz; r++) {
                    const uint tid1 = get_c(i, r, c_);
                    const uint tid2 = get_c(i, ((r + 1) % cnt_sz), c_);
                    const uint cid1 = tcon_[tid1];
                    const uint cid2 = tcon_[tid2];
                    // compute index distance
                    const int dst = distanceRingIntsModulo(cid1, cid2);
                    switch (dst) {
                        case 0: {
                            Triangle tr;
                            tr.v[0] = vertices[tid1].g_idx;
                            tr.v[1] = vertices[tid2].g_idx;
                            tr.v[2] = tg_idx[cid1];
                            m_triangles.push_back(tr);
                            lt_.push_back(tr);
                        } break;
                        case 1: {
                            // measure diagonals
                            // triangulate along shortest diagonal
                            double u_edge = e_vert(tid1, 0);
                            double v_edge = e_vert(tid1, 1);
                            double w_edge = e_vert(tid1, 2);
                            const double l1 = (u_edge - hvt[cid2][0]) * (u_edge - hvt[cid2][0]) +
                                              (v_edge - hvt[cid2][1]) * (v_edge - hvt[cid2][1]) +
                                              (w_edge - hvt[cid2][2]) * (w_edge - hvt[cid2][2]);
                            u_edge = e_vert(tid2, 0);
                            v_edge = e_vert(tid2, 1);
                            w_edge = e_vert(tid2, 2);
                            const double l2 = (u_edge - hvt[cid1][0]) * (u_edge - hvt[cid1][0]) +
                                              (v_edge - hvt[cid1][1]) * (v_edge - hvt[cid1][1]) +
                                              (w_edge - hvt[cid1][2]) * (w_edge - hvt[cid1][2]);
                            Triangle t1;
                            Triangle t2;
                            if (l1 < l2) {
                                t1.v[0] = vertices[tid1].g_idx;
                                t1.v[1] = vertices[tid2].g_idx;
                                t1.v[2] = tg_idx[cid2];
                                t2.v[0] = vertices[tid1].g_idx;
                                t2.v[1] = tg_idx[cid2];
                                t2.v[2] = tg_idx[cid1];
                            } else {
                                t1.v[0] = vertices[tid1].g_idx;
                                t1.v[1] = vertices[tid2].g_idx;
                                t1.v[2] = tg_idx[cid1];
                                t2.v[0] = vertices[tid2].g_idx;
                                t2.v[1] = tg_idx[cid2];
                                t2.v[2] = tg_idx[cid1];
                            }
                            m_triangles.push_back(t1);
                            m_triangles.push_back(t2);
                            lt_.push_back(t1);
                            lt_.push_back(t2);
                        } break;
                        case 2: {
                            const int cidm = midpointRingIntModulo(cid1, cid2);
                            Triangle t1;
                            Triangle t2;
                            Triangle t3;
                            t1.v[0] = vertices[tid1].g_idx;
                            t1.v[1] = vertices[tid2].g_idx;
                            t1.v[2] = tg_idx[cidm];

                            t2.v[0] = vertices[tid1].g_idx;
                            t2.v[1] = tg_idx[cidm];
                            t2.v[2] = tg_idx[cid1];

                            t3.v[0] = vertices[tid2].g_idx;
                            t3.v[1] = tg_idx[cid2];
                            t3.v[2] = tg_idx[cidm];

                            m_triangles.push_back(t1);
                            m_triangles.push_back(t2);
                            m_triangles.push_back(t3);
                            lt_.push_back(t1);
                            lt_.push_back(t2);
                            lt_.push_back(t3);
                        } break;
                    }  // switch
                }      // for loop over the vertices of the contour
            }          // if (_not_tunnel)
        }              // for loop over contours
        if (cnt_ == 1) {
            m_ccases_tunnel = m_ccases_tunnel - 1;
            m_ccases_12cont++;
            // there is a single contour
            // triangulate and close inner hexagon
            // triangle must have the correct orientation
            // use asymptotic_decider() to see if positive vertices
            // are separated, in thic case orientation must be changed
            const bool s_ = (asymptotic_decider(F[0], F[1], F[2], F[3]) <= i0);
            const bool of_ = (wi[1] < wi[0]) ? s_ : !s_;
            Triangle t1;
            Triangle t2;
            Triangle t3;
            Triangle t4;
            if (!of_) {
                t1.v[0] = tg_idx[0];
                t1.v[1] = tg_idx[2];
                t1.v[2] = tg_idx[1];
                t2.v[0] = tg_idx[2];
                t2.v[1] = tg_idx[4];
                t2.v[2] = tg_idx[3];
                t3.v[0] = tg_idx[0];
                t3.v[1] = tg_idx[5];
                t3.v[2] = tg_idx[4];
                t4.v[0] = tg_idx[0];
                t4.v[1] = tg_idx[4];
                t4.v[2] = tg_idx[2];
            } else {
                t1.v[0] = tg_idx[0];
                t1.v[1] = tg_idx[1];
                t1.v[2] = tg_idx[2];
                t2.v[0] = tg_idx[2];
                t2.v[1] = tg_idx[3];
                t2.v[2] = tg_idx[4];
                t3.v[0] = tg_idx[0];
                t3.v[1] = tg_idx[4];
                t3.v[2] = tg_idx[5];
                t4.v[0] = tg_idx[0];
                t4.v[1] = tg_idx[2];
                t4.v[2] = tg_idx[4];
            }
            m_triangles.push_back(t1);
            m_triangles.push_back(t2);
            m_triangles.push_back(t3);
            m_triangles.push_back(t4);
        }
    } else {
        // there is no tunnel
        // handle case with no saddle point as simple polygons with 3, 4, 5 or six vertices
        const unsigned char nr_u{(unsigned char)nrQSolFace(0, q_sol)};
        const unsigned char nr_v{(unsigned char)nrQSolFace(1, q_sol)};
        const unsigned char nr_w{(unsigned char)nrQSolFace(2, q_sol)};
        const unsigned char nr_t{(unsigned char)(nr_u + nr_v + nr_w)};
        if (nr_t == nr_u || nr_t == nr_v || nr_t == nr_w) {
            // loop over all contours
            for (int i = 0; i < (int)cnt_; i++) {
                switch (get_cnt_size(i, c_)) {
                    case 3: {
                        Triangle t1;
                        t1.v[0] = vertices[get_c(i, 0, c_)].g_idx;
                        t1.v[1] = vertices[get_c(i, 1, c_)].g_idx;
                        t1.v[2] = vertices[get_c(i, 2, c_)].g_idx;
                        m_triangles.push_back(t1);
                        m_ccases_3++;
                    } break;
                    case 4: {
                        Triangle t1;
                        Triangle t2;
                        t1.v[0] = vertices[get_c(i, 0, c_)].g_idx;
                        t1.v[1] = vertices[get_c(i, 1, c_)].g_idx;
                        t1.v[2] = vertices[get_c(i, 2, c_)].g_idx;
                        t2.v[0] = vertices[get_c(i, 0, c_)].g_idx;
                        t2.v[1] = vertices[get_c(i, 2, c_)].g_idx;
                        t2.v[2] = vertices[get_c(i, 3, c_)].g_idx;
                        m_triangles.push_back(t1);
                        m_triangles.push_back(t2);
                        m_ccases_4++;
                    } break;
                    case 5: {
                        Triangle t1;
                        Triangle t2;
                        Triangle t3;
                        t1.v[0] = vertices[get_c(i, 0, c_)].g_idx;
                        t1.v[1] = vertices[get_c(i, 1, c_)].g_idx;
                        t1.v[2] = vertices[get_c(i, 2, c_)].g_idx;
                        t2.v[0] = vertices[get_c(i, 0, c_)].g_idx;
                        t2.v[1] = vertices[get_c(i, 2, c_)].g_idx;
                        t2.v[2] = vertices[get_c(i, 3, c_)].g_idx;
                        t3.v[0] = vertices[get_c(i, 0, c_)].g_idx;
                        t3.v[1] = vertices[get_c(i, 3, c_)].g_idx;
                        t3.v[2] = vertices[get_c(i, 4, c_)].g_idx;
                        m_triangles.push_back(t1);
                        m_triangles.push_back(t2);
                        m_triangles.push_back(t3);
                        m_ccases_5++;
                    } break;
                    case 6: {
                        Triangle t1;
                        Triangle t2;
                        Triangle t3;
                        Triangle t4;
                        t1.v[0] = vertices[get_c(i, 0, c_)].g_idx;
                        t1.v[1] = vertices[get_c(i, 1, c_)].g_idx;
                        t1.v[2] = vertices[get_c(i, 3, c_)].g_idx;
                        t2.v[0] = vertices[get_c(i, 1, c_)].g_idx;
                        t2.v[1] = vertices[get_c(i, 2, c_)].g_idx;
                        t2.v[2] = vertices[get_c(i, 3, c_)].g_idx;
                        t3.v[0] = vertices[get_c(i, 0, c_)].g_idx;
                        t3.v[1] = vertices[get_c(i, 3, c_)].g_idx;
                        t3.v[2] = vertices[get_c(i, 4, c_)].g_idx;
                        t4.v[0] = vertices[get_c(i, 0, c_)].g_idx;
                        t4.v[1] = vertices[get_c(i, 4, c_)].g_idx;
                        t4.v[2] = vertices[get_c(i, 5, c_)].g_idx;
                        m_triangles.push_back(t1);
                        m_triangles.push_back(t2);
                        m_triangles.push_back(t3);
                        m_triangles.push_back(t4);
                        m_ccases_6++;
                    } break;
                }  // switch over size of contour
            }      // loop over contorus
        }          // thre are no saddle points
        else {
            // there are saddle points
            // fc1 = fs(1, 1)*fs(2, 1) + fs(1, 2)*fs(2, 2);
            // fc2 = fs(1, 1)*fs(3, 1) + fs(1, 2)*fs(3, 2);
            // fc3 = fs(2, 1)*fs(3, 2) + fs(2, 2)*fs(3, 1);
            unsigned char fs[3][2]{{(uchar)(q_sol & 1), (uchar)((q_sol >> 1) & 1)},
                                   {(uchar)((q_sol >> 2) & 1), (uchar)((q_sol >> 3) & 1)},
                                   {(uchar)((q_sol >> 4) & 1), (uchar)((q_sol >> 5) & 1)}};

            const unsigned char fc1 = fs[0][0] * fs[1][0] + fs[0][1] * fs[1][1];
            const unsigned char fc2 = fs[0][0] * fs[2][0] + fs[0][1] * fs[2][1];
            const unsigned char fc3 = fs[1][0] * fs[2][1] + fs[1][1] * fs[2][0];
            const unsigned char c_faces = fc1 + fc2 + fc3;
            double ucoord{};
            double vcoord{};
            double wcoord{};
            switch (c_faces) {
                case 2: {
                    if (fc1 == 0) {
                        ucoord = fs[0][0] * ui[0] + fs[0][1] * ui[1];
                        vcoord = fs[1][0] * vi[0] + fs[1][1] * vi[1];
                        wcoord = fs[1][0] * wi[1] + fs[1][1] * wi[0];
                    } else if (fc2 == 0) {
                        ucoord = fs[0][0] * ui[0] + fs[0][1] * ui[1];
                        vcoord = fs[0][0] * vi[0] + fs[0][1] * vi[1];
                        wcoord = fs[0][0] * wi[1] + fs[0][1] * wi[0];
                    } else if (fc3 == 0) {
                        ucoord = fs[1][0] * ui[0] + fs[1][1] * ui[1];
                        vcoord = fs[1][0] * vi[0] + fs[1][1] * vi[1];
                        wcoord = fs[1][0] * wi[0] + fs[1][1] * wi[1];
                    }
                } break;
                case 3: {
                    ucoord = (fs[0][0] * ui[0] + fs[0][1] * ui[1]) / (fs[0][0] + fs[0][1]);
                    vcoord = (fs[1][0] * vi[0] + fs[1][1] * vi[1]) / (fs[1][0] + fs[1][1]);
                    wcoord = (fs[2][0] * wi[0] + fs[2][1] * wi[1]) / (fs[2][0] + fs[2][1]);
                } break;
                case 4: {
                    const unsigned char nr_u = fs[0][0] + fs[0][1];
                    const unsigned char nr_v = fs[1][0] + fs[1][1];
                    const unsigned char nr_w = fs[2][0] + fs[2][1];
                    if (nr_w == 1) {
                        ucoord = fs[2][0] * ui[0] + fs[2][1] * ui[1];
                        vcoord = fs[2][1] * vi[0] + fs[2][0] * vi[1];
                        wcoord = fs[2][0] * wi[0] + fs[2][1] * wi[1];
                    } else if (nr_v == 1) {
                        ucoord = fs[1][0] * ui[0] + fs[1][1] * ui[1];
                        vcoord = fs[1][0] * vi[0] + fs[1][1] * vi[1];
                        wcoord = fs[1][1] * wi[0] + fs[1][0] * wi[1];
                    } else if (nr_u == 1) {
                        ucoord = fs[0][0] * ui[0] + fs[0][1] * ui[1];
                        vcoord = fs[0][0] * vi[0] + fs[0][1] * vi[1];
                        wcoord = fs[0][0] * wi[0] + fs[0][1] * wi[1];
                    }
                } break;
            }  // switch(c_faces)

            // create inner vertex
            Point ip;
            Normal in;
            ip[0] = (1 - wcoord) * ((1 - vcoord) * (p[0][0] + ucoord * (p[1][0] - p[0][0])) +
                                    vcoord * (p[2][0] + ucoord * (p[3][0] - p[2][0]))) +
                    wcoord * ((1 - vcoord) * (p[4][0] + ucoord * (p[5][0] - p[4][0])) +
                              vcoord * (p[6][0] + ucoord * (p[7][0] - p[6][0])));
            ip[1] = (1 - wcoord) * ((1 - vcoord) * (p[0][1] + ucoord * (p[1][1] - p[0][1])) +
                                    vcoord * (p[2][1] + ucoord * (p[3][1] - p[2][1]))) +
                    wcoord * ((1 - vcoord) * (p[4][1] + ucoord * (p[5][1] - p[4][1])) +
                              vcoord * (p[6][1] + ucoord * (p[7][1] - p[6][1])));
            ip[2] = (1 - wcoord) * ((1 - vcoord) * (p[0][2] + ucoord * (p[1][2] - p[0][2])) +
                                    vcoord * (p[2][2] + ucoord * (p[3][2] - p[2][2]))) +
                    wcoord * ((1 - vcoord) * (p[4][2] + ucoord * (p[5][2] - p[4][2])) +
                              vcoord * (p[6][2] + ucoord * (p[7][2] - p[6][2])));
            in[0] = (1 - wcoord) * ((1 - vcoord) * (n[0][0] + ucoord * (n[1][0] - n[0][0])) +
                                    vcoord * (n[2][0] + ucoord * (n[3][0] - n[2][0]))) +
                    wcoord * ((1 - vcoord) * (n[4][0] + ucoord * (n[5][0] - n[4][0])) +
                              vcoord * (n[6][0] + ucoord * (n[7][0] - n[6][0])));
            in[1] = (1 - wcoord) * ((1 - vcoord) * (n[0][1] + ucoord * (n[1][1] - n[0][1])) +
                                    vcoord * (n[2][1] + ucoord * (n[3][1] - n[2][1]))) +
                    wcoord * ((1 - vcoord) * (n[4][1] + ucoord * (n[5][1] - n[4][1])) +
                              vcoord * (n[6][1] + ucoord * (n[7][1] - n[6][1])));
            in[2] = (1 - wcoord) * ((1 - vcoord) * (n[0][2] + ucoord * (n[1][2] - n[0][2])) +
                                    vcoord * (n[2][2] + ucoord * (n[3][2] - n[2][2]))) +
                    wcoord * ((1 - vcoord) * (n[4][2] + ucoord * (n[5][2] - n[4][2])) +
                              vcoord * (n[6][2] + ucoord * (n[7][2] - n[6][2])));
            // normalize normal
            const double factor = std::sqrt(in[0] * in[0] + in[1] * in[1] + in[2] * in[2]);
            in[0] = in[0] / factor;
            in[1] = in[1] / factor;
            in[2] = in[2] / factor;
            const uint g_index = (uint)m_points.size();
            // loop over the contorus
            bool pt_used = false;
            for (int i = 0; i < (int)cnt_; i++) {
                const unsigned char cnt_sz = (unsigned char)get_cnt_size(i, c_);
                if (cnt_sz == 3) {
                    Triangle t1;
                    t1.v[0] = vertices[get_c(i, 0, c_)].g_idx;
                    t1.v[1] = vertices[get_c(i, 1, c_)].g_idx;
                    t1.v[2] = vertices[get_c(i, 2, c_)].g_idx;
                    m_triangles.push_back(t1);
                } else {
                    pt_used = true;
                    for (int t = 0; t < cnt_sz; t++) {
                        Triangle ts;
                        ts.v[0] = vertices[get_c(i, t, c_)].g_idx;
                        ts.v[1] = vertices[get_c(i, (t + 1) % cnt_sz, c_)].g_idx;
                        ts.v[2] = g_index;
                        m_triangles.push_back(ts);
                    }
                }

                switch (cnt_sz) {
                    case 3:
                        m_ccases_3++;
                        break;
                    case 4:
                        m_ccases_4++;
                        break;
                    case 5:
                        m_ccases_5++;
                        break;
                    case 6:
                        m_ccases_6a++;
                        break;
                    case 7:
                        m_ccases_7++;
                        break;
                    case 8:
                        m_ccases_8++;
                        break;
                    case 9:
                        m_ccases_9++;
                        break;
                    default:
                        break;
                }
            }
            if (pt_used) {
                m_points.push_back(ip);
                m_pnorms.push_back(in);
            }
        }  // else - there are saddle points
    }

}

class TMC {

void t_mc(const double i0) {
    // edges are uniquely characterized by the two end vertices, which have a unique vertex id
    // the end vertices of the edge are computed in the cell by giving the indices (i,j,k).
    // These indices are obtained from the cell index by adding 0 or 1 to i, j or k respectively
    // Example: edge 0: (i,j,k) - (i+1,j,k)
    //          edge 1: (i+1,j,k) - (i+1,j+1,k)
    // The first 3 indices are for the first vertex and the second 3 for the second vertex.
    // there are 12 edges, assign to each vertex three edges, the global edge numbering
    // consist of 3*global_vertex_id + edge_offset.
    const int global_edge_id[][4] = {{0, 0, 0, 0}, {1, 0, 0, 1}, {0, 1, 0, 0}, {0, 0, 0, 1},
                                     {0, 0, 1, 0}, {1, 0, 1, 1}, {0, 1, 1, 0}, {0, 0, 1, 1},
                                     {0, 0, 0, 2}, {1, 0, 0, 2}, {1, 1, 0, 2}, {0, 1, 0, 2}};
    // the end vertices of an edge
    int l_edges[12][2] = {{0, 1}, {1, 3}, {2, 3}, {0, 2}, {4, 5}, {5, 7},
                          {6, 7}, {4, 6}, {0, 4}, {1, 5}, {3, 7}, {2, 6}};
    // compute sizes
    const int nx = m_ugrid.x_size();
    const int ny = m_ugrid.y_size();
    const int nz = m_ugrid.z_size();

    // we need to compute up to 3 vertices at the interior of a cell, therefore
    // the cell shift factor is set to 3+3 = 6, i.e. 3 edges assigned to a cell for global numberig
    // and 3 vertices in the interior of the cell
    m_cell_shift_factor = 6;

    // there can be at most 12 intersections
    std::vector<Vertex> vertices(12);
    std::vector<Point> ip(12);
    std::vector<Normal> in(12);

    int i_case_count = 0;
    // marching cubes
    for (int k = 0; k < (nz - 1); k++) {
        m_kindex = k;
        for (int j = 0; j < (ny - 1); j++) {
            m_jindex = j;
            for (int i = 0; i < (nx - 1); i++) {
                m_iindex = i;
                // slice hex
                // collect function values and build index
                double u[8];
                Point p[8];
                Normal n[8];
                int vi{0};
                std::bitset<8> index = 0;
                for (int kl = 0; kl <= 1; kl++) {
                    for (int jl = 0; jl <= 1; jl++) {
                        for (int il = 0; il <= 1; il++) {
                            // collect scalar values and computex index
                            p[vi] = m_ugrid.point(i + il, j + jl, k + kl);
                            u[vi] = m_ugrid.scalar(i + il, j + jl, k + kl);
                            if (u[vi] >= i0) {
                                // index.set(VertexMapping[vi]);
                                index.set(vi);
                            }
                            // probably better get normals here
                            n[vi] = m_ugrid.normal(i + il, j + jl, k + kl);
                            // next cell vertex
                            vi++;
                        }
                    }
                }

                // collect edges from table and
                // interpolate triangle vertex positon
                int i_case = int(index.to_ullong());
                // int tcm = (int)t_entry[i_case];
                int tcm = (int)t_ambig[i_case];
                if (tcm == 105) {
                    i_case_count++;
                    // t_slice(i,j,k,i0,u,p,n,i_case);
                    p_slice(i, j, k, i0, u, p, n, i_case);
                } else {
                    // compute for this case the vertices
                    ushort flag = 1;
                    for (int eg = 0; eg < 12; eg++) {
                        if (flag & e_pattern[i_case]) {
                            // the edge global index is given by the vertex global index + the edge offset
                            const int ix = i + global_edge_id[eg][0];
                            const int iy = j + global_edge_id[eg][1];
                            const int iz = k + global_edge_id[eg][2];
                            vertices[eg].g_edg =
                                uint(m_cell_shift_factor * m_ugrid.global_index(ix, iy, iz) + global_edge_id[eg][3]);
                            // generate vertex here, do not care at this point if vertex already exist
                            int* vert = l_edges[eg];
                            // interpolation weight
                            const int v0 = vert[0];
                            const int v1 = vert[1];
                            double l = (i0 - u[v0]) / (u[v1] - u[v0]);
                            // interpolate vertex
                            ip[eg][0] = (1 - l) * p[v0][0] + l * p[v1][0];
                            ip[eg][1] = (1 - l) * p[v0][1] + l * p[v1][1];
                            ip[eg][2] = (1 - l) * p[v0][2] + l * p[v1][2];

                            // interpolate normal
                            in[eg][0] = (1 - l) * n[v0][0] + l * n[v1][0];
                            in[eg][1] = (1 - l) * n[v0][1] + l * n[v1][1];
                            in[eg][2] = (1 - l) * n[v0][2] + l * n[v1][2];
                            const double nlength =
                                std::sqrt(in[eg][0] * in[eg][0] + in[eg][1] * in[eg][1] + in[eg][2] * in[eg][2]);
                            in[eg][0] = in[eg][0] / nlength;
                            in[eg][1] = in[eg][1] / nlength;
                            in[eg][2] = in[eg][2] / nlength;

                            // set vertex index
                            auto s_index = m_vertices.find(vertices[eg].g_edg);
                            if (s_index == m_vertices.end()) {
                                // index not found!
                                const int g_idx = (int)m_points.size();
                                m_vertices[vertices[eg].g_edg] = g_idx;
                                vertices[eg].g_idx = g_idx;
                                m_points.push_back(ip[eg]);
                                m_pnorms.push_back(in[eg]);
                            } else {
                                vertices[eg].g_idx = s_index->second;  // this is vertex global index g_idx
                            }
                        }
                        flag <<= 1;
                    }

                    // construct triangles
                    for (int t = 0; t < 16; t += 3) {
                        const int t_index = i_case * 16 + t;
                        // if (e_tris_list[t_index] == 0x7f)
                        if (t_pattern[t_index] == -1) break;
                        Triangle tri;
                        //                        const int eg0 = e_tris_list[t_index];
                        //                        const int eg1 = e_tris_list[t_index+1];
                        //                        const int eg2 = e_tris_list[t_index+2];
                        const int eg0 = t_pattern[t_index];
                        const int eg1 = t_pattern[t_index + 1];
                        const int eg2 = t_pattern[t_index + 2];
                        // change MC triangle orientation
                        // positive vertices must be outside the surface
                        tri.v[0] = (int)vertices[eg2].g_idx;
                        tri.v[1] = (int)vertices[eg1].g_idx;
                        tri.v[2] = (int)vertices[eg0].g_idx;

                        // insert new triangle in list
                        m_triangles.push_back(tri);
                    }
                }
            }
        }
    }
}
};

}  // namespace internal
}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_TMC_INTERNAL_TMC_H