#ifndef CGAL_CARTESIAN_GRID_3_H
#define CGAL_CARTESIAN_GRID_3_H

#include <type_traits>

#include <CGAL/Image_3.h>

namespace CGAL {

// TODO: not sure if anything other than float works
template<typename T>
class Cartesian_grid_3 {
public:
    typedef T value_type;

public:
    Cartesian_grid_3(const std::size_t xdim, const std::size_t ydim, const std::size_t zdim) {
        create_image(xdim, ydim, zdim);
    }

    value_type value(const std::size_t i, const std::size_t j, const std::size_t k) const {
        return grid.value(i, j, k);
    }

    value_type& value(const std::size_t i, const std::size_t j, const std::size_t k) {
        value_type* data = (value_type*) grid.image()->data;
        return data[(k * ydim() + j) * xdim() + i];
    }

    std::size_t xdim() const { return grid.xdim(); }
    std::size_t ydim() const { return grid.ydim(); }
    std::size_t zdim() const { return grid.zdim(); }

    // TODO: better return types
    double voxel_x() const { return grid.vx(); }
    double voxel_y() const { return grid.vy(); }
    double voxel_z() const { return grid.vz(); }

    float offset_x() const { return grid.tx(); }
    float offset_y() const { return grid.ty(); }
    float offset_z() const { return grid.tz(); }

private:
    void create_image(const std::size_t xdim, const std::size_t ydim, const std::size_t zdim);

private:
    Image_3 grid;
};


template<typename T>
void Cartesian_grid_3<T>::create_image(const std::size_t xdim, const std::size_t ydim, const std::size_t zdim) {

    WORD_KIND wordkind;
    if (std::is_floating_point<value_type>::value)
        wordkind = WK_FLOAT;
    else
        wordkind = WK_FIXED;

    SIGN sign;
    if (std::is_signed<value_type>::value)
        sign = SGN_SIGNED;
    else
        sign = SGN_UNSIGNED;
    
    _image *im = _createImage(xdim, ydim, zdim,
                    1,                      //vectorial dimension
                    1, 1, 1,                //voxel size
                    sizeof(value_type),     //image word size in bytes
                    wordkind,               //image word kind WK_FIXED, WK_FLOAT, WK_UNKNOWN
                    sign);                  //image word sign
    
    if (im == nullptr || im->data == nullptr) {
        throw std::bad_alloc();  // TODO: idk?
    }
    
    grid = Image_3(im, Image_3::OWN_THE_DATA);
}

} // end namespace CGAL

#endif // CGAL_CARTESIAN_GRID_3_H