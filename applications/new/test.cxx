#include <Eigen/CXX11/Tensor>
#include "numerics.hxx"
#include "utility.hxx"

int
main(int argc, char** argv)
{
    Eigen::TensorFixedSize<double, Eigen::Sizes<4, 4, 2>> tmp;
    Eigen::TensorFixedSize<double, Eigen::Sizes<3>> pnt;
    Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3>> mat1;
    Eigen::Matrix<double, 3, 3> mat;
    auto res = mat1 * pnt;
}