//
// Created by dsiedel on 01/04/2022.
//

#ifndef FETA_FETA_CORE_GEOMETRY_GEOMETRY_HXX
#define FETA_FETA_CORE_GEOMETRY_GEOMETRY_HXX

#include "new/_feta_matrix.hxx"

namespace feta::core::geometry
{

    template<typename M>
    inline
    auto
    getBarycenter(
            M const &
            point_args
    )
    {
        auto const constexpr rows = M::CompileTimeTraits::RowsAtCompileTime;
        auto const constexpr cols = M::CompileTimeTraits::ColsAtCompileTime;
        auto barycenter = Matrix<Real, rows>();
        for (Indx i = 0; i < cols; ++i) {
            barycenter += point_args.col(i);
        }
        barycenter /= Real(cols);
        return barycenter;
    }

    template<typename M>
    inline
    auto
    getNormalVector(
            M const &
            vector_args
    )
    {
        auto const constexpr rows = M::CompileTimeTraits::RowsAtCompileTime;
        auto const constexpr cols = M::CompileTimeTraits::ColsAtCompileTime;
        if constexpr (rows == 1 && cols == 1) {
            return vector_args;
        }
        else if constexpr (rows == 2 && cols == 1) {
            auto norm = vector_args.norm();
            return Matrix<Real, 2>{vector_args(1)/norm, -vector_args(0)/norm};
        }
        else if constexpr (rows == 3 && cols == 2) {
            auto v0 = vector_args.col(0);
            auto v1 = vector_args.col(1);
            v0 /= v0.norm();
            v1 /= v1.norm();
            return v0.cross(v1);
        }
        else {
            assert(false);
        }
    }


    template<typename M>
    inline
    auto
    getRotationMatrix(
            M const &
            point_args
    )
    {
        auto const constexpr rows = M::CompileTimeTraits::RowsAtCompileTime;
        auto const constexpr cols = M::CompileTimeTraits::ColsAtCompileTime;
        if constexpr (rows == 1 && cols == 1) {
            return point_args;
        }
        else if constexpr (rows == 2 && cols == 2) {
            auto rotation_matrix = Matrix<Real, 2, 2>();
            auto edge = point_args.col(1) - point_args.col(0);
            edge /= edge.norm();
            rotation_matrix(0, 0) = edge(0);
            rotation_matrix(0, 1) = edge(1);
            rotation_matrix(1, 0) = edge(1);
            rotation_matrix(1, 1) = -edge(0);
            return rotation_matrix;
        }
        else if constexpr (rows == 3 && cols == 3) {
            auto rotation_matrix = Matrix<Real, 3, 3>();
            auto x_axis = point_args.col(2) - point_args.col(0);
            auto vector = point_args.col(1) - point_args.col(0);
            auto z_axis = x_axis.cross(vector);
            auto y_axis = z_axis.cross(x_axis);
            x_axis /= x_axis.norm();
            y_axis /= y_axis.norm();
            z_axis /= z_axis.norm();
            rotation_matrix.row(0) = x_axis;
            rotation_matrix.row(1) = y_axis;
            rotation_matrix.row(2) = z_axis;
            return rotation_matrix;
        }
        else {
            assert(false);
        }

    }

}

#endif //FETA_FETA_CORE_GEOMETRY_GEOMETRY_HXX
