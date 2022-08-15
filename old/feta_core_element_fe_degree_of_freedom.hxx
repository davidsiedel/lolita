//
// Created by dsiedel on 06/04/2022.
//

#ifndef FETA_FETA_CORE_ELEMENT_FE_DEGREE_OF_FREEDOM_HXX
#define FETA_FETA_CORE_ELEMENT_FE_DEGREE_OF_FREEDOM_HXX

#include "new/feta_core_element_element_description.hxx"
#include "new/feta_core_field_description.hxx"
#include "new/feta_core_element_basis.hxx"
#include "new/unique_pointer.hxx"

namespace feta::core::element::finite_element
{

    template<ElementDescription E, BasisDescription B>
    struct Binding
    {

        auto const static constexpr dim_unknowns = FiniteElementBasis<E, B>::dim_basis;

        Binding()
        :
        indices(),
        values(Matrix<Real, dim_unknowns>::Zero())
        {}

        explicit
        Binding(
                auto &
                index_arg
        )
        :
        indices(setIndex(index_arg)),
        values(Matrix<Real, dim_unknowns>::Zero())
        {}

        Matrix<Real, dim_unknowns> values;

        Array<Indx, dim_unknowns> indices;

    private:

        static
        auto
        setIndex(
                auto &
                binding_index_arg
        )
        {
            auto indices_arg = Array<Indx, dim_unknowns>();
            for (int i = 0; i < dim_unknowns; ++i) {
                indices_arg.get(i) = binding_index_arg;
                binding_index_arg += 1;
            }
            return indices_arg;
        }

    };

    template<ElementDescription E, BasisDescription B, FieldDescription F>
    struct Unknowns
    {

        auto const static constexpr dim_unknowns = FiniteElementBasis<E, B>::dim_basis;

        auto const static constexpr num_unknowns = FiniteElementBasis<E, B>::dim_basis * F.tensor_description.size;

        Unknowns()
        :
        indices(),
        values(Matrix<Real, num_unknowns>::Zero())
        {}

        explicit
        Unknowns(
                auto &
                unknown_index_arg
        )
        :
        indices(setIndex(unknown_index_arg)),
        values(Matrix<Real, num_unknowns>::Zero())
        {}

        auto
        getComponentValues(
                auto
                row_arg,
                auto
                col_arg
        )
        const
        {
            auto p = Indx(row_arg * F.tensor_description.cols + col_arg) * dim_unknowns;
            auto c = Matrix<Matrix<Real, dim_unknowns> const>(values.data() + p);
            return c;
        }

        auto
        getComponentValues(
                auto
                row_arg,
                auto
                col_arg
        )
        {
            auto p = Indx(row_arg * F.tensor_description.cols + col_arg) * dim_unknowns;
            auto c = Matrix<Matrix<Real, dim_unknowns>>(values.data() + p);
            return c;
        }

        auto
        getComponentIndices(
                auto
                row_arg,
                auto
                col_arg
        )
        const
        {
            auto p = Indx(row_arg * F.tensor_description.cols + col_arg) * dim_unknowns;
            auto c = Matrix<Matrix<Indx, dim_unknowns> const>(indices.data.data() + p);
            return c;
        }

        auto
        getComponentIndices(
                auto
                row_arg,
                auto
                col_arg
        )
        {
            auto p = Indx(row_arg * F.tensor_description.cols + col_arg) * dim_unknowns;
            auto c = Matrix<Matrix<Indx, dim_unknowns>>(indices.data.data() + p);
            return c;
        }

        Matrix<Real, num_unknowns> values;

        Array<Indx, num_unknowns> indices;

    private:

        static
        auto
        setIndex(
                auto &
                unknown_index_arg
        )
        {
            auto indices_arg = Array<Indx, num_unknowns>();
            for (int i = 0; i < num_unknowns; ++i) {
                indices_arg.get(i) = unknown_index_arg;
                unknown_index_arg += 1;
            }
            return indices_arg;
        }

    };

    template<ElementDescription E, FieldDescription F, BasisDescription B>
    struct DegreesOfFreedom
    {

    private:

        auto const static constexpr rows = F.tensor_description.rows;

        auto const static constexpr cols = F.tensor_description.cols;

    public:

        auto const static constexpr num_unknowns = Unknowns<E, B, F>::num_unknowns;

        auto const static constexpr dim_unknowns = Unknowns<E, B, F>::dim_unknowns;

        DegreesOfFreedom()
        :
        unknowns(),
        bindings()
        {}

        explicit
        DegreesOfFreedom(
                auto &
                unknown_index_arg
        )
        :
        unknowns(unknown_index_arg),
        bindings()
        {}

        DegreesOfFreedom(
                auto &
                unknown_index_arg,
                auto &
                binding_index_arg,
                auto const &
                directions_arg
        )
        :
        unknowns(unknown_index_arg),
        bindings(setBindings(binding_index_arg, directions_arg))
        {}

    private:

        static
        auto
        setBindings(
                auto &
                binding_index_arg,
                auto const &
                directions_arg
        )
        {
            auto bindings_arg = Array<UniquePointer<Binding<E, B>>, rows, cols>();
            for (auto const & dir: directions_arg) {
                bindings_arg.get(dir.row, dir.col) = UniquePointer<Binding<E, B>>(binding_index_arg);
            }
            return bindings_arg;
        }

    public:

        Unknowns<E, B, F> unknowns;

        Array<UniquePointer<Binding<E, B>>, rows, cols> bindings;

    };

//    template<ElementDescription E, BasisDescription B>
//    struct UnknownComponent
//    {
//
//        auto const static constexpr num_unknowns = FiniteElementBasis<E, B>::dim_basis;
//
//        UnknownComponent()
//                :
//                ptr_unknown_coefficients(),
//                ptr_binding_coefficients(),
//                ptr_unknown_indices(),
//                ptr_binding_indices()
//        {}
//
//        explicit
//        UnknownComponent(
//                auto &
//                unknown_index
//        )
//                :
//                ptr_unknown_coefficients(Matrix<Real, FiniteElementBasis<E, B>::dim_basis>::Zero()),
//                ptr_unknown_indices(setIndicesArray(unknown_index)),
//                ptr_binding_coefficients(),
//                ptr_binding_indices()
//        {}
//
//        UnknownComponent(
//                auto &
//                unknown_index,
//                auto &
//                binding_index
//        )
//                :
//                ptr_unknown_coefficients(Matrix<Real, FiniteElementBasis<E, B>::dim_basis>::Zero()),
//                ptr_unknown_indices(setIndicesArray(unknown_index)),
//                ptr_binding_coefficients(Matrix<Real, FiniteElementBasis<E, B>::dim_basis>::Zero()),
//                ptr_binding_indices(setIndicesArray(binding_index))
//        {}
//
//        Bool
//        operator==(
//                UnknownComponent const &
//                other
//        )
//        const
//        {
//            auto eq_0 = ptr_unknown_coefficients == other.ptr_unknown_coefficients;
//            auto eq_1 = ptr_unknown_indices == other.ptr_unknown_indices;
//            auto eq_2 = ptr_binding_coefficients == other.ptr_binding_coefficients;
//            auto eq_3 = ptr_binding_indices == other.ptr_binding_indices;
//            return eq_0 && eq_1 && eq_2 && eq_3;
//        }
//
//        Bool
//        operator!=(
//                UnknownComponent const &
//                other
//        )
//        const
//        {
//            return !(other == * this);
//        }
//
//        auto
//        isBound()
//        const
//        {
//            if (ptr_binding_indices.data == nullptr) {
//                return false;
//            }
//            else {
//                return true;
//            }
//        }
//
//    private:
//
//        static
//        Array<Indx, FiniteElementBasis<E, B>::dim_basis>
//        setIndicesArray(
//                auto &
//                index_arg
//        )
//        {
//            auto indices = Array<Indx, FiniteElementBasis<E, B>::dim_basis>();
//            for (Indx i = 0; i < FiniteElementBasis<E, B>::dim_basis; ++i) {
//                indices(i) = static_cast<Indx>(index_arg);
//                index_arg += 1;
//            }
//            return indices;
//        }
//
//    public:
//
//        UniquePointer<Matrix<Real, FiniteElementBasis<E, B>::dim_basis>> ptr_unknown_coefficients;
//
//        UniquePointer<Matrix<Real, FiniteElementBasis<E, B>::dim_basis>> ptr_binding_coefficients;
//
//        UniquePointer<Array<Indx, FiniteElementBasis<E, B>::dim_basis>> ptr_unknown_indices;
//
//        UniquePointer<Array<Indx, FiniteElementBasis<E, B>::dim_basis>> ptr_binding_indices;
//
//    };

//    template<FieldDescription F, ElementDescription E, BasisDescription B>
//    struct Unknowns : public Array<UnknownComponent<E, B>, F.tensor_description.rows, F.tensor_description.cols>
//    {
//
//    private:
//
//        auto const static constexpr rows = F.tensor_description.rows;
//
//        auto const static constexpr cols = F.tensor_description.cols;
//
//        using Base = Array<UnknownComponent<E, B>, rows, cols>;
//
//        using UC = UnknownComponent<E, B>;
//
//    public:
//
//        auto const static constexpr num_unknowns = UnknownComponent<E, B>::num_unknowns * Base::size();
//
//        Unknowns()
//        :
//        Base()
//        {}
//
//        explicit
//        Unknowns(
//                auto &
//                unknown_index
//        )
//        :
//        Base(setComponents(unknown_index))
//        {}
//
//        Unknowns(
//                auto &
//                unknown_index_arg,
//                auto &
//                binding_index_arg
//        )
//        :
//        Base(setComponents(unknown_index_arg, binding_index_arg))
//        {}
//
//        auto
//        getUnknowns()
//        const
//        {
//            auto u = Matrix<Real, num_unknowns>();
//            auto c = Indx(0);
//            for (int i = 0; i < rows; ++i) {
//                for (int j = 0; j < cols; ++j) {
//                    u.template segment<UC::num_unknowns>(c) = this->get(i, j).ptr_unknown_coefficients.get();
//                    c += UC::num_unknowns;
//                }
//            }
//            return u;
//        }
//
//    private:
//
//        static
//        auto
//        setComponents(
//                auto &
//                unknown_index_arg,
//                auto &
//                binding_index_arg
//        )
//        {
//            auto data = Array<UnknownComponent<E, B>, rows, cols>();
//            for (int i = 0; i < rows; ++i) {
//                for (int j = 0; j < cols; ++j) {
//                    data.get(i, j) = UnknownComponent<E, B>(unknown_index_arg, binding_index_arg);
//                }
//            }
//            return data;
//        }
//
//        static
//        auto
//        setComponents(
//                auto &
//                unknown_index_arg
//        )
//        {
//            auto data = Array<UnknownComponent<E, B>, rows, cols>();
//            for (int i = 0; i < rows; ++i) {
//                for (int j = 0; j < cols; ++j) {
//                    data.get(i, j) = UnknownComponent<E, B>(unknown_index_arg);
//                }
//            }
//            return data;
//        }
//
//    };

}

#endif //FETA_FETA_CORE_ELEMENT_FE_DEGREE_OF_FREEDOM_HXX
