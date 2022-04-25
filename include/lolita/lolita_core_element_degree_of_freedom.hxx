//
// Created by dsiedel on 16/04/22.
//

#ifndef LOLITA_LOLITA_CORE_ELEMENT_DEGREE_OF_FREEDOM_HXX
#define LOLITA_LOLITA_CORE_ELEMENT_DEGREE_OF_FREEDOM_HXX

#include "lolita/lolita_core_element_geometry.hxx"
#include "lolita/lolita_core_element_basis.hxx"
#include "lolita/lolita_pointers.hxx"

namespace lolita::core::element
{

    struct DOF
    {

        DOF() = default;

        DOF(
                Indx
                a,
                Indx &
                b
        )
        :
        coefficients(a),
        index(b)
        {
            b += a;
        }

        Vector<Real> coefficients;

        Indx index;

    };

    template<FField F>
    struct DegOfFre
    {

//        DegOfFre() = default;
//
//        DegOfFre(
//                Indx
//                a,
//                Indx
//                b
//        )
//        :
//        unknowns(a, b)
//        {}
//
//        template<Indx I>
//        auto
//        getVec()
//        const
//        {
//            auto vec = Vector<Real, F.size() * I>();
//            auto count = Indx(0);
//            for (int i = 0; i < F.rows(); ++i) {
//                for (int j = 0; j < F.cols(); ++j) {
//                    vec.template segment<I>(count) = unknowns.get(i, j).get().template segment<I>(0);
//                    count += I;
//                }
//            }
//            return vec;
//        }
//
//        template<Indx I>
//        auto
//        getVec2()
//        const
//        {
//            return Vector<Vector<Real, F.size() * I> const>(us.data());
//        }
//
//        template<Indx I>
//        auto
//        getVec2(Indx a, Indx b)
//        const
//        {
//            return Vector<Vector<Real, F.size() * I> const>(us.data());
//        }

//        Array<UniquePointer<DOF>, F.rows(), F.cols()> unknowns;
//        Array<UniquePointer<DOF>, F.rows(), F.cols()> bindings;

        Vector<Real> us;
        Vector<Real> bs;
        Array<Pair<Intg>, F.rows(), F.cols()> usi;
        Array<Pair<Intg>, F.rows(), F.cols()> bsi;

    };

    template<FField F>
    struct DeOfFr
    {

        DOF unknowns;
        Array<DOF, F.rows(), F.cols()> bindings;


    };

    template<Element E, BasisDescription B>
    struct Binding
    {

        auto const static constexpr dim_unknowns = FiniteElementBasis<E, B>::dim_basis;

        Binding()
        :
        indices(),
        values(Vector<Real, dim_unknowns>::Zero())
        {}

        explicit
        Binding(
                auto &
                index_arg
        )
        :
        indices(setIndex(index_arg)),
        values(Vector<Real, dim_unknowns>::Zero())
        {}

        Vector<Real, dim_unknowns> values;

        Vector<Indx, dim_unknowns> indices;

    private:

        static
        auto
        setIndex(
                auto &
                binding_index_arg
        )
        {
            auto indices_arg = Vector<Indx, dim_unknowns>();
            for (int i = 0; i < dim_unknowns; ++i) {
                indices_arg(i) = binding_index_arg;
                binding_index_arg += 1;
            }
            return indices_arg;
        }

    };

    template<Element E, BasisDescription B, FField F>
    struct Unknowns
    {

        auto const static constexpr dim_unknowns = FiniteElementBasis<E, B>::dim_basis;

        auto const static constexpr num_unknowns = FiniteElementBasis<E, B>::dim_basis * F.size();

        Unknowns()
        :
        indices(),
        values(Vector<Real, num_unknowns>::Zero())
        {}

        explicit
        Unknowns(
                auto &
                unknown_index_arg
        )
        :
        indices(setIndex(unknown_index_arg)),
        values(Vector<Real, num_unknowns>::Zero())
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
            auto p = Indx(row_arg * F.cols() + col_arg) * dim_unknowns;
            auto c = Vector<Vector<Real, dim_unknowns> const>(values.data() + p);
            return c;
        }

        auto
        getComponentValues(
                auto
                row_arg,
                auto
                col_arg
        ) {
            auto p = Indx(row_arg * F.cols() + col_arg) * dim_unknowns;
            auto c = Vector<Vector<Real, dim_unknowns>>(values.data() + p);
            return c;
        }

        auto
        getComponentIndices(
                auto
                row_arg,
                auto
                col_arg
        )
        const {
            auto p = Indx(row_arg * F.cols() + col_arg) * dim_unknowns;
            auto c = Vector<Vector<Indx, dim_unknowns> const>(indices.data() + p);
            return c;
        }

        auto
        getComponentIndices(
                auto
                row_arg,
                auto
                col_arg
        ) {
            auto p = Indx(row_arg * F.cols() + col_arg) * dim_unknowns;
            auto c = Vector<Vector<Indx, dim_unknowns>>(indices.data() + p);
            return c;
        }

        Vector<Real, num_unknowns> values;

        Vector<Indx, num_unknowns> indices;

    private:

        static
        auto
        setIndex(
                auto &
                unknown_index_arg
        ) {
            auto indices_arg = Vector<Indx, num_unknowns>();
            for (int i = 0; i < num_unknowns; ++i) {
                indices_arg(i) = unknown_index_arg;
                unknown_index_arg += 1;
            }
            return indices_arg;
        }

    };

    template<Element E, FField F, BasisDescription B>
    struct DegreesOfFreedom
    {

    private:

        auto const static constexpr rows = F.rows();

        auto const static constexpr cols = F.cols();

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

}

#endif //LOLITA_LOLITA_CORE_ELEMENT_DEGREE_OF_FREEDOM_HXX
