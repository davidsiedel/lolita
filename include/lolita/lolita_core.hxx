//
// Created by dsiedel on 14/04/2022.
//

#ifndef LOLITA_LOLITA_CORE_HXX
#define LOLITA_LOLITA_CORE_HXX

#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/Behaviour.h>
#include <MGIS/Behaviour/MaterialStateManager.h>
#include <MGIS/Behaviour/MaterialStateManager.hxx>
#include <MGIS/Behaviour/Integrate.hxx>
#include <MGIS/Behaviour/Integrate.h>

#include "lolita.hxx"
#include "lolita_matrix.hxx"
#include "lolita_pointers.hxx"

namespace lolita::core
{

    enum struct HypothesisType
    {

        PlaneStrain,
        PlaneStress,
        AxySymmetric,
        Volumetric

    };

    enum struct StrainType
    {

        LargeStrain,
        SmallStrain

    };

    enum struct OutputResult
    {

        BehaviourLawIntegrationFailure,
        BehaviourLawIntegrationSuccess

    };

    enum struct MeshFormatType
    {

        Gmsh,

    };

    enum struct LoadType
    {

        Natural,
        Constraint,

    };

    enum struct EuclideanFrame
    {

        Cartesian,
        AxiSymmetric

    };

    enum struct Model
    {

        Solid,
        Shell,
        Beam,

    };

    enum struct QuadratureRule
    {

        Gauss

    };

    enum struct BasisName
    {

        Monomial,
        Lagrange

    };

    enum struct MappingOperator
    {

        Gradient,
        SymmetricGradient,
        Divergence,
        SmallStrain,
        LargeStrain,
        Identity,

    };

    enum struct FiniteElementMethod
    {

        Lagrange,
        HybridHighOrder

    };

    auto const static constexpr fem_hho = FiniteElementMethod::HybridHighOrder;
    auto const static constexpr fem_lag = FiniteElementMethod::Lagrange;

    struct Basis
    {

        constexpr
        Basis()
                :
                basis(),
                ord()
        {}

        constexpr
        Basis(
                BasisName
                basis_arg,
                Indx
                order_arg
        )
                :
                basis(basis_arg),
                ord(order_arg)
        {}

        constexpr
        Bool
        operator==(
                Basis const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                Basis const &
                other
        )
        const = default;

        BasisName basis;

        Indx ord;

    };

    struct Quadrature
    {

        constexpr
        Quadrature(
                QuadratureRule
                quadrature_rule_arg,
                Indx
                ord_arg
        )
        :
        quadrature(quadrature_rule_arg),
        ord(ord_arg)
        {}

        constexpr
        Bool
        operator==(
                Quadrature const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                Quadrature const &
                other
        )
        const = default;

        QuadratureRule quadrature;

        Indx ord;

    };

    struct Domain
    {

        constexpr
        Domain() = default;

        constexpr
        Domain(
                Indx
                dim_arg,
                EuclideanFrame
                euclidean_frame_arg
        )
        :
        dim(dim_arg),
        frame(euclidean_frame_arg)
        {}

        constexpr
        Bool
        operator==(
                Domain const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                Domain const &
                other
        )
        const = default;

        constexpr
        auto
        ordIntegration(
                auto
                ord_arg
        )
        const
        {
            return frame == EuclideanFrame::AxiSymmetric ? Indx(2 * ord_arg + 1) : Indx(2 * ord_arg);
        }

        Indx dim;

        EuclideanFrame frame;

    };

    struct Field
    {

        constexpr
        Field() = default;

        constexpr
        Field(
                Indx
                ord_arg,
                Indx
                dim_arg
        )
        :
        ord(ord_arg),
        dim(dim_arg)
        {}

        static constexpr
        auto
        fromMapping(
                Field const &
                field_description_arg,
                MappingOperator
                operator_type_arg
        )
        {
            if (isIn(operator_type_arg, MappingOperator::Gradient, MappingOperator::SymmetricGradient)) {
                return Field(numerics::pow(2, field_description_arg.ord), field_description_arg.dim);
            }
            else if (isIn(operator_type_arg, MappingOperator::LargeStrain, MappingOperator::SmallStrain)) {
                return Field(numerics::pow(2, field_description_arg.ord), 3);
            }
            else {
                return field_description_arg;
            }
        }

        constexpr
        Bool
        operator==(
                Field const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                Field const &
                other
        )
        const = default;

        constexpr
        auto
        rows()
        const
        {
            if (isIn(ord, 0)) {
                return Indx(1);
            }
            else if (isIn(ord, 1, 2)) {
                return dim;
            }
            else if (isIn(ord, 3, 4)) {
                return dim * dim;
            }
            else {
                assert(false);
            }
        }

        constexpr
        auto
        cols()
        const
        {
            if (isIn(ord, 0, 1)) {
                return Indx(1);
            }
            else if (isIn(ord, 2, 3)) {
                return dim;
            }
            else if (isIn(ord, 4)) {
                return dim * dim;
            }
            else {
                assert(false);
            }
        }

        constexpr
        auto
        size()
        const
        {
            return rows() * cols();
        }

        Indx ord;

        Indx dim;

    };

    template<FiniteElementMethod Fem>
    struct Discretization;

    template<>
    struct Discretization<FiniteElementMethod::HybridHighOrder>
    {

        constexpr
        Discretization(
                Indx
                ord_cell_arg,
                Indx
                ord_face_arg
        )
        :
        ord_cell(ord_cell_arg),
        ord_face(ord_face_arg)
        {}

        constexpr
        Bool
        operator==(
                Discretization const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                Discretization const &
                other
        )
        const = default;

        constexpr
        auto
        ordMapping(
                MappingOperator
                mapping_arg
        )
        const
        {
            return mapping_arg == MappingOperator::Identity ? ord_cell : ord_face;
        }

        Indx ord_cell;

        Indx ord_face;

    };

    template<>
    struct Discretization<FiniteElementMethod::Lagrange>
    {

        constexpr
        Discretization(
                Indx
                ord_cell_arg
        )
        :
        ord_cell(ord_cell_arg)
        {}

        constexpr
        Bool
        operator==(
                Discretization const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                Discretization const &
                other
        )
        const = default;

        constexpr
        auto
        ordMapping(
                MappingOperator
                mapping_arg
        )
        const
        {
            return mapping_arg == MappingOperator::Identity ? ord_cell : ord_cell - 1;
        }

        Indx ord_cell;

    };

    using HHO = Discretization<FiniteElementMethod::HybridHighOrder>;

    template<Domain D>
    struct LoadComponent
    {

        auto const static constexpr zero = [] (auto const &, auto const &) constexpr { return Real(0); };

        LoadComponent()
        :
        function(zero),
        load_type(LoadType::Natural)
        {}

        LoadComponent(
                LoadType
                load_type_arg
        )
        :
        function(zero),
        load_type(load_type_arg)
        {}

        LoadComponent(
                auto &&
                function_arg,
                LoadType
                load_type_arg
        )
        :
        function(function_arg),
        load_type(load_type_arg)
        {}

        auto
        getImposedValue(
                auto const &
                point_arg,
                auto const &
                time_arg
        )
        const
        {
            return function(point_arg, time_arg);
        }

        std::function<Real(Vector<Real, D.dim> const &, Real const &)> function;

        LoadType load_type;

    };

    template<Domain D>
    struct LoadL
    {

        LoadL(
                auto &&
                felabel_arg,
                auto &&
                label_arg,
                auto
                i_arg,
                auto
                j_arg,
                auto &&
                fun
        )
        :
        fe_label(felabel_arg),
        label(label_arg),
        components(i_arg, j_arg),
        load(SharedPointer<LoadComponent<D>>(fun))
        {}

        Strg label;

        Char fe_label;

        Pair<Indx> components;

        SharedPointer<LoadComponent<D>> load;

    };

    template<Domain D, auto F>
    static constexpr
    auto
    field()
    {
        return Field(F.ord_field, D.dim);
    }

//    template<Domain D, auto F>
//    struct Load : public Array<LoadComponent<D>, field<D, F>().rows(), field<D, F>().cols()>
//    {
//
//        using Base = Array<LoadComponent<D>, field<D, F>().rows(), field<D, F>().cols()>;
//
//        auto const static constexpr finite_element = F;
//
//        Load()
////        :
////        Base(setDefault())
//        {}
//
//        Load(
//                auto &&
//                label_arg,
//                auto &&...
//                load_component_args
//        )
//        requires(sizeof...(load_component_args) == Base::size())
//        :
//        Base{load_component_args...},
//        label(label_arg)
//        {}
//
////        static
////        auto
////        setDefault()
////        {
////            auto base = Array<SharedPointer<LoadComponent<D>>, field<D, F>().rows(), field<D, F>().cols()>();
////            for (int i = 0; i < field<D, F>().rows(); ++i) {
////                for (int j = 0; j < field<D, F>().cols(); ++j) {
////                    base.get(i, j) = SharedPointer<LoadComponent<D>>(LoadComponent<D>());
////                }
////            }
////            return base;
////        }
//
//        Strg label;
//
//    };

//    template<Domain D, auto... F>
//    struct MixedElementLoads : public Collection<Load<D, F>...>
//    {
//
//    };
//
//    template<Domain D, auto... M>
//    struct CoupledElementLoads : public Collection<Load<D, F>...>
//    {
//
//    };

    template<FiniteElementMethod Fem, typename... OperatorType>
    struct FiniteElement
    {

        auto const static constexpr method = Fem;

        static_assert(IsSameAs<MappingOperator, OperatorType...> && sizeof...(OperatorType) > 0);

        constexpr
        FiniteElement(
                Char
                tag_arg,
                Indx
                ord_field_arg,
                Discretization<Fem>
                discretization_arg,
                OperatorType...
                operator_args
        )
        :
        tag(tag_arg),
        ord_field(ord_field_arg),
        discretization(discretization_arg),
        mappings({operator_args...})
        {}

        constexpr
        Bool
        operator==(
                FiniteElement const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                FiniteElement const &
                other
        )
        const = default;

        constexpr
        auto
        ordIntegration(
                Domain const &
                domain_arg
        )
        const
        {
            auto ord_integration_arg = Indx(0);
            auto set_ord_integration = [& ord_integration_arg, this](auto const & m) constexpr mutable {
                ord_integration_arg = numerics::max(ord_integration_arg, discretization.ordMapping(m));
            };
            std::for_each(mappings.data.begin(), mappings.data.end(), set_ord_integration);
            ord_integration_arg *= 2;
            if (domain_arg.frame == EuclideanFrame::AxiSymmetric) {
                ord_integration_arg += 1;
            }
            return ord_integration_arg;
        }

        Char tag;

        Indx ord_field;

        Discretization<Fem> discretization;

        Array<MappingOperator, sizeof...(OperatorType)> mappings;

    };

    template<auto... F>
    struct MixedElement : public Aggregate<std::remove_cvref_t<decltype(F)>...>
    {

        template<template<auto> typename T>
        using Elements = Collection<T<F>...>;

        template<template<auto...> typename T, auto E, Domain D, auto M, auto C>
        using Elements2 = Collection<T<E, D, F, M, C>...>;

    private:

        using FiniteElementAggregate = Aggregate<std::remove_cvref_t<decltype(F)>...>;

    public:

        constexpr
        MixedElement(
                Model
                model_arg
        )
        :
        FiniteElementAggregate({F...}),
        model(model_arg)
        {}

        constexpr
        Bool
        operator==(
                MixedElement const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                MixedElement const &
                other
        )
        const = default;

        constexpr
        auto
        ordIntegration()
        const
        {
            auto res = Indx(0);
            auto set_ord_integration = [& res](auto const & x) constexpr mutable {
                for (int i = 0; i < x.mappings.size(); ++i) {
                    res = numerics::max(res, x.discretization.ordMapping(x.mappings.get(i)));
                }
            };
            aggregate::apply(set_ord_integration, * this);
            return res;
        }

        constexpr
        auto
        ordIntegration(
                Domain const &
                domain_arg
        )
        const
        {
            return numerics::max(F.ordIntegration(domain_arg)...);
        }

        Model model;

        // une loi de comportement

    };

    template<auto... M>
    struct CoupledElement : public Aggregate<std::remove_cvref_t<decltype(M)>...>
    {

        template<template<auto> typename T>
        using Elements = Collection<T<M>...>;

        template<template<auto...> typename T, auto E, Domain D, auto C>
        using Elements2 = Collection<T<E, D, M, C>...>;

    private:

        using MixedElementAggregate = Aggregate<std::remove_cvref_t<decltype(M)>...>;

    public:

        constexpr
        CoupledElement(
                QuadratureRule
                quadrature_rule_arg
        )
        :
        MixedElementAggregate({M...}),
        quadrature_rule(quadrature_rule_arg)
        {}

        constexpr
        Bool
        operator==(
                CoupledElement const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                CoupledElement const &
                other
        )
        const = default;

        constexpr
        auto
        ordIntegration()
        const
        {
            return numerics::max(M.ordIntegration()...);
        }

        constexpr
        auto
        ordIntegration(
                Domain const &
                domain_arg
        )
        const
        {
            return numerics::max(M.ordIntegration(domain_arg)...);
        }

        QuadratureRule quadrature_rule;

        // un point de gauss, plusiuers lois de comportment

    };

    struct Behaviour
    {

        using BehaviourPointer = SharedPointer<mgis::behaviour::Behaviour>;

        Behaviour()
        :
        ptr_behaviour(),
        strain_type(StrainType::SmallStrain)
        {}

        Behaviour(
                Strg const &
                path,
                Strg const &
                name,
                HypothesisType const &
                hypothesis_type,
                StrainType const &
                strain_type_arg = StrainType::SmallStrain
        )
        :
        strain_type(strain_type_arg),
        ptr_behaviour(
        getMgisBehaviour(
                path,
                name,
                hypothesis_type,
                strain_type_arg
        )
        )
        {}

    private:

        inline static
        mgis::behaviour::Hypothesis
        getMgisHypothesis(
                HypothesisType
                hypothesis_type_arg
        )
        {
            mgis::behaviour::Hypothesis hypothesis;
            switch (hypothesis_type_arg) {
                case HypothesisType::PlaneStrain:
                    hypothesis = mgis::behaviour::Hypothesis::PLANESTRAIN;
                    break;
                case HypothesisType::PlaneStress:
                    hypothesis = mgis::behaviour::Hypothesis::PLANESTRESS;
                    break;
                case HypothesisType::AxySymmetric:
                    hypothesis = mgis::behaviour::Hypothesis::AXISYMMETRICAL;
                    break;
                case HypothesisType::Volumetric:
                    hypothesis = mgis::behaviour::Hypothesis::TRIDIMENSIONAL;
                    break;
                default:
                    assert(false);
            }
            return hypothesis;
        }

        static inline
        BehaviourPointer
        getMgisBehaviour(
                Strg const &
                path_arg,
                Strg const &
                name_arg,
                HypothesisType
                hypothesis_type_arg,
                StrainType
                strain_type_arg
        )
        {
            BehaviourPointer a;
            mgis::behaviour::Hypothesis hyp = getMgisHypothesis(hypothesis_type_arg);
            if (strain_type_arg == StrainType::SmallStrain) {
                a = BehaviourPointer(mgis::behaviour::load(path_arg, name_arg, hyp));
                return a;
            }
            else {
                mgis::behaviour::FiniteStrainBehaviourOptions opt;
                opt.stress_measure = mgis::behaviour::FiniteStrainBehaviourOptions::PK1;
                opt.tangent_operator = mgis::behaviour::FiniteStrainBehaviourOptions::TangentOperator::DPK1_DF;
                a = BehaviourPointer(mgis::behaviour::load(opt, path_arg, name_arg, hyp));
                return a;
            }
        }

    public:

        BehaviourPointer ptr_behaviour;

        StrainType strain_type;

    };

    template<auto M>
    struct MeshInteriorDomain
    {

        MeshInteriorDomain() = default;

    };

//    template<auto A>
//    struct MeshInteriorDomain
//    {
//
//
//
//    };

    template<typename T>
    struct DDomain
    {

        //

    };

    template<auto... A>
    struct Mesh
    {



    };

}

#endif //LOLITA_LOLITA_CORE_HXX
