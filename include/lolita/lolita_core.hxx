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
#include <ostream>

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
                EuclideanFrame
                euclidean_frame_arg,
                Indx
                dim_arg
        )
        :
        frame(euclidean_frame_arg),
        dim(dim_arg)
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

        EuclideanFrame frame;

        Indx dim;

    };

    static constexpr
    auto
    ordIntegration(
            EuclideanFrame
            euclidean_frame,
            Indx
            ord_polynomial
    )
    {
        return euclidean_frame == EuclideanFrame::AxiSymmetric ? 2 * (ord_polynomial) + 1 : 2 * (ord_polynomial);
    }

    struct Field
    {

        constexpr
        Field()
        :
        ord(),
        dim()
        {}

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

    struct Label2
    {

        constexpr
        Label2(
                auto inpt,
                Indx i1,
                Indx i2
        )
        :
        sv(inpt),
        i(i1),
        j(i2)
        {}

        std::string_view sv;

        Indx i;
        Indx j;

    };

    template<typename T>
    struct
    StrView2
    {

        auto const static constexpr max_size_ = Indx(30);

    private:

        static constexpr
        auto
        setData(
                std::basic_string_view<T> &&
                str
        )
        {
            auto data = std::array<T, max_size_>();
            for (int i = 0; i < max_size_; i++) {
                if (i < str.size()) {
                    data[i] = str[i];
                }
                else {
                    data[i] = '#';
                }
            }
            return data;
        }

    public:

        constexpr
        StrView2(
                std::basic_string_view<T> &&
                str
        )
        noexcept
        :
        size_(str.size()),
        data_(setData(std::forward<std::basic_string_view<T>>(str)))
        {}

        constexpr
        Bool
        operator==(
                StrView2 const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                StrView2 const &
                other
        )
        const = default;

        auto
        view()
        const
        {
            return std::basic_string_view<T>(data_.data(), size_);
        }

        auto
        view()
        {
            return std::basic_string_view<T>(data_.data(), size_);
        }

        friend
        std::ostream &
        operator<<(
                std::ostream &
                os,
                StrView2 const &
                view_2
        )
        {
            os << view_2.view();
            return os;
        }

        Indx size_;

        std::array<T, max_size_> data_;

    };

    using StrView2Char = StrView2<Char>;

    template<auto A>
    struct StrTest
    {

        auto const static constexpr a = A;

        StrTest(std::basic_string_view<Char> && va, Indx i)
        :
        v(std::forward<std::basic_string_view<Char>>(va))
        {}

        StrView2Char v;

    };

    struct
    StrView
    {

        constexpr
        StrView(const Char * __str) noexcept
                : _M_len{std::char_traits<Char>::length(__str)},
                  _M_str{__str}
        { }

        size_t	     _M_len;
        const Char * _M_str;

        constexpr
        Bool
        operator==(
                StrView const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                StrView const &
                other
        )
        const = default;

    };

    struct Label
    {

        constexpr
        Label()
        :
        label_()
        {}

        constexpr
        Label(
                auto &&
                label_arg
        )
        :
        label_(label_arg)
        {}

        constexpr
        Bool
        operator==(
                Label const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                Label const &
                other
        )
        const = default;

        auto
        toString()
        const
        {
            return Strg(label_);
        }

        Char const * label_;

    };

    template<Domain D>
    struct LoadComponent
    {

        auto const static constexpr zero = [] (auto const &, auto const &) constexpr { return Real(0); };

        constexpr
        LoadComponent()
        :
        function(zero),
        load_type(LoadType::Natural)
        {}

        constexpr explicit
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
    struct Load
    {

        Load(
                auto &&
                finite_element_label_arg,
                auto &&
                domain_label_arg,
                auto
                i_arg,
                auto
                j_arg,
                auto &&
                fun
        )
        :
        finite_element_label(finite_element_label_arg),
        domain_label(domain_label_arg),
        components(i_arg, j_arg),
        load(SharedPointer<LoadComponent<D>>(fun))
        {}

        Char finite_element_label;

        Strg domain_label;

        Pair<Indx> components;

        SharedPointer<LoadComponent<D>> load;

    };

    enum struct
    FLDTyp
    {

        Primal,
        Internal,
        External,
        Material

    };

    template<typename... OperatorType>
    requires(std::same_as<MappingOperator, OperatorType> && ...)
    struct StructuralField
    {

        auto const static constexpr fld_type = FLDTyp::Primal;

        constexpr
        StructuralField(
                Char &&
                tag_arg,
                Indx &&
                ord_field_arg,
                OperatorType &&...
                mpo
        )
        :
        tag(tag_arg),
        ord_field(ord_field_arg),
        mappings({std::forward<MappingOperator>(mpo)...})
        {}

        constexpr
        Bool
        operator==(
                StructuralField const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                StructuralField const &
                other
        )
        const = default;

        Char
        tag;

        Indx
        ord_field;

        Array<MappingOperator, sizeof...(OperatorType)>
        mappings;

    };

    struct AuxiliaryField
    {

        constexpr
        AuxiliaryField(
                Char &&
                tag_arg,
                Indx &&
                ord_field_arg,
                FLDTyp &&
                fld_type_arg
        )
        :
        tag(tag_arg),
        ord_field(ord_field_arg),
        fld_type(fld_type_arg)
        {}

        constexpr
        Bool
        operator==(
                AuxiliaryField const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                AuxiliaryField const &
                other
        )
        const = default;

        Char
        tag;

        Indx
        ord_field;

        FLDTyp
        fld_type;

    };

    template<typename T>
    concept FieldType = requires(
            std::remove_reference_t<T> const
            t
    )
    {
        { t.tag } -> std::convertible_to<Indx>;
        { t.ord_field } -> std::convertible_to<Indx>;
        { t.fld_type } -> std::convertible_to<FLDTyp>;
//        // C++ 23 features for gcc-12
//        { auto(t.tag) } -> std::same_as<Indx>;
//        { auto(t.ord_field) } -> std::same_as<Indx>;
//        { auto(t.fld_type) } -> std::same_as<FLDTyp>;
    };

    namespace detail
    {

        template<typename T>
        struct FLDTraits : public std::false_type
        {};

        template<>
        struct FLDTraits<AuxiliaryField> : public std::true_type
        {};

        template<typename... OperatorType>
        struct FLDTraits<StructuralField<OperatorType...>> : public std::true_type
        {};

    }

    template<typename T>
    concept FLDType = detail::FLDTraits<T>::value;

    namespace detail
    {

        struct MappingData
        {

            constexpr
            MappingData()
                    :
                    coefficient_(1)
            {}

            constexpr
            MappingData(
                    Real &&
                    coefficient_arg
            )
                    :
                    coefficient_(coefficient_arg)
            {}

            Real coefficient_;

        };

        template<Domain, FLDType auto, MappingOperator>
        struct MappingPolicy;

        template<Domain Dmn, FLDType auto Fld, MappingOperator Mpg>
        requires(Mpg == MappingOperator::Identity)
        struct MappingPolicy<Dmn, Fld, Mpg>
        {

            auto const static constexpr rows_ = Indx(1);

            auto const static constexpr cols_ = Fld.ord_;

            auto const static constexpr data_ = Array<MappingData, rows_, cols_>();

        };

        template<Domain Dmn, FLDType auto Fld, MappingOperator Mpg>
        requires(Mpg == MappingOperator::Gradient)
        struct MappingPolicy<Dmn, Fld, Mpg>
        {

            auto const static constexpr rows_ = Dmn.dim;

            auto const static constexpr cols_ = Dmn.dim;

            auto const static constexpr data_ = Array<MappingData, rows_, cols_>();

        };

    }

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

        constexpr
        auto
        makeHigher()
        const
        {
            return Discretization(ord_cell + 1, ord_face + 1);
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

        constexpr
        auto
        makeHigher()
        const
        {
            return Discretization(ord_cell + 1);
        }

        Indx ord_cell;

    };

    namespace detail
    {

        template<typename T>
        struct DiscretizationTrait
        {

            auto const static constexpr value = false;

        };

        template<FiniteElementMethod Fem>
        struct DiscretizationTrait<Discretization<Fem>>
        {

            auto const static constexpr value = true;

        };

    }

    template<typename T>
    concept DiscretizationType = detail::DiscretizationTrait<T>::value;

    using HHO = Discretization<FiniteElementMethod::HybridHighOrder>;

    template<typename... FLDt>
    requires(FieldType<FLDt> && ...)
    struct BHV2 : public Collection<FLDt...>
    {

        constexpr
        BHV2(
                FLDt...
                fl_dt
        )
        :
        Collection<FLDt...>(fl_dt...)
//        num_internal_variables(num_internal_variables_arg),
//        num_external_variables(num_external_variables_arg),
//        num_material_variables(num_material_variables_arg)
        {}

//        Indx
//        num_internal_variables;
//        Indx
//        num_external_variables;
//        Indx
//        num_material_variables;

    };

    template<FiniteElementMethod Fem, typename FLDt>
    requires(FLDType<FLDt>)
    struct FiniteElement2
    {

        auto const static constexpr method = Fem;

        constexpr
        FiniteElement2(
                FLDt
                fl_dt,
                Discretization<Fem>
                discretization_arg,
                QuadratureRule
                quadrature_rule_arg,
                Indx
                ord_quadrature_arg
        )
        :
        field_(fl_dt),
        discretization_(discretization_arg),
        quadrature_rule_(quadrature_rule_arg),
        ord_quadrature_(ord_quadrature_arg)
        {}

        constexpr
        FiniteElement2(
                FLDt
                fl_dt,
                Discretization<Fem>
                discretization_arg,
                QuadratureRule
                quadrature_rule_arg
        )
        :
        field_(fl_dt),
        discretization_(discretization_arg),
        quadrature_rule_(quadrature_rule_arg),
        ord_quadrature_(ordOperator())
        {}

        constexpr
        Bool
        operator==(
                FiniteElement2 const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                FiniteElement2 const &
                other
        )
        const = default;

        constexpr
        auto
        ordOperator()
        const
        {
            auto ord_operator_arg = Indx(0);
            for (auto const & mapping : field_.mappings.data) {
                ord_operator_arg = numerics::max(ord_operator_arg, discretization_.ordMapping(mapping));
            }
            return ord_operator_arg;
        }

        constexpr
        auto
        field()
        const
        {
            return field_;
        }

        constexpr
        auto
        discretization()
        const
        {
            return discretization_;
        }

        constexpr
        auto const &
        quadratureRule()
        const
        {
            return quadrature_rule_;
        }

        FLDt field_;

        Discretization<Fem> discretization_;

        QuadratureRule quadrature_rule_;

        Indx ord_quadrature_;

    };

    namespace detail
    {

        template<typename T>
        struct IsFEM
        {

            auto const static constexpr value = false;

        };

        template<FiniteElementMethod Fem, typename FLDt>
        struct IsFEM<FiniteElement2<Fem, FLDt>>
        {

            auto const static constexpr value = true;

        };

    }

    template<typename T>
    concept FEMType = detail::IsFEM<T>::value;

    template<Domain D, auto F>
    static constexpr
    auto
    field()
    {
        return Field(F.field().ord_field, D.dim);
    }

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
                QuadratureRule
                quadrature_rule_arg,
                Indx
                ord_quadrature_arg,
                OperatorType...
                operator_args
        )
        :
        mappings({operator_args...}),
        tag(tag_arg),
        ord_field(ord_field_arg),
        discretization(discretization_arg),
        quadrature_rule(quadrature_rule_arg),
        ord_quadrature(ord_quadrature_arg)
        {}

        constexpr
        FiniteElement(
                Char
                tag_arg,
                Indx
                ord_field_arg,
                Discretization<Fem>
                discretization_arg,
                QuadratureRule
                quadrature_rule_arg,
                OperatorType...
                operator_args
        )
        :
        mappings({operator_args...}),
        tag(tag_arg),
        ord_field(ord_field_arg),
        discretization(discretization_arg),
        quadrature_rule(quadrature_rule_arg),
        ord_quadrature(ordOperator())
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
        ordOperator()
        const
        {
            auto ord_integration_arg = Indx(0);
            for (auto const & mapping : mappings.data) {
                ord_integration_arg = numerics::max(ord_integration_arg, discretization.ordMapping(mapping));
            }
            return ord_integration_arg;
        }

        constexpr
        auto
        ordIntegration(
                Domain const &
                domain_arg
        )
        const
        {
            return domain_arg.template ordIntegration(ordOperator());
        }

        Array<MappingOperator, sizeof...(OperatorType)> mappings;

        Char tag;

        Indx ord_field;

        Discretization<Fem> discretization;

        QuadratureRule quadrature_rule;

        Indx ord_quadrature;

    };

    template<auto... F>
    struct Elements : public Aggregate<std::remove_cvref_t<decltype(F)>...>
    {

        template<template<auto, Domain, auto> typename T, auto E, Domain D>
        using MixedElementType = Collection<SharedPointer<T<E, D, F>>...>;
//        using Type = Collection<UniquePointer<T<E, D, F>>...>;
//        using Type = Collection<T<E, D, F>...>;

        constexpr
        Elements()
        :
        Aggregate<std::remove_cvref_t<decltype(F)>...>({F...})
        {}

        constexpr
        Bool
        operator==(
                Elements const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                Elements const &
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
            return numerics::max(F.ordIntegration(domain_arg)...);
        }

    };

    template<auto... F>
    struct MixedElement : public Aggregate<std::remove_cvref_t<decltype(F)>...>
    {

        template<template<auto> typename T>
        using Elements = Collection<T<F>...>;

        template<template<auto...> typename T, auto E, Domain D>
        using Elements2 = Collection<T<E, D, F>...>;

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
        ptr_behaviour()
        {}

        Behaviour(
                Char &&
                finite_element_tag_arg,
                Indx &&
                finite_element_index_arg,
                SharedPointer<mgis::behaviour::Behaviour> const &
                bhv_arg,
                auto &&...
                domains_args
        )
        :
        finite_element_tag(finite_element_tag_arg),
        finite_element_index(finite_element_index_arg),
        ptr_behaviour(bhv_arg),
        domains({domains_args...})
        {}

        Char finite_element_tag;

        Indx finite_element_index;

        BehaviourPointer ptr_behaviour;

        Array<Strg> domains;

    };

}

#endif //LOLITA_LOLITA_CORE_HXX
