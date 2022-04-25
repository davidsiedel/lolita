//
// Created by dsiedel on 14/04/2022.
//

#ifndef LOLITA_LOLITA_CORE_HXX
#define LOLITA_LOLITA_CORE_HXX

#include "lolita.hxx"
#include "lolita_lolita.hxx"
#include "lolita_matrix.hxx"

namespace lolita::core
{

    enum struct MeshFormatType
    {

        Gmsh,

    };

    enum struct LoadType
    {

        Volumetric,
        Dirichlet,
        Neumann

    };

    enum struct EuclideanFrame
    {

        Cartesian,
        AxiSymmetric

    };

//    enum struct Body
//    {
//
//        Cell,
//        Face,
//        Edge,
//        Node
//
//    };

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

    enum struct Basis
    {

        Monomial,
        Lagrange

    };

    struct BasisDescription
    {

        constexpr
        BasisDescription()
                :
                basis(),
                ord()
        {}

        constexpr
        BasisDescription(
                Basis
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
                BasisDescription const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                BasisDescription const &
                other
        )
        const = default;

        Basis basis;

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

    enum struct MappingOperator
    {

        Gradient,
        SymmetricGradient,
        Divergence,
        SmallStrain,
        LargeStrain,
        Identity,

    };

//    struct Mapping
//    {
//
//        constexpr
//        Bool
//        operator==(
//                Mapping const &
//                other
//        )
//        const = default;
//
//        constexpr
//        Bool
//        operator!=(
//                Mapping const &
//                other
//        )
//        const = default;
//
//        MappingOperator mapping;
//
//        Indx ord;
//
//    };

    enum struct FiniteElementMethod
    {

        Lagrange,
        HybridHighOrder

    };

    auto const static constexpr fem_hho = FiniteElementMethod::HybridHighOrder;
    auto const static constexpr fem_lag = FiniteElementMethod::Lagrange;

    template<auto M>
    struct MeshInteriorDomain
    {

        MeshInteriorDomain() = default;

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

        Indx dim;

        EuclideanFrame frame;

    };

    struct FField
    {

        constexpr
        FField() = default;

        constexpr
        FField(
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
        FField
        fromMapping(
                FField const &
                field_description_arg,
                MappingOperator
                operator_type_arg
        )
        {
            if (isIn(operator_type_arg, MappingOperator::Gradient, MappingOperator::SymmetricGradient)) {
                return FField(numerics::pow(2, field_description_arg.ord), field_description_arg.dim);
            }
            else if (isIn(operator_type_arg, MappingOperator::LargeStrain, MappingOperator::SmallStrain)) {
                return FField(numerics::pow(2, field_description_arg.ord), 3);
            }
            else {
                return field_description_arg;
            }
        }

        constexpr
        Bool
        operator==(
                FField const &
                other
        )
        const = default;

        constexpr
        Bool
        operator!=(
                FField const &
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

    template<FiniteElementMethod Fem, typename... OperatorType>
    struct FiniteElement
    {

        auto const static constexpr method = Fem;

        static_assert(areSame<MappingOperator, OperatorType...>() && sizeof...(OperatorType) > 0);

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
            auto find_index = [& res](auto const & x) constexpr mutable {
                for (int i = 0; i < x.mappings.size(); ++i) {
                    res = numerics::max(res, x.discretization.ordMapping(x.mappings.get(i)));
                }
            };
            aggregate::apply([&](auto const &... x){(..., find_index(x));}, * this);
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

//        constexpr
//        auto
//        ordIntegration()
//        const
//        {
//            return numerics::max(M.ordIntegration()...);
//        }

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

    namespace detail
    {

        template<FiniteElementMethod, auto>
        struct FiniteElementDescriptionPolicy;

        template<auto discrete_field_arg>
        struct FiniteElementDescriptionPolicy<FiniteElementMethod::HybridHighOrder, discrete_field_arg>
        {

            using DiscreteField = std::remove_cvref_t<decltype(discrete_field_arg)>;

        private:

            auto const static constexpr num_mappings = discrete_field_arg.mappings.size();

        public:

            constexpr
            FiniteElementDescriptionPolicy(
                    Model
                    model_arg,
                    Indx
                    cell_order_arg,
                    Indx
                    face_order_arg,
                    auto...
                    mapping_order_args
            )
            :
                    cell_order(cell_order_arg),
                    face_order(face_order_arg),
                    discrete_field(setDiscreteField(Array<Indx, num_mappings>{static_cast<Indx>(mapping_order_args)...}))
            {
                static_assert(sizeof...(mapping_order_args) == num_mappings);
            }

            constexpr
            FiniteElementDescriptionPolicy(
                    Model &&
                    model_arg,
                    Indx &&
                    cell_order_arg,
                    Indx &&
                    face_order_arg
            )
            :
                    cell_order(cell_order_arg),
                    face_order(face_order_arg),
                    discrete_field(setDiscreteField(cell_order_arg, face_order_arg))
            {}

//            constexpr bool operator==(FiniteElementPolicy const &) const = default;
//
//            constexpr bool operator!=(FiniteElementPolicy const &) const = default;

            constexpr
            Bool
            operator==(
                    FiniteElementDescriptionPolicy const &
                    other
            )
            const
            {
                auto eq_0 = cell_order == other.cell_order;
                auto eq_1 = face_order == other.face_order;
                auto eq_2 = discrete_field == other.discrete_field;
                return eq_0 && eq_1 && eq_2;
            }

            constexpr
            Bool
            operator!=(
                    FiniteElementDescriptionPolicy const &
                    other
            )
            const
            {
                return !(other == * this);
            }

            constexpr
            auto
            getIntegrationOrder()
            const
            {
                auto val = Indx(0);
                for (int i = 0; i < num_mappings; ++i) {
                    auto const & mapping_description = discrete_field.mappings.get(i);
                    if (mapping_description.ord_mapping > val) {
                        val = mapping_description.ord_mapping;
                    }
                }
                return val;
            }

        private:

            static constexpr
            auto
            setDiscreteField(
                    auto const &
                    mapping_order_array_arg
            )
            {
                auto other = discrete_field_arg;
                for (int i = 0; i < num_mappings; ++i) {
                    other.mappings.get(i).ord_mapping = mapping_order_array_arg.get(i);
                }
                return other;
            }

            static constexpr
            auto
            setDiscreteField(
                    Indx
                    cell_order_arg,
                    Indx
                    face_order_arg
            )
            {
                auto other = discrete_field_arg;
                for (int i = 0; i < num_mappings; ++i) {
                    auto & mapping_description = other.mappings.get(i);
                    if (mapping_description.mapping == MappingOperator::Identity) {
                        mapping_description.ord_mapping = cell_order_arg;
                    } else {
                        mapping_description.ord_mapping = face_order_arg;
                    }
                }
                return other;
            }

        public:

            Indx cell_order;

            Indx face_order;

            DiscreteField discrete_field;

        };

    }

    template<auto discrete_field_arg>
    using FiniteElementDescription = detail::FiniteElementDescriptionPolicy<discrete_field_arg.finite_element_method, discrete_field_arg>;

//    auto const static constexpr fld_scalar = Field::Scalar;
//    auto const static constexpr fld_vector = Field::Vector;

//    struct TensorDescription
//    {
//
//        constexpr
//        TensorDescription()
//        :
//        rows(-1),
//        cols(-1),
//        size(-1),
//        storage_option(matrix::StorageOption::RowMajor)
//        {}
//
//        constexpr explicit
//        TensorDescription(
//                Intg
//                rows_arg,
//                matrix::StorageOption
//                storage_option_arg = matrix::StorageOption::RowMajor
//        )
//        :
//        rows(rows_arg),
//        cols(1),
//        size(rows_arg),
//        storage_option(storage_option_arg)
//        {
//            assert(rows_arg >= -1);
//        }
//
//        constexpr
//        TensorDescription(
//                Intg
//                rows_arg,
//                Intg
//                cols_arg,
//                matrix::StorageOption
//                storage_option_arg = matrix::StorageOption::RowMajor
//        )
//        :
//        rows(rows_arg),
//        cols(cols_arg),
//        size(rows_arg * cols_arg),
//        storage_option(storage_option_arg)
//        {
//            assert(rows_arg >= -1);
//            assert(cols_arg >= -1);
//        }
//
//        constexpr
//        Bool
//        operator==(
//                TensorDescription const &
//                other
//        )
//        const = default;
//
//        constexpr
//        Bool
//        operator!=(
//                TensorDescription const &
//                other
//        )
//        const = default;
//
//        Intg rows;
//
//        Intg cols;
//
//        Intg size;
//
//        matrix::StorageOption storage_option;
//
//    };

//    struct FieldDescription
//    {
//
//        constexpr
//        FieldDescription()
//        :
//        field(Field::Scalar),
//        dim_field(0),
//        tensor_description()
//        {}
//
//        constexpr inline
//        FieldDescription(
//                Field
//                field_arg,
//                Indx
//                dim_field_arg
//        )
//        :
//        field(field_arg),
//        dim_field(dim_field_arg),
//        tensor_description(setTensorDescription(field_arg, dim_field_arg))
//        {}
//
//        constexpr
//        Bool
//        operator==(
//                FieldDescription const &
//                other
//        )
//        const
//        {
//            auto eq_0 = field == other.field;
//            auto eq_1 = dim_field == other.dim_field;
//            auto eq_2 = tensor_description == other.tensor_description;
//            return eq_0 && eq_1 && eq_2;
//        }
//
//        constexpr
//        Bool
//        operator!=(
//                FieldDescription const &
//                other
//        )
//        const
//        {
//            return !(other == * this);
//        }
//
//    private:
//
//        static constexpr inline
//        TensorDescription
//        setTensorDescription(
//                Field
//                field_arg,
//                Indx
//                dim_field_arg
//        )
//        {
//            if (field_arg == Field::Scalar) {
//                return TensorDescription(1);
//            }
//            else if (field_arg == Field::Vector) {
//                return TensorDescription(dim_field_arg);
//            }
//            else if (field_arg == Field::Tensor) {
//                return TensorDescription(dim_field_arg, dim_field_arg);
//            }
//            else {
//                return TensorDescription();
//            }
//        }
//
//    public:
//
//        static constexpr inline
//        FieldDescription
//        fromMapping(
//                FieldDescription const &
//                field_description_arg,
//                Mapping
//                operator_type_arg
//        )
//        {
//            if (operator_type_arg == Mapping::Gradient || operator_type_arg == Mapping::SymmetricGradient) {
//                if (field_description_arg.field == Field::Scalar) {
//                    return FieldDescription(Field::Vector, field_description_arg.dim_field);
//                }
//                else if (field_description_arg.field == Field::Vector) {
//                    return FieldDescription(Field::Tensor, field_description_arg.dim_field);
//                }
//                else {
//                    assert(false);
//                }
//            }
//            else if (operator_type_arg == Mapping::LargeStrain || operator_type_arg == Mapping::SmallStrain) {
//                if (field_description_arg.field == Field::Scalar) {
//                    return FieldDescription(Field::Vector, 3);
//                }
//                else if (field_description_arg.field == Field::Vector) {
//                    return FieldDescription(Field::Tensor, 3);
//                }
//                else {
//                    assert(false);
//                }
//            }
//            else {
//                return field_description_arg;
//            }
//        }
//
//        Field field;
//
//        Indx dim_field;
//
//        TensorDescription tensor_description;
//
//    };
//
//
//    namespace detail
//    {
//
//        template<FiniteElementMethod, auto>
//        struct FiniteElementDescriptionPolicy;
//
//        template<auto discrete_field_arg>
//        struct FiniteElementDescriptionPolicy<FiniteElementMethod::HybridHighOrder, discrete_field_arg>
//        {
//
//            using DiscreteField = std::remove_cvref_t<decltype(discrete_field_arg)>;
//
//        private:
//
//            auto const static constexpr num_mappings = discrete_field_arg.mappings.size();
//
//        public:
//
//            constexpr
//            FiniteElementDescriptionPolicy(
//                    elt::Model &&
//                    model_arg,
//                    Indx &&
//                    cell_order_arg,
//                    Indx &&
//                    face_order_arg,
//                    auto &&...
//                    mapping_order_args
//            )
//            :
//                    cell_order(cell_order_arg),
//                    face_order(face_order_arg),
//                    discrete_field(setDiscreteField(Array<Indx, num_mappings>{static_cast<Indx>(mapping_order_args)...}))
//            {
//                static_assert(sizeof...(mapping_order_args) == num_mappings);
//            }
//
//            constexpr
//            FiniteElementDescriptionPolicy(
//                    elt::Model &&
//                    model_arg,
//                    Indx &&
//                    cell_order_arg,
//                    Indx &&
//                    face_order_arg
//            )
//            :
//                    cell_order(cell_order_arg),
//                    face_order(face_order_arg),
//                    discrete_field(setDiscreteField(cell_order_arg, face_order_arg))
//            {}
//
////            constexpr bool operator==(FiniteElementPolicy const &) const = default;
////
////            constexpr bool operator!=(FiniteElementPolicy const &) const = default;
//
//            constexpr
//            Bool
//            operator==(
//                    FiniteElementDescriptionPolicy const &
//                    other
//            )
//            const
//            {
//                auto eq_0 = cell_order == other.cell_order;
//                auto eq_1 = face_order == other.face_order;
//                auto eq_2 = discrete_field == other.discrete_field;
//                return eq_0 && eq_1 && eq_2;
//            }
//
//            constexpr
//            Bool
//            operator!=(
//                    FiniteElementDescriptionPolicy const &
//                    other
//            )
//            const
//            {
//                return !(other == * this);
//            }
//
//            constexpr
//            auto
//            getIntegrationOrder()
//            const
//            {
//                auto val = Indx(0);
//                for (int i = 0; i < num_mappings; ++i) {
//                    auto const & mapping_description = discrete_field.mappings.get(i);
//                    if (mapping_description.ord_mapping > val) {
//                        val = mapping_description.ord_mapping;
//                    }
//                }
//                return val;
//            }
//
//        private:
//
//            static constexpr
//            auto
//            setDiscreteField(
//                    auto const &
//                    mapping_order_array_arg
//            )
//            {
//                auto other = discrete_field_arg;
//                for (int i = 0; i < num_mappings; ++i) {
//                    other.mappings.get(i).ord_mapping = mapping_order_array_arg.get(i);
//                }
//                return other;
//            }
//
//            static constexpr
//            auto
//            setDiscreteField(
//                    Indx
//                    cell_order_arg,
//                    Indx
//                    face_order_arg
//            )
//            {
//                auto other = discrete_field_arg;
//                for (int i = 0; i < num_mappings; ++i) {
//                    auto & mapping_description = other.mappings.get(i);
//                    if (mapping_description.mapping == Mapping::Identity) {
//                        mapping_description.ord_mapping = cell_order_arg;
//                    } else {
//                        mapping_description.ord_mapping = face_order_arg;
//                    }
//                }
//                return other;
//            }
//
//        public:
//
//            Indx cell_order;
//
//            Indx face_order;
//
//            DiscreteField discrete_field;
//
//        };
//
//    }
//
//    template<auto discrete_field_arg>
//    using FiniteElementDescription = detail::FiniteElementDescriptionPolicy<discrete_field_arg.finite_element_method, discrete_field_arg>;
//
//    template<auto... finite_element_args>
//    struct MixedElementDescription : public Aggregate<std::remove_cvref_t<decltype(finite_element_args)>...>
//    {
//
//        template<template<auto> typename T>
//        using Elements = Collection<T<finite_element_args>...>;
//
//    private:
//
//        using FiniteElementAggregate = Aggregate<std::remove_cvref_t<decltype(finite_element_args)>...>;
//
//    public:
//
////        auto const static constexpr finite_elements = FiniteElementAggregate{finite_element_args...};
//
//        constexpr
//        MixedElementDescription(
//                elt::Model &&
//                model_arg,
//                Indx &&
//                dim_euclidean_arg,
//                EuclideanFrame &&
//                frame_type_arg,
//                Quadrature &&
//                quadrature_type_arg
//        )
//        :
//        FiniteElementAggregate({finite_element_args...}),
//        model(model_arg),
//        dim_euclidean(dim_euclidean_arg),
//        frame_type(frame_type_arg),
//        quadrature_description(quadrature_type_arg, setIntegrationOrder(frame_type_arg))
//        {}
//
//        constexpr
//        Bool
//        operator==(
//                MixedElementDescription const &
//                other
//        )
//        const
//        {
//            auto eq_0 = model == other.model;
//            auto eq_1 = dim_euclidean == other.dim_euclidean;
//            auto eq_2 = frame_type == other.frame_type;
//            auto eq_3 = quadrature_description == other.quadrature_description;
//            return eq_0 && eq_1 && eq_2 && eq_3;
//        }
//
//        constexpr
//        Bool
//        operator!=(
//                MixedElementDescription const &
//                other
//        )
//        const
//        {
//            return !(other == * this);
//        }
//
//        static constexpr
//        auto
//        setIntegrationOrder(
//                EuclideanFrame const &
//                frame_arg
//        )
//        {
//            auto val = max(finite_element_args.getIntegrationOrder()...);
//            if (frame_arg == EuclideanFrame::AxiSymmetric) {
//                val += 1;
//            }
//            return val;
//        }
//
//
//        QuadratureDescription quadrature_description;
//
//        Indx dim_euclidean;
//
//        EuclideanFrame frame_type;
//
//        elt::Model model;
//
//        // une loi de comportement
//
//    };
//
//    template<auto... A>
//    struct CoupledElement
//    {
//
//        // un point de gauss, plusiuers lois de comportment
//
//    };
//
//    template<auto A>
//    struct MeshInteriorDomain
//    {
//
//
//
//    };
//
//    template<typename T>
//    struct DDomain
//    {
//
//        //
//
//    };
//
//    template<auto... A>
//    struct Mesh
//    {
//
//
//
//    };

}

#endif //LOLITA_LOLITA_CORE_HXX
