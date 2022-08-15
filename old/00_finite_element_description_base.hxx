//
// Created by dsiedel on 23/03/2022.
//

#ifndef FETA_00_FINITE_ELEMENT_DESCRIPTION_BASE_HXX
#define FETA_00_FINITE_ELEMENT_DESCRIPTION_BASE_HXX

//#include "new/finite_element_discrete_fields.hxx"
//#include "new/operator_output.hxx"
#include "new/feta_core_model_description.hxx"
#include "new/00_finite_element_type.hxx"
#include "new/feta_core_field_description.hxx"
#include "new/feta_core_model_description.hxx"
#include "new/feta_core_field_description.hxx"

namespace feta::core::element
{

    struct DiscreteFieldComponent
    {

        constexpr
        DiscreteFieldComponent()
        :
        structure(Body::Cell),
        polynomial_order(0)
        {}

        constexpr
        DiscreteFieldComponent(
                Body
                structure_arg,
                Indx
                polynomial_order_arg
        )
        :
        structure(structure_arg),
        polynomial_order(polynomial_order_arg)
        {}

        Body structure;

        Indx polynomial_order;

    };

    template<typename... C>
    struct DiscreteFieldDescription
    {

        static_assert(UniqueType<DiscreteFieldComponent, C...>::value);

        constexpr
        DiscreteFieldDescription()
        :
        field(Field::Scalar),
        discrete_field_components()
        {}

        constexpr
        DiscreteFieldDescription(
                Field
                field_arg,
                C...
                discrete_field_component_args
        )
        :
        field(field_arg),
        discrete_field_components({std::forward<DiscreteFieldComponent>(discrete_field_component_args)...})
        {}

        Field field;

        Array<DiscreteFieldComponent, sizeof...(C)> discrete_field_components;

    };

    struct MappingDescription
    {

        constexpr
        MappingDescription()
        :
        mapping(Mapping::Identity),
        polynomial_order(-1)
        {}

        constexpr
        MappingDescription(
                Mapping
                operator_type_arg
        )
        :
        mapping(operator_type_arg),
        polynomial_order(-1)
        {}

        constexpr
        MappingDescription(
                Mapping
                operator_type_arg,
                Indx
                polynomial_order_arg
        )
        :
        mapping(operator_type_arg),
        polynomial_order(polynomial_order_arg)
        {}

        Mapping mapping;

        Intg polynomial_order;

    };

    template<typename F, typename... M>
    struct FiniteElementDescription2
    {

        static_assert(UniqueType<MappingDescription, M...>::value);

        using DiscreteFieldDescription = F;

        constexpr
        FiniteElementDescription2(
                Char
                label_arg,
                FiniteElementType
                finite_element_type_arg,
                F
                discrete_field_description_arg,
                M...
                mapping_args
        )
        :
        label(label_arg),
        finite_element_type(finite_element_type_arg),
        discrete_field_description(discrete_field_description_arg),
//        mappings({static_cast<MappingDescription>(mapping_args)...})
        mappings({std::forward<MappingDescription>(mapping_args)...})
        {}

        Char label;

        FiniteElementType finite_element_type;

        DiscreteFieldDescription discrete_field_description;

        Array<MappingDescription, sizeof...(M)> mappings;

        constexpr
        Indx
        getIntegrationOrder()
        const
        {
            auto integration_order_arg = Indx(0);
            for (Indx i = 0; i < sizeof...(M); ++i) {
                auto const operator_polynomial_order = mappings(i).polynomial_order;
                if (integration_order_arg < operator_polynomial_order) {
                    integration_order_arg = operator_polynomial_order;
                }
            }
            return integration_order_arg;
        }

    };

    template<FiniteElementType F, typename... M>
    struct Setter;

    template<typename... M>
    struct Setter<FiniteElementType::HybridHighOrder, M...>
    {

        static constexpr
        auto
        set(
                Field
                field_arg,
                Char
                label_arg,
                Indx
                cell_order_arg,
                Indx
                face_order_arg,
                M ...
                operator_order_args
        )
        {
            auto cell_field = DiscreteFieldComponent(Body::Cell, cell_order_arg);
            auto face_field = DiscreteFieldComponent(Body::Face, face_order_arg);
            auto fff = DiscreteFieldDescription(field_arg, cell_field, face_field);
//            auto mmm = Array<MappingDescription, sizeof...(M)>{MappingDescription(operator_order_args, 1)...};
            auto gg = FiniteElementDescription2(label_arg, FiniteElementType::HybridHighOrder, fff, operator_order_args...);
        }

        Setter()
        {}

        Setter(
                Indx
                cell_order_arg,
                Indx
                face_order_arg,
                auto ...
                operator_order_args
        )
        {}

    };

    template<auto A, typename F, typename... O>
    struct FiniteElementDescription3;

    template<typename F, typename... O>
    struct FiniteElementDescription3<FiniteElementType::HybridHighOrder, F, O...>
    {

        constexpr
        FiniteElementDescription3() = default;

        constexpr
        FiniteElementDescription3(
                F
                f,
                O...
                o
        )
        {}



    };

    template<typename F, typename... O>
    struct FiniteElementDescription3<FiniteElementType::Lagrange, F, O...>
    {



    };

    template<
            FiniteElementType FiniteElementTypeArg,
            Field FieldTypeArg,
            Mapping ...OperatorTypeArgs
    >
    struct FiniteElementDescription;



    template<
            FiniteElementType FiniteElementTypeArg,
            Field FieldTypeArg,
            Mapping ...OperatorTypeArgs
    >
    struct FiniteElementDescriptionBase
    {

//        Array<FieldDescription, sizeof...(OperatorTypeArgs)> mapped_fields;
        FieldDescription primal_field;

        using OperatorTypesT = Array<Mapping, sizeof...(OperatorTypeArgs)>;

        using OperatorDescriptionsT = Array<MappingDescription, sizeof...(OperatorTypeArgs)>;

        constexpr
        FiniteElementDescriptionBase(
                Char
                label_arg,
                Indx
                dim_euclidean_arg,
                auto ...
                operator_polynomial_order_args
        )
        :
        primal_field(FieldTypeArg, dim_euclidean_arg),
        label(label_arg),
        dim_euclidean(dim_euclidean_arg),
        operator_descriptions(setOperatorDescriptions(operator_polynomial_order_args...))
        {}

        static constexpr
        Indx
        getNumOperators()
        {
            return sizeof...(OperatorTypeArgs);
        }

        static constexpr
        OperatorTypesT
        getOperatorTypes()
        {
            return {OperatorTypeArgs...};
        }

        static constexpr
        FiniteElementType
        getFiniteElementType()
        {
            return FiniteElementTypeArg;
        }

        constexpr
        FieldDescription
        getFieldType()
        const
        {
            return primal_field;
        }

        constexpr
        OperatorDescriptionsT const &
        getOperatorDescriptions()
        const
        {
            return operator_descriptions;
        }

        constexpr
        Indx
        getDimEuclidean()
        const
        {
            return dim_euclidean;
        }

        constexpr
        Indx
        getIntegrationOrder()
        const
        {
            Indx integration_order_arg = 0;
            for (Indx i = 0; i < getNumOperators(); ++i) {
                auto const operator_polynomial_order = operator_descriptions(i).getPolynomialOrder();
                if (integration_order_arg < operator_polynomial_order) {
                    integration_order_arg = operator_polynomial_order;
                }
            }
            return integration_order_arg;
        }

        constexpr
        Char
        getLabel()
        const
        {
            return label;
        }

    private:

        using OperatorOrdersT = Array<Indx, getNumOperators()>;

        static constexpr
        OperatorOrdersT
        getOperatorOrders(
                auto ...
                operator_polynomial_order_args
        )
        {
            return OperatorOrdersT{static_cast<Indx>(operator_polynomial_order_args)...};
        }

        static constexpr
        OperatorDescriptionsT
        setOperatorDescriptions(
                auto ...
                operator_polynomial_order_args
        )
        {
            static_assert(sizeof...(operator_polynomial_order_args) == getNumOperators());
            OperatorDescriptionsT operator_descriptions_arg;
            auto const operator_orders = getOperatorOrders(operator_polynomial_order_args...);
            auto const operator_types = getOperatorTypes();
            for (Indx i = 0; i < getNumOperators(); ++i) {
                operator_descriptions_arg(i) = MappingDescription(operator_types(i), operator_orders(i));
            }
            return operator_descriptions_arg;
        }

    public:

        Char label;

        Indx dim_euclidean;

        OperatorDescriptionsT operator_descriptions;

    };

}

#endif //FETA_00_FINITE_ELEMENT_DESCRIPTION_BASE_HXX
