//
// Created by dsiedel on 26/02/2022.
//

#ifndef FETA_00_FINITE_ELEMENT_DESCRIPTION_HHO_HXX
#define FETA_00_FINITE_ELEMENT_DESCRIPTION_HHO_HXX

//#include "feta/core/finite_element_discrete_fields.hxx"
//#include "feta/core/operator_output.hxx"
//#include "feta/core/finite_element_type.hxx"
//#include "feta/core/field_type.hxx"
//#include "new/00_finite_element_description_base.hxx"
#include "new/feta_core_field_description.hxx"

namespace feta::core::element
{

    template<
            Field FieldTypeArg,
            Mapping ...OperatorTypeArgs
    >
    struct FiniteElementDescription<
            FiniteElementType::HybridHighOrder,
            FieldTypeArg,
            OperatorTypeArgs...
    >
            :
                    public FiniteElementDescriptionBase<
                            FiniteElementType::HybridHighOrder,
                            FieldTypeArg,
                            OperatorTypeArgs...
                    >
    {

        using Base =  FiniteElementDescriptionBase<
                FiniteElementType::HybridHighOrder,
                FieldTypeArg,
                OperatorTypeArgs...
        >;

        constexpr
        FiniteElementDescription(
                Char
                label_arg,
                Indx
                dim_euclidean_arg,
                Indx
                face_polynomial_order_arg,
                Indx
                cell_polynomial_order_arg,
                auto ...
                operator_polynomial_order_args
        )
                :
                Base(label_arg, dim_euclidean_arg, operator_polynomial_order_args...),
                face_polynomial_order(face_polynomial_order_arg),
                cell_polynomial_order(cell_polynomial_order_arg)
        {}

        constexpr
        Indx
        getFacePolynomialOrder()
        const
        {
            return face_polynomial_order;
        }

        constexpr
        Indx
        getCellPolynomialOrder()
        const
        {
            return cell_polynomial_order;
        }

        Indx face_polynomial_order;

        Indx cell_polynomial_order;

    };

//    struct OperatorDescription
//    {
//
//        constexpr
//        OperatorDescription() = default;
//
//        constexpr
//        OperatorDescription(
//                OperatorType
//                operator_type_arg,
//                Indx
//                polynomial_order_arg
//        )
//        :
//        operator_type(operator_type_arg),
//        polynomial_order(polynomial_order_arg)
//        {}
//
//        constexpr inline
//        Indx
//        getPolynomialOrder()
//        const
//        {
//            return polynomial_order;
//        }
//
//        constexpr inline
//        OperatorType
//        getOperatorType()
//        const
//        {
//            return operator_type;
//        }
//
//        OperatorType operator_type;
//
//        Indx polynomial_order;
//
//    };
//
//    template<
//            FiniteElementType FiniteElementTypeArg,
//            FieldType FieldTypeArg,
//            OperatorType ...OperatorTypeArgs
//    >
//    struct FiniteElementDescriptionBase
//    {
//
//        using OperatorTypesT = StaticArray<
//                OperatorType,
//                sizeof...(OperatorTypeArgs)
//        >;
//
//        using OperatorDescriptionsT = StaticArray<
//                OperatorDescription,
//                sizeof...(OperatorTypeArgs)
//        >;
//
//        constexpr
//        FiniteElementDescriptionBase(
//                Char
//                label_arg,
//                Indx
//                dim_euclidean_arg,
//                auto ...
//                operator_polynomial_order_args
//        )
//        :
//                label(label_arg),
//                dim_euclidean(dim_euclidean_arg),
//                operator_descriptions(setOperatorDescriptions(operator_polynomial_order_args...))
//        {}
//
//        static constexpr
//        Indx
//        getNumOperators()
//        {
//            return sizeof...(OperatorTypeArgs);
//        }
//
//        static constexpr
//        OperatorTypesT
//        getOperatorTypes()
//        {
//            return {OperatorTypeArgs...};
//        }
//
//        static constexpr
//        FiniteElementType
//        getFiniteElementType()
//        {
//            return FiniteElementTypeArg;
//        }
//
//        static constexpr
//        FieldType
//        getFieldType()
//        {
//            return FieldTypeArg;
//        }
//
//        constexpr
//        OperatorDescriptionsT const &
//        getOperatorDescriptions()
//        const
//        {
//            return operator_descriptions;
//        }
//
//        constexpr
//        Indx
//        getDimEuclidean()
//        const
//        {
//            return dim_euclidean;
//        }
//
//        constexpr
//        Indx
//        getIntegrationOrder()
//        const
//        {
//            Indx integration_order_arg = 0;
//            for (Indx i = 0; i < getNumOperators(); ++i) {
//                auto const operator_polynomial_order = operator_descriptions(i).getPolynomialOrder();
//                if (integration_order_arg < operator_polynomial_order) {
//                    integration_order_arg = operator_polynomial_order;
//                }
//            }
//            return integration_order_arg;
//        }
//
//        constexpr
//        Char
//        getLabel()
//        const
//        {
//            return label;
//        }
//
//    private:
//
//        using OperatorOrdersT = StaticArray<
//                Indx,
//                getNumOperators()
//        >;
//
//        static constexpr
//        OperatorOrdersT
//        getOperatorOrders(
//                auto ...
//                operator_polynomial_order_args
//        )
//        {
//            static_assert(sizeof...(operator_polynomial_order_args) == getNumOperators());
//            return OperatorOrdersT{static_cast<Indx>(operator_polynomial_order_args)...};
//        }
//
//        static constexpr
//        OperatorDescriptionsT
//        setOperatorDescriptions(
//                auto ...
//                operator_polynomial_order_args
//        )
//        {
//            OperatorDescriptionsT operator_descriptions_arg;
//            auto const operator_orders = getOperatorOrders(operator_polynomial_order_args...);
//            auto const operator_types = getOperatorTypes();
//            for (Indx i = 0; i < getNumOperators(); ++i) {
//                operator_descriptions_arg(i) = OperatorDescription(operator_types(i), operator_orders(i));
//            }
//            return operator_descriptions_arg;
//        }
//
//    public:
//
//        Char label;
//
//        Indx dim_euclidean;
//
//        OperatorDescriptionsT operator_descriptions;
//
//    };
//
//    template<
//            FiniteElementType FiniteElementTypeArg,
//            FieldType FieldTypeArg,
//            OperatorType ...OperatorTypeArgs
//    >
//    struct FiniteElementDescription;

//    template<
//            Indx DimEuclideanArg,
//            FiniteElementType FiniteElementTypeArg,
//            FieldType FieldTypeArg,
//            Indx NumMeshCoreUnknownArg,
//            Indx NumCellCoreUnknownArg,
//            Indx NumCellPlugUnknownArg,
//            OperatorType... OperatorTypeArg
//    >
//    struct FiniteElementDescriptionBase
//    {
//
//        using OperatorTypesT = StaticArray<
//                OperatorType,
//                sizeof...(OperatorTypeArg)
//        >;
//
//        using OperatorOrdersT = StaticArray<
//                Indx,
//                sizeof...(OperatorTypeArg)
//        >;
//
//        using MeshCoreDiscreteFieldDescriptionT = DiscreteFieldDescription<
//                DiscreteFieldSupportType::Mesh,
//                DiscreteFieldUnknownType::Core
//        >;
//
//        using CellCoreDiscreteFieldDescriptionT = DiscreteFieldDescription<
//                DiscreteFieldSupportType::Cell,
//                DiscreteFieldUnknownType::Core
//        >;
//
//        using CellPlugDiscreteFieldDescriptionT = DiscreteFieldDescription<
//                DiscreteFieldSupportType::Cell,
//                DiscreteFieldUnknownType::Plug
//        >;
//
//        using FiniteElementDiscreteFieldsT = FiniteElementDiscreteFields<
//                NumMeshCoreUnknownArg,
//                NumCellCoreUnknownArg,
//                NumCellPlugUnknownArg,
//                sizeof...(OperatorTypeArg)
//        >;
//
//        constexpr
//        FiniteElementDescriptionBase()
//        :
//        finite_element_discrete_fields(FiniteElementDiscreteFieldsT{}),
//        label('N')
//        {}
//
//        constexpr
//        FiniteElementDescriptionBase(
//                FiniteElementDiscreteFieldsT const &
//                finite_element_discrete_fields_arg,
//                Char
//                label_arg
//        )
//        :
//        finite_element_discrete_fields(finite_element_discrete_fields_arg),
//        label(label_arg)
//        {}
//
//        constexpr
//        FiniteElementDescriptionBase(
//                FiniteElementDiscreteFieldsT &&
//                finite_element_discrete_fields_arg,
//                Char
//                label_arg
//        )
//                :
//                finite_element_discrete_fields(finite_element_discrete_fields_arg),
//                label(label_arg)
//        {}
//
//        static constexpr
//        Indx
//        getDimEuclidean()
//        {
//            return DimEuclideanArg;
//        }
//
//        static constexpr
//        Indx
//        getNumMeshCoreUnknown()
//        {
//            return NumMeshCoreUnknownArg;
//        }
//
//        static constexpr
//        Indx
//        getNumCellCoreUnknown()
//        {
//            return NumCellCoreUnknownArg;
//        }
//
//        static constexpr
//        Indx
//        getNumCellPlugUnknown()
//        {
//            return NumCellPlugUnknownArg;
//        }
//
//        /**
//         * @brief
//         * @return
//         */
//        static constexpr
//        Indx
//        getNumOperators()
//        {
//            return sizeof...(OperatorTypeArg);
//        }
//
//        /**
//         * @brief
//         * @return
//         */
//        static constexpr
//        OperatorTypesT
//        getOperatorTypes()
//        {
//            return {OperatorTypeArg...};
//        }
//
//        /**
//         * @brief
//         * @return
//         */
//        static constexpr
//        FiniteElementType
//        getFiniteElementType()
//        {
//            return FiniteElementTypeArg;
//        }
//
//        /**
//         * @brief
//         * @return
//         */
//        static constexpr
//        FieldType
//        getFieldType()
//        {
//            return FieldTypeArg;
//        }
//
//        FiniteElementDiscreteFieldsT finite_element_discrete_fields;
//
//        Char label;
//
//    };
//
//    template<
//            Indx DimEuclideanArg,
//            FiniteElementType FiniteElementTypeArg,
//            FieldType FieldTypeArg,
//            OperatorType... OperatorTypeArg
//    >
//    struct FiniteElementDescription;
//
//    template<
//            Indx DimEuclideanArg,
//            FieldType FieldTypeArg,
//            OperatorType... OperatorTypeArg
//    >
//    struct FiniteElementDescription<
//            DimEuclideanArg,
//            FiniteElementType::HybridHighOrder,
//            FieldTypeArg,
//            OperatorTypeArg...
//    >
//    :
//    public
//    FiniteElementDescriptionBase<
//            DimEuclideanArg,
//            FiniteElementType::HybridHighOrder,
//            FieldTypeArg,
//            1,
//            1,
//            sizeof...(OperatorTypeArg) + 1,
//            OperatorTypeArg...
//    >
//    {
//
//        using BaseT = FiniteElementDescriptionBase<
//                DimEuclideanArg,
//                FiniteElementType::HybridHighOrder,
//                FieldTypeArg,
//                1,
//                1,
//                sizeof...(OperatorTypeArg) + 1,
//                OperatorTypeArg...
//        >;
//
//        constexpr
//        FiniteElementDescription(
//                Char
//                label_arg,
//                Indx
//                face_order_arg,
//                Indx
//                cell_order_arg,
//                auto ...
//                operator_orders_arg
//        )
//        :
//        BaseT(getFiniteElementDescription(face_order_arg, cell_order_arg, operator_orders_arg...), label_arg)
//        {}
//
//    private:
//
//        template<
//                Indx I = 0
//        >
//        static constexpr
//        void
//        setFiniteElementDescription(
//                typename BaseT::FiniteElementDiscreteFieldsT &
//                finite_element_discrete_fields_arg,
//                auto ...
//                operator_orders_arg
//        )
//        {
//            typename BaseT::OperatorTypesT constexpr operator_types = {OperatorTypeArg...};
//            typename BaseT::OperatorOrdersT operator_orders = {static_cast<Indx>(operator_orders_arg)...};
//            using OperatorFieldOutputT = OperatorOutput<
//                    DimEuclideanArg,
//                    FieldTypeArg,
//                    operator_types(I)
//            >;
//            finite_element_discrete_fields_arg.cell_plug_unknown_discrete_field_descriptions(I + 1) = typename BaseT::CellPlugDiscreteFieldDescriptionT{
//                    .dim_euclidean = OperatorFieldOutputT::getDimField(),
//                    .field_type = OperatorFieldOutputT::getFieldType(),
//                    .polynomial_order = operator_orders(I),
//                    .basis_type = BasisType::Monomial,
//                    .cell_layer = 0,
//            };
//            finite_element_discrete_fields_arg.operator_descriptions(I) = DiscreteOperatorDescription{
//                    .operator_type = operator_types(I),
//                    .polynomial_order = operator_orders(I),
//            };
//            if constexpr (I < BaseT::getNumOperators() - 1) {
//                setFiniteElementDescription<I + 1>(finite_element_discrete_fields_arg, operator_orders_arg...);
//            }
//        }
//
//        static constexpr
//        typename BaseT::FiniteElementDiscreteFieldsT
//        getFiniteElementDescription(
//                Indx
//                face_order_arg,
//                Indx
//                cell_order_arg,
//                auto ...
//                operator_orders_arg
//        )
//        {
//            static_assert(sizeof...(operator_orders_arg) == BaseT::getNumOperators());
//            typename BaseT::FiniteElementDiscreteFieldsT finite_element_discrete_fields_arg;
//            finite_element_discrete_fields_arg.mesh_core_unknown_discrete_field_descriptions(0) = typename BaseT::MeshCoreDiscreteFieldDescriptionT{
//                    .dim_euclidean = DimEuclideanArg,
//                    .field_type = FieldTypeArg,
//                    .polynomial_order = face_order_arg,
//                    .basis_type = BasisType::Monomial,
//                    .cell_layer = 1,
//                    .face_layer = 0,
//                    .mesh_layer = DimEuclideanArg - 1,
//            };
//            finite_element_discrete_fields_arg.cell_core_unknown_discrete_field_descriptions(0) = typename BaseT::CellCoreDiscreteFieldDescriptionT{
//                    .dim_euclidean = DimEuclideanArg,
//                    .field_type = FieldTypeArg,
//                    .polynomial_order = cell_order_arg,
//                    .basis_type = BasisType::Monomial,
//                    .cell_layer = 0,
//            };
//            finite_element_discrete_fields_arg.cell_plug_unknown_discrete_field_descriptions(0) = typename BaseT::CellPlugDiscreteFieldDescriptionT{
//                    .dim_euclidean = DimEuclideanArg,
//                    .field_type = FieldTypeArg,
//                    .polynomial_order = face_order_arg + 1,
//                    .basis_type = BasisType::Monomial,
//                    .cell_layer = 0,
//            };
//            setFiniteElementDescription(finite_element_discrete_fields_arg, operator_orders_arg...);
//            return finite_element_discrete_fields_arg;
//        }
//
//    };

}

#endif //FETA_00_FINITE_ELEMENT_DESCRIPTION_HHO_HXX
