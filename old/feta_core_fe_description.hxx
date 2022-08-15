//
// Created by dsiedel on 03/04/2022.
//

#ifndef FETA_FETA_CORE_FE_DESCRIPTION_HXX
#define FETA_FETA_CORE_FE_DESCRIPTION_HXX

#include "new/feta_core_model_description.hxx"
//#include "new/finite_element_type.hxx"
#include "new/feta_core_field_description.hxx"
#include "new/feta_core_element_quadrature_description.hxx"
#include "new/feta_core_element_basis_description.hxx"
#include "new/feta_core_mapping_description.hxx"
#include "new/feta_core_model_description.hxx"

namespace feta::core
{

    enum struct FiniteElementMethod
    {

        Lagrange,
        HybridHighOrder

    };

    auto const static constexpr fem_hho = FiniteElementMethod::HybridHighOrder;
    auto const static constexpr fem_lag = FiniteElementMethod::Lagrange;

    enum struct EuclideanFrame
    {

        Cartesian,
        AxiSymmetric

    };

    template<typename... T>
    struct DiscretizationDescription
    {

        static_assert(UniqueType<Mapping, T...>::value && sizeof...(T) > 0);

        constexpr
        DiscretizationDescription(
                Char &&
                label_arg,
                Field &&
                field_arg,
                FiniteElementMethod &&
                finite_element_method_arg,
                T &&...
                mapping_args
        )
        :
                label(label_arg),
                field(field_arg),
                finite_element_method(finite_element_method_arg),
                mappings({MappingDescription(mapping_args)...})
        {}

        constexpr
        Bool
        operator==(
                DiscretizationDescription const &
                other
        )
        const
        {
            auto eq_0 = label == other.label;
            auto eq_1 = field == other.field;
            auto eq_2 = finite_element_method == other.finite_element_method;
            auto eq_3 = mappings == other.mappings;
            return eq_0 && eq_1 && eq_2 && eq_3;
        }

        constexpr
        Bool
        operator!=(
                DiscretizationDescription const &
                other
        )
        const
        {
            return !(other == * this);
        }

        Char label;

        Field field;

        FiniteElementMethod finite_element_method;

        Array<MappingDescription, sizeof...(T)> mappings;

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
                    elt::Model &&
                    model_arg,
                    Indx &&
                    cell_order_arg,
                    Indx &&
                    face_order_arg,
                    auto &&...
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
                    elt::Model &&
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
                    if (mapping_description.mapping == Mapping::Identity) {
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

    template<auto... finite_element_args>
    struct MixedElementDescription : public Aggregate<std::remove_cvref_t<decltype(finite_element_args)>...>
    {

        template<template<auto> typename T>
        using Elements = Collection<T<finite_element_args>...>;

    private:

        using FiniteElementAggregate = Aggregate<std::remove_cvref_t<decltype(finite_element_args)>...>;

    public:

//        auto const static constexpr finite_elements = FiniteElementAggregate{finite_element_args...};

        constexpr
        MixedElementDescription(
                elt::Model &&
                model_arg,
                Indx &&
                dim_euclidean_arg,
                EuclideanFrame &&
                frame_type_arg,
                Quadrature &&
                quadrature_type_arg
        )
        :
        FiniteElementAggregate({finite_element_args...}),
        model(model_arg),
        dim_euclidean(dim_euclidean_arg),
        frame_type(frame_type_arg),
        quadrature_description(quadrature_type_arg, setIntegrationOrder(frame_type_arg))
        {}

        constexpr
        Bool
        operator==(
                MixedElementDescription const &
                other
        )
        const
        {
            auto eq_0 = model == other.model;
            auto eq_1 = dim_euclidean == other.dim_euclidean;
            auto eq_2 = frame_type == other.frame_type;
            auto eq_3 = quadrature_description == other.quadrature_description;
            return eq_0 && eq_1 && eq_2 && eq_3;
        }

        constexpr
        Bool
        operator!=(
                MixedElementDescription const &
                other
        )
        const
        {
            return !(other == * this);
        }

        static constexpr
        auto
        setIntegrationOrder(
                EuclideanFrame const &
                frame_arg
        )
        {
            auto val = max(finite_element_args.getIntegrationOrder()...);
            if (frame_arg == EuclideanFrame::AxiSymmetric) {
                val += 1;
            }
            return val;
        }


        QuadratureDescription quadrature_description;

        Indx dim_euclidean;

        EuclideanFrame frame_type;

        elt::Model model;

        // une loi de comportement

    };

    template<auto... A>
    struct CoupledElement
    {

        // un point de gauss, plusiuers lois de comportment

    };

    template<auto A>
    struct MeshInteriorDomain
    {



    };

    template<typename T>
    struct DDomain
    {

        //

    };

    template<auto... A>
    struct Mesh
    {



    };
//
//    template<typename... T>
//    struct DiscreteElementDescription
//    {
//
//        static_assert(UniqueType<Mapping, T...>::value && sizeof...(T) > 0);
//
//        auto const static constexpr num_mappings = sizeof...(T);
//
//        constexpr
//        DiscreteElementDescription()
//                :
//                label(),
//                field(),
//                discretization(),
//                mappings()
//        {}
//
//        constexpr
//        DiscreteElementDescription(
//                Char
//                label_arg,
//                Field
//                field_arg,
//                Discretization
//                discretization_arg,
//                T...
//                mapping_args
//        )
//                :
//                label(label_arg),
//                field(field_arg),
//                discretization(discretization_arg),
//                mappings(setMapDesc(mapping_args...))
//        {}
//
//        static constexpr
//        auto
//        setMapDesc(
//                T...
//                mapping_args
//        )
//        {
//            auto mappingslabs = Array<Mapping, num_mappings>{mapping_args...};
//            auto mappingss = Array<MappingDescription, num_mappings>();
//            for (Indx i = 0; i < num_mappings; ++i) {
//                mappingss.get(i) = MappingDescription(mappingslabs.get(i));
//            }
//            return mappingss;
//        }
//
//        Char label;
//
//        Field field;
//
//        Discretization discretization;
//
//        Array<MappingDescription, num_mappings> mappings;
//
//    };
//
//    namespace detail
//    {
//
//        template<Discretization, auto>
//        struct DiscreteFieldPolicy;
//
//        template<auto A>
//        struct DiscreteFieldPolicy<Discretization::HybridHighOrder, A> : public RawType<decltype(A)>
//        {
//
//        private:
//
//            using Self = DiscreteFieldPolicy;
//
//        public:
//
//            constexpr
//            DiscreteFieldPolicy()
//                    :
//                    RawType<decltype(A)>(setBase(0, 0)),
//                    cell_order(),
//                    face_order()
////                    mapping_orders()
//                    {}
//
//            constexpr
//            DiscreteFieldPolicy(
//                    Indx
//                    cell_order_arg,
//                    Indx
//                    face_order_arg,
//                    auto &&...
//                    mapping_order_args
//            )
//                    :
//                    RawType<decltype(A)>(setBase(mapping_order_args...)),
//                    cell_order(cell_order_arg),
//                    face_order(face_order_arg)
//            {}
////                    mapping_orders(setMappingOrders(mapping_order_args...)) {}
//
//            constexpr
//            DiscreteFieldPolicy(
//                    Indx
//                    cell_order_arg,
//                    Indx
//                    face_order_arg
//            )
//                    :
//                    RawType<decltype(A)>(setBase(cell_order_arg, face_order_arg)),
//                    cell_order(cell_order_arg),
//                    face_order(face_order_arg)
//            {}
////                    mapping_orders(setMappingOrders(cell_order_arg, face_order_arg)) {}
//
//        private:
//
//            static constexpr
//            auto
//            setBase(
//                    Indx
//                    cell_order_arg,
//                    Indx
//                    face_order_arg
//            )
//            {
//                auto aobj = A;
//                for (Indx i = 0; i < Self::num_mappings; ++i) {
//                    auto mapping = A.mappings.get(i).mapping;
//                    if (mapping == Mapping::Identity) {
//                        aobj.mappings.get(i) = MappingDescription(mapping, cell_order_arg);
//                    } else {
//                        aobj.mappings.get(i) = MappingDescription(mapping, face_order_arg);
//                    }
//                }
//                return aobj;
//            }
//
//            static constexpr
//            auto
//            setBase(
//                    auto &&...
//                    mapping_order_args
//            )
//            {
//                static_assert(sizeof...(mapping_order_args) == Self::num_mappings);
//                auto valss = Array <Indx, Self::num_mappings > {static_cast<Indx>(mapping_order_args)...};
//                auto aobj = A;
//                for (Indx i = 0; i < Self::num_mappings; ++i) {
//                    aobj.mappings.get(i) = MappingDescription(A.mappings.get(i).mapping, valss.get(i));
//                }
//                return aobj;
//            }
//
////            static constexpr
////            auto
////            setMappingOrders(
////                    auto...
////                    mapping_orders_args
////            )
////            {
////                static_assert(sizeof...(mapping_orders_args) == Self::num_mappings);
////                return Array < Indx, Self::num_mappings > {static_cast<Indx>(mapping_orders_args)...};
////            }
////
////            static constexpr
////            auto
////            setMappingOrders(
////                    Indx
////                    cell_order_arg,
////                    Indx
////                    face_order_arg
////            )
////            {
////                auto mapping_orders_val = Array<Indx, Self::num_mappings>();
////                for (int i = 0; i < Self::num_mappings; ++i) {
////                    auto const mapping = A.mappings.get(i).mapping;
////                    if (mapping == Mapping::Identity) {
////                        mapping_orders_val.get(i) = cell_order_arg;
////                    } else {
////                        mapping_orders_val.get(i) = face_order_arg;
////                    }
////                }
////                return mapping_orders_val;
////            }
//
//        public:
//
//            Indx cell_order;
//
//            Indx face_order;
//
////            Array<Indx, Self::num_mappings> mapping_orders;
//
//        };
//
//    }

//    struct MappingDescription
//    {
//
//        constexpr
//        MappingDescription()
//        :
//        mapping(),
//        order(-1)
//        {}
//
//        constexpr
//        MappingDescription(
//                Mapping
//                mapping_arg,
//                Indx
//                order_arg
//        )
//        :
//        mapping(mapping_arg),
//        order(order_arg)
//        {}
//
//        Mapping mapping;
//
//        Intg order;
//
//    };

    namespace detail
    {

        template<auto A>
        struct DiscreteFieldBase
        {

//        private:
//
//            auto const static constexpr num_mappings = A.mappings.getSize();
//
//        public:

            auto const static constexpr field = A.discrete_field;

            auto const static constexpr discretization = A.finite_element_method_v;

            auto const static constexpr mappings = A.mappings;

//            constexpr
//            DiscreteFieldBase()
//            :
//            mapping_descriptions()
//            {}
//
//            constexpr
//            DiscreteFieldBase(
//                    auto &&...
//                    mapping_order_args
//            )
//            :
//            mapping_descriptions(setMappings(mapping_order_args...))
//            {}
//
//            constexpr
//            auto
//            getMaxMappingOrder()
//            const
//            {
//                auto val = Indx(0);
//                for (Indx i = 0; i < num_mappings; ++i) {
//                    auto mapping_order = mapping_descriptions.get(i).order;
//                    if (mapping_order > val) {
//                        val = mapping_order;
//                    }
//                }
//                return val;
//            }
//
//        private:
//
//            static constexpr
//            auto
//            setMappings(
//                    auto const &...
//                    mapping_ord_args
//            )
//            {
//                static_assert(sizeof...(mapping_ord_args) == num_mappings);
//                auto mapping_ords = Array<Indx, num_mappings>{static_cast<Indx>(mapping_ord_args)...};
//                auto mapping_vals = Array<MappingDescription, num_mappings>();
//                for (Indx i = 0; i < num_mappings; ++i) {
//                    mapping_vals.get(i) = MappingDescription(A.mappings.get(i), mapping_ords.get(i));
//                }
//                return mapping_vals;
//            }
//
//        public:
//
//            Array<Mapping, num_mappings> mapping_descriptions;

        };

    }


//    enum struct FiniteElementType
//    {
//
//        HybridHighOrder,
//        Lagrange
//
//    };
//
//    struct UnknownComponentDescription
//    {
//
//        constexpr
//        UnknownComponentDescription()
//        :
//        model_description(elt::Body::Cell, elt::Model::Solid),
//        basis_description(elt::Basis::Monomial, 0)
//        {}
//
//        constexpr
//        UnknownComponentDescription(
//                elt::Body
//                body_arg,
//                elt::Model
//                model_arg,
//                elt::Basis
//                basis_arg,
//                Indx
//                polynomial_order_arg
//        )
//        :
//        model_description(body_arg, model_arg),
//        basis_description(basis_arg, polynomial_order_arg)
//        {}
//
//        elt::ModelDescription model_description;
//
//        elt::BasisDescription basis_description;
//
//    };
//
//    template<typename... C>
//    struct DiscreteFieldDescription
//    {
//
//        static_assert(UniqueType<UnknownComponentDescription, C...>::value);
//
//        constexpr
//        DiscreteFieldDescription()
//        :
//        field(Field::Scalar),
//        unknowns()
//        {}
//
//        constexpr
//        DiscreteFieldDescription(
//                Field
//                field_arg,
//                C...
//                discrete_field_component_args
//        )
//        :
//        field(field_arg),
//        unknowns({std::forward<UnknownComponentDescription>(discrete_field_component_args)...})
//        {}
//
//        Field field;
//
//        Array<UnknownComponentDescription, sizeof...(C)> unknowns;
//
//    };
//
//    struct MappingDescription2
//    {
//
//        constexpr
//        MappingDescription2()
//        :
//        mapping(Mapping::Identity),
//        polynomial_order(-1)
//        {}
//
//        constexpr
//        MappingDescription2(
//                Mapping
//                operator_type_arg
//        )
//        :
//        mapping(operator_type_arg),
//        polynomial_order(-1)
//        {}
//
//        constexpr
//        MappingDescription2(
//                Mapping
//                operator_type_arg,
//                Indx
//                polynomial_order_arg
//        )
//        :
//        mapping(operator_type_arg),
//        polynomial_order(polynomial_order_arg)
//        {}
//
//        Mapping mapping;
//
//        Intg polynomial_order;
//
//    };

//    template<Discretization>
//    struct DiscreteField;
//
//    template<>
//    struct DiscreteField<Discretization::HybridHighOrder>
//    {
//
//        Field f;
//
//        Indx a;
//
//        Indx b;
//
//    };
//
//    template<Mapping>
//    struct DiscreteOperator
//    {
//
//        constexpr
//        DiscreteOperator()
//        :
//        order(-1)
//        {}
//
//        constexpr
//        DiscreteOperator(
//                Indx
//                order_arg
//        )
//        :
//        order(order_arg)
//        {}
//
//        Intg order;
//
//    };
//
//    template<Discretization D, typename ...T>
//    struct Hello
//    {
//
//
//
//    };

//    template<Field F, Discretization D, Indx N>
//    struct DiscretizationDescriptionBase
//    {
//
//        auto const static constexpr field = F;
//
//        auto const static constexpr discretization = D;
//
//        constexpr
//        DiscretizationDescriptionBase()
//        {}
//
//        constexpr
//        DiscretizationDescriptionBase(
//                auto...
//                mapping_order_args
//        )
//        {}
//
//        struct MappingDescription
//        {
//
//            constexpr
//            MappingDescription()
//            :
//            mapping(),
//            order()
//            {}
//
//            constexpr
//            MappingDescription(
//                    Mapping
//                    mapping_arg,
//                    Indx
//                    order_arg
//            )
//            :
//            mapping(mapping_arg),
//            order(order_arg)
//            {}
//
//            Mapping mapping;
//
//            Indx order;
//
//        };
//
//        Array<MappingDescription, N> mappings;
//
//    };

//    template<Field, Discretization, Indx N>
//    struct DiscretizationDescription;
//
//    template<Field F, Indx N>
//    struct DiscretizationDescription<F, hho, N> : public DiscretizationDescriptionBase<F, hho, N>
//    {
//
//        using Base = DiscretizationDescriptionBase<F, hho, N>;
//
//        using MDesc = typename Base::MappingDescription;
//
////        DiscretizationDescription(
////                Indx
////                cell_order_arg,
////                Indx
////                face_order_arg,
////                auto...
////                mapping_order_args
////
////        )
////        :
////        DiscretizationDescriptionBase<F, hho, N>()
////        {}
//
//        Indx cell_order;
//
//        Indx face_order;
//
//        static constexpr
//        auto
//        setG(
//                auto...
//                mapping_order_args
//        )
//        {
//            static_assert(sizeof...(mapping_order_args) == N);
//            auto mappings = Array<MDesc, N>();
//            for (int i = 0; i < N; ++i) {
//                mappings.get(i) = MDesc();
//            }
//        }
//
//    };
//
//    template<auto M>
//    using DiscretizationD = DiscretizationDescription<M.field, M.discretization, M.num_mappings>;



//
//    template<typename F, typename... M>
//    struct FiniteElementDescription2
//    {
//
//        static_assert(UniqueType<MappingDescription, M...>::value);
//
//        constexpr
//        FiniteElementDescription2(
//                Char
//                label_arg,
//                FiniteElementType
//                finite_element_type_arg,
//                F
//                discrete_field_description_arg,
//                M...
//                mapping_args
//        )
//                :
//                label(label_arg),
//                finite_element_type(finite_element_type_arg),
//                discrete_field_description(discrete_field_description_arg),
//                mappings({std::forward<MappingDescription>(mapping_args)...}) {}
//
//        Char label;
//
//        FiniteElementType finite_element_type;
//
//        F discrete_field_description;
//
//        Array<MappingDescription, sizeof...(M)> mappings;
//
//        constexpr
//        Indx
//        getIntegrationOrder()
//        const {
//            auto integration_order_arg = Indx(0);
//            for (Indx i = 0; i < sizeof...(M); ++i) {
//                auto const operator_polynomial_order = mappings(i).polynomial_order;
//                if (integration_order_arg < operator_polynomial_order) {
//                    integration_order_arg = operator_polynomial_order;
//                }
//            }
//            return integration_order_arg;
//        }
//
//    };

//    template<typename... M>
//    struct Mappings : public Array<MappingDescription2, sizeof...(M)>
//    {
//
//        using Data = typename UniqueType<MappingDescription2, M...>::Type;
//
////        static_assert(UniqueType<Mapping, M...>::value);
//
////        constexpr
////        Mappings()
////        :
////        Array<Mapping, sizeof...(M)>()
////        {}
//
//        constexpr
//        Mappings(
//                M &&...
//                mapping_args
//        )
//                :
//                Array<Data, sizeof...(M)>({mapping_args...})
//        {}
//
//        constexpr
//        Mappings(
//                M &&...
//                mapping_args,
//                Indx
//                a
//        )
//                :
//                Array<Data, sizeof...(M)>({mapping_args...})
//        {}
//
//
//    };

//    template<FiniteElementType F, typename... M>
//    struct Setter;
//
//    template<typename... M>
//    struct Setter<FiniteElementType::HybridHighOrder, M...>
//    {
//
//        Indx cell_order;
//
//        Indx face_order;
//
//        Field field;
//
//        Array<MappingDescription, sizeof...(M)> mappings;
//
//        static constexpr
//        auto
//        set(
//                Field
//                field_arg,
//                Char
//                label_arg,
//                Indx
//                cell_order_arg,
//                Indx
//                face_order_arg,
//                M ...
//                operator_order_args
//        )
//        {
////            auto cell_field = DiscreteFieldComponent(elt::Body::Cell, cell_order_arg);
////            auto face_field = DiscreteFieldComponent(elt::Body::Face, face_order_arg);
////            auto fff = DiscreteFieldDescription(field_arg, cell_field, face_field);
//////            auto mmm = Array<MappingDescription, sizeof...(M)>{MappingDescription(operator_order_args, 1)...};
////            auto gg = FiniteElementDescription2(label_arg, FiniteElementType::HybridHighOrder, fff,
////                                                operator_order_args...);
//        }
//
//        Setter() {}
//
//        Setter(
//                Indx
//                cell_order_arg,
//                Indx
//                face_order_arg,
//                auto ...
//                operator_order_args
//        ) {}
//
//    };

}

#endif //FETA_FETA_CORE_FE_DESCRIPTION_HXX
