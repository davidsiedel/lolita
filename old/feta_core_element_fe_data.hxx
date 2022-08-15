//
// Created by dsiedel on 06/04/2022.
//

#ifndef FETA_FETA_CORE_ELEMENT_FE_DATA_HXX
#define FETA_FETA_CORE_ELEMENT_FE_DATA_HXX

#include "new/feta_core_element_element_connectivity.hxx"
#include "new/feta_core_element_quadrature_gauss_segment.hxx"
#include "new/feta_core_geometry_geometry.hxx"
#include "new/feta_core_element_basis.hxx"
#include "new/feta_core_model_description.hxx"
#include "new/feta_core_element_fe_degree_of_freedom.hxx"

namespace feta::core::element::finite_element
{

    template<auto F, auto M>
    constexpr inline
    auto
    fieldDescription()
    {
        return FieldDescription(F.discrete_field.field, M.dim_euclidean);
    }

    template<auto F, auto M>
    constexpr inline
    auto
    fieldDescription(
            Mapping
            mapping_arg
    )
    {
        return FieldDescription::fromMapping(fieldDescription<F, M>(), mapping_arg);
    }

    constexpr inline
    auto
    fieldSize(
            FieldDescription
            field_description_arg
    )
    {
        return field_description_arg.tensor_description.size;
    }

    template<ElementDescription E, auto M>
    constexpr inline
    auto
    modelDescription()
    {
        return ModelDescription(E, M.dim_euclidean, M.model);
    }

    template<auto F>
    constexpr inline
    auto
    finiteElementMethod()
    {
        return F.discrete_field.finite_element_method;
    }

    namespace detail
    {

        template<ElementDescription E, ModelDescription Md, FiniteElementMethod D, auto F, auto M>
        struct FiniteElement;

        template<ElementDescription E, FiniteElementMethod D, auto F, auto M>
        struct FiniteElementPolicy
        {

        private:

            using Self = FiniteElement<E, modelDescription<E, M>(), D, F, M>;

            template<ElementDescription Ec>
            using ComponentT = FiniteElement<Ec, modelDescription<Ec, M>(), D, F, M>;

        public:

            using Components = typename ElementReference<E>::template Components<ComponentT>;

            template<Indx I, Indx J>
            using Component = typename Components::template Type<I>::template Type<J>::Type;

            static constexpr
            auto
            mappingArray()
            {
                auto mappings_arg = Array<Mapping, F.discrete_field.mappings.size()>();
                for (int i = 0; i < F.discrete_field.mappings.size(); ++i) {
                    mappings_arg.get(i) = F.discrete_field.mappings.get(i).mapping;
                }
                return mappings_arg;
            }

            static constexpr
            auto
            mappingMatrixSizes()
            {
                auto s = Indx(0);
                for (int i = 0; i < F.discrete_field.mappings.size(); ++i) {
                    auto m = F.discrete_field.mappings.get(i).mapping;
                    s += FieldDescription::fromMapping(fieldDescription<F, M>(), m).tensor_description.size;
                }
                return s;
            }

            auto const static constexpr mappings = mappingArray();

            auto const static constexpr dim_mapping = mappingMatrixSizes();

            struct ImplementationBase : public Element<E, M>
            {

                using Self = ImplementationBase;

                using FiniteElement = typename Self::template Type<M.index(F)>;

                auto
                fetch()
                {
                    return this->template get<M.index(F)>();
                }

                auto
                fetch()
                const
                {
                    return this->template get<M.index(F)>();
                }

                template<Indx I, Indx J>
                auto
                fetch(
                        auto
                        index_arg
                )
                {
                    auto & p = this->components.template get<I>().template get<J>().template get(index_arg).get();
                    return p.template get<M.index(F)>();
                }

                template<Indx I, Indx J>
                auto
                fetch(
                        auto
                        index_arg
                )
                const
                {
                    auto const & p = this->components.template get<I>().template get<J>().template get(index_arg).get();
                    return p.template get<M.index(F)>();
                }

            };

        };

        /*
         * CELL
         */

        template<ElementDescription E, auto F, auto M>
        struct FiniteElement<E, cell_solid, fem_hho, F, M>
        {

        private:

            using Self = FiniteElement;

            auto const static constexpr basis_description = BasisDescription(Basis::Monomial, F.cell_order);

        public:

            using DegreesOfFreedom = finite_element::DegreesOfFreedom<E, fieldDescription<F, M>(), basis_description>;

        private:

            using FiniteElementPolicy = detail::FiniteElementPolicy<E, fem_hho, F, M>;

            template<Indx I>
            static constexpr
            auto
            numUnknowns()
            {
                return FiniteElementPolicy::template Component<0, I>::DegreesOfFreedom::num_unknowns;
            }

            template<Indx I = 0>
            static constexpr
            auto
            numUnknowns(
                    auto &
                    value
            )
            {
                value += numUnknowns<I>() * numComponents<E, 0, I>();
                if constexpr (I < numComponents<E, 0>() - 1) {
                    numUnknowns<I + 1>(value);
                }
            }

            static constexpr
            auto
            numUnknowns()
            {
                auto value = DegreesOfFreedom::num_unknowns;
                numUnknowns(value);
                return value;
            }

            auto const static constexpr num_unknowns = numUnknowns();

        public:

            using Stabilization = Matrix<Real, num_unknowns, num_unknowns>;

//            template<Mapping O>
//            using MappingMatrix = Matrix<Real, mappingMatrixSize<O, F, M>(), num_unknowns>;

            template<Mapping O>
            using MappingMatrix = Matrix<Real, fieldSize(fieldDescription<F, M>(O)), num_unknowns>;

            using Mappings = typename ArrayCollectionWrapper<FiniteElementPolicy::mappings>::template Wrapper<MappingMatrix>;

            FiniteElement()
            :
            degrees_of_freedom(),
            stabilization(Stabilization::Zero())
            {}

            FiniteElement(
                    DegreesOfFreedom const &
                    degrees_of_freedom_arg,
                    Stabilization const &
                    stabilization_arg
            )
            :
            degrees_of_freedom(degrees_of_freedom_arg),
            stabilization(stabilization_arg)
            {}

            FiniteElement(
                    DegreesOfFreedom &&
                    degrees_of_freedom_arg,
                    Stabilization &&
                    stabilization_arg
            )
            :
            degrees_of_freedom(degrees_of_freedom_arg),
            stabilization(stabilization_arg)
            {}

            DegreesOfFreedom degrees_of_freedom;

            Stabilization stabilization;

            struct Implementation : public FiniteElementPolicy::ImplementationBase
            {

                using Self = Implementation;

                template<Mapping O>
                auto
                getMappingMatrix(
//                        auto const &
//                        point_arg
                )
                const
                {
                    return MappingMatrix<O>().setZero();
                }

                auto
                setStabilization()
                const
                {
                    return Stabilization::Identity();
                }

                template<Indx I>
                auto
                getUnknowns(
                        auto &
                        v,
                        auto &
                        in
                )
                const
                {
                    for (int i = 0; i < numComponents<E, 0>(); ++i) {
                        auto const & n = this->template fetch<0, I>(i).degrees_of_freedom.unknowns.coefficients;
                        v.template segment<numUnknowns<I>()>(in) = n;
                        in += numUnknowns<I>();
                    }
                    if constexpr (I < numComponents<E, 0>() - 1) {
                        getUnknowns<I + 1>(v, in);
                    }
                }

                auto
                getUnknowns()
                const
                {
                    auto v = Matrix<Real, num_unknowns>();
                    auto in = Indx(0);
                    auto const & n = this->template fetch().degrees_of_freedom.unknowns.coefficients;
                    v.template segment<DegreesOfFreedom::num_unknowns>(in) = n;
                    in += DegreesOfFreedom::num_unknowns;
                    getUnknowns(v, in);
                    return v;
                }

                auto
                make()
                const
                {
                    auto unknown_index = Indx(0);
                    auto binding_index = Indx(0);
//                    auto u = DegreesOfFreedom(unknown_index, binding_index, {{1, 2}});
                    auto u = DegreesOfFreedom(unknown_index);
                    auto s = setStabilization();
                    return FiniteElement(u, s);
                }

                void
                hey()
                const
                {
                    print("hey a :", num_unknowns);
                }

            };

        };

        /*
         * FACE
         */

        template<ElementDescription E, auto F, auto M>
        struct FiniteElement<E, face_solid, FiniteElementMethod::HybridHighOrder, F, M>
        {

        private:

            auto const static constexpr basis_description = BasisDescription(Basis::Monomial, F.face_order);

        public:

            using DegreesOfFreedom = finite_element::DegreesOfFreedom<E, fieldDescription<F, M>(), basis_description>;

            DegreesOfFreedom degrees_of_freedom;

            struct Implementation : public Element<E, M>
            {

                using Self = Implementation;

                static
                FiniteElement
                make()
                {
                    return FiniteElement();
                }

            };

        };

        /*
         * EDGE
         */

        template<ElementDescription E, auto F, auto M>
        struct FiniteElement<E, edge_solid, FiniteElementMethod::HybridHighOrder, F, M>
        {

            using DegreesOfFreedom = Void;

            struct Implementation : public Element<E, M>
            {

                using Self = Implementation;

                static
                FiniteElement
                make()
                {
                    return FiniteElement();
                }

            };

        };

        /*
         * NODE
         */

        template<ElementDescription E, auto F, auto M>
        struct FiniteElement<E, node_solid, FiniteElementMethod::HybridHighOrder, F, M>
        {

            using DegreesOfFreedom = Void;

            struct Implementation : public Element<E, M>
            {

                using Self = Implementation;

                static
                FiniteElement
                make()
                {
                    return FiniteElement();
                }

            };

        };

    }

    template<ElementDescription E, auto F, auto M>
    using FiniteElement = detail::FiniteElement<E, modelDescription<E, M>(), finiteElementMethod<F>(), F, M>;

}

//template<ElementDescription E, FiniteElementMethod D, auto F, auto M>
//struct FiniteElementPolicy
//{
//
//private:
//
//    using Self = FiniteElement<E, modelDescription<E, M>(), D, F, M>;
//
//    template<ElementDescription Ec>
//    using Other = FiniteElement<Ec, modelDescription<Ec, M>(), D, F, M>;
//
//    using Components = typename ElementReference<E>::template Components<Other>;
//
//    template<Indx I, Indx J>
//    using Component = typename Components::template Type<I>::template Type<J>::Type;
//
//public:
//
//    template<Indx I, Indx J>
//    static constexpr
//    auto
//    numComponentUnknowns()
//    {
//        using ComponentDegreesOfFreedom = typename Component<I, J>::DegreesOfFreedom;
//        if constexpr (!is_same_v<ComponentDegreesOfFreedom, Void>) {
//            return ComponentDegreesOfFreedom::num_unknowns;
//        }
//        else {
//            return Indx(0);
//        }
//    }
//
//private:
//
//    template<Indx I = 0, Indx J = 0>
//    static constexpr
//    auto
//    setNumComponentUnknowns(
//            auto &
//            value
//    )
//    {
//        value += numComponentUnknowns<I, J>() * numComponents<E, I, J>();
//        if constexpr (J < numComponents<E, I>() - 1) {
//            setNumComponentUnknowns<I, J + 1>(value);
//        }
//        if constexpr (I < numComponents<E>() - 1) {
//            setNumComponentUnknowns<I + 1, 0>(value);
//        }
//    }
//
//public:
//
//    static constexpr
//    auto
//    numComponentUnknowns()
//    {
//        auto a = Indx(0);
//        setNumComponentUnknowns(a);
//        return a;
//    }
//
//    static constexpr
//    auto
//    numElementUnknowns()
//    {
//        using ElementDegreesOfFreedom = typename Self::DegreesOfFreedom;
//        if constexpr (!is_same_v<ElementDegreesOfFreedom, Void>) {
//            return ElementDegreesOfFreedom::num_unknowns;
//        }
//        else {
//            return Indx(0);
//        }
//    }
//
//    static constexpr
//    auto
//    mappingArray()
//    {
//        auto mappings_arg = Array<Mapping, F.discrete_field.mappings.size()>();
//        for (int i = 0; i < F.discrete_field.mappings.size(); ++i) {
//            mappings_arg.get(i) = F.discrete_field.mappings.get(i).mapping;
//        }
//        return mappings_arg;
//    }
//
//    auto const static constexpr num_unknowns = numElementUnknowns() + numComponentUnknowns();
//
//    struct ImplementationBase : public Element<E, M>
//    {
//
//        using Self = ImplementationBase;
//
//        using Elem = typename Self::template Type<M.index(F)>;
//
//        auto
//        getIt()
//        {
//            return this->template get<M.index(F)>();
//        }
//
//        auto
//        getIt()
//        const
//        {
//            return this->template get<M.index(F)>();
//        }
//
//        template<Indx I, Indx J>
//        auto
//        getIt(
//                auto
//                index_arg
//        )
//        {
//            auto & p = this->components.template get<I>().template get<J>().template get(index_arg).get();
//            return p.template get<M.index(F)>();
//        }
//
//        template<Indx I, Indx J>
//        auto
//        getIt(
//                auto
//                index_arg
//        )
//        const
//        {
//            auto const & p = this->components.template get<I>().template get<J>().template get(index_arg).get();
//            return p.template get<M.index(F)>();
//        }
//
//    private:
//
//        template<Indx I = 0, Indx J = 0>
//        void
//        setUnknowns(
//                auto &
//                u,
//                auto &
//                i
//        )
//        const
//        {
//            if constexpr (numComponentUnknowns<I, J>() > 0) {
//                for (int j = 0; j < numComponents<E, I, J>(); ++j) {
//                    auto const & fil = getIt<I, J>(j).degrees_of_freedom.unknowns.values;
//                    u.template segment<numComponentUnknowns<I, J>()>(i) = fil;
//                    i += numComponentUnknowns<I, J>();
//                }
//            }
//            if constexpr (J < numComponents<E, I>() - 1) {
//                setUnknowns<I, J + 1>(u, i);
//            }
//            if constexpr (I < numComponents<E>() - 1) {
//                setUnknowns<I + 1, J>(u, i);
//            }
//        }
//
//    public:
//
//        auto
//        getUnknowns()
//        const
//        {
//            auto u = Matrix<Real, num_unknowns>();
//            auto i = Indx(0);
//            if constexpr (numElementUnknowns() > 0) {
//                u.template segment<numElementUnknowns()>(i) = getIt().degrees_of_freedom.unknowns.values;
//                i += numElementUnknowns();
//            }
//            setUnknowns(u, i);
//            return u;
//        }
//
//    };
//
//};

        /*
         *
         */

//        template<ElementDescription E, ModelDescription Md, FiniteElementMethod D, auto F, auto M>
//        struct ElementUnknownsPolicy;
//
//        template<ElementDescription E, auto F, auto M>
//        struct ElementUnknownsPolicy<E, cell_solid, FiniteElementMethod::HybridHighOrder, F, M>
//        {
//
//        private:
//
//            auto const static constexpr field_description = FieldDescription(F.discrete_field.field, M.dim_euclidean);
//
//            auto const static constexpr basis_description = BasisDescription(Basis::Monomial, F.cell_order);
//
//        public:
//
//            using Type = element::Unknowns<field_description, E, basis_description>;
//
//            using Unknowns = element::Unknowns<field_description, E, basis_description>;
//
//            auto const static constexpr num_u = Unknowns::num_unknowns;
//
//        };
//
//        template<ElementDescription E, auto F, auto M>
//        struct ElementUnknownsPolicy<E, face_solid, FiniteElementMethod::HybridHighOrder, F, M>
//        {
//
//        private:
//
//            auto const static constexpr field_description = FieldDescription(F.discrete_field.field, M.dim_euclidean);
//
//            auto const static constexpr basis_description = BasisDescription(Basis::Monomial, F.face_order);
//
//        public:
//
//            using Type = element::Unknowns<field_description, E, basis_description>;
//
//            using Unknowns = element::Unknowns<field_description, E, basis_description>;
//
//        };
//
//        template<ElementDescription E, auto F, auto M>
//        struct ElementUnknownsPolicy<E, edge_solid, FiniteElementMethod::HybridHighOrder, F, M>
//        {
//
//            using Type = Void;
//
//            using Unknowns = Void;
//
//        };
//
//        template<ElementDescription E, auto F, auto M>
//        struct ElementUnknownsPolicy<E, node_solid, FiniteElementMethod::HybridHighOrder, F, M>
//        {
//
//            using Type = Void;
//
//            using Unknowns = Void;
//
//        };
//
//        /*
//         *
//         */
//
//        template<ElementDescription E, ModelDescription Md, FiniteElementMethod D, auto F, auto M>
//        struct ElementImplementation;
//
//        template<ElementDescription E, auto F, auto M>
//        struct ElementImplementation<E, cell_solid, FiniteElementMethod::HybridHighOrder, F, M> : public Element<E, M>
//        {
//
//        private:
//
//            static_assert(M.has(F));
//
//            using Cell = ElementImplementation;
//
//            using CellUnknowns = typename Cell::Unknowns::template Type<M.index(F)>;
//
//            using Faces = typename Cell::Components::template Type<0>;
//
//            template<Indx I>
//            using Face = typename Cell::Components::template Type<0>::template Type<I>::Type::Type;
//
//            template<Indx I>
//            using FaceUnknowns = typename Face<I>::Unknowns::template Type<M.index(F)>;
//
//            template<Indx I = 0>
//            static constexpr
//            void
//            setNumUnknowns(
//                    auto &
//                    num
//            )
//            {
//                num += FaceUnknowns<I>::dim_unknowns * FaceUnknowns<I>::num_unknowns;
//                if constexpr (I < Faces::size() - 1) {
//                    setNumUnknowns<I + 1, 0>(num);
//                }
//            }
//
//            static constexpr
//            auto
//            setNumUnknowns()
//            {
//
//                auto num = Indx(CellUnknowns::dim_unknowns * CellUnknowns::num_unknowns);
//                setNumUnknowns(num);
//                return num;
//            }
//
//        public:
//
//            static constexpr
//            auto
//            getNumUnknowns()
//            {
//                return setNumUnknowns();
//            }
//
//        };
//
//    }
//
//    template<ElementDescription E, auto F, auto M>
//    using ElementUnknowns = typename detail::ElementUnknownsPolicy<E, ModelDescription(E, M.dim_euclidean, M.model), F.discrete_field.finite_element_method, F, M>::Type;
//
//    template<ElementDescription E, auto F, auto M>
//    using ElementImplementation = detail::ElementImplementation<E, ModelDescription(E, M.dim_euclidean, M.model), F.discrete_field.finite_element_method, F, M>;

//    template<ElementDescription E, auto F, auto M>
//    struct ElementData
//    {
//
//    private:
//
//        auto const static constexpr model_description = ModelDescription(E, M.dim_euclidean, M.model);
//
//        auto const static constexpr discretization = F.discrete_field.discretization;
//
//        using ElementDataPolicy = detail::ElementDataPolicy<E, model_description, discretization, F, M>;
//
//        static constexpr
//        auto
//        setNumElementUnknowns()
//        {
//            if constexpr(std::is_same_v<typename ElementDataPolicy::Unknowns, Void>) {
//                return 0;
//            } else {
//                return ElementDataPolicy::Unknowns::size() * ElementDataPolicy::Unknowns::Type::Basis::dim_basis;
////                return 0;
//            }
//        }
//
//    public:
//
//        auto const static constexpr dim_unknowns = setNumElementUnknowns();
//
//        using Unknowns = typename ElementDataPolicy::Unknowns;
//
//        Unknowns unknowns;
//
//    };

//        template<ElementDescription E, auto F, auto M>
//        struct FiniteElementImplementationBase
//        {
//
//            using Element = element::Element<E, M>;
//
//            using FiniteElement = typename Element::template Type<M.index(F)>;
//
//            template<Indx I, Indx J>
//            using ElementComponent = typename Element::Components::template Type<I>::template Type<J>::Type::Type;
//
//            template<Indx I, Indx J>
//            using FiniteElementComponent = typename ElementComponent<I, J>::template Type<M.index(F)>;
//
//        private:
//
//            template<Indx I = 0, Indx J = 0>
//            static constexpr
//            auto
//            setStructuralUnknownSizes(
//                    auto &
//                    structural_unknown_sizes_arg
//            )
//            {
//                if constexpr (I == 0) {
//                    using FU = typename FiniteElement::StructuralUnknowns;
//                    if constexpr (!is_same_v<FU, Void>) {
////                    if constexpr (!TypeTraits<FU>::template are_same<Void>) {
//                        structural_unknown_sizes_arg.get(0) = FU::dim_unknowns;
//                    }
//                    else {
//                        structural_unknown_sizes_arg.get(0) = 0;
//                    }
//                }
//                else {
//                    using FU = typename FiniteElementComponent<I - 1, J>::StructuralUnknowns;
//                    if (J == 0) {
//                        structural_unknown_sizes_arg.get(I) = 0;
//                    }
////                    if constexpr (!TypeTraits<FU>::template are_same<Void>) {
//                    if constexpr (!is_same_v<FU, Void>) {
//                        auto num_faces = Element::Components::template Type<I - 1>::template Type<J>::size();
//                        structural_unknown_sizes_arg.get(I) += num_faces * FU::dim_unknowns;
//                    }
//                    if constexpr (J < Element::Components::template Type<I - 1>::size() - 1) {
//                        setStructuralUnknownSizes<I, J + 1>(structural_unknown_sizes_arg);
//                    }
//                }
//                if constexpr (I < Element::Components::size()) {
//                    setStructuralUnknownSizes<I + 1, 0>(structural_unknown_sizes_arg);
//                }
//            }
//
//            template<Indx I = 0, Indx J = 0>
//            static constexpr
//            auto
//            setSubsidiaryUnknownSizes(
//                    auto &
//                    subsidiary_unknown_sizes_arg
//            )
//            {
//                if constexpr (I == 0) {
//                    using FU = typename FiniteElement::SubsidiaryUnknowns;
//                    if constexpr (!TypeTraits<FU>::template are_same<Void>) {
//                        subsidiary_unknown_sizes_arg.get(0) = FU::dim_unknowns;
//                    }
//                    else {
//                        subsidiary_unknown_sizes_arg.get(0) = 0;
//                    }
//                }
//                else {
//                    using FU = typename FiniteElementComponent<I - 1, J>::SubsidiaryUnknowns;
//                    if (J == 0) {
//                        subsidiary_unknown_sizes_arg.get(I) = 0;
//                    }
//                    if constexpr (!TypeTraits<FU>::template are_same<Void>) {
//                        auto num_faces = Element::Components::template Type<I - 1>::template Type<J>::size();
//                        subsidiary_unknown_sizes_arg.get(I) += num_faces * FU::dim_unknowns;
//                    }
//                    if constexpr (J < Element::Components::template Type<I - 1>::size() - 1) {
//                        setSubsidiaryUnknownSizes<I, J + 1>(subsidiary_unknown_sizes_arg);
//                    }
//                }
//                if constexpr (I < Element::Components::size()) {
//                    setSubsidiaryUnknownSizes<I + 1, 0>(subsidiary_unknown_sizes_arg);
//                }
//            }
//
//            static constexpr
//            auto
//            setStructuralUnknownSizes()
//            {
//                auto structural_unknown_sizes_arg = Array<Indx, E.shape_description.ord_shape + 1>();
//                setStructuralUnknownSizes(structural_unknown_sizes_arg);
//                return structural_unknown_sizes_arg;
//            }
//
//            static constexpr
//            auto
//            setSubsidiaryUnknownSizes()
//            {
//                auto subsidiary_unknown_sizes_arg = Array<Indx, E.shape_description.ord_shape + 1>();
//                setSubsidiaryUnknownSizes(subsidiary_unknown_sizes_arg);
//                return subsidiary_unknown_sizes_arg;
//            }
//
//            static constexpr
//            auto
//            setSize()
//            {
//                auto sum = Indx(0);
//                for (int i = 0; i < E.shape_description.ord_shape + 1; ++i) {
//                    sum += structural_unknown_sizes.get(i) + subsidiary_unknown_sizes.get(i);
//                }
//                return sum;
//            }
//
//        public:
//
//            auto const static constexpr structural_unknown_sizes = setStructuralUnknownSizes();
//
//            auto const static constexpr subsidiary_unknown_sizes = setSubsidiaryUnknownSizes();
//
//            auto const static constexpr sum = setSize();
//
//        };

#endif //FETA_FETA_CORE_ELEMENT_FE_DATA_HXX
