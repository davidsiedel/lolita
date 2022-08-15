//
// Created by dsiedel on 23/03/2022.
//

#ifndef FETA_FETA_CORE_ELEMENT_ELEMENT_CONNECTIVITY_HXX
#define FETA_FETA_CORE_ELEMENT_ELEMENT_CONNECTIVITY_HXX

//#include "new/00_finite_element_description_hho.hxx"
//#include "new/00_frame_type.hxx"
#include "new/feta_core_model_description.hxx"
//#include "new/feta_core_mesh_interior_domain.hxx"
#include "new/feta_core_element_element_reference.hxx"
#include "new/feta_core_fe_description.hxx"
#include "new/_feta_shared_pointer.hxx"
#include "new/_feta_matrix.hxx"

namespace feta::core::element
{

    template<ElementDescription E, auto M>
    struct ElementConnectivity : public ElementReference<E>
    {

        using ElementReference = element::ElementReference<E>;

        auto const static constexpr element_description = E;

        auto const static constexpr dim_euclidean = M.dim_euclidean;

        auto const static constexpr ord_shape = element_description.shape_description.ord_shape;

        auto const static constexpr dim_shape = element_description.shape_description.dim_shape;

    private:

        using Self = ElementConnectivity;

        template<ElementDescription C>
        using ElementPointer = SharedPointer<Element<C, M>>;

    public:

        template<ElementDescription C>
        struct Component : public ElementPointer<C>
        {

            using ElementPointer = Self::ElementPointer<C>;

            Component()
            :
            ElementPointer(),
            orientation(1)
            {}

            Component(
                    ElementPointer const &
                    base_arg,
                    Intg
                    orientation_arg
            )
            :
            ElementPointer(base_arg),
            orientation(orientation_arg)
            {}

            Component(
                    typename ElementPointer::Type const &
                    element_arg,
                    Intg
                    orientation_arg
            )
            :
            ElementPointer(element_arg),
            orientation(orientation_arg)
            {}

            Bool
            operator==(
                    Component const &
                    other
            )
            const
            {
                auto eq_0 = ElementPointer::operator==(other);
                auto eq_1 = orientation == other.orientation;
                return eq_0 && eq_1;
            }

            Bool
            operator!=(
                    ElementDescription const &
                    other
            )
            const
            {
                return !(other == * this);
            }

            Intg orientation;

        };

        using Components = typename ElementReference::template Components<Component>;

        using Neighbours = ElementNeighbourArray<E, dim_euclidean, ElementPointer>;

        using DomainPointer = SharedPointer<core::MeshInteriorDomain<M>>;

        ElementConnectivity()
        :
        tag(),
        components(Components()),
        neighbours(Neighbours()),
        domain(DomainPointer())
        {}

        Bool
        operator==(
                ElementConnectivity const &
                other
        )
        const
        {
            auto eq_0 = tag == other.tag;
            return eq_0;
        }

        Bool
        operator!=(
                ElementConnectivity const &
                other
        )
        const
        {
            return !(other == * this);
        }

        auto
        getCurrentCoordinates()
        const
        {
            auto current_nodes_coordinates = Matrix<Real, Self::dim_euclidean, Self::num_nodes>();
            auto const & nodes = components.template get<ord_shape - 1>().template get<0>();
            for (Indx i = 0; i < Self::num_nodes; ++i) {
                current_nodes_coordinates.col(i) = nodes.get(i).get().getCurrentCoordinates();
            }
            return current_nodes_coordinates;
        }

        static
        auto
        getReferenceCoordinates()
        {
            using ReferenceCoordinates = Matrix<Matrix<Real, dim_shape, Self::num_nodes, matrix::col_major> const>;
            return ReferenceCoordinates(Self::reference_nodes.data.data());
        }

        auto
        getHash()
        const
        {
            StrgStream ss;
            auto const & nodes = components.template get<ord_shape - 1>().template get<0>();
            for (Indx i = 0; i < Self::num_nodes; ++i) {
                ss << nodes(i).get().getHash();
            }
            return ss.str();
        }

        Indx tag;

        Components components;

        Neighbours neighbours;

        DomainPointer domain;

    };

    template<auto M>
    struct ElementConnectivity<pnt_0, M> : public ElementReference<pnt_0>
    {

        using ElementReference = element::ElementReference<pnt_0>;

        auto const static constexpr element_description = pnt_0;

        auto const static constexpr dim_euclidean = M.dim_euclidean;

        auto const static constexpr ord_shape = element_description.shape_description.ord_shape;

        auto const static constexpr dim_shape = element_description.shape_description.dim_shape;

    private:

        using Self = ElementConnectivity;

        template<ElementDescription N>
        using ElementPointer = SharedPointer<Element<N, M>>;

    public:

        using Coordinates = Matrix<Real, dim_euclidean>;

        using Neighbours = ElementNeighbourArray<element_description, dim_euclidean, ElementPointer>;

        using DomainPointer = SharedPointer<core::MeshInteriorDomain<M>>;

        ElementConnectivity()
        :
        tag(),
        coordinates(Coordinates::Zero()),
        neighbours(Neighbours()),
        domain(DomainPointer())
        {}

        Bool
        operator==(
                ElementConnectivity const &
                other
        )
        const
        {
            auto eq_0 = tag == other.tag;
            return eq_0;
        }

        Bool
        operator!=(
                ElementConnectivity const &
                other
        )
        const
        {
            return !(other == * this);
        }

        auto &
        getCurrentCoordinates()
        {
            return coordinates;
        }

        auto const &
        getCurrentCoordinates()
        const
        {
            return coordinates;
        }

        static
        auto
        getReferenceCoordinates()
        {
            using ReferenceCoordinates = Matrix<Matrix<Real, dim_euclidean, Self::num_nodes, matrix::col_major> const>;
            return ReferenceCoordinates(Self::reference_nodes.data.data());
        }

        auto
        getHash()
        const
        {
            StrgStream ss;
            ss << tag + 1;
            return ss.str();
        }

        Indx tag;

        Coordinates coordinates;

        Neighbours neighbours;

        DomainPointer domain;

    };

}

#endif //FETA_FETA_CORE_ELEMENT_ELEMENT_CONNECTIVITY_HXX
