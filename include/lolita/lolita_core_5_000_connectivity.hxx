//
// Created by dsiedel on 14/06/22.
//

#ifndef LOLITA_LOLITA_CORE_5_000_CONNECTIVITY_HXX
#define LOLITA_LOLITA_CORE_5_000_CONNECTIVITY_HXX

#include <execution>

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"
#include "lolita/lolita_core_1.hxx"
#include "lolita/lolita_core_2.hxx"
#include "lolita/lolita_core_4.hxx"

namespace lolita::core::finite_element
{

    namespace core_fem = lolita::core::finite_element;

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElementTraits;

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core::geometry::Element t_element, lolita::domain::Domain t_domain, auto t_finite_element>
    struct FiniteElement;

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core::geometry::Element t_element, lolita::domain::Domain t_domain, auto t_finite_element>
    struct FiniteElementConnectivity
    {

    private:

        /**
         * @brief
         */
        using t_ElementTraits = lolita::core::geometry::ElementTraits<t_element, t_domain>;

        /**
         * @brief
         */
        template<lolita::core::geometry::Element t__element, lolita::domain::Domain t__domain, auto... t__finite_element>
        using t_ElementPointer = std::shared_ptr<FiniteElement<t__element, t__domain, t__finite_element...>>;

    public:

        /**
         * @brief
         */
        using Components = typename t_ElementTraits::template InnerConnectivity<t_ElementPointer, t_finite_element>;

        /**
         * @brief
         */
        using Neighbours = typename t_ElementTraits::template OuterConnectivity<t_ElementPointer, t_finite_element>;

        /**
         * @brief Basic equality comparison operator
         * @param other
         * @return
         */
        lolita::boolean
        operator==(
                FiniteElementConnectivity const & other
        )
        const = default;

        /**
         * @brief Basic equality comparison operator
         * @param other
         * @return
         */
        lolita::boolean
        operator!=(
                FiniteElementConnectivity const & other
        )
        const = default;

        /**
         * @brief
         * @return
         */
        std::basic_string<lolita::character>
        hash()
        const
        requires(t_element.isPoint())
        {
            return std::to_string(this->tag_);
        }

        /**
         * @brief
         * @return
         */
        std::basic_string<lolita::character>
        hash()
        const
        {
            std::basic_stringstream<lolita::character> hash;
            auto const & nodes = getComponents<t_element.dim_ - 1, 0>();
            for (auto const & node : nodes) {
                hash << node->hash();
            }
            return hash.str();
        }

        /**
         * @brief
         * @tparam t_i
         * @tparam t_j
         * @param i
         * @param j
         * @return
         */
        template<lolita::index t_i, lolita::index t_j>
        static constexpr
        lolita::index
        getComponentNodeConnection(
                lolita::index i,
                lolita::index j
        )
        requires(!t_element.isPoint())
        {
            return std::get<t_j>(std::get<t_i>(lolita::core::geometry::ElementTraits<t_element, t_domain>::node_connectivity_))[i][j];
        }

        /**
         * @brief
         * @tparam t_i
         * @tparam t_j
         * @return
         */
        template<lolita::index t_i, lolita::index t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, Neighbours>> &
        getNeighbours()
        {
            return std::get<t_j>(std::get<t_i>(neighbours_));
        }

        /**
         * @brief
         * @tparam t_i
         * @tparam t_j
         * @return
         */
        template<lolita::index t_i, lolita::index t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, Neighbours>> const &
        getNeighbours()
        const
        {
            return std::get<t_j>(std::get<t_i>(neighbours_));
        }

        /**
         * @brief
         * @tparam t_i
         * @tparam t_j
         * @return
         */
        template<lolita::index t_i, lolita::index t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, Components>> &
        getComponents()
        requires(!t_element.isPoint())
        {
            return std::get<t_j>(std::get<t_i>(components_));
        }

        /**
         * @brief
         * @tparam t_i
         * @tparam t_j
         * @return
         */
        template<lolita::index t_i, lolita::index t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, Components>> const &
        getComponents()
        const
        requires(!t_element.isPoint())
        {
            return std::get<t_j>(std::get<t_i>(components_));
        }

        /**
         * @brief
         * @tparam t_i
         * @tparam t_j
         * @param i
         * @return
         */
        template<lolita::index t_i, lolita::index t_j>
        lolita::index
        getComponentIndex(
                lolita::index i
        )
        const
        requires(!t_element.isPoint())
        {
            auto constexpr t_component = t_ElementTraits::template getComponent<t_i, t_j>();
            using t_NeighbourTraits = lolita::core::geometry::ElementTraits<t_component, t_domain>;
            auto constexpr t_coordinates = t_NeighbourTraits::template getNeighbourCoordinates<t_element>();
            auto const & items = getComponents<t_i, t_j>()[i]->template getNeighbours<t_coordinates.dim_, t_coordinates.tag_>();
            auto is_equal = [&] (std::shared_ptr<FiniteElementConnectivity> const & ptr_element) {
                return * ptr_element == * this;
            };
            return std::distance(items.begin(), std::find_if(items.begin(), items.end(), is_equal));
        }

        /**
         * @brief
         * @tparam t_i
         * @tparam t_j
         * @param i
         * @return
         */
        template<lolita::index t_i, lolita::index t_j>
        lolita::integer
        getComponentOrientation(
                lolita::index i
        )
        const
        requires(!t_element.isPoint())
        {
            return getComponentIndex<t_i, t_j>(i) == 0 ? 1 : -1;
        }

        /**
         * @brief
         */
        lolita::natural tag_;

        /**
         * @brief
         */
        Neighbours neighbours_;

        /**
         * @brief
         */
        Components components_;

        /**
         * @brief
         */
        std::vector<std::shared_ptr<std::basic_string<lolita::character>>> domains_;

        /**
         * @brief
         */
        std::shared_ptr<lolita::domain::Point> coordinates_;

    };

}

#endif //LOLITA_LOLITA_CORE_5_000_CONNECTIVITY_HXX
