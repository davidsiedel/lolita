//
// Created by dsiedel on 04/06/22.
//

#ifndef LOLITA_LOLITA_CORE_3_HXX
#define LOLITA_LOLITA_CORE_3_HXX

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"
#include "lolita/lolita_core_1.hxx"
#include "lolita/lolita_core_2.hxx"

namespace lolita::core::geometry
{

    /**
     * @brief
     * @tparam t_element
     * @tparam t_quadrature
     * @tparam t_ord
     */
    template<lolita::core::geometry::Element t_element, lolita::finite_element::Quadrature t_quadrature, lolita::index t_ord>
    struct ElementQuadratureRuleTraits;

    /**
     * @brief
     * @tparam t_element
     * @tparam t_quadrature
     * @tparam t_ord
     */
    template<lolita::core::geometry::Element t_element, lolita::finite_element::Quadrature t_quadrature, lolita::index t_ord>
    requires(t_element.isPoint() || !t_element.isPoint())
    struct ElementQuadratureRuleTraits<t_element, t_quadrature, t_ord>
    {

        /**
         * @brief
         */
        lolita::index const static constexpr dim_ = 1;

        /**
         * @brief
         */
        std::array<std::array<lolita::real, 3>, dim_> const static constexpr reference_points_ = {
                +0.0000000000000000, +0.0000000000000000, +0.0000000000000000
        };

        /**
         * @brief
         */
        std::array<lolita::real, dim_> const static constexpr reference_weights_ = {
                +1.0000000000000000
        };

    };

    /**
     * @brief
     * @tparam t_quadrature
     * @tparam t_ord
     */
    template<lolita::core::geometry::Element t_element, lolita::finite_element::Quadrature t_quadrature, lolita::index t_ord>
    struct ElementQuadratureTraits : public ElementQuadratureRuleTraits<t_element, t_quadrature, t_ord>
    {

        /**
         * @brief
         */
        lolita::finite_element::Quadrature const static constexpr quadrature_ = t_quadrature;

        /**
         * @brief
         */
        lolita::index const static constexpr ord_ = t_ord;

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::finite_element::Quadrature
        getQuadrature()
        {
            return t_quadrature;
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::integer
        getOrd()
        {
            return t_ord;
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::integer
        getDim()
        {
            return ElementQuadratureRuleTraits<t_element, t_quadrature, t_ord>::dim_;
        }

    };

}

#endif //LOLITA_LOLITA_CORE_3_HXX
