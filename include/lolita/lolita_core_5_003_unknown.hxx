//
// Created by dsiedel on 14/06/22.
//

#ifndef LOLITA_LOLITA_CORE_5_003_UNKNOWN_HXX
#define LOLITA_LOLITA_CORE_5_003_UNKNOWN_HXX

#include <execution>

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"
#include "lolita/lolita_core_1.hxx"
#include "lolita/lolita_core_2.hxx"
#include "lolita/lolita_core_4.hxx"
#include "lolita/lolita_core_5_000_connectivity.hxx"
#include "lolita/lolita_core_5_001_base.hxx"
#include "lolita/lolita_core_5_002_basis.hxx"

namespace lolita::core::finite_element
{

    namespace lolita_fem = lolita::finite_element;
    namespace lolita_dom = lolita::domain;
    namespace lolita_fld = lolita::field;
    namespace core_fld = lolita::core::field;
    namespace core_geo = lolita::core::geometry;
    namespace core_fem = lolita::core::finite_element;

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<core_geo::Element t_element, lolita_dom::Domain t_domain, lolita_fem::FiniteElementConcept auto t_finite_element>
    struct FiniteElementFieldUnknowns;

    namespace unknown
    {
        
        /**
         * @brief 
         * 
         * @tparam t_dim 
         */
        template<lolita::integer t_dim>
        struct ScalarUnknown
        {

            lolita::integer unknown_index_;

            lolita::integer binding_index_;

            std::shared_ptr<lolita::core::system::SparseSystem> system_;

            /**
             * @brief The binding coefficient vector type
             */
            using CoordinateVector = lolita::matrix::Vector<lolita::integer, t_dim>;

            /**
             * @brief The binding coefficient vector type
             */
            using CoefficientVector = lolita::matrix::Vector<lolita::real, t_dim>;

            /**
             * @brief
             * @return
             */
            lolita::boolean
            isBound()
            const
            {
                return !(bindings_coefficients_ == nullptr);
            }

            /**
             * @brief
             * @return
             */
            lolita::boolean
            isStructural()
            const
            {
                return !(unknowns_coordinates_ == nullptr);
            }

            /**
             * @brief
             * @return
             */
            CoefficientVector &
            getUnknownsCoefficients()
            {
                return * unknowns_coefficients_;
            }

            /**
             * @brief
             * @return
             */
            CoefficientVector const &
            getUnknownsCoefficients()
            const
            {
                return * unknowns_coefficients_;
            }

            /**
             * @brief
             * @return
             */
            CoefficientVector &
            getBindingsCoefficients()
            {
                return * bindings_coefficients_;
            }

            /**
             * @brief
             * @return
             */
            CoefficientVector const &
            getBindingsCoefficients()
            const
            {
                return * bindings_coefficients_;
            }

            /**
             * @brief
             * @return
             */
            CoordinateVector &
            getUnknownsCoordinates()
            {
                return * unknowns_coordinates_;
            }

            /**
             * @brief
             * @return
             */
            CoordinateVector const &
            getUnknownsCoordinates()
            const
            {
                return * unknowns_coordinates_;
            }

            /**
             * @brief
             * @return
             */
            CoordinateVector &
            getBindingsCoordinates()
            {
                return * bindings_coordinates_;
            }

            /**
             * @brief
             * @return
             */
            CoordinateVector const &
            getBindingsCoordinates()
            const
            {
                return * bindings_coordinates_;
            }

            /**
             * @brief
             */
            std::unique_ptr<CoefficientVector> unknowns_coefficients_;

            /**
             * @brief
             */
            std::unique_ptr<CoefficientVector> bindings_coefficients_;

            /**
             * @brief
             */
            std::unique_ptr<CoordinateVector> unknowns_coordinates_;

            /**
             * @brief
             */
            std::unique_ptr<CoordinateVector> bindings_coordinates_;

        };

    };

}

#endif //LOLITA_LOLITA_CORE_5_003_UNKNOWN_HXX
