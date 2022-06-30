//
// Created by dsiedel on 14/06/22.
//

#ifndef LOLITA_LOLITA_CORE_5_002_BASIS_HXX
#define LOLITA_LOLITA_CORE_5_002_BASIS_HXX

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

namespace lolita::core2::finite_element
{

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElementBasis;

    namespace basis
    {

        struct Basis2
        {

        private:

            /**
             * @brief 
             * 
             */
            enum Type
            {

                Monomial,
                Lagrange,

            };

            /**
             * @brief 
             * 
             */
            constexpr
            Basis2(
                core_fem::basis::Basis2::Type type,
                lolita::integer ord
            )
            :
            type_(type),
            ord_(ord)
            {}
            
            /**
             * @brief 
             * 
             * @param other 
             * @return constexpr lolita::boolean 
             */
            constexpr
            lolita::boolean
            operator==(
                    Basis2 const & other
            )
            const = default;
            
            /**
             * @brief 
             * 
             * @param other 
             * @return constexpr lolita::boolean 
             */
            constexpr
            lolita::boolean
            operator!=(
                    Basis2 const & other
            )
            const = default;

        public:

            /**
             * @brief 
             * 
             * @param ord 
             * @return constexpr core_fem::basis::Basis2 
             */
            static constexpr
            core_fem::basis::Basis2
            monomial(
                lolita::integer ord
            )
            {
                return Basis2(core_fem::basis::Basis2::Type::Monomial, ord);
            }

            /**
             * @brief 
             * 
             * @return constexpr lolita::boolean 
             */
            constexpr
            lolita::boolean
            isMonomial()
            const
            {
                return type_ == core_fem::basis::Basis2::Type::Monomial;
            }

            /**
             * @brief 
             * 
             */
            core_fem::basis::Basis2::Type type_;

            /**
             * @brief 
             * 
             */
            lolita::integer ord_;

        };

        /**
         * @brief 
         * 
         * @tparam t_element 
         * @tparam t_basis 
         */
        template<lolita::core2::geometry::Element t_element, lolita::core2::finite_element::basis::Basis2 t_basis>
        struct FiniteElementBasisTraits2;

        template<lolita::core2::geometry::Element t_element, lolita::core2::finite_element::basis::Basis2 t_basis>
        requires(t_basis.isMonomial())
        struct FiniteElementBasisTraits2<t_element, t_basis>
        {

            /**
             * @brief The basis dimension or cardinality
             * 
             */
            lolita::index const static constexpr dim_ = numerics::binomial(t_element.dim_ + t_basis.ord_, t_element.dim_);
            
            /**
             * @brief 
             * 
             * @tparam t_domain 
             * @tparam t_finite_element 
             */
            template<lolita::domain::Domain t_domain, auto t_finite_element>
            struct Implementation : FiniteElementBasis<t_element, t_domain, t_finite_element>
            {

            private:
            
                /**
                 * @brief 
                 * 
                 * @return constexpr std::array<std::array<lolita::index, 3>, dim_> 
                 */
                static constexpr
                std::array<std::array<lolita::index, 3>, dim_>
                make_exponents()
                {
                    auto exponents = std::array<std::array<lolita::index, 3>, dim_>();
                    auto row = lolita::index(0);
                    if constexpr (t_element.dim_ == 0) {
                        exponents[row][0] = 0;
                        exponents[row][1] = 0;
                        exponents[row][2] = 0;
                    }
                    else if constexpr (t_element.dim_ == 1) {
                        for (lolita::index i = 0; i < t_basis.ord_ + 1; ++i) {
                            exponents[row][0] = i;
                            exponents[row][1] = 0;
                            exponents[row][2] = 0;
                            row += 1;
                        }
                    }
                    else if constexpr (t_element.dim_ == 2) {
                        for (lolita::index i = 0; i < t_basis.ord_ + 1; ++i) {
                            for (lolita::index j = 0; j < i + 1; ++j) {
                                exponents[row][0] = i - j;
                                exponents[row][1] = j;
                                exponents[row][2] = 0;
                                row += 1;
                            }
                        }
                    }
                    else if constexpr (t_element.dim_ == 3) {
                        for (lolita::index i = 0; i < t_basis.ord_ + 1; ++i) {
                            for (lolita::index j = 0; j < i + 1; ++j) {
                                for (lolita::index k = 0; k < i + 1; ++k) {
                                    if (j + k < i + 1) {
                                        exponents[row][0] = i - (j + k);
                                        exponents[row][1] = k;
                                        exponents[row][2] = j;
                                        row += 1;
                                    }
                                }
                            }
                        }
                    }
                    return exponents;
                }
                
                /**
                 * @brief 
                 * 
                 */
                std::array<std::array<lolita::index, 3>, dim_> static constexpr exponents_ = make_exponents();

            public:

                /**
                 * @brief evaluate the polynomial basis function at some point
                 * 
                 * @param point the point where to evaluate the polynomial basis function
                 * @return a scalar-valued vector of size the dimension of the polynomial basis
                 */
                lolita::matrix::Vector<lolita::real, dim_>
                getBasisEvaluation(
                        lolita::domain::Point const & point
                )
                const
                {
                    auto basis_vector_values = lolita::matrix::Vector<lolita::real, dim_>();
                    auto const centroid = this->getReferenceCentroid();
                    auto const diameters = this->getCurrentDiameters();
                    for (lolita::index i = 0; i < dim_; ++i) {
                        auto value = lolita::real(1);
                        for (lolita::index j = 0; j < t_element.dim_; ++j) {
                            auto dist = this->getRiemannianDistance(centroid, point, j);
                            value *= std::pow(2.0 * dist / diameters(j), exponents_[i][j]);
                        }
                        basis_vector_values(i) = value;
                    }
                    return basis_vector_values;
                }

                /**
                 * @brief evaluate the polynomial basis derivative at some point
                 * 
                 * @param point the point where to evaluate the polynomial basis function
                 * @param derivative_direction the derivative direction, that must be less than or equal to the element dimension
                 * @return a scalar-valued vector of size the dimension of the polynomial basis
                 */
                lolita::matrix::Vector<lolita::real, dim_>
                getBasisDerivative(
                        lolita::domain::Point const & point,
                        lolita::index derivative_direction
                )
                const
                {
                    auto basis_vector_values = lolita::matrix::Vector<lolita::real, dim_>();
                    auto const centroid = this->getReferenceCentroid();
                    auto const diameters = this->getCurrentDiameters();
                    for (lolita::index i = 0; i < dim_; ++i) {
                        auto value = lolita::real(1);
                        for (lolita::index j = 0; j < t_element.dim_; ++j) {
                            if (j != derivative_direction) {
                                auto dist = this->getRiemannianDistance(centroid, point, j);
                                value *= std::pow(2.0 * (dist) / diameters(j), exponents_[i][j]);
                            }
                            else {
                                if (exponents_[i][j] > 0) {
                                    auto c = 2.0 * exponents_[i][j] / diameters(j);
                                    auto dist = this->getRiemannianDistance(centroid, point, j);
                                    value *= c * std::pow(2.0 * (dist) / diameters(j), exponents_[i][j] - 1);
                                }
                                else {
                                    value *= 0.0;
                                }
                            }
                        }
                        basis_vector_values(i) = value;
                    }
                    return basis_vector_values;
                }

            };

        };

        /**
         * @brief
         */
        struct Basis : public lolita::utility::Enumeration<Basis>
        {

            /**
             * @brief
             * @param tag
             */
            constexpr
            Basis(
                    std::basic_string_view<lolita::character> && tag
            )
            :
            lolita::utility::Enumeration<Basis>(std::forward<std::basic_string_view<lolita::character>>(tag))
            {}

            /**
             * @brief
             * @return
             */
            static constexpr
            Basis
            Monomial()
            {
                return Basis("Monomial");
            }

            /**
             * @brief
             * @return
             */
            constexpr
            lolita::boolean
            isMonomial()
            const
            {
                return * this == Monomial();
            }

            /**
             * @brief
             * @return
             */
            static constexpr
            Basis
            Lagrange()
            {
                return Basis("Lagrange");
            }

            /**
             * @brief
             * @return
             */
            constexpr
            lolita::boolean
            isLagrange()
            const
            {
                return * this == Lagrange();
            }

        };

        /**
         * @brief Implementation object for the definition of a polynomial basis function
         * @tparam t_element the element on which the monomial is defined
         * @tparam t_basis the basis enum
         * @tparam t_ord the polynomial order of the basis
         */
        template<lolita::core2::geometry::Element t_element, lolita::core2::finite_element::basis::Basis t_basis, lolita::index t_ord>
        struct FiniteElementBasisTraits;

        /**
         * @brief Implementation object for the definition of a monomial basis function
         * @tparam t_element the element on which the monomial is defined
         * @tparam t_basis the monomial basis enum
         * @tparam t_ord the polynomial order of the basis
         */
        template<lolita::core2::geometry::Element t_element, lolita::core2::finite_element::basis::Basis t_basis, lolita::index t_ord>
        requires(t_basis.isMonomial())
        struct FiniteElementBasisTraits<t_element, t_basis, t_ord>
        {

            /**
             * @brief
             * @return
             */
            static constexpr
            lolita::integer
            dim()
            {
                return numerics::binomial(t_element.dim_ + t_ord, t_element.dim_);
            }

            /**
             * @brief The basis dimension or cardinality
             */
            lolita::index const static constexpr dim_ = numerics::binomial(t_element.dim_ + t_ord, t_element.dim_);

            /**
             * @brief Implementation of the evaluation and derivative function for the given polynomial basis
             * @tparam t_domain
             * @tparam t_finite_element
             */
            template<lolita::domain::Domain t_domain, auto t_finite_element>
            struct Implementation : FiniteElementBasis<t_element, t_domain, t_finite_element>
            {

            private:

                /**
                 * @brief
                 * @return
                 */
                static constexpr
                std::array<std::array<lolita::index, 3>, dim_>
                make_exponents()
                {
                    auto exponents = std::array<std::array<lolita::index, 3>, dim_>();
                    auto row = lolita::index(0);
                    if constexpr (t_element.dim_ == 0) {
                        exponents[row][0] = 0;
                        exponents[row][1] = 0;
                        exponents[row][2] = 0;
                    }
                    else if constexpr (t_element.dim_ == 1) {
                        for (lolita::index i = 0; i < t_ord + 1; ++i) {
                            exponents[row][0] = i;
                            exponents[row][1] = 0;
                            exponents[row][2] = 0;
                            row += 1;
                        }
                    }
                    else if constexpr (t_element.dim_ == 2) {
                        for (lolita::index i = 0; i < t_ord + 1; ++i) {
                            for (lolita::index j = 0; j < i + 1; ++j) {
                                exponents[row][0] = i - j;
                                exponents[row][1] = j;
                                exponents[row][2] = 0;
                                row += 1;
                            }
                        }
                    }
                    else if constexpr (t_element.dim_ == 3) {
                        for (lolita::index i = 0; i < t_ord + 1; ++i) {
                            for (lolita::index j = 0; j < i + 1; ++j) {
                                for (lolita::index k = 0; k < i + 1; ++k) {
                                    if (j + k < i + 1) {
                                        exponents[row][0] = i - (j + k);
                                        exponents[row][1] = k;
                                        exponents[row][2] = j;
                                        row += 1;
                                    }
                                }
                            }
                        }
                    }
                    return exponents;
                }

                /**
                 * @brief
                 */
                std::array<std::array<lolita::index, 3>, dim_> const static constexpr exponents_ = make_exponents();

            public:

                /**
                 * @brief evaluate the polynomial basis function at some point
                 * @param point the point where to evaluate the polynomial basis function
                 * @return a scalar-valued vector of size the dimension of the polynomial basis
                 */
                lolita::matrix::Vector<lolita::real, dim_>
                getBasisEvaluation(
                        lolita::domain::Point const & point
                )
                const
                {
                    auto basis_vector_values = lolita::matrix::Vector<lolita::real, dim_>();
                    auto const centroid = this->getReferenceCentroid();
                    auto const diameters = this->getCurrentDiameters();
                    for (lolita::index i = 0; i < dim_; ++i) {
                        auto value = lolita::real(1);
                        for (lolita::index j = 0; j < t_element.dim_; ++j) {
                            auto dist = this->getRiemannianDistance(centroid, point, j);
                            //value *= numerics::pow(2.0 * dist / diameters(j), exponents.get(i, j));
                            value *= std::pow(2.0 * dist / diameters(j), exponents_[i][j]);
                        }
                        basis_vector_values(i) = value;
                    }
                    return basis_vector_values;
                }

                /**
                 * @brief evaluate the polynomial basis derivative at some point
                 * @param point the point where to evaluate the polynomial basis function
                 * @param derivative_direction the derivative direction, that must be less than or equal to the element dimension
                 * @return a scalar-valued vector of size the dimension of the polynomial basis
                 */
                lolita::matrix::Vector<lolita::real, dim_>
                getBasisDerivative(
                        lolita::domain::Point const & point,
                        lolita::index derivative_direction
                )
                const
                {
                    auto basis_vector_values = lolita::matrix::Vector<lolita::real, dim_>();
                    auto const centroid = this->getReferenceCentroid();
                    auto const diameters = this->getCurrentDiameters();
                    for (lolita::index i = 0; i < dim_; ++i) {
                        auto value = lolita::real(1);
                        for (lolita::index j = 0; j < t_element.dim_; ++j) {
                            if (j != derivative_direction) {
                                auto dist = this->getRiemannianDistance(centroid, point, j);
                                value *= std::pow(2.0 * (dist) / diameters(j), exponents_[i][j]);
                            }
                            else {
                                if (exponents_[i][j] > 0) {
                                    auto c = 2.0 * exponents_[i][j] / diameters(j);
                                    auto dist = this->getRiemannianDistance(centroid, point, j);
                                    value *= c * std::pow(2.0 * (dist) / diameters(j), exponents_[i][j] - 1);
//                                value *= c * numerics::pow(2.0 * (dist) / diameters(j), exponents.get(i, j) - 1);
                                }
                                else {
                                    value *= 0.0;
                                }
                            }
                        }
                        basis_vector_values(i) = value;
                    }
                    return basis_vector_values;
                }

            };

        };

    }

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElementBasis : virtual FiniteElementBase<t_element, t_domain, t_finite_element>
    {

        /**
         * @brief
         * @tparam t_basis
         * @tparam t_ord
         * @param point
         * @return
         */
        template<lolita::core2::finite_element::basis::Basis t_basis, lolita::index t_ord>
        lolita::matrix::Vector<lolita::real, lolita::core2::finite_element::basis::FiniteElementBasisTraits<t_element, t_basis, t_ord>::dim_>
        getBasisEvaluation(
                lolita::domain::Point const & point
        )
        const
        {
            using t_Implementation = typename basis::FiniteElementBasisTraits<t_element, t_basis, t_ord>::template Implementation<t_domain, t_finite_element>;
            return static_cast<t_Implementation const *>(this)->getBasisEvaluation(point);
        }

        /**
         * @brief
         * @tparam t_basis
         * @tparam t_ord
         * @param point
         * @param derivative_direction
         * @return
         */
        template<lolita::core2::finite_element::basis::Basis t_basis, lolita::index t_ord>
        lolita::matrix::Vector<lolita::real, lolita::core2::finite_element::basis::FiniteElementBasisTraits<t_element, t_basis, t_ord>::dim_>
        getBasisDerivative(
                lolita::domain::Point const & point,
                lolita::index derivative_direction
        )
        const
        {
            using t_Implementation = typename basis::FiniteElementBasisTraits<t_element, t_basis, t_ord>::template Implementation<t_domain, t_finite_element>;
            return static_cast<t_Implementation const *>(this)->getBasisDerivative(point, derivative_direction);
        }

        /**
         * @brief
         * @tparam t_basis
         * @tparam t_ord
         * @param point
         * @return
         */
        template<lolita::core2::finite_element::basis::Basis2 t_basis>
        lolita::matrix::Vector<lolita::real, lolita::core2::finite_element::basis::FiniteElementBasisTraits2<t_element, t_basis>::dim_>
        getBasisEvaluation(
                lolita::domain::Point const & point
        )
        const
        {
            using t_Implementation = typename core_fem::basis::FiniteElementBasisTraits2<t_element, t_basis>::template Implementation<t_domain, t_finite_element>;
            return static_cast<t_Implementation const *>(this)->getBasisEvaluation(point);
        }

        /**
         * @brief
         * @tparam t_basis
         * @tparam t_ord
         * @param point
         * @param derivative_direction
         * @return
         */
        template<lolita::core2::finite_element::basis::Basis2 t_basis>
        lolita::matrix::Vector<lolita::real, lolita::core2::finite_element::basis::FiniteElementBasisTraits2<t_element, t_basis>::dim_>
        getBasisDerivative(
                lolita::domain::Point const & point,
                lolita::index derivative_direction
        )
        const
        {
            using t_Implementation = typename core_fem::basis::FiniteElementBasisTraits2<t_element, t_basis>::template Implementation<t_domain, t_finite_element>;
            return static_cast<t_Implementation const *>(this)->getBasisDerivative(point, derivative_direction);
        }

    };

}

#endif //LOLITA_LOLITA_CORE_5_002_BASIS_HXX
