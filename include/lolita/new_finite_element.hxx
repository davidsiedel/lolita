//
// Created by dsiedel on 28/05/22.
//

#ifndef LOLITA_NEW_FINITE_ELEMENT_HXX
#define LOLITA_NEW_FINITE_ELEMENT_HXX

#include <execution>

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"
#include "lolita/new.hxx"
#include "lolita/new_element.hxx"

namespace lolita::core::finite_element
{

    namespace basis
    {

        /**
         * @brief Owning object for the definition of a polynomial basis function
         * @tparam _element the element on which the monomial is defined
         * @tparam _basis the basis enum
         * @tparam _ord the polynomial order of the basis
         */
        template<lolita::core::Element _element, lolita::finite_element::Basis _basis, lolita::index _ord>
        struct FiniteElementBasisBase
        {

            /**
             * @brief The basis object
             */
            lolita::finite_element::Basis const static constexpr basis_ = _basis;

            /**
             * @brief The polynomial order
             */
            lolita::index const static constexpr ord_ = _ord;

        };

        /**
         * @brief Implementation object for the definition of a polynomial basis function
         * @tparam _element the element on which the monomial is defined
         * @tparam _basis the basis enum
         * @tparam _ord the polynomial order of the basis
         */
        template<lolita::core::Element _element, lolita::finite_element::Basis _basis, lolita::index _ord>
        struct FiniteElementBasis;

        /**
         * @brief Implementation object for the definition of a monomial basis function
         * @tparam _element the element on which the monomial is defined
         * @tparam _basis the monomial basis enum
         * @tparam _ord the polynomial order of the basis
         */
        template<lolita::core::Element _element, lolita::finite_element::Basis _basis, lolita::index _ord>
        requires(_basis == lolita::finite_element::Basis::Monomial)
        struct FiniteElementBasis<_element, _basis, _ord> : public lolita::core::finite_element::basis::FiniteElementBasisBase<_element, _basis, _ord>
        {

        private:

            /**
             * @brief
             * @return
             */
            static constexpr
            lolita::integer
            _dim()
            {
                return numerics::binomial(_element.dim_ + _ord, _element.dim_);
            }

        public:

            /**
             * @brief The basis dimension or cardinality
             */
            lolita::index const static constexpr dim_ = _dim();

            /**
             * @brief Implementation of the evaluation and derivative function for the given polynomial basis
             * @tparam _T
             * @tparam _domain
             * @tparam _arg
             */
            template<
                    template<lolita::core::Element, lolita::geometry::Domain, auto...> typename _T,
                    lolita::geometry::Domain _domain,
                    lolita::finite_element::FiniteElementConcept auto... _arg
            >
            struct Implementation : public lolita::core::finite_element::FiniteElementGeometry<_T, _element, _domain, _arg...>
            {

            private:

                /**
                 * @brief
                 * @return
                 */
                static constexpr
                std::array<std::array<lolita::index, 3>, dim_>
                _exponents()
                {
                    auto exponents = std::array<std::array<lolita::index, 3>, dim_>();
                    auto row = lolita::index(0);
                    if constexpr (_element.dim_ == 0) {
                        exponents[row][0] = 0;
                        exponents[row][1] = 0;
                        exponents[row][2] = 0;
                    }
                    else if constexpr (_element.dim_ == 1) {
                        for (lolita::index i = 0; i < _ord + 1; ++i) {
                            exponents[row][0] = i;
                            exponents[row][1] = 0;
                            exponents[row][2] = 0;
                            row += 1;
                        }
                    }
                    else if constexpr (_element.dim_ == 2) {
                        for (lolita::index i = 0; i < _ord + 1; ++i) {
                            for (lolita::index j = 0; j < i + 1; ++j) {
                                exponents[row][0] = i - j;
                                exponents[row][1] = j;
                                exponents[row][2] = 0;
                                row += 1;
                            }
                        }
                    }
                    else if constexpr (_element.dim_ == 3) {
                        for (lolita::index i = 0; i < _ord + 1; ++i) {
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
                std::array<std::array<lolita::index, 3>, dim_> const static constexpr exponents_ = _exponents();

            public:

                /**
                 * @brief evaluate the polynomial basis function at some point
                 * @param point the point where to evaluate the polynomial basis function
                 * @return a scalar-valued vector of size the dimension of the polynomial basis
                 */
                lolita::matrix::Vector<lolita::real, dim_>
                evaluate(
                        lolita::geometry::Point const & point
                )
                const
                {
                    auto basis_vector_values = lolita::matrix::Vector<lolita::real, dim_>();
                    auto const centroid = this->getReferenceCentroid();
                    auto const diameters = this->getCurrentDiameters();
                    for (lolita::index i = 0; i < dim_; ++i) {
                        auto value = lolita::real(1);
                        for (lolita::index j = 0; j < _element.dim_; ++j) {
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
                evaluate(
                        lolita::geometry::Point const & point,
                        lolita::index derivative_direction
                )
                const
                {
                    auto basis_vector_values = lolita::matrix::Vector<lolita::real, dim_>();
                    auto const centroid = this->getReferenceCentroid();
                    auto const diameters = this->getCurrentDiameters();
                    for (lolita::index i = 0; i < dim_; ++i) {
                        auto value = lolita::real(1);
                        for (lolita::index j = 0; j < _element.dim_; ++j) {
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

    namespace unknown
    {

        /**
         * @brief
         * @tparam _domain
         * @tparam _field
         * @tparam _dim
         */
        template<
                lolita::geometry::Domain _domain,
                lolita::field::Tensor _field,
                lolita::integer _dim,
                lolita::core::finite_element::unknown::UnknownType _unknown_type
        >
        struct Unknown
        {

            /**
             * @brief The unknown type object
             */
            lolita::core::finite_element::unknown::UnknownType const static constexpr unknown_type_ = _unknown_type;

            /**
             * @brief The field object
             */
            lolita::field::Tensor const static constexpr field_ = _field;

            /**
             * @brief The dimension or cardinality of the unknown
             */
            lolita::integer const static constexpr dim_ = _dim;

        private:

            /**
             * @brief The binding coefficient vector type
             */
            using _CoordinatesVector = lolita::matrix::Vector<lolita::integer, _dim>;

            /**
             * @brief The binding coefficient vector type
             */
            using _CoefficientVector = lolita::matrix::Vector<lolita::real, _dim>;

            /**
             * @brief The field type
             */
            using _Field = lolita::core::field::TensorPolicy<_field, _domain.dim_>;

        public:

            template<lolita::core::finite_element::unknown::UnknownType>
            struct DegreeOfFreedom;

            /**
             * @brief
             * @tparam __type
             */
            template<lolita::core::finite_element::unknown::UnknownType __type>
            requires(__type == lolita::core::finite_element::unknown::UnknownType::Structural())
            struct DegreeOfFreedom<__type>
            {

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
                 */
                std::unique_ptr<_CoefficientVector> unknowns_coefficients_;

                /**
                 * @brief
                 */
                std::unique_ptr<_CoefficientVector> bindings_coefficients_;

                /**
                 * @brief
                 */
                std::unique_ptr<_CoordinatesVector> unknowns_coordinates_;

                /**
                 * @brief
                 */
                std::unique_ptr<_CoordinatesVector> bindings_coordinates_;

            };

            /**
             * @brief
             * @tparam __type
             */
            template<lolita::core::finite_element::unknown::UnknownType __type>
            requires(__type == lolita::core::finite_element::unknown::UnknownType::Subsidiary())
            struct DegreeOfFreedom<__type>
            {

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
                 */
                std::unique_ptr<_CoefficientVector> unknowns_coefficients_;

                /**
                 * @brief
                 */
                std::unique_ptr<_CoefficientVector> bindings_coefficients_;

            };

            /**
             * @brief
             * @param i
             * @param j
             */
            void
            setUnknownCoefficients(
                    lolita::integer i,
                    lolita::integer j
            )
            requires(_dim > 0 && _unknown_type == lolita::core::finite_element::unknown::UnknownType::Subsidiary())
            {
                degrees_of_freedom_[i][j].unknowns_coefficients = std::make_unique<_CoefficientVector>(_CoefficientVector::Zero());
            }

            /**
             * @brief
             * @param i
             * @param j
             * @param dim_unknown
             */
            void
            setUnknownCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::index dim_unknown
            )
            requires(_dim == -1 && _unknown_type == lolita::core::finite_element::unknown::UnknownType::Subsidiary())
            {
                degrees_of_freedom_[i][j].unknowns_coefficients_ = std::make_unique<_CoefficientVector>(_CoefficientVector::Zero(dim_unknown));
            }

            /**
             * @brief
             * @tparam _finite_element
             * @tparam __finite_element
             * @param i
             * @param j
             * @param mesh
             */
            template<
                    lolita::core::Element _element,
                    lolita::finite_element::FiniteElementConcept auto _finite_element,
                    lolita::finite_element::FiniteElementConcept auto... __finite_element
            >
            void
            setUnknownCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::core::mesh::Mesh<_domain, __finite_element...> & mesh
            )
            requires(_dim > 0 && _unknown_type == lolita::core::finite_element::unknown::UnknownType::Structural())
            {
                using _MixedElementDescription = lolita::core::MixedElementDescription<_element, _domain, __finite_element...>;
                auto constexpr _finite_element_index = _MixedElementDescription::template getFiniteElementIndex<_finite_element>();
                degrees_of_freedom_[i][j].unknowns_coefficients_ = std::make_unique<_CoefficientVector>(_CoefficientVector::Zero());
                auto & k = mesh.systems_[_finite_element_index].num_unknowns_;
                degrees_of_freedom_[i][j].unknowns_coordinates_ = std::make_unique<_CoordinatesVector>(_CoordinatesVector::LinSpaced(k, k + _dim));
                k += _dim;
            }

            /**
             * @brief
             * @tparam _finite_element
             * @tparam __finite_element
             * @param i
             * @param j
             * @param dim_unknown
             * @param mesh
             */
            template<
                    lolita::core::Element _element,
                    lolita::finite_element::FiniteElementConcept auto _finite_element,
                    lolita::finite_element::FiniteElementConcept auto... __finite_element
            >
            void
            setUnknownCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::index dim_unknown,
                    lolita::core::mesh::Mesh<_domain, __finite_element...> & mesh
            )
            requires(_dim == -1 && _unknown_type == lolita::core::finite_element::unknown::UnknownType::Structural())
            {
                using _MixedElementDescription = lolita::core::MixedElementDescription<_element, _domain, __finite_element...>;
                auto constexpr _finite_element_index = _MixedElementDescription::template getFiniteElementIndex<_finite_element>();
                degrees_of_freedom_[i][j].unknowns_coefficients_ = std::make_unique<_CoefficientVector>(_CoefficientVector::Zero(dim_unknown));
                auto & k = mesh.systems_[_finite_element_index].num_unknowns_;
                degrees_of_freedom_[i][j].unknowns_coordinates_ = std::make_unique<_CoordinatesVector>(_CoordinatesVector::LinSpaced(dim_unknown, k, k + dim_unknown));
                k += dim_unknown;
            }

            /**
             * @brief
             * @param i
             * @param j
             */
            void
            setBindingCoefficients(
                    lolita::integer i,
                    lolita::integer j
            )
            requires(_dim > 0 && _unknown_type == lolita::core::finite_element::unknown::UnknownType::Subsidiary())
            {
                degrees_of_freedom_[i][j].bindings_coefficients_ = std::make_unique<_CoefficientVector>(_CoefficientVector::Zero());
            }

            /**
             * @brief
             * @param i
             * @param j
             * @param dim_unknown
             */
            void
            setBindingCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::index dim_unknown
            )
            requires(_dim == -1 && _unknown_type == lolita::core::finite_element::unknown::UnknownType::Subsidiary())
            {
                degrees_of_freedom_[i][j].bindings_coefficients_ = std::make_unique<_CoefficientVector>(_CoefficientVector::Zero(dim_unknown));
            }

            /**
             * @brief
             * @tparam _finite_element
             * @tparam __finite_element
             * @param i
             * @param j
             * @param mesh
             */
            template<
                    lolita::core::Element _element,
                    lolita::finite_element::FiniteElementConcept auto _finite_element,
                    lolita::finite_element::FiniteElementConcept auto... __finite_element
            >
            void
            setBindingCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::core::mesh::Mesh<_domain, __finite_element...> & mesh
            )
            requires(_dim > 0 && _unknown_type == lolita::core::finite_element::unknown::UnknownType::Structural())
            {
                using _MixedElementDescription = lolita::core::MixedElementDescription<_element, _domain, __finite_element...>;
                auto constexpr _finite_element_index = _MixedElementDescription::template getFiniteElementIndex<_finite_element>();
                degrees_of_freedom_[i][j].bindings_coefficients_ = std::make_unique<_CoefficientVector>(_CoefficientVector::Zero());
                auto & k = mesh.systems_[_finite_element_index].num_bindings_;
                degrees_of_freedom_[i][j].bindings_coordinates_ = std::make_unique<_CoordinatesVector>(_CoordinatesVector::LinSpaced(k, k + _dim));
                k += _dim;
            }

            /**
             * @brief
             * @tparam _finite_element
             * @tparam __finite_element
             * @param i
             * @param j
             * @param dim_unknown
             * @param mesh
             */
            template<
                    lolita::core::Element _element,
                    lolita::finite_element::FiniteElementConcept auto _finite_element,
                    lolita::finite_element::FiniteElementConcept auto... __finite_element
            >
            void
            setBindingCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::index dim_unknown,
                    lolita::core::mesh::Mesh<_domain, __finite_element...> & mesh
            )
            requires(_dim == -1 && _unknown_type == lolita::core::finite_element::unknown::UnknownType::Structural())
            {
                using _MixedElementDescription = lolita::core::MixedElementDescription<_element, _domain, __finite_element...>;
                auto constexpr _finite_element_index = _MixedElementDescription::template getFiniteElementIndex<_finite_element>();
                degrees_of_freedom_[i][j].bindings_coefficients_ = std::make_unique<_CoefficientVector>(_CoefficientVector::Zero(dim_unknown));
                auto & k = mesh.systems_[_finite_element_index].num_bindings_;
                degrees_of_freedom_[i][j].bindings_coordinates_ = std::make_unique<_CoordinatesVector>(_CoordinatesVector::LinSpaced(dim_unknown, k, k + dim_unknown));
                k += dim_unknown;
            }

            /**
             * @brief
             */
            std::array<std::array<DegreeOfFreedom<_unknown_type>, _Field::shape_.rows_>, _Field::shape_.cols_> degrees_of_freedom_;

        };

        /**
         * @brief
         * @tparam _T
         */
        template<lolita::core::finite_element::unknown::UnknownConcept... _T>
        struct Unknowns
        {

        private:

            /**
             * @brief
             */
            using _Unknowns = std::tuple<_T...>;

        public:

            /**
             * @brief
             */
            std::array<lolita::integer, sizeof...(_T)> const static constexpr dim_unknowns_ = {_T::dim_...};

            /**
             * @brief
             */
            lolita::integer const static constexpr num_unknowns_ = sizeof...(_T);

            /**
             * @brief
             * @return
             */
            static constexpr
            lolita::boolean
            hasStructuralUnknowns()
            {
                return ((_T::unknown_type_ == lolita::core::finite_element::unknown::UnknownType::Structural) || ...);
            }

            /**
             * @brief
             * @tparam _domain
             * @tparam _unknown_type
             * @return
             */
            template<lolita::geometry::Domain _domain, lolita::core::finite_element::unknown::UnknownType _unknown_type>
            static constexpr
            lolita::integer
            getNumUnknowns()
            {
                auto num_unknowns = lolita::integer(0);
                auto set_num_unknowns  = [&] <lolita::integer _i = 0> (auto & self) constexpr mutable {
                    if constexpr (std::tuple_size_v<_Unknowns> > 0) {
                        using _Unknown = std::tuple_element_t<_i, _Unknowns>;
                        if constexpr (_Unknown::unknown_type_ == _unknown_type) {
                            num_unknowns += _Unknown::dim_ * lolita::core::field::TensorPolicy<_Unknown::field_, _domain.dim_>::shape_.size_;
                        }
                        if constexpr (_i < std::tuple_size_v<_Unknowns> - 1) {
                            self.template operator()<_i + 1>(self);
                        }
                    }
                };
                set_num_unknowns(set_num_unknowns);
                return num_unknowns;
            }

            /**
             * @brief
             * @tparam _i
             * @return
             */
            template<lolita::integer _i>
            std::tuple_element_t<_i, _Unknowns> const &
            getUnknown()
            const
            {
                return std::get<_i>(unknowns_);
            }

            /**
             * @brief
             * @tparam _i
             * @return
             */
            template<lolita::integer _i>
            std::tuple_element_t<_i, _Unknowns> &
            getUnknown()
            {
                return std::get<_i>(unknowns_);
            }

            /**
             * @brief
             */
            _Unknowns unknowns_;

        };

        /**
         * @brief Default initialization is zero unknowns
         * @tparam _element
         * @tparam _domain
         * @tparam _finite_element
         */
        template<lolita::core::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
        struct FiniteElementUnknowns
        {

            /**
             * @brief
             */
            using Unknowns = lolita::core::finite_element::unknown::Unknowns<>;

        };

        /**
         * @brief
         * @tparam _element
         * @tparam _domain
         * @tparam _finite_element
         */
        template<lolita::core::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
        requires(_finite_element.hasMethod(lolita::finite_element::Method::HHO) && lolita::core::ElementDescription<_element, _domain>().isCell())
        struct FiniteElementUnknowns<_element, _domain, _finite_element>
        {

        private:

            /**
             * @brief
             */
            using _Basis = lolita::core::finite_element::basis::FiniteElementBasis<
                    _element,
                    lolita::finite_element::Basis::Monomial,
                    _finite_element.discretization_.ord_cell_
            >;

        public:

            /**
             * @brief
             */
            using Unknowns = lolita::core::finite_element::unknown::Unknowns<
                    lolita::core::finite_element::unknown::Unknown<
                            _domain,
                            _finite_element.unknown_.tensor_,
                            _Basis::dim_,
                            lolita::core::finite_element::unknown::UnknownType::Subsidiary()
                    >
            >;

        };

        /**
         * @brief
         * @tparam _element
         * @tparam _domain
         * @tparam _finite_element
         */
        template<lolita::core::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
        requires(_finite_element.hasMethod(lolita::finite_element::Method::HHO) && lolita::core::ElementDescription<_element, _domain>().isFace())
        struct FiniteElementUnknowns<_element, _domain, _finite_element>
        {

        private:

            /**
             * @brief
             */
            using _Basis = lolita::core::finite_element::basis::FiniteElementBasis<
                    _element,
                    lolita::finite_element::Basis::Monomial,
                    _finite_element.discretization_.ord_face_
            >;

        public:

            /**
             * @brief
             */
            using Unknowns = lolita::core::finite_element::unknown::Unknowns<
                    lolita::core::finite_element::unknown::Unknown<
                            _domain,
                            _finite_element.unknown_.tensor_,
                            _Basis::dim_,
                            lolita::core::finite_element::unknown::UnknownType::Structural()
                    >
            >;

        };

    }

    namespace data
    {

        /**
         * @brief
         * @tparam _element
         * @tparam _domain
         * @tparam _finite_element
         */
        template<lolita::core::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
        struct FiniteElementData;

        template<lolita::core::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
        requires(_finite_element.hasMethd(lolita::finite_element::Method::HHO))
        struct FiniteElementData<_element, _domain, _finite_element>
        {

            struct Module
            {



            };

            struct Implementation : public lolita::core::finite_element::FEObject<_element, _domain, _finite_element>
            {



            };

        };

    }

    /**
     * @brief
     * @tparam _T
     * @tparam _element
     * @tparam _domain
     * @tparam _finite_element
     */
    template<
            template<lolita::core::Element, lolita::geometry::Domain, lolita::finite_element::FiniteElementConcept auto...> typename _T,
            lolita::core::Element _element,
            lolita::geometry::Domain _domain,
            lolita::finite_element::FiniteElementConcept auto... _finite_element
    >
    struct FiniteElementGeometry
    {

    private:

        /**
         * @brief
         */
        using _ElementDescription = lolita::core::ElementDescription<_element, _domain>;

        /**
         * @brief
         */
        using _ElementGeometry = lolita::core::element::ElementGeometry<_element, _domain>;

        /**
         * @brief
         */
        template<lolita::core::Element __element, lolita::geometry::Domain __domain, auto... __finite_element>
        using _ElementPointer = std::shared_ptr<_T<__element, __domain, __finite_element...>>;

        /**
         * @brief
         */
        template<lolita::core::Element>
        struct _ConnexionPolicy;

        /**
         * @brief
         * @tparam __element
         */
        template<lolita::core::Element __element>
        requires(!__element.isPoint())
        struct _ConnexionPolicy<__element>
        {

            /**
             * @brief
             */
            using _Components = typename lolita::core::element::ElementGeometry<__element, _domain>::template Components<_ElementPointer, _finite_element...>;

            /**
             * @brief
             */
            using _Neighbours = typename lolita::core::element::ElementGeometry<__element, _domain>::template Neighbours<_ElementPointer, _finite_element...>;

        };

        /**
         * @brief
         * @tparam __element
         */
        template<lolita::core::Element __element>
        requires(__element.isPoint())
        struct _ConnexionPolicy<__element>
        {

            /**
             * @brief
             */
            using _Components = std::shared_ptr<typename lolita::core::element::ElementGeometry<__element, _domain>::Components>;

            /**
             * @brief
             */
            using _Neighbours = typename lolita::core::element::ElementGeometry<__element, _domain>::template Neighbours<_ElementPointer, _finite_element...>;

        };

    public:

        /**
         * @brief The element object
         */
        lolita::core::Element const static constexpr element_ = _element;

        /**
         * @brief Element components as shared pointers
         */
        using Components = typename _ConnexionPolicy<_element>::_Components;

        /**
         * @brief Element neighbours as shared pointers
         */
        using Neighbours = typename _ConnexionPolicy<_element>::_Neighbours;

        /**
         * @brief Basic equality comparison operator
         * @param other
         * @return
         */
        lolita::boolean
        operator==(
                FiniteElementGeometry const & other
        )
        const = default;

        /**
         * @brief Basic equality comparison operator
         * @param other
         * @return
         */
        lolita::boolean
        operator!=(
                FiniteElementGeometry const & other
        )
        const = default;

        /**
         * @brief
         * @return
         */
        std::basic_string<lolita::character>
        hash()
        const
        {
            if constexpr (_element.isPoint()) {
                return std::to_string(this->tag_);
            }
            else {
                std::basic_stringstream<lolita::character> hash;
                auto const & nodes = this->template getComponents<_element.dim_ - 1, 0>();
                for (auto const & node : nodes) {
                    hash << node->hash();
                }
                return hash.str();
            }
        }

        /**
         * @brief Get the current coordinates of the point/node in the mesh
         * @return the vector with x, y and z coordinates
         */
        lolita::geometry::Point &
        getCurrentCoordinates()
        requires(_element.isPoint())
        {
            return * components_;
        }

        /**
         * @brief Get the current coordinates of the point/node in the mesh
         * @return the vector with x, y and z coordinates
         */
        lolita::geometry::Point const &
        getCurrentCoordinates()
        const
        requires(_element.isPoint())
        {
            return * components_;
        }

        /**
         * @brief Get the current coordinates of the element in the mesh
         * @return the matrix composed of the column vectors with x, y and z coordinates of each node
         */
        lolita::matrix::Matrix<lolita::real, 3, _element.num_nodes_>
        getCurrentCoordinates()
        const
        requires(!_element.isPoint())
        {
            auto current_nodes_coordinates = lolita::matrix::Matrix<lolita::real, 3, _element.num_nodes_>();
            auto count = lolita::index(0);
            for (auto const & node : getComponents<_element.dim_ - 1, 0>()) {
                current_nodes_coordinates.col(count) = node->getCurrentCoordinates();
                count ++;
            }
            return current_nodes_coordinates;
        }

        /**
         * @brief
         * @return
         */
        static
        lolita::matrix::Span<lolita::matrix::Matrix<lolita::real, 3, _element.num_nodes_, matrix::col_major> const>
        getReferenceCoordinates()
        {
            using _ReferenceCoordinates = lolita::matrix::Span<lolita::matrix::Matrix<lolita::real, 3, _element.num_nodes_, matrix::col_major> const>;
            return _ReferenceCoordinates(lolita::core::element::ElementGeometry<_element, _domain>::reference_nodes_.begin()->begin());
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @param i
         * @param j
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        static constexpr
        lolita::index
        getComponentNodeConnection(
                lolita::index i,
                lolita::index j
        )
        requires(!_element.isPoint())
        {
            return std::get<_j>(std::get<_i>(lolita::core::element::ElementGeometry<_element, _domain>::node_connectivity_))[i][j];
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, Neighbours>> &
        getNeighbours()
        {
            return std::get<_j>(std::get<_i>(neighbours_));
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, Neighbours>> const &
        getNeighbours()
        const
        {
            return std::get<_j>(std::get<_i>(neighbours_));
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, Components>> &
        getComponents()
        requires(!_element.isPoint())
        {
            return std::get<_j>(std::get<_i>(components_));
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, Components>> const &
        getComponents()
        const
        requires(!_element.isPoint())
        {
            return std::get<_j>(std::get<_i>(components_));
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @param i
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        lolita::index
        getComponentIndex(
                lolita::index i
        )
        const
        requires(!_element.isPoint())
        {
//            auto const constexpr _component = lolita::core::component<_element, _domain, _i, _j>();
//            auto const constexpr _position = lolita::core::neighbourPosition<_component, _domain, _element>();
            auto const constexpr _component = _ElementDescription::template getComponent<_i, _j>();
            auto const constexpr _position = lolita::core::ElementDescription<_component, _domain>::template getNeighbourCoordinates<_element>();
            auto const & items = getComponents<_i, _j>()[i]->template getNeighbours<_position.dim_, _position.tag_>();
            auto is_equal = [&] (
                    std::shared_ptr<FiniteElementGeometry> const & ptr_element
            )
            {
                return * ptr_element == * this;
            };
            return std::distance(items.begin(), std::find_if(items.begin(), items.end(), is_equal));
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @param i
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        lolita::integer
        getComponentOrientation(
                lolita::index i
        )
        const
        requires(!_element.isPoint())
        {
            return getComponentIndex<_i, _j>(i) == 0 ? 1 : -1;
        }

        /**
         * @brief
         * @tparam _basis
         * @tparam _ord
         * @param point
         * @return
         */
        template<lolita::finite_element::Basis _basis, lolita::index _ord>
        lolita::matrix::Vector<lolita::real, lolita::core::finite_element::basis::FiniteElementBasis<_element, _basis, _ord>::dim_>
        getBasisEvaluation(
                lolita::geometry::Point const & point
        )
        const
        {
            using _FiniteElementBasis = lolita::core::finite_element::basis::FiniteElementBasis<_element, _basis, _ord>;
            using _FiniteElementBasisImplementation = typename _FiniteElementBasis::template Implementation<_T, _domain, _finite_element...>;
            return static_cast<_FiniteElementBasisImplementation const *>(this)->evaluate(point);
        }

        /**
         * @brief
         * @tparam _basis
         * @tparam _ord
         * @param point
         * @param derivative_direction
         * @return
         */
        template<lolita::finite_element::Basis _basis, lolita::index _ord>
        lolita::matrix::Vector<lolita::real, lolita::core::finite_element::basis::FiniteElementBasis<_element, _basis, _ord>::dim_>
        getBasisDerivative(
                lolita::geometry::Point const & point,
                lolita::index derivative_direction
        )
        const
        {
            using _FiniteElementBasis = lolita::core::finite_element::basis::FiniteElementBasis<_element, _basis, _ord>;
            using _FiniteElementBasisImplementation = typename _FiniteElementBasis::template Implementation<_T, _domain, _finite_element...>;
            return static_cast<_FiniteElementBasisImplementation const *>(this)->evaluate(point, derivative_direction);
        }

        /**
         * @brief
         * @param nodal_field_values
         * @param reference_point
         * @return
         */
        static
        lolita::real
        getShapeMappingEvaluation(
                lolita::matrix::Vector<lolita::real, _element.num_nodes_> const & nodal_field_values,
                lolita::geometry::Point const & reference_point
        )
        {
            return lolita::core::element::ElementGeometry<_element, _domain>::getShapeMappingEvaluation(nodal_field_values, reference_point);
        }

        /**
         * @brief
         * @param nodal_field_values
         * @param reference_point
         * @param derivative_direction
         * @return
         */
        static
        lolita::real
        getShapeMappingDerivative(
                lolita::matrix::Vector<lolita::real, _element.num_nodes_> const & nodal_field_values,
                lolita::geometry::Point const & reference_point,
                lolita::index derivative_direction
        )
        {
            return lolita::core::element::ElementGeometry<_element, _domain>::getShapeMappingDerivative(nodal_field_values, reference_point, derivative_direction);
        }

        /*
         *
         */

        /**
         * @brief
         * @param point
         * @return
         */
        lolita::real
        getShapeMappingDifferential(
                lolita::geometry::Point const & point
        )
        const
        {
            auto const current_coordinates = this->getCurrentCoordinates();
            auto ru = lolita::matrix::Matrix<lolita::real, 3, 3>();
            auto du = lolita::real(0);
            ru.setZero();
            for (lolita::index i = 0; i < _domain.dim_; ++i) {
                for (lolita::index j = 0; j < _element.dim_; ++j) {
                    ru(i, j) = FiniteElementGeometry::getShapeMappingDerivative(current_coordinates.row(i), point, j);
                }
            }
            if constexpr (_element.dim_ == 3) {
                du = numerics::abs((ru.col(0).template cross(ru.col(1))).template dot(ru.col(2)));
            }
            else if constexpr (_element.dim_ == 2) {
                du = numerics::abs((ru.col(0).template cross(ru.col(1))).norm());
            }
            else {
                du = numerics::abs(ru.col(0).norm());
            }
            if constexpr (_domain.frame_ == lolita::geometry::Frame::AxiSymmetric) {
                lolita::real r0 = getShapeMappingEvaluation(current_coordinates.row(0), point);
                if (r0 < 1.e-10) {
                    r0 = 1.e-10;
                }
                du *= 2.0 * numerics::pi * r0;
            }
            return du;
        }

        /**
         * @brief
         * @param first_point
         * @param second_point
         * @param direction
         * @return
         */
        lolita::real
        getRiemannianDistance(
                lolita::geometry::Point const & first_point,
                lolita::geometry::Point const & second_point,
                lolita::integer direction = -1
        )
        const
        {
            if constexpr (_ElementDescription::isCell()) {
                auto const & current_coordinates = this->getCurrentCoordinates();
                auto distance = lolita::real();
                auto mp0 = lolita::geometry::Point();
                auto mp1 = lolita::geometry::Point();
                for (lolita::index i = 0; i < _element.dim_; ++i) {
                    mp0(i) = FiniteElementGeometry::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
                    mp1(i) = FiniteElementGeometry::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
                }
                direction == -1 ? distance = (mp1 - mp0).norm() : distance = (mp1 - mp0)(direction);
                return distance;
            }
            else {
                auto const constexpr _quadrature = lolita::finite_element::Quadrature::Gauss();
                using SegmentQuadrature = lolita::core::element::ElementQuadrature<lolita::core::Element::LinearSegment(), _quadrature, 4>;
                auto distance = lolita::real(0);
                auto dt = lolita::real();
                auto const current_nodes_coordinates = this->getCurrentCoordinates();
                for (lolita::index q = 0; q < SegmentQuadrature::dim_; ++q) {
                    auto pq = SegmentQuadrature::reference_points_[q][0];
                    auto wq = SegmentQuadrature::reference_weights_[q];
                    auto ru = lolita::matrix::Matrix<lolita::real, 3, 3>().setZero();
                    auto difference = second_point - first_point;
                    auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
                    for (lolita::index i = 0; i < _domain.dim_; ++i) {
                        for (lolita::index j = 0; j < _element.dim_; ++j) {
                            if (direction == -1 || i == static_cast<lolita::index>(direction)) {
                                auto du = (1.0 / 2.0) * (second_point(j) - first_point(j));
                                auto dx = FiniteElementGeometry::getShapeMappingDerivative(current_nodes_coordinates.row(i), uq, j);
                                ru(i, j) = dx * du;
                            }
                        }
                    }
                    if constexpr (_element.isFacet()) {
                        auto Eff = ru.col(0).template dot(ru.col(0));
                        auto Fff = ru.col(0).template dot(ru.col(1));
                        auto Gff = ru.col(1).template dot(ru.col(1));
                        dt = std::sqrt(Eff + 2.0 * Fff + Gff);
                    }
                    else if constexpr (_element.isCurve()) {
                        auto Eff = ru.col(0).template dot(ru.col(0));
                        dt = std::sqrt(Eff);
                    }
                    else {
                        dt = 0;
                    }
                    distance += wq * dt;
                }
                return distance;
            }
        }

        /*
         * GEOMETRY
         */

        /**
         * @brief
         * @return
         */
        lolita::geometry::Point
        getCurrentDiameters()
        const
        {
            auto reference_coordinates = FiniteElementGeometry::getReferenceCoordinates();
            auto current_diameters = lolita::geometry::Point().setZero();
            for (lolita::index i = 0; i < _element.num_nodes_; ++i) {
                for (lolita::index j = i + 1; j < _element.num_nodes_; ++j) {
                    auto const & pt0 = reference_coordinates.col(i);
                    auto const & pt1 = reference_coordinates.col(j);
                    for (lolita::index k = 0; k < 3; ++k) {
                        auto new_value = numerics::abs(getRiemannianDistance(pt0, pt1, k));
                        auto & current_value = current_diameters(k);
                        if (new_value > current_value) {
                            current_value = new_value;
                        }
                    }
                }
            }
            return current_diameters;
        }

        /**
         * @brief
         * @return
         */
        lolita::geometry::Point
        getCurrentCentroid()
        const
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            return lolita::geometry::getBarycenter(current_nodes_coordinates);
        }

        /**
         * @brief
         * @return
         */
        static
        lolita::geometry::Point
        getReferenceCentroid()
        {
            auto nds = FiniteElementGeometry::getReferenceCoordinates();
            return lolita::geometry::getBarycenter(nds);
        }

        /**
         * @brief
         * @return
         */
        static
        lolita::geometry::Point
        getReferenceDiameters()
        {
            auto dts = lolita::geometry::Point();
            auto nds = FiniteElementGeometry::getReferenceCoordinates();
            dts.setZero();
            for (lolita::index i = 0; i < _element.num_nodes_; ++i) {
                for (lolita::index j = i + 1; j < _element.num_nodes_; ++j) {
                    for (lolita::index k = 0; k < 3; ++k) {
                        auto new_value = numerics::abs(nds(k, i) - nds(k, j));
                        auto & current_value = dts(k);
                        if (new_value > current_value) {
                            current_value = new_value;
                        }
                    }
                }
            }
            return dts;
        }

        /*
         * QUADRATURE
         */

        /**
         * @brief
         * @tparam _quadrature
         * @tparam _ord
         * @param index
         * @return
         */
        template<lolita::finite_element::Quadrature _quadrature, lolita::index _ord>
        static
        lolita::real
        getReferenceQuadratureWeight(
                lolita::index index
        )
        {
            return lolita::core::element::ElementQuadrature<_element, _quadrature, _ord>::reference_weights_[index];
        }

        /**
         * @brief
         * @tparam _quadrature
         * @tparam _ord
         * @param index
         * @return
         */
        template<lolita::finite_element::Quadrature _quadrature, lolita::index _ord>
        static
        lolita::matrix::Span<lolita::geometry::Point const>
        getReferenceQuadraturePoint(
                lolita::index index
        )
        {
            return lolita::matrix::Span<lolita::geometry::Point const>(
                    lolita::core::element::ElementQuadrature<_element, _quadrature, _ord>::reference_points_[index].begin()
            );
        }

        /**
         * @brief
         * @tparam _quadrature
         * @tparam _ord
         * @param index
         * @return
         */
        template<lolita::finite_element::Quadrature _quadrature, lolita::index _ord>
        lolita::real
        getCurrentQuadratureWeight(
                lolita::index index
        )
        const
        {
            auto w = getReferenceQuadratureWeight<_quadrature, _ord>(index);
            return w * getShapeMappingDifferential(getReferenceQuadraturePoint<_quadrature, _ord>(index));
        }

        /**
         * @brief
         * @tparam _quadrature
         * @tparam _ord
         * @param index
         * @return
         */
        template<lolita::finite_element::Quadrature _quadrature, lolita::index _ord>
        lolita::geometry::Point
        getCurrentQuadraturePoint(
                lolita::index index
        )
        const
        {
            auto p = lolita::geometry::Point();
            auto const nds = this->getCurrentCoordinates();
            for (lolita::index j = 0; j < 3; ++j) {
                p(j) = FiniteElementGeometry::getShapeMappingEvaluation(nds.row(j), getReferenceQuadraturePoint<_quadrature, _ord>(index));
            }
            return p;
        }

        /**
         * @brief
         * @tparam _quadrature
         * @tparam _ord
         * @tparam _i
         * @tparam _j
         * @param index
         * @return
         */
        template<lolita::finite_element::Quadrature _quadrature, lolita::index _ord, lolita::index _i, lolita::index _j>
        static
        lolita::real
        getComponentReferenceQuadratureWeight(
                lolita::index index
        )
        requires(!_element.isPoint())
        {
            auto const constexpr _component = _ElementDescription::template getComponent<_i, _j>();
            using ComponentGeometry = FiniteElementGeometry<_T, _component, _domain, _finite_element...>;
            return ComponentGeometry::template getReferenceQuadratureWeight<_quadrature, _ord>(index);
        }

        /**
         * @brief
         * @tparam _quadrature
         * @tparam _ord
         * @tparam _i
         * @tparam _j
         * @param component_index
         * @param index
         * @return
         */
        template<lolita::finite_element::Quadrature _quadrature, lolita::index _ord, lolita::index _i, lolita::index _j>
        static
        lolita::geometry::Point
        getComponentReferenceQuadraturePoint(
                lolita::index component_index,
                lolita::index index
        )
        requires(!_element.isPoint())
        {
//            auto const constexpr _component = lolita::core::component<_element, _domain, _i, _j>();
            auto const constexpr _component = _ElementDescription::template getComponent<_i, _j>();
            auto p = lolita::geometry::Point();
            using ComponentGeometry = FiniteElementGeometry<_T, _component, _domain, _finite_element...>;
            auto const & elt_reference_nodes = lolita::core::element::ElementGeometry<_element, _domain>::reference_nodes_;
            for (lolita::index i = 0; i < 3; ++i) {
                auto cpt_coordinates = lolita::matrix::Vector<lolita::real, _component.num_nodes_>();
                for (lolita::index j = 0; j < _component.num_nodes_; ++j) {
                    auto const node_tag = getComponentNodeConnection<_i, _j>(component_index, j);//.get(component_index).get(j);
                    cpt_coordinates(j) = elt_reference_nodes[node_tag][i];
                }
                auto cpt_reference_point = ComponentGeometry::template getReferenceQuadraturePoint<_quadrature, _ord>(index);
                p(i) = ComponentGeometry::getShapeMappingEvaluation(cpt_coordinates, cpt_reference_point);
            }
            return p;
        }

        /**
         * @brief
         * @tparam _quadrature
         * @tparam _ord
         * @tparam _i
         * @tparam _j
         * @param component_index
         * @param index
         * @return
         */
        template<lolita::finite_element::Quadrature _quadrature, lolita::index _ord, lolita::index _i, lolita::index _j>
        lolita::real
        getComponentCurrentQuadratureWeight(
                lolita::index component_index,
                lolita::index index
        )
        const
        requires(!_element.isPoint())
        {
            auto const & cmp =  this->template getComponents<_i, _j>()[component_index];//.template get<I>().template get<J>().get(component_index).get();
            return cmp->template getCurrentQuadratureWeight<_quadrature, _ord>(index);
        }

        /**
         * @brief
         * @tparam _quadrature
         * @tparam _ord
         * @tparam _i
         * @tparam _j
         * @param component_index
         * @param index
         * @return
         */
        template<lolita::finite_element::Quadrature _quadrature, lolita::index _ord, lolita::index _i, lolita::index _j>
        lolita::geometry::Point
        getComponentCurrentQuadraturePoint(
                lolita::index component_index,
                lolita::index index
        )
        const
        requires(!_element.isPoint())
        {
            auto p = lolita::geometry::Point();
            auto const cpt_ref_pnt = getComponentReferenceQuadraturePoint<_quadrature, _ord, _i, _j>(component_index, index);
            auto const nds = this->getCurrentCoordinates();
            for (lolita::index j = 0; j < 3; ++j) {
                p(j) = FiniteElementGeometry::getShapeMappingEvaluation(nds.row(j), cpt_ref_pnt);
            }
            return p;
        }

        /**
         * @brief
         * @param point
         * @return
         */
        lolita::geometry::Point
        getNormalVector(
                lolita::geometry::Point const & point
        )
        const
        requires(_ElementDescription::isFace())
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto ru = lolita::matrix::Matrix<lolita::real, 3, _element.dim_>();
            ru.setZero();
            for (lolita::index i = 0; i < 3; ++i) {
                for (lolita::index j = 0; j < _element.dim_; ++j) {
                    ru(i, j) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, j);
                }
            }
            if constexpr (_element.isPoint()) {
                return lolita::geometry::Point{0, 0, 0};
            }
            else if constexpr (_element.isCurve()) {
                return lolita::geometry::Point{ru(1)/ru.norm(), -ru(0)/ru.norm(), 0};
            }
            else {
                return (ru.col(0) / ru.col(0).norm()).cross((ru.col(1) / ru.col(1).norm()));
            }
        }

        /**
         * @brief
         * @param point
         * @param direction
         * @return
         */
        lolita::geometry::Point
        getTangentVector(
                lolita::geometry::Point const & point,
                lolita::index direction
        )
        const
        requires(_element.isCurve())
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto tangent_vector = lolita::geometry::Point();
            for (lolita::index i = 0; i < 3; ++i) {
                tangent_vector(i) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, direction);
            }
            return tangent_vector;
        }

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
        lolita::natural tag_;

    };


    template<lolita::core::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    struct IntegrationPoint
    {

    private:

        /**
         * @brief
         */
        using _FiniteElementDescription = lolita::core::FiniteElementDescription<_element, _domain, _finite_element>;

        /**
         * @brief
         */
        using _Quadrature = typename _FiniteElementDescription::Quadrature;

        template<lolita::core::Element __element>
        struct _OperatorsPolicy;

        template<lolita::core::Element __element>
        requires(ElementDescription<__element, _domain>::isCell())
        struct _OperatorsPolicy<__element>
        {

            using _LhsOperator = lolita::matrix::Matrix<lolita::real, _FiniteElementDescription::getOperatorNumRows(), _FiniteElementDescription::getNumUnknowns()>;

        };

        template<lolita::core::Element __element>
        requires(ElementDescription<__element, _domain>::isFace())
        struct _OperatorsPolicy<__element>
        {



        };

        /**
         * @brief
         */
        using _LhsOperator = lolita::matrix::Matrix<lolita::real, _FiniteElementDescription::getOperatorNumRows(), _FiniteElementDescription::getNumUnknowns()>;

        /**
         * @brief
         */
        using _RhsOperator = lolita::matrix::Vector<lolita::real, _FiniteElementDescription::getNumOwnedUnknowns()>;

    public:

    };

    template<lolita::core::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    struct FEObject;

    /**
     * @brief
     */
    auto const static null_load_ptr = std::make_shared<lolita::finite_element::LoadComponent>(lolita::finite_element::LoadComponent());

    /**
     * @brief
     * @tparam _element
     * @tparam _domain
     * @tparam _finite_element
     */
    template<lolita::core::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    struct FiniteElementBase1 : public lolita::core::finite_element::FiniteElementGeometry<FEObject, _element, _domain, _finite_element>
    {

    private:

        /**
         * @brief
         */
        using _FiniteElementDescription = lolita::core::FiniteElementDescription<_element, _domain, _finite_element>;

        /**
         * @brief
         */
        using _Quadrature = typename _FiniteElementDescription::Quadrature;

    public:

        /**
         * @brief
         */
        using Quadrature = lolita::core::element::ElementQuadrature<_element, _finite_element.quadrature_, _finite_element.ord_quadrature_>;

        /**
         * @brief
         */
        using Field = lolita::core::field::TensorPolicy<_finite_element.unknown_.tensor_, _domain.dim_>;

    private:

        /**
         * @brief loads initializer
         * @return loads member initialized with the zero natural load
         */
        static
        std::array<std::array<std::shared_ptr<lolita::finite_element::LoadComponent>, Field::shape_.rows_>, Field::shape_.cols_>
        _loads()
        {
            auto loads = std::array<std::array<std::shared_ptr<lolita::finite_element::LoadComponent>, Field::shape_.rows_>, Field::shape_.cols_>();
            for (int i = 0; i < Field::shape_.rows_; ++i) {
                for (int j = 0; j < Field::shape_.cols_; ++j) {
                    loads[j][i] = null_load_ptr;
                }
            }
            return loads;
        }

    public:

        /**
         * @brief
         * @param row
         * @param col
         * @return
         */
        std::shared_ptr<lolita::finite_element::LoadComponent> const &
        getLoad(
                lolita::integer row,
                lolita::integer col
        )
        const
        {
            return loads_[col][row];
        }

        /**
         * @brief
         * @param row
         * @param col
         * @return
         */
        std::shared_ptr<lolita::finite_element::LoadComponent> &
        getLoad(
                lolita::integer row,
                lolita::integer col
        )
        {
            return loads_[col][row];
        }

        /**
         * @brief
         * @tparam __finite_element
         * @param mesh
         */
        template<lolita::finite_element::FiniteElementConcept auto... __finite_element>
        void
        setLoads(
                lolita::core::mesh::Mesh<_domain, __finite_element...> & mesh
        )
        {
            for (auto const & load : mesh.loads_) {
                auto is_equal = [&] (
                        std::shared_ptr<std::basic_string<lolita::character>> const & domain
                )
                {
                    return * domain == load.domain_tag_;
                };
                auto const has_domain = std::find_if(this->domains_.begin(), this->domains_.end(), is_equal) != this->domains_.end();
                auto const has_unknown = load.unknown_tag_ == lolita::utility::readLabel(_finite_element.unknown_.tensor_.tag_);
                if (has_domain && has_unknown) {
                    std::cout << "setting load : " << load.domain_tag_ << std::endl;
                    auto const i = load.components_.row_;
                    auto const j = load.components_.col_;
                    loads_[i][j] = load.load_;
                }
            }
        }

        /**
         * @brief
         * @tparam __finite_element
         * @param mesh
         */
        template<lolita::finite_element::FiniteElementConcept auto... __finite_element>
        void
        setBehaviour(
                lolita::core::mesh::Mesh<_domain, __finite_element...> & mesh
        )
        {
            for (auto const & behaviour : mesh.behaviours_) {
                auto is_equal = [&] (
                        std::shared_ptr<std::basic_string<lolita::character>> const & domain
                )
                {
                    return * domain == behaviour.domain_tag_;
                };
                auto const has_domain = std::find_if(this->domains_.begin(), this->domains_.end(), is_equal) != this->domains_.end();
                auto const has_unknown = behaviour.unknown_tag_ == lolita::utility::readLabel(_finite_element.unknown_.tensor_.tag_);
                if (has_domain && has_unknown) {
                    std::cout << "setting behaviour : " << behaviour.domain_tag_ << std::endl;
                    behaviour_ = behaviour.behaviour_data_;
                }
            }
        }

        /**
         * @brief
         */
        std::array<std::array<std::shared_ptr<lolita::finite_element::LoadComponent>, Field::shape_.rows_>, Field::shape_.cols_> loads_ = _loads();

        /**
         * @brief
         */
        std::shared_ptr<lolita::behaviour::MgisBehaviourData> behaviour_;

    };


    template<lolita::core::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    struct FEObject : public lolita::core::finite_element::FiniteElementBase1<_element, _domain, _finite_element>
    {

        template<auto... __args>
        void
        initialize(
                lolita::core::mesh::Mesh<_domain, __args...> & mesh
        )
        {

        }

    };

    template<
            template<lolita::core::Element, lolita::geometry::Domain, auto...> typename _Base,
            template<lolita::core::Element, lolita::geometry::Domain, auto...> typename _Data
    >
    struct ElementObjectLayerPolicy
    {

        template<lolita::core::Element _element, lolita::geometry::Domain _domain, auto... _args>
        struct ElementObjectLayer : public lolita::core::finite_element::FiniteElementGeometry<_Base, _element, _domain, _args...>
        {

            /**
             * @brief
             */
            using ElementPointers = std::tuple<std::shared_ptr<_Data<_element, _domain, _args>>...>;

            /**
             * @brief
             */
            using Elements = std::tuple<_Data<_element, _domain, _args>...>;

            /**
             * @brief
             * @tparam _i
             * @return
             */
            template<lolita::index _i>
            std::tuple_element_t<_i, ElementPointers> const &
            getElement()
            const
            {
                return std::get<_i>(elements_);
            }

            /**
             * @brief
             * @tparam _i
             * @return
             */
            template<lolita::index _i>
            std::tuple_element_t<_i, ElementPointers> &
            getElement()
            {
                return std::get<_i>(elements_);
            }

            template<lolita::core::Element __element, auto... __args>
            void
            initialize(
                    lolita::core::mesh::Mesh<_domain, _args...> & mesh_data
            );

            template<lolita::core::Element __element, auto... __args>
            void
            initialize(
                    lolita::core::mesh::Mesh<_domain, __args...> & mesh_data
            )
            requires(!__element.isPoint())
            {
                /*
                 *
                 */
                auto initialize_element = [&] <lolita::index _k = 0u> (
                        auto & initialize_element_imp
                )
                        mutable
                {
                    /*
                     *
                     */
                    auto initialize_components = [&] <lolita::index __i = 0u, lolita::index __j = 0u> (
                            auto & initialize_components_imp
                    )
                            mutable
                    {
                        for (int i = 0; i < lolita::core::ElementDescription<_element, _domain>::template getNumComponents<__i, __j>(); ++i) {
                            auto & rhs = this->template getComponents<__i, __j>()[i]->template getElement<_k>();
                            auto & lhs = this->template getElement<_k>()->template getComponents<__i, __j>()[i];
                            lhs = rhs;
                        }
                        if constexpr (__j < lolita::core::ElementDescription<_element, _domain>::template getNumComponents<__i>() - 1) {
                            initialize_components_imp.template operator()<__i, __j + 1u>(initialize_components_imp);
                        }
                        else if constexpr (__i < lolita::core::ElementDescription<_element, _domain>::getNumComponents() - 1) {
                            initialize_components_imp.template operator()<__i + 1u, 0u>(initialize_components_imp);
                        }
                    };
                    /*
                     *
                     */
                    auto initialize_neighbours = [&] <lolita::index __i = 0u, lolita::index __j = 0u> (
                            auto & initialize_neighbours_imp
                    )
                            mutable
                    {
                        for (int i = 0; i < this->template getNeighbours<__i, __j>().size(); ++i) {
                            auto & rhs = this->template getNeighbours<__i, __j>()[i]->template getElement<_k>();
                            auto & lhs = this->template getElement<_k>()->template getNeighbours<__i, __j>();
                            lhs.push_back(rhs);
                        }
                        if constexpr (__j < lolita::core::ElementDescription<_element, _domain>::template getNumNeighbours<__i>() - 1) {
                            initialize_neighbours_imp.template operator()<__i, __j + 1u>(initialize_neighbours_imp);
                        }
                        else if constexpr (__i < lolita::core::ElementDescription<_element, _domain>::getNumNeighbours() - 1) {
                            initialize_neighbours_imp.template operator()<__i + 1u, 0u>(initialize_neighbours_imp);
                        }
                    };
                    /*
                     *
                     */
                    initialize_components(initialize_components);
                    initialize_neighbours(initialize_neighbours);
                    this->template getElement<_k>()->tag_ = this->tag_;
                    this->template getElement<_k>()->domains_ = this->domains_;
                    this->template getElement<_k>()->initialize(mesh_data);
                    if constexpr (_k < sizeof...(_args) - 1) {
                        initialize_element_imp.template operator()<_k + 1u>(initialize_element_imp);
                    }
                };
                /*
                 *
                 */
                initialize_element(initialize_element);
            }

            template<lolita::core::Element __element, auto... __args>
            void
            initialize(
                    lolita::core::mesh::Mesh<_domain, __args...> & mesh_data
            )
            requires(__element.isPoint())
            {
                /*
                  *
                  */
                auto initialize_element = [&] <lolita::index _k = 0u> (
                        auto & initialize_element_imp
                )
                        constexpr mutable
                {
                    /*
                     *
                     */
                    auto initialize_coordinates = [&] ()
                    {
                        this->template getElement<_k>()->components_ = this->components_;
                    };
                    /*
                     *
                     */
                    auto initialize_neighbours = [&] <lolita::index __i = 0u, lolita::index __j = 0u> (
                            auto & initialize_neighbours_imp
                    )
                            mutable
                    {
                        for (int i = 0; i < this->template getNeighbours<__i, __j>().size(); ++i) {
                            auto & rhs = this->template getNeighbours<__i, __j>()[i]->template getElement<_k>();
                            auto & lhs = this->template getElement<_k>()->template getNeighbours<__i, __j>();
                            lhs.push_back(rhs);
                        }
                        if constexpr (__j < lolita::core::ElementDescription<_element, _domain>::template getNumNeighbours<__i>() - 1) {
                            initialize_neighbours_imp.template operator()<__i, __j + 1u>(initialize_neighbours_imp);
                        }
                        else if constexpr (__i < lolita::core::ElementDescription<_element, _domain>::getNumNeighbours() - 1) {
                            initialize_neighbours_imp.template operator()<__i + 1u, 0u>(initialize_neighbours_imp);
                        }
                    };
                    /*
                     *
                     */
                    initialize_coordinates();
                    initialize_neighbours(initialize_neighbours);
                    this->template getElement<_k>()->tag_ = this->tag_;
                    this->template getElement<_k>()->domains_ = this->domains_;
                    this->template getElement<_k>()->initialize(mesh_data);
                    if constexpr (_k < sizeof...(_args) - 1) {
                        initialize_element_imp.template operator()<_k + 1u>(initialize_element_imp);
                    }
                };
                /*
                 *
                 */
                initialize_element(initialize_element);
            }

            ElementPointers elements_;

        };

    };

    template<
            template<lolita::core::Element, lolita::geometry::Domain> typename _Base,
            template<lolita::core::Element, lolita::geometry::Domain> typename _Data
    >
    struct ElementObjectLayerPolicy2
    {

        template<lolita::core::Element _element, lolita::geometry::Domain _domain>
        struct ElementObjectLayer : public lolita::core::finite_element::FiniteElementGeometry<_Base, _element, _domain>
        {

            using DataPointer = _Data<_element, _domain>; // _Data = template<_element, _domain> std::tuple<std::shared_ptr<Type, _element, _domain, _args...>>

            /**
             * @brief
             * @tparam _i
             * @return
             */
            template<lolita::index _i>
            std::tuple_element_t<_i, DataPointer> const &
            getElement()
            const
            {
                return std::get<_i>(elements_);
            }

            /**
             * @brief
             * @tparam _i
             * @return
             */
            template<lolita::index _i>
            std::tuple_element_t<_i, DataPointer> &
            getElement()
            {
                return std::get<_i>(elements_);
            }

            template<lolita::core::Element __element, auto... __args>
            void
            initialize(
                    lolita::core::mesh::Mesh<_domain, __args...> & mesh_data
            );

            template<lolita::core::Element __element, auto... __args>
            void
            initialize(
                    lolita::core::mesh::Mesh<_domain, __args...> & mesh_data
            )
            requires(!__element.isPoint())
            {
                /*
                 *
                 */
                auto initialize_element = [&] <lolita::index _k = 0u> (
                        auto & initialize_element_imp
                )
                        mutable
                {
                    /*
                     *
                     */
                    auto initialize_components = [&] <lolita::index __i = 0u, lolita::index __j = 0u> (
                            auto & initialize_components_imp
                    )
                            mutable
                    {
                        for (int i = 0; i < lolita::core::ElementDescription<_element, _domain>::template getNumComponents<__i, __j>(); ++i) {
                            auto & rhs = this->template getComponents<__i, __j>()[i]->template getElement<_k>();
                            auto & lhs = this->template getElement<_k>()->template getComponents<__i, __j>()[i];
                            lhs = rhs;
                        }
                        if constexpr (__j < lolita::core::ElementDescription<_element, _domain>::template getNumComponents<__i>() - 1) {
                            initialize_components_imp.template operator()<__i, __j + 1u>(initialize_components_imp);
                        }
                        else if constexpr (__i < lolita::core::ElementDescription<_element, _domain>::getNumComponents() - 1) {
                            initialize_components_imp.template operator()<__i + 1u, 0u>(initialize_components_imp);
                        }
                    };
                    /*
                     *
                     */
                    auto initialize_neighbours = [&] <lolita::index __i = 0u, lolita::index __j = 0u> (
                            auto & initialize_neighbours_imp
                    )
                            mutable
                    {
                        for (int i = 0; i < this->template getNeighbours<__i, __j>().size(); ++i) {
                            auto & rhs = this->template getNeighbours<__i, __j>()[i]->template getElement<_k>();
                            auto & lhs = this->template getElement<_k>()->template getNeighbours<__i, __j>();
                            lhs.push_back(rhs);
                        }
                        if constexpr (__j < lolita::core::ElementDescription<_element, _domain>::template getNumNeighbours<__i>() - 1) {
                            initialize_neighbours_imp.template operator()<__i, __j + 1u>(initialize_neighbours_imp);
                        }
                        else if constexpr (__i < lolita::core::ElementDescription<_element, _domain>::getNumNeighbours() - 1) {
                            initialize_neighbours_imp.template operator()<__i + 1u, 0u>(initialize_neighbours_imp);
                        }
                    };
                    /*
                     *
                     */
                    initialize_components(initialize_components);
                    initialize_neighbours(initialize_neighbours);
                    this->template getElement<_k>()->tag_ = this->tag_;
                    this->template getElement<_k>()->domains_ = this->domains_;
                    this->template getElement<_k>()->initialize(mesh_data);
                    if constexpr (_k < std::tuple_size_v<DataPointer> - 1) {
                        initialize_element_imp.template operator()<_k + 1u>(initialize_element_imp);
                    }
                };
                /*
                 *
                 */
                initialize_element(initialize_element);
            }

            template<lolita::core::Element __element, auto... __args>
            void
            initialize(
                    lolita::core::mesh::Mesh<_domain, __args...> & mesh_data
            )
            requires(__element.isPoint())
            {
                /*
                  *
                  */
                auto initialize_element = [&] <lolita::index _k = 0u> (
                        auto & initialize_element_imp
                )
                        constexpr mutable
                {
                    /*
                     *
                     */
                    auto initialize_coordinates = [&] ()
                    {
                        this->template getElement<_k>()->components_ = this->components_;
                    };
                    /*
                     *
                     */
                    auto initialize_neighbours = [&] <lolita::index __i = 0u, lolita::index __j = 0u> (
                            auto & initialize_neighbours_imp
                    )
                            mutable
                    {
                        for (int i = 0; i < this->template getNeighbours<__i, __j>().size(); ++i) {
                            auto & rhs = this->template getNeighbours<__i, __j>()[i]->template getElement<_k>();
                            auto & lhs = this->template getElement<_k>()->template getNeighbours<__i, __j>();
                            lhs.push_back(rhs);
                        }
                        if constexpr (__j < lolita::core::ElementDescription<_element, _domain>::template getNumNeighbours<__i>() - 1) {
                            initialize_neighbours_imp.template operator()<__i, __j + 1u>(initialize_neighbours_imp);
                        }
                        else if constexpr (__i < lolita::core::ElementDescription<_element, _domain>::getNumNeighbours() - 1) {
                            initialize_neighbours_imp.template operator()<__i + 1u, 0u>(initialize_neighbours_imp);
                        }
                    };
                    /*
                     *
                     */
                    initialize_coordinates();
                    initialize_neighbours(initialize_neighbours);
                    this->template getElement<_k>()->tag_ = this->tag_;
                    this->template getElement<_k>()->domains_ = this->domains_;
                    this->template getElement<_k>()->initialize(mesh_data);
                    if constexpr (_k < std::tuple_size_v<DataPointer> - 1) {
                        initialize_element_imp.template operator()<_k + 1u>(initialize_element_imp);
                    }
                };
                /*
                 *
                 */
                initialize_element(initialize_element);
            }

            DataPointer elements_;

        };

    };

    template<auto _coupled_element, lolita::integer _n>
    struct ElementGroupPolicy
    {

        template<lolita::core::Element _element, lolita::geometry::Domain _domain, auto _mixed_element>
        using HHH = typename ElementGroupPolicy<_mixed_element, _n - 1>::template Object<_element, _domain>;

        template<lolita::core::Element _element, lolita::geometry::Domain _domain>
        using Data = typename std::remove_cvref_t<decltype(_coupled_element)>::template Type<HHH, _element, _domain>;

        template<lolita::core::Element _element, lolita::geometry::Domain _domain>
        struct Object : public ElementObjectLayerPolicy2<Object, Data> {};

    };

    template<auto _coupled_element, lolita::integer _n>
    requires(_n == 0)
    struct ElementGroupPolicy<_coupled_element, _n>
    {

        template<lolita::core::Element _element, lolita::geometry::Domain _domain, auto _finite_element>
        using HHH = lolita::core::finite_element::FEObject<_element, _domain, _finite_element>;

        template<lolita::core::Element _element, lolita::geometry::Domain _domain>
        using Data = typename std::remove_cvref_t<decltype(_coupled_element)>::template Type<HHH, _element, _domain>;

        template<lolita::core::Element _element, lolita::geometry::Domain _domain>
        struct Object : public ElementObjectLayerPolicy2<Object, Data>
        {};

    };

    template<lolita::core::Element _element, lolita::geometry::Domain _domain, auto _arg>
    using HHHHH = typename ElementGroupPolicy<_arg, 1>::template Object<_element, _domain>;

//    template<lolita::core::Element _element, lolita::geometry::Domain _domain, auto... _finite_element>
//    struct FiniteElementFinal : public ElementObjectLayerPolicy<
//            lolita::core::finite_element::FiniteElementFinal,
//            lolita::core::finite_element::FEObject>::template ElementObjectLayer<_element, _domain, _finite_element...>
//    {
//
//        using HHH = ElementObjectLayerPolicy<lolita::core::finite_element::FiniteElementFinal, lolita::core::finite_element::FEObject>;
//
//        using _Elements = typename HHH::template ElementObjectLayer<_element, _domain, _finite_element...>::Elements;
//
//        /**
//         * @brief
//         * @tparam __element
//         * @param initialization_data
//         * @param mesh_data
//         */
//        template<lolita::core::Element __element>
//        static
//        void
//        makeElement(
//                lolita::core::mesh::ElementInitializationData<_element> const & initialization_data,
//                lolita::core::mesh::Mesh<_domain, _finite_element...> & mesh_data
//        );
//
//        /**
//         * @brief
//         * @tparam ___element
//         * @param initialization_data
//         * @param mesh_data
//         */
//        template<lolita::core::Element ___element>
//        static
//        void
//        makeElement(
//                lolita::core::mesh::ElementInitializationData<_element> const & initialization_data,
//                lolita::core::mesh::Mesh<_domain, _finite_element...> & mesh_data
//        )
//        requires(!___element.isPoint())
//        {
//            /*
//             *
//             */
//            auto get_element_hash = [&] <lolita::core::Element __element> (
//                    std::array<lolita::index, __element.num_nodes_> node_tags
//            )
//            {
//                auto element_hash = std::basic_stringstream<lolita::character>();
//                if constexpr (node_tags.size() > 0) {
//                    std::sort(std::execution::par_unseq, node_tags.begin(), node_tags.end());
//                }
//                for (int i = 0; i < node_tags.size(); ++i) {
//                    element_hash << node_tags[i];
//                }
//                return element_hash.str();
//            };
//            /*
//             * make element
//             */
//            auto make_element = [&] <lolita::core::Element __element = _element, lolita::index _i = 0, lolita::index _j = 0> (
//                    std::shared_ptr<lolita::core::finite_element::FiniteElementFinal<__element, _domain, _finite_element...>> & ptr_element,
//                    lolita::core::mesh::ElementInitializationData<__element> const & element_initialization_data,
//                    auto & make_element_imp
//            )
//                    mutable
//            {
//                if constexpr (!__element.isPoint()) {
//                    using __ElementDescription = lolita::core::ElementDescription<__element, _domain>;
//                    using __ComponentDescription = lolita::core::ElementDescription<__ElementDescription::template getComponent<_i, _j>(), _domain>;
//                    using __MeshDescription = lolita::core::MeshDescription<_domain>;
//                    auto const constexpr _is_initialized = _i == 0 && _j == 0;
//                    auto const constexpr _component = __ElementDescription::template getComponent<_i, _j>();
//                    auto const constexpr _component_coordinates = __MeshDescription::template getElementCoordinates<_component>();
//                    auto const constexpr _neighbour_coordinates = __ComponentDescription::template getNeighbourCoordinates<__element>();
//                    auto const constexpr _element_coordinates = __MeshDescription::template getElementCoordinates<__element>();
//                    auto const constexpr _node_coordinates = __ElementDescription::template getComponentCoordinates<lolita::core::Element::Node()>();
//                    using __Component = lolita::core::finite_element::FiniteElementFinal<_component, _domain, _finite_element...>;
//                    using __Self = lolita::core::finite_element::FiniteElementFinal<__element, _domain, _finite_element...>;
//                    using __Elements = std::tuple<lolita::core::finite_element::FEObject<__element, _domain, _finite_element>...>;
//                    auto & components = mesh_data.elements_.template getElements<_component_coordinates.dim_, _component_coordinates.tag_>();
//                    auto & element_component_array = ptr_element->template getComponents<_i, _j>();
//                    for (auto i = 0; i < element_component_array.size(); ++i) {
//                        auto component_hash = std::basic_string<lolita::character>();
//                        if constexpr(!_component.isPoint()) {
//                            auto component_initialization_data = lolita::core::mesh::ElementInitializationData<_component>();
//                            for (int j = 0; j < _component.num_nodes_; ++j) {
//                                auto const k = __Self::template getComponentNodeConnection<_i, _j>(i, j);
//                                component_initialization_data.node_tags_[j] = element_initialization_data.node_tags_[k];
//                            }
//                            component_initialization_data.tag_ = components.size();
//                            component_hash = get_element_hash.template operator ()<_component>(component_initialization_data.node_tags_);
//                            if (!components.contains(component_hash)) {
//                                auto ptr_component = std::make_shared<__Component>(__Component());
//                                make_element_imp.template operator ()<_component, 0, 0>(ptr_component, component_initialization_data, make_element_imp);
//                            }
//                        }
//                        else {
//                            component_hash = std::to_string(element_initialization_data.node_tags_[__Self::template getComponentNodeConnection<_i, _j>(i, 0)]);
//                        }
//                        element_component_array[i] = components[component_hash];
//                        components[component_hash]->template getNeighbours<_neighbour_coordinates.dim_, _neighbour_coordinates.tag_>().push_back(ptr_element);
//                    }
//                    if constexpr (_j < __ElementDescription::template getNumComponents<_i>() - 1) {
//                        make_element_imp.template operator()<__element, _i, _j + 1u>(ptr_element, element_initialization_data, make_element_imp);
//                    }
//                    else if constexpr (_i < __ElementDescription::getNumComponents() - 1) {
//                        make_element_imp.template operator()<__element, _i + 1u, 0u>(ptr_element, element_initialization_data, make_element_imp);
//                    }
//                    if constexpr (_is_initialized) {
//                        auto & elements = mesh_data.elements_.template getElements<_element_coordinates.dim_, _element_coordinates.tag_>();
//                        auto const & nodes = ptr_element->template getComponents<_node_coordinates.dim_, _node_coordinates.tag_>();
//                        auto domains = std::unordered_set<std::shared_ptr<std::basic_string<lolita::character>>>();
//                        for (auto const & domain : nodes[0]->domains_) {
//                            auto has_domain = true;
//                            for (int j = 1; j < __element.num_nodes_; ++j) {
//                                has_domain = std::find(nodes[j]->domains_.begin(), nodes[j]->domains_.end(), domain) != nodes[j]->domains_.end();
//                            }
//                            if (has_domain) {
//                                domains.insert(domain);
//                            }
//                        }
//                        ptr_element->tag_ = element_initialization_data.tag_;
//                        ptr_element->domains_.assign(domains.begin(), domains.end());
//                        /*
//                         *
//                         */
//                        auto make_elements = [&] <lolita::index _k = 0> (
//                                auto & make_elements_imp
//                        )
//                        {
//                            using __Element = std::tuple_element_t<_k, __Elements>;
//                            ptr_element->template getElement<_k>() = std::make_shared<__Element>(__Element());
//                            if constexpr (_k < sizeof...(_finite_element) - 1) {
//                                make_elements_imp.template operator()<_k + 1u>(make_elements_imp);
//                            }
//                        };
//                        /*
//                         *
//                         */
//                        make_elements(make_elements);
//                        auto element_hash = get_element_hash.template operator ()<__element>(element_initialization_data.node_tags_);
//                        elements.insert({element_hash, ptr_element});
//                    }
//                }
//            };
//            /*
//             *
//             */
//            auto ptr_element = std::make_shared<FiniteElementFinal>(FiniteElementFinal());
//            make_element(ptr_element, initialization_data, make_element);
//        }
//
//        template<lolita::core::Element ___element>
//        static
//        void
//        makeElement(
//                lolita::core::mesh::ElementInitializationData<_element> const & initialization_data,
//                lolita::core::mesh::Mesh<_domain, _finite_element...> & mesh_data
//        )
//        requires(___element.isPoint())
//        {
//            /*
//             *
//             */
//            auto make_element = [&] (
//                    std::shared_ptr<lolita::core::finite_element::FiniteElementFinal<_element, _domain, _finite_element...>> & ptr_element
//            )
//                    mutable
//            {
//                auto const constexpr _element_coordinates = lolita::core::MeshDescription<_domain>::template getElementCoordinates<_element>();
//                ptr_element->tag_ = initialization_data.tag_;
//                ptr_element->domains_ = initialization_data.domains_;
//                ptr_element->components_ = std::make_shared<lolita::geometry::Point>(initialization_data.coordinates_);
//                /*
//                 *
//                 */
//                auto make_elements = [&] <lolita::index _k = 0> (
//                        auto & make_elements_imp
//                )
//                {
//                    ptr_element->template getElement<_k>() = std::make_shared<std::tuple_element_t<_k, _Elements>>(std::tuple_element_t<_k, _Elements>());
//                    if constexpr (_k < sizeof...(_finite_element) - 1) {
//                        make_elements_imp.template operator()<_k + 1u>(make_elements_imp);
//                    }
//                };
//                /*
//                 *
//                 */
//                make_elements(make_elements);
//                auto elem_hash = std::to_string(initialization_data.tag_);
//                mesh_data.elements_.template getElements<_element_coordinates.dim_, _element_coordinates.tag_>().insert({elem_hash, ptr_element});
//            };
//            /*
//             *
//             */
//            auto ptr_element = std::make_shared<FiniteElementFinal>(FiniteElementFinal());
//            make_element(ptr_element);
//        }
//
//    };

    /**
     * @brief
     * @tparam _element
     * @tparam _domain
     * @tparam _finite_element
     */
    template<lolita::core::Element _element, lolita::geometry::Domain _domain, auto... _finite_element>
    struct FiniteElementFinal
    :
    public lolita::core::finite_element::FiniteElementGeometry<lolita::core::finite_element::FiniteElementFinal, _element, _domain, _finite_element...>
    {

    private:

        /**
         * @brief
         */
        using _ElementsPtr = std::tuple<std::shared_ptr<lolita::core::finite_element::FEObject<_element, _domain, _finite_element>>...>;

        /**
         * @brief
         */
        using _Elements = std::tuple<lolita::core::finite_element::FEObject<_element, _domain, _finite_element>...>;

    public:

        /**
         * @brief
         * @tparam _i
         * @return
         */
        template<lolita::index _i>
        std::tuple_element_t<_i, _ElementsPtr> const &
        getElement()
        const
        {
            return std::get<_i>(elements_);
        }

        /**
         * @brief
         * @tparam _i
         * @return
         */
        template<lolita::index _i>
        std::tuple_element_t<_i, _ElementsPtr> &
        getElement()
        {
            return std::get<_i>(elements_);
        }

        template<lolita::core::Element __element>
        static
        void
        makeElement(
                lolita::core::mesh::ElementInitializationData<_element> const & initialization_data,
                lolita::core::mesh::Mesh<_domain, _finite_element...> & mesh_data
        );

        template<lolita::core::Element ___element>
        static
        void
        makeElement(
                lolita::core::mesh::ElementInitializationData<_element> const & initialization_data,
                lolita::core::mesh::Mesh<_domain, _finite_element...> & mesh_data
        )
        requires(!___element.isPoint())
        {
            /*
             *
             */
            auto get_element_hash = [&] <lolita::core::Element __element> (
                    std::array<lolita::index, __element.num_nodes_> node_tags
            )
            {
                auto element_hash = std::basic_stringstream<lolita::character>();
                if constexpr (node_tags.size() > 0) {
                    std::sort(std::execution::par_unseq, node_tags.begin(), node_tags.end());
                }
                for (int i = 0; i < node_tags.size(); ++i) {
                    element_hash << node_tags[i];
                }
                return element_hash.str();
            };
            /*
             * make element
             */
            auto make_element = [&] <lolita::core::Element __element = _element, lolita::index _i = 0, lolita::index _j = 0> (
                    std::shared_ptr<lolita::core::finite_element::FiniteElementFinal<__element, _domain, _finite_element...>> & ptr_element,
                    lolita::core::mesh::ElementInitializationData<__element> const & element_initialization_data,
                    auto & make_element_imp
            )
            mutable
            {
                if constexpr (!__element.isPoint()) {
                    using __ElementDescription = lolita::core::ElementDescription<__element, _domain>;
                    using __ComponentDescription = lolita::core::ElementDescription<__ElementDescription::template getComponent<_i, _j>(), _domain>;
                    using __MeshDescription = lolita::core::MeshDescription<_domain>;
                    auto const constexpr _is_initialized = _i == 0 && _j == 0;
                    auto const constexpr _component = __ElementDescription::template getComponent<_i, _j>();
                    auto const constexpr _component_coordinates = __MeshDescription::template getElementCoordinates<_component>();
                    auto const constexpr _neighbour_coordinates = __ComponentDescription::template getNeighbourCoordinates<__element>();
                    auto const constexpr _element_coordinates = __MeshDescription::template getElementCoordinates<__element>();
                    auto const constexpr _node_coordinates = __ElementDescription::template getComponentCoordinates<lolita::core::Element::Node()>();
                    using __Component = lolita::core::finite_element::FiniteElementFinal<_component, _domain, _finite_element...>;
                    using __Self = lolita::core::finite_element::FiniteElementFinal<__element, _domain, _finite_element...>;
                    using __Elements = std::tuple<lolita::core::finite_element::FEObject<__element, _domain, _finite_element>...>;
                    auto & components = mesh_data.elements_.template getElements<_component_coordinates.dim_, _component_coordinates.tag_>();
                    auto & element_component_array = ptr_element->template getComponents<_i, _j>();
                    for (auto i = 0; i < element_component_array.size(); ++i) {
                        auto component_hash = std::basic_string<lolita::character>();
                        if constexpr(!_component.isPoint()) {
                            auto component_initialization_data = lolita::core::mesh::ElementInitializationData<_component>();
                            for (int j = 0; j < _component.num_nodes_; ++j) {
                                auto const k = __Self::template getComponentNodeConnection<_i, _j>(i, j);
                                component_initialization_data.node_tags_[j] = element_initialization_data.node_tags_[k];
                            }
                            component_initialization_data.tag_ = components.size();
                            component_hash = get_element_hash.template operator ()<_component>(component_initialization_data.node_tags_);
                            if (!components.contains(component_hash)) {
                                auto ptr_component = std::make_shared<__Component>(__Component());
                                make_element_imp.template operator ()<_component, 0, 0>(ptr_component, component_initialization_data, make_element_imp);
                            }
                        }
                        else {
                            component_hash = std::to_string(element_initialization_data.node_tags_[__Self::template getComponentNodeConnection<_i, _j>(i, 0)]);
                        }
                        element_component_array[i] = components[component_hash];
                        components[component_hash]->template getNeighbours<_neighbour_coordinates.dim_, _neighbour_coordinates.tag_>().push_back(ptr_element);
                    }
                    if constexpr (_j < __ElementDescription::template getNumComponents<_i>() - 1) {
                        make_element_imp.template operator()<__element, _i, _j + 1u>(ptr_element, element_initialization_data, make_element_imp);
                    }
                    else if constexpr (_i < __ElementDescription::getNumComponents() - 1) {
                        make_element_imp.template operator()<__element, _i + 1u, 0u>(ptr_element, element_initialization_data, make_element_imp);
                    }
                    if constexpr (_is_initialized) {
                        auto & elements = mesh_data.elements_.template getElements<_element_coordinates.dim_, _element_coordinates.tag_>();
                        auto const & nodes = ptr_element->template getComponents<_node_coordinates.dim_, _node_coordinates.tag_>();
                        auto domains = std::unordered_set<std::shared_ptr<std::basic_string<lolita::character>>>();
                        for (auto const & domain : nodes[0]->domains_) {
                            auto has_domain = true;
                            for (int j = 1; j < __element.num_nodes_; ++j) {
                                has_domain = std::find(nodes[j]->domains_.begin(), nodes[j]->domains_.end(), domain) != nodes[j]->domains_.end();
                            }
                            if (has_domain) {
                                domains.insert(domain);
                            }
                        }
                        ptr_element->tag_ = element_initialization_data.tag_;
                        ptr_element->domains_.assign(domains.begin(), domains.end());
                        /*
                         *
                         */
                        auto make_elements = [&] <lolita::index _k = 0> (
                                auto & make_elements_imp
                        )
                        {
                            using __Element = std::tuple_element_t<_k, __Elements>;
                            ptr_element->template getElement<_k>() = std::make_shared<__Element>(__Element());
                            if constexpr (_k < sizeof...(_finite_element) - 1) {
                                make_elements_imp.template operator()<_k + 1u>(make_elements_imp);
                            }
                        };
                        /*
                         *
                         */
                        make_elements(make_elements);
                        auto element_hash = get_element_hash.template operator ()<__element>(element_initialization_data.node_tags_);
                        elements.insert({element_hash, ptr_element});
                    }
                }
            };
            /*
             *
             */
            auto ptr_element = std::make_shared<FiniteElementFinal>(FiniteElementFinal());
            make_element(ptr_element, initialization_data, make_element);
        }

        template<lolita::core::Element ___element>
        static
        void
        makeElement(
                lolita::core::mesh::ElementInitializationData<_element> const & initialization_data,
                lolita::core::mesh::Mesh<_domain, _finite_element...> & mesh_data
        )
        requires(___element.isPoint())
        {
            /*
             *
             */
            auto make_element = [&] (
                    std::shared_ptr<lolita::core::finite_element::FiniteElementFinal<_element, _domain, _finite_element...>> & ptr_element
            )
                    mutable
            {
                auto const constexpr _element_coordinates = lolita::core::MeshDescription<_domain>::template getElementCoordinates<_element>();
                ptr_element->tag_ = initialization_data.tag_;
                ptr_element->domains_ = initialization_data.domains_;
                ptr_element->components_ = std::make_shared<lolita::geometry::Point>(initialization_data.coordinates_);
                /*
                 *
                 */
                auto make_elements = [&] <lolita::index _k = 0> (
                        auto & make_elements_imp
                )
                {
                    ptr_element->template getElement<_k>() = std::make_shared<std::tuple_element_t<_k, _Elements>>(std::tuple_element_t<_k, _Elements>());
                    if constexpr (_k < sizeof...(_finite_element) - 1) {
                        make_elements_imp.template operator()<_k + 1u>(make_elements_imp);
                    }
                };
                /*
                 *
                 */
                make_elements(make_elements);
                auto elem_hash = std::to_string(initialization_data.tag_);
                mesh_data.elements_.template getElements<_element_coordinates.dim_, _element_coordinates.tag_>().insert({elem_hash, ptr_element});
            };
            /*
             *
             */
            auto ptr_element = std::make_shared<FiniteElementFinal>(FiniteElementFinal());
            make_element(ptr_element);
        }

        template<lolita::core::Element __element>
        void
        initialize(
                lolita::core::mesh::Mesh<_domain, _finite_element...> & mesh_data
        );

        template<lolita::core::Element __element>
        void
        initialize(
                lolita::core::mesh::Mesh<_domain, _finite_element...> & mesh_data
        )
        requires(!__element.isPoint())
        {
            /*
             *
             */
            auto initialize_element = [&] <lolita::index _k = 0u> (
                    auto & initialize_element_imp
            )
            mutable
            {
                /*
                 *
                 */
                auto initialize_components = [&] <lolita::index __i = 0u, lolita::index __j = 0u> (
                        auto & initialize_components_imp
                )
                        mutable
                {
                    for (int i = 0; i < lolita::core::ElementDescription<_element, _domain>::template getNumComponents<__i, __j>(); ++i) {
                        auto & rhs = this->template getComponents<__i, __j>()[i]->template getElement<_k>();
                        auto & lhs = this->template getElement<_k>()->template getComponents<__i, __j>()[i];
                        lhs = rhs;
                    }
                    if constexpr (__j < lolita::core::ElementDescription<_element, _domain>::template getNumComponents<__i>() - 1) {
                        initialize_components_imp.template operator()<__i, __j + 1u>(initialize_components_imp);
                    }
                    else if constexpr (__i < lolita::core::ElementDescription<_element, _domain>::getNumComponents() - 1) {
                        initialize_components_imp.template operator()<__i + 1u, 0u>(initialize_components_imp);
                    }
                };
                /*
                 *
                 */
                auto initialize_neighbours = [&] <lolita::index __i = 0u, lolita::index __j = 0u> (
                        auto & initialize_neighbours_imp
                )
                        mutable
                {
                    for (int i = 0; i < this->template getNeighbours<__i, __j>().size(); ++i) {
                        auto & rhs = this->template getNeighbours<__i, __j>()[i]->template getElement<_k>();
                        auto & lhs = this->template getElement<_k>()->template getNeighbours<__i, __j>();
                        lhs.push_back(rhs);
                    }
                    if constexpr (__j < lolita::core::ElementDescription<_element, _domain>::template getNumNeighbours<__i>() - 1) {
                        initialize_neighbours_imp.template operator()<__i, __j + 1u>(initialize_neighbours_imp);
                    }
                    else if constexpr (__i < lolita::core::ElementDescription<_element, _domain>::getNumNeighbours() - 1) {
                        initialize_neighbours_imp.template operator()<__i + 1u, 0u>(initialize_neighbours_imp);
                    }
                };
                /*
                 *
                 */
                initialize_components(initialize_components);
                initialize_neighbours(initialize_neighbours);
                this->template getElement<_k>()->tag_ = this->tag_;
                this->template getElement<_k>()->domains_ = this->domains_;
                this->template getElement<_k>()->initialize(mesh_data);
                if constexpr (_k < sizeof...(_finite_element) - 1) {
                    initialize_element_imp.template operator()<_k + 1u>(initialize_element_imp);
                }
            };
            /*
             *
             */
            initialize_element(initialize_element);
        }

        template<lolita::core::Element __element>
        void
        initialize(
                lolita::core::mesh::Mesh<_domain, _finite_element...> & mesh_data
        )
        requires(__element.isPoint())
        {
            /*
              *
              */
            auto initialize_element = [&] <lolita::index _k = 0u> (
                    auto & initialize_element_imp
            )
                    constexpr mutable
            {
                /*
                 *
                 */
                auto initialize_coordinates = [&] ()
                {
                    this->template getElement<_k>()->components_ = this->components_;
                };
                /*
                 *
                 */
                auto initialize_neighbours = [&] <lolita::index __i = 0u, lolita::index __j = 0u> (
                        auto & initialize_neighbours_imp
                )
                        mutable
                {
                    for (int i = 0; i < this->template getNeighbours<__i, __j>().size(); ++i) {
                        auto & rhs = this->template getNeighbours<__i, __j>()[i]->template getElement<_k>();
                        auto & lhs = this->template getElement<_k>()->template getNeighbours<__i, __j>();
                        lhs.push_back(rhs);
                    }
                    if constexpr (__j < lolita::core::ElementDescription<_element, _domain>::template getNumNeighbours<__i>() - 1) {
                        initialize_neighbours_imp.template operator()<__i, __j + 1u>(initialize_neighbours_imp);
                    }
                    else if constexpr (__i < lolita::core::ElementDescription<_element, _domain>::getNumNeighbours() - 1) {
                        initialize_neighbours_imp.template operator()<__i + 1u, 0u>(initialize_neighbours_imp);
                    }
                };
                /*
                 *
                 */
                initialize_coordinates();
                initialize_neighbours(initialize_neighbours);
                this->template getElement<_k>()->tag_ = this->tag_;
                this->template getElement<_k>()->domains_ = this->domains_;
                this->template getElement<_k>()->initialize(mesh_data);
                if constexpr (_k < sizeof...(_finite_element) - 1) {
                    initialize_element_imp.template operator()<_k + 1u>(initialize_element_imp);
                }
            };
            /*
             *
             */
            initialize_element(initialize_element);
        }

        _ElementsPtr elements_;

    };

//    template<lolita::core::Element _element, lolita::geometry::Domain _domain, auto... _finite_element>
//    struct MixedElementObject : public lolita::core::finite_element::FiniteElementGeometry<MixedElementObject, _element, _domain, _finite_element...>
//    {
//
//    private:
//
//        /**
//         * @brief
//         */
//        using _ElementsPtr = std::tuple<std::shared_ptr<lolita::core::finite_element::FEObject<_element, _domain, _finite_element>>...>;
//
//        /**
//         * @brief
//         */
//        using _Elements = std::tuple<lolita::core::finite_element::FEObject<_element, _domain, _finite_element>...>;
//
//    public:
//
//        /**
//         * @brief
//         * @tparam _i
//         * @return
//         */
//        template<lolita::index _i>
//        std::tuple_element_t<_i, _ElementsPtr> const &
//        getElement()
//        const
//        {
//            return std::get<_i>(elements_);
//        }
//
//        /**
//         * @brief
//         * @tparam _i
//         * @return
//         */
//        template<lolita::index _i>
//        std::tuple_element_t<_i, _ElementsPtr> &
//        getElement()
//        {
//            return std::get<_i>(elements_);
//        }
//
//        /**
//         * @brief
//         * @tparam __finite_element
//         * @param mesh
//         */
//        template<auto... __finite_element>
//        void
//        initialize(
//                lolita::core::mesh::Mesh<_domain, __finite_element...> & mesh
//        )
//        {
//            auto aaa = [&] <lolita::integer _i = 0> (auto & self) constexpr mutable {
//                getElement<_i>()->initialize(mesh);
//                if constexpr (_i < std::tuple_size_v<_ElementsPtr> - 1) {
//                    self.template operator()<_i + 1>(self);
//                }
//            };
//            aaa(aaa);
//        }
//
//        _ElementsPtr elements_;
//
//    };

//    template<auto _mixed_element>
//    struct MixedElementObjectPolicy
//    {
//
//        template<lolita::core::Element _element, lolita::geometry::Domain _domain, auto _finite_element>
//        using HHH = lolita::core::finite_element::FEObject<_element, _domain, _finite_element>;
//
//        template<lolita::core::Element __element, lolita::geometry::Domain __domain>
//        using Data = typename std::remove_cvref_t<decltype(_mixed_element)>::template Type<HHH, __element, __domain>;
//
//        template<lolita::core::Element __element, lolita::geometry::Domain __domain>
//        struct MixedElementObject : public ElementObjectLayerPolicy2<MixedElementObject, Data>
//        {};
//
//    };
//
//    template<auto _coupled_element>
//    struct CoupledElementObjectPolicy2
//    {
//
//        template<lolita::core::Element _element, lolita::geometry::Domain _domain, auto _mixed_element>
//        using HHH = typename MixedElementObjectPolicy<_mixed_element>::template MixedElementObject<_element, _domain>;
//
//        template<lolita::core::Element _element, lolita::geometry::Domain _domain>
//        using Data = typename std::remove_cvref_t<decltype(_coupled_element)>::template Type<HHH, _element, _domain>;
//
//        template<lolita::core::Element _element, lolita::geometry::Domain _domain>
//        struct CoupledElementObject : public ElementObjectLayerPolicy2<CoupledElementObject, Data>
//        {};
//
//    };
//
//    template<lolita::core::Element _element, lolita::geometry::Domain _domain, auto... _finite_elements>
//    struct MixedElementObject : public ElementObjectLayerPolicy<
//            lolita::core::finite_element::MixedElementObject,
//            lolita::core::finite_element::FEObject>::template ElementObjectLayer<_element, _domain, _finite_elements...>
//    {
//
//    };
//
//    template<lolita::core::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::MixedElementConcept auto... _mixed_elements>
//    struct CoupledElementObjectPolicy
//    {
//
//        using MixedElements = std::tuple<typename std::remove_cvref_t<decltype(_mixed_elements)>::template Type<MixedElementObject, _element, _domain>...>;
//
//    };
//
//    template<lolita::core::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::MixedElementConcept auto... _mixed_elements>
//    struct CoupledElementObject : public ElementObjectLayerPolicy<
//            lolita::core::finite_element::CoupledElementObject,
//            lolita::core::finite_element::MixedElementObject>::template ElementObjectLayer<_element, _domain, _mixed_elements...>
////    : public lolita::core::finite_element::FiniteElementGeometry<MixedElementObject, _element, _domain, _mixed_elements...>
//    {
//
//        using MixedElements = std::tuple<typename std::remove_cvref_t<decltype(_mixed_elements)>::template Type<MixedElementObject, _element, _domain>...>;
//
//    };

}

#endif //LOLITA_NEW_FINITE_ELEMENT_HXX
