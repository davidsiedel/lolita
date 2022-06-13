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

namespace lolita::core2::finite_element
{

    /**
     * @brief
     * @tparam t_element
     * @tparam t_domain
     * @tparam t_finite_element
     */
    template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
    struct FiniteElementFieldUnknowns;

    namespace unknown
    {

        /**
         * @brief
         */
        struct Unknown : public lolita::utility::Enumeration<Unknown>
        {

            /**
             * @brief
             * @param tag
             */
            constexpr
            Unknown(
                    std::basic_string_view<lolita::character> && tag
            )
            :
            lolita::utility::Enumeration<Unknown>(std::forward<std::basic_string_view<lolita::character>>(tag))
            {}

            /**
             * @brief
             * @return
             */
            static constexpr
            Unknown
            Structural()
            {
                return Unknown("Structural");
            }

            /**
             * @brief
             * @return
             */
            constexpr
            lolita::boolean
            isStructural()
            const
            {
                return * this == Structural();
            }

            /**
             * @brief
             * @return
             */
            static constexpr
            Unknown
            Subsidiary()
            {
                return Unknown("Subsidiary");
            }

            /**
             * @brief
             * @return
             */
            constexpr
            lolita::boolean
            isSubsidiary()
            const
            {
                return * this == Subsidiary();
            }

            /**
             * @brief
             * @return
             */
            static constexpr
            std::array<Unknown, 2>
            Unknowns()
            {
                return std::array<Unknown, 2>{
                        Unknown::Subsidiary(),
                        Unknown::Structural(),
                };
            }

        };

        /**
         * @brief
         * @tparam t_unknown
         * @tparam t_dim
         */
        template<lolita::core2::finite_element::unknown::Unknown t_unknown, lolita::integer t_dim>
        struct ScalarUnknown;

        /**
         * @brief
         * @tparam t_unknown
         * @tparam t_dim
         */
        template<lolita::core2::finite_element::unknown::Unknown t_unknown, lolita::integer t_dim>
        requires(t_unknown.isStructural())
        struct ScalarUnknown<t_unknown, t_dim>
        {

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

        /**
         * @brief
         * @tparam t_unknown
         * @tparam t_dim
         */
        template<lolita::core2::finite_element::unknown::Unknown t_unknown, lolita::integer t_dim>
        requires(t_unknown.isSubsidiary())
        struct ScalarUnknown<t_unknown, t_dim>
        {

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
             */
            std::unique_ptr<CoefficientVector> unknowns_coefficients_;

            /**
             * @brief
             */
            std::unique_ptr<CoefficientVector> bindings_coefficients_;

        };

        /**
         * @brief
         * @tparam t_domain
         * @tparam t_field
         * @tparam t_dim
         */
        template<
                lolita::domain::Domain t_domain,
                lolita::field::Field t_field,
                lolita::integer t_dim,
                lolita::core2::finite_element::unknown::Unknown t_unknown
        >
        struct FieldUnknown
        {

            /**
             * @brief The field type
             */
            using Field = lolita::core2::field::TensorPolicy<t_field, t_domain.dim_>;

            /**
             * @brief The unknown type object
             */
            lolita::core2::finite_element::unknown::Unknown const static constexpr unknown_ = t_unknown;

            /**
             * @brief The field object
             */
            lolita::field::Field const static constexpr field_ = t_field;

            /**
             * @brief The dimension or cardinality of the unknown
             */
            lolita::integer const static constexpr dim_ = t_dim;

            /**
             * @brief
             * @param row
             * @param col
             * @return
             */
            ScalarUnknown<t_unknown, t_dim> &
            getUnknownComponent(
                    lolita::integer row,
                    lolita::integer col
            )
            {
                return field_unknown_[col][row];
            }

            /**
             * @brief
             * @param row
             * @param col
             * @return
             */
            ScalarUnknown<t_unknown, t_dim> const &
            getUnknownComponent(
                    lolita::integer row,
                    lolita::integer col
            )
            const
            {
                return field_unknown_[col][row];
            }

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
            requires(t_dim > 0 && t_unknown.isSubsidiary())
            {
                using t_CoefficientVector = typename ScalarUnknown<t_unknown, t_dim>::CoefficientVector;
                getUnknownComponent(i, j).unknowns_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::Zero());
                //getUnknownComponent(i, j).unknowns_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::LinSpaced(0, 0 + t_dim - 1));
            }

            /**
             * @brief
             * @param i
             * @param j
             * @param dim
             */
            void
            setUnknownCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::index dim
            )
            requires(t_dim == -1 && t_unknown.isSubsidiary())
            {
                using t_CoefficientVector = typename ScalarUnknown<t_unknown, t_dim>::CoefficientVector;
                getUnknownComponent(i, j).unknowns_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::Zero(dim));
            }

            /**
             * @brief
             * @tparam t_finite_element
             * @tparam t_element_group
             * @param i
             * @param j
             * @param mesh
             */
            template<auto t_finite_element, auto t_element_group>
            void
            setUnknownCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
            )
            requires(t_dim > 0 && t_unknown.isStructural())
            {
                using t_CoefficientVector = typename ScalarUnknown<t_unknown, t_dim>::CoefficientVector;
                using t_CoordinateVector = typename ScalarUnknown<t_unknown, t_dim>::CoordinateVector;
                auto constexpr _finite_element_index = t_element_group.template getFiniteElementIndex<t_finite_element>();
                getUnknownComponent(i, j).unknowns_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::Zero());
                //getUnknownComponent(i, j).unknowns_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::LinSpaced(0, 0 + t_dim - 1));
                auto & k = mesh.systems_[_finite_element_index].num_unknowns_;
                getUnknownComponent(i, j).unknowns_coordinates_ = std::make_unique<t_CoordinateVector>(t_CoordinateVector::LinSpaced(k, k + t_dim - 1));
                k += t_dim;
            }

            /**
             * @brief
             * @tparam t_finite_element
             * @tparam t_element_group
             * @param i
             * @param j
             * @param dim
             * @param mesh
             */
            template<auto t_finite_element, auto t_element_group>
            void
            setUnknownCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::index dim,
                    lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
            )
            requires(t_dim == -1 && t_unknown.isStructural())
            {
                using t_CoefficientVector = typename ScalarUnknown<t_unknown, t_dim>::CoefficientVector;
                using t_CoordinateVector = typename ScalarUnknown<t_unknown, t_dim>::CoordinateVector;
                auto constexpr _finite_element_index = t_element_group.template getFiniteElementIndex<t_finite_element>();
                getUnknownComponent(i, j).unknowns_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::Zero(dim));
                auto & k = mesh.systems_[_finite_element_index].num_unknowns_;
                getUnknownComponent(i, j).unknowns_coordinates_ = std::make_unique<t_CoordinateVector>(t_CoordinateVector::LinSpaced(dim, k, k + dim - 1));
                k += dim;
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
            requires(t_dim > 0 && t_unknown.isSubsidiary())
            {
                using t_CoefficientVector = typename ScalarUnknown<t_unknown, t_dim>::CoefficientVector;
                getUnknownComponent(i, j).bindings_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::Zero());
            }

            /**
             * @brief
             * @param i
             * @param j
             * @param dim
             */
            void
            setBindingCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::index dim
            )
            requires(t_dim == -1 && t_unknown.isSubsidiary())
            {
                using t_CoefficientVector = typename ScalarUnknown<t_unknown, t_dim>::CoefficientVector;
                getUnknownComponent(i, j).bindings_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::Zero(dim));
            }

            /**
             * @brief
             * @tparam t_finite_element
             * @tparam t_element_group
             * @param i
             * @param j
             * @param mesh
             */
            template<auto t_finite_element, auto t_element_group>
            void
            setBindingCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
            )
            requires(t_dim > 0 && t_unknown.isStructural())
            {
                using t_CoefficientVector = typename ScalarUnknown<t_unknown, t_dim>::CoefficientVector;
                using t_CoordinateVector = typename ScalarUnknown<t_unknown, t_dim>::CoordinateVector;
                auto constexpr _finite_element_index = t_element_group.template getFiniteElementIndex<t_finite_element>();
                getUnknownComponent(i, j).bindings_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::Zero());
                auto & k = mesh.systems_[_finite_element_index].num_bindings_;
                getUnknownComponent(i, j).bindings_coordinates_ = std::make_unique<t_CoordinateVector>(t_CoordinateVector::LinSpaced(k, k + t_dim - 1));
                k += t_dim;
            }

            /**
             * @brief
             * @tparam t_finite_element
             * @tparam t_element_group
             * @param i
             * @param j
             * @param dim
             * @param mesh
             */
            template<auto t_finite_element, auto t_element_group>
            void
            setBindingCoefficients(
                    lolita::integer i,
                    lolita::integer j,
                    lolita::index dim,
                    lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
            )
            requires(t_dim == -1 && t_unknown.isStructural())
            {
                using t_CoefficientVector = typename ScalarUnknown<t_unknown, t_dim>::CoefficientVector;
                using t_CoordinateVector = typename ScalarUnknown<t_unknown, t_dim>::CoordinateVector;
                auto constexpr _finite_element_index = t_element_group.template getFiniteElementIndex<t_finite_element>();
                getUnknownComponent(i, j).bindings_coefficients_ = std::make_unique<t_CoefficientVector>(t_CoefficientVector::Zero(dim));
                auto & k = mesh.systems_[_finite_element_index].num_bindings_;
                getUnknownComponent(i, j).bindings_coordinates_ = std::make_unique<t_CoordinateVector>(t_CoordinateVector::LinSpaced(dim, k, k + dim - 1));
                k += dim;
            }

            /**
             * @brief
             */
            std::array<std::array<ScalarUnknown<t_unknown, t_dim>, Field::shape_.rows_>, Field::shape_.cols_> field_unknown_;

        };

        namespace detail
        {

            template<typename t_T>
            struct is_unknown : public std::false_type {};

            template<
                    lolita::domain::Domain t_domain,
                    lolita::field::Field t_field,
                    lolita::integer t_dim,
                    lolita::core2::finite_element::unknown::Unknown t_unknown
            >
            struct is_unknown<lolita::core2::finite_element::unknown::FieldUnknown<t_domain, t_field, t_dim, t_unknown>> : public std::true_type {};

        }

        /**
         * @brief
         * @tparam t_T
         */
        template<typename t_T>
        concept FieldUnknownConcept = lolita::core2::finite_element::unknown::detail::is_unknown<t_T>::value;

        /**
         * @brief
         * @tparam t_FieldUnknown
         */
        template<lolita::core2::finite_element::unknown::FieldUnknownConcept... t_FieldUnknown>
        struct FieldUnknownCollection
        {

            /**
             * @brief
             */
            using Unknowns = std::tuple<t_FieldUnknown...>;

            /**
             * @brief
             */
            std::array<lolita::integer, sizeof...(t_FieldUnknown)> const static constexpr dim_unknowns_ = {t_FieldUnknown::dim_...};

            /**
             * @brief
             */
            lolita::integer const static constexpr num_unknowns_ = sizeof...(t_FieldUnknown);

            /**
             * @brief
             * @return
             */
            static constexpr
            lolita::boolean
            hasStructuralUnknowns()
            {
                return ((t_FieldUnknown::unknown_.isStructural()) || ...);
            }

            /**
             * @brief
             * @tparam t_unknown
             * @return
             */
            template<lolita::core2::finite_element::unknown::Unknown... t_unknown>
            static constexpr
            lolita::integer
            getDimUnknowns()
            {
                auto constexpr unknowns = std::array<lolita::core2::finite_element::unknown::Unknown, sizeof...(t_unknown)>{t_unknown...};
                auto dim_unknowns = lolita::integer(0);
                auto set_dim_unknowns  = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (auto & self) constexpr mutable {
                    if constexpr (std::tuple_size_v<Unknowns> > 0) {
                        using t_Unknown = std::tuple_element_t<t_j, Unknowns>;
                        if constexpr (t_Unknown::unknown_ == unknowns[t_i]) {
                            dim_unknowns += t_Unknown::dim_;
                        }
                        if constexpr (t_j < std::tuple_size_v<Unknowns> - 1) {
                            self.template operator()<t_i, t_j + 1>(self);
                        }
                        else if constexpr (t_i < sizeof...(t_unknown) - 1) {
                            self.template operator()<t_i + 1, 0>(self);
                        }
                    }
                };
                set_dim_unknowns(set_dim_unknowns);
                return dim_unknowns;
            }

            /**
             * @brief
             * @tparam t_domain
             * @tparam t_unknown
             * @return
             */
            template<lolita::domain::Domain t_domain, lolita::core2::finite_element::unknown::Unknown... t_unknown>
            static constexpr
            lolita::integer
            getNumUnknowns()
            {
//                auto num_unknowns = lolita::integer(0);
//                auto set_num_unknowns  = [&] <lolita::integer t_i = 0> (auto & self) constexpr mutable {
//                    if constexpr (std::tuple_size_v<Unknowns> > 0) {
//                        using t_Unknown = std::tuple_element_t<t_i, Unknowns>;
//                        if constexpr (t_Unknown::unknown_ == t_unknown) {
//                            num_unknowns += t_Unknown::dim_ * t_Unknown::Field::shape_.size_;
//                        }
//                        if constexpr (t_i < std::tuple_size_v<Unknowns> - 1) {
//                            self.template operator()<t_i + 1>(self);
//                        }
//                    }
//                };
//                set_num_unknowns(set_num_unknowns);
//                return num_unknowns;
                auto constexpr unknowns = std::array<lolita::core2::finite_element::unknown::Unknown, sizeof...(t_unknown)>{t_unknown...};
                auto num_unknowns = lolita::integer(0);
                auto set_num_unknowns  = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (auto & self) constexpr mutable {
                    if constexpr (std::tuple_size_v<Unknowns> > 0) {
                        using t_Unknown = std::tuple_element_t<t_j, Unknowns>;
                        if constexpr (t_Unknown::unknown_ == unknowns[t_i]) {
                            num_unknowns += t_Unknown::dim_ * t_Unknown::Field::shape_.size_;
                        }
                        if constexpr (t_j < std::tuple_size_v<Unknowns> - 1) {
                            self.template operator()<t_i, t_j + 1>(self);
                        }
                        else if constexpr (t_i < sizeof...(t_unknown) - 1) {
                            self.template operator()<t_i + 1, 0>(self);
                        }
                    }
                };
                set_num_unknowns(set_num_unknowns);
                return num_unknowns;
            }

            /**
             * @brief
             * @tparam t_i
             * @return
             */
            template<lolita::integer t_i>
            std::tuple_element_t<t_i, Unknowns> const &
            getUnknown()
            const
            {
                return std::get<t_i>(unknowns_);
            }

            /**
             * @brief
             * @tparam t_i
             * @return
             */
            template<lolita::integer t_i>
            std::tuple_element_t<t_i, Unknowns> &
            getUnknown()
            {
                return std::get<t_i>(unknowns_);
            }

            /**
             * @brief
             */
            Unknowns unknowns_;

        };

        /**
         * @brief Default initialization is zero unknowns
         * @tparam t_element
         * @tparam t_domain
         * @tparam t_finite_element
         */
        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
        struct FiniteElementFieldUnknownsTraits
        {

            /**
             * @brief
             */
            using FieldUnknowns = lolita::core2::finite_element::unknown::FieldUnknownCollection<>;

            /**
             * @brief
             */
            struct Implementation : FiniteElementFieldUnknowns<t_element, t_domain, t_finite_element>
            {

                /**
                 * @brief
                 * @tparam t_element_group
                 * @param mesh
                 */
                template<auto t_element_group>
                void
                setUnknowns(
                        lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
                )
                {
                    std::cout << "----> setting unknown empty" << std::endl;
                }

            };

        };

        /**
         * @brief
         * @tparam t_element
         * @tparam t_domain
         * @tparam t_finite_element
         */
        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
        requires(t_finite_element.discretization_.isHHO() && t_element.isSub(t_domain, 0))
        struct FiniteElementFieldUnknownsTraits<t_element, t_domain, t_finite_element>
        {

        private:

            /**
             * @brief
             */
            using t_Field = typename lolita::core2::finite_element::FiniteElementTraits<t_element, t_domain, t_finite_element>::Field;

        public:

            /**
             * @brief
             */
            using Basis = lolita::core2::finite_element::basis::FiniteElementBasisTraits<
                    t_element,
                    lolita::core2::finite_element::basis::Basis::Monomial(),
                    t_finite_element.discretization_.ord_cell_
            >;

            /**
             * @brief
             */
            using FieldUnknowns = lolita::core2::finite_element::unknown::FieldUnknownCollection<
                    lolita::core2::finite_element::unknown::FieldUnknown<
                            t_domain,
                            t_finite_element.unknown_.tensor_,
                            Basis::dim_,
                            lolita::core2::finite_element::unknown::Unknown::Subsidiary()
                    >
            >;

            /**
             * @brief
             */
            struct Implementation : FiniteElementFieldUnknowns<t_element, t_domain, t_finite_element>
            {

                /**
                 * @brief
                 * @tparam t_element_group
                 * @param mesh
                 */
                template<auto t_element_group>
                void
                setUnknowns(
                        lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
                )
                {
                    for (int i = 0; i < t_Field::shape_.rows_; ++i) {
                        for (int j = 0; j < t_Field::shape_.cols_; ++j) {
                            this->field_unknowns_.template getUnknown<0>().setUnknownCoefficients(i, j);
                            if (this->getLoadComponent(i, j)->load_.isConstraint()) {
                                this->field_unknowns_.template getUnknown<0>().setBindingCoefficients(i, j);
                            }
                        }
                    }
                }

            };

        };

        /**
         * @brief
         * @tparam t_element
         * @tparam t_domain
         * @tparam t_finite_element
         */
        template<lolita::core2::geometry::Element t_element, lolita::domain::Domain t_domain, lolita::finite_element::FiniteElementConcept auto t_finite_element>
        requires(t_finite_element.discretization_.isHHO() && t_element.isSub(t_domain, 1))
        struct FiniteElementFieldUnknownsTraits<t_element, t_domain, t_finite_element>
        {

        private:

            /**
             * @brief
             */
            using t_Field = typename lolita::core2::finite_element::FiniteElementTraits<t_element, t_domain, t_finite_element>::Field;

        public:

            /**
             * @brief
             */
            using Basis = lolita::core2::finite_element::basis::FiniteElementBasisTraits<
                    t_element,
                    lolita::core2::finite_element::basis::Basis::Monomial(),
                    t_finite_element.discretization_.ord_face_
            >;

            /**
             * @brief
             */
            using FieldUnknowns = lolita::core2::finite_element::unknown::FieldUnknownCollection<
                    lolita::core2::finite_element::unknown::FieldUnknown<
                            t_domain,
                            t_finite_element.unknown_.tensor_,
                            Basis::dim_,
                            lolita::core2::finite_element::unknown::Unknown::Structural()
                    >
            >;

            /**
             * @brief
             */
            struct Implementation : FiniteElementFieldUnknowns<t_element, t_domain, t_finite_element>
            {

                /**
                 * @brief
                 * @tparam t_element_group
                 * @param mesh
                 */
                template<auto t_element_group>
                void
                setUnknowns(
                        lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
                )
                {
                    for (int i = 0; i < t_Field::shape_.rows_; ++i) {
                        for (int j = 0; j < t_Field::shape_.cols_; ++j) {
                            this->field_unknowns_.template getUnknown<0>().template setUnknownCoefficients<t_finite_element>(i, j, mesh);
                            if (this->getLoadComponent(i, j)->load_.isConstraint()) {
                                std::cout << "setting const " << std::endl;
                                this->field_unknowns_.template getUnknown<0>().template setBindingCoefficients<t_finite_element>(i, j, mesh);
                            }
                        }
                    }
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
    struct FiniteElementFieldUnknowns : virtual FiniteElementBase<t_element, t_domain, t_finite_element>
    {

        /**
         * @brief
         */
        using t_FiniteElementFieldUnknownsTraits = typename unknown::FiniteElementFieldUnknownsTraits<t_element, t_domain, t_finite_element>;

        /**
         * @brief
         */
        using t_FieldUnknowns = typename t_FiniteElementFieldUnknownsTraits::FieldUnknowns;

        /**
         * @brief
         */
        using t_FiniteElementTraits = FiniteElementTraits<t_element, t_domain, t_finite_element>;

        /**
         * @brief
         * @tparam t_unknown
         * @return
         */
        template<lolita::core2::finite_element::unknown::Unknown... t_unknown>
        static constexpr
        lolita::integer
        getNumUnknowns()
        {
            return t_FieldUnknowns::template getNumUnknowns<t_domain, t_unknown...>();
        }

        /**
         * @brief
         * @tparam t_unknown
         * @return
         */
        template<lolita::core2::finite_element::unknown::Unknown... t_unknown>
        static constexpr
        lolita::integer
        getDimUnknowns()
        {
            return t_FieldUnknowns::template getDimUnknowns<t_unknown...>();
        }

        /**
         * @brief
         * @tparam t_unknown
         * @return
         */
        template<lolita::core2::finite_element::unknown::Unknown t_unknown>
        lolita::matrix::Vector<lolita::real, t_FiniteElementTraits::template getNumUnknowns<t_unknown>()>
        getUnknowns()
        const
        {
            using t_UnknownVector = lolita::matrix::Vector<lolita::real, t_FiniteElementTraits::template getNumUnknowns<t_unknown>()>;
            auto unknowns = lolita::matrix::Zero<t_UnknownVector>();
            auto count = 0;
            auto set_element_unknowns = [&] <lolita::integer t_k = 0> (auto & self) constexpr mutable {
                using t_ElementUnknowns = typename unknown::FiniteElementFieldUnknownsTraits<t_element, t_domain, t_finite_element>::FieldUnknowns;
                if constexpr (std::tuple_size_v<typename t_ElementUnknowns::Unknowns> > 0) {
                    using t_ElementUnknown = std::tuple_element_t<t_k, typename t_ElementUnknowns::Unknowns>;
                    if constexpr (t_ElementUnknown::unknown_ == t_unknown) {
                        for (int i = 0; i < t_ElementUnknown::Field::shape_.rows_; ++i) {
                            for (int j = 0; j < t_ElementUnknown::Field::shape_.cols_; ++j) {
                                auto const & rhs = this->field_unknowns_.template getUnknown<t_k>().getUnknownComponent(i, j).getUnknownsCoefficients();
                                unknowns.template segment<t_ElementUnknown::dim_>(count) = rhs;
                                count += t_ElementUnknown::dim_;
                            }
                        }
                    }
                }
                if constexpr (t_k < lolita::integer(std::tuple_size_v<typename t_ElementUnknowns::Unknowns>) - 1) {
                    self.template operator()<t_k + 1>(self);
                }
            };
            set_element_unknowns(set_element_unknowns);
            auto set_neighbours_unknowns = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0, lolita::integer t_k = 0> (auto & self) constexpr mutable {
                //std::cout << t_element << " " << "t_i : " << t_i << " " << "t_j : " << t_j << " " << "t_k : " << t_k << std::endl;
                auto constexpr t_face = lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getComponent<t_i, t_j>();
                using t_FaceUnknowns = typename unknown::FiniteElementFieldUnknownsTraits<t_face, t_domain, t_finite_element>::FieldUnknowns;
                if constexpr (std::tuple_size_v<typename t_FaceUnknowns::Unknowns> > 0) {
                    using t_FaceUnknown = std::tuple_element_t<t_k, typename t_FaceUnknowns::Unknowns>;
                    if constexpr (t_FaceUnknown::unknown_ == t_unknown) {
                        auto const & faces = this->template getComponents<t_i, t_j>();
                        for (auto const & face : faces) {
                            for (int i = 0; i < t_FaceUnknown::Field::shape_.rows_; ++i) {
                                for (int j = 0; j < t_FaceUnknown::Field::shape_.cols_; ++j) {
                                    auto const & rhs = face->field_unknowns_.template getUnknown<t_k>().getUnknownComponent(i, j).getUnknownsCoefficients();
                                    unknowns.template segment<t_FaceUnknown::dim_>(count) = rhs;
                                    count += t_FaceUnknown::dim_;
                                }
                            }
                        }
                    }
                }
                if constexpr (t_k < lolita::integer(std::tuple_size_v<typename t_FaceUnknowns::Unknowns>) - 1) {
                    self.template operator()<t_i, t_j, t_k + 1>(self);
                }
                else if constexpr (t_j < lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getNumComponents<t_i>() - 1) {
                    self.template operator()<t_i, t_j + 1, 0>(self);
                }
                else if constexpr (t_i < lolita::core2::geometry::ElementTraits<t_element, t_domain>::template getNumComponents<>() - 1) {
                    self.template operator()<t_i + 1, 0, 0>(self);
                }
            };
            if constexpr (!t_element.isPoint()) {
                set_neighbours_unknowns(set_neighbours_unknowns);
            }
            //std::cout << t_element << " " << t_unknown << " unknowns :" << std::endl;
            //std::cout << unknowns << std::endl;
            return unknowns;
        }

        /**
         * @brief
         * @return
         */
        lolita::matrix::Vector<lolita::real, t_FiniteElementTraits::template getNumUnknowns()>
        getUnknowns()
        const
        {
            auto unknowns = lolita::matrix::Vector<lolita::real, t_FiniteElementTraits::template getNumUnknowns()>();
            auto count = lolita::integer(0);
            auto set_unknown_block = [&] <lolita::integer t_i = 0> (auto & self) mutable {
                auto constexpr t_unknown = lolita::core2::finite_element::unknown::Unknown::Unknowns()[t_i];
                auto constexpr num_unknowns = t_FiniteElementTraits::template getNumUnknowns<t_unknown>();
                unknowns.template segment<num_unknowns>(count) = getUnknowns<t_unknown>();
                count += num_unknowns;
                if constexpr (t_i < lolita::core2::finite_element::unknown::Unknown::Unknowns().size() - 1) {
                    self.template operator()<t_i + 1>(self);
                }
            };
            set_unknown_block(set_unknown_block);
            return unknowns;
        }

        /**
         * @brief
         * @tparam t_element_group
         * @param mesh
         */
        template<auto t_element_group>
        void
        setUnknowns(
                lolita::core2::mesh::Mesh<t_domain, t_element_group> & mesh
        )
        {
            static_cast<typename t_FiniteElementFieldUnknownsTraits::Implementation *>(this)->setUnknowns(mesh);
        }

        /**
         * @brief
         */
        t_FieldUnknowns field_unknowns_;

    };

}

#endif //LOLITA_LOLITA_CORE_5_003_UNKNOWN_HXX
