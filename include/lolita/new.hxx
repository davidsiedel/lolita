//
// Created by dsiedel on 28/05/22.
//

#ifndef LOLITA_NEW_HXX
#define LOLITA_NEW_HXX

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"

namespace lolita::core::field
{

    template<lolita::field::Tensor _tensor, lolita::index _dim_euclidean>
    struct TensorPolicy;

    template<lolita::field::Tensor _tensor, lolita::index _dim_euclidean>
    requires(_tensor.ord_ == 0)
    struct TensorPolicy<_tensor, _dim_euclidean>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{
                lolita::numerics::pow(_dim_euclidean, 0),
                lolita::numerics::pow(_dim_euclidean, 0)
        };

        lolita::index const static constexpr dim_ = shape_.size_;

    };

    template<lolita::field::Tensor _tensor, lolita::index _dim_euclidean>
    requires(_tensor.ord_ == 1)
    struct TensorPolicy<_tensor, _dim_euclidean>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{
                lolita::numerics::pow(_dim_euclidean, 0),
                lolita::numerics::pow(_dim_euclidean, 1)
        };

        lolita::index const static constexpr dim_ = shape_.size_;

    };

    template<lolita::field::Tensor _tensor, lolita::index _dim_euclidean>
    requires(_tensor.ord_ == 2)
    struct TensorPolicy<_tensor, _dim_euclidean>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{
                lolita::numerics::pow(_dim_euclidean, 1),
                lolita::numerics::pow(_dim_euclidean, 1)
        };

        lolita::index const static constexpr dim_ = shape_.size_;

    };

    template<lolita::field::Tensor _tensor, lolita::index _dim_euclidean>
    requires(_tensor.ord_ == 3)
    struct TensorPolicy<_tensor, _dim_euclidean>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{
                lolita::numerics::pow(_dim_euclidean, 1),
                lolita::numerics::pow(_dim_euclidean, 2)
        };

        lolita::index const static constexpr dim_ = shape_.size_;

    };

    template<lolita::field::Tensor _tensor, lolita::index _dim_euclidean>
    requires(_tensor.ord_ == 4)
    struct TensorPolicy<_tensor, _dim_euclidean>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{
                lolita::numerics::pow(_dim_euclidean, 2),
                lolita::numerics::pow(_dim_euclidean, 2)
        };

        lolita::index const static constexpr dim_ = shape_.size_;

    };

    struct MappingValues
    {

        lolita::index row_;

        lolita::index col_;

        lolita::index position_;

        lolita::real coefficient_;

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    struct MappingPolicy;

    /*
     * SCALAR
     */

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 0 && _mapping == lolita::field::Mapping::Identity)
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{1, 1};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
        };

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {}

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 0 && _domain.dim_ == 2 && _mapping == lolita::field::Mapping::Gradient)
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{2, 1};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{0, 1, 1, 1},
        };

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {}

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 0 && _domain.dim_ == 3 && _mapping == lolita::field::Mapping::Gradient)
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{3, 1};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{0, 1, 1, 1},
                lolita::core::field::MappingValues{0, 2, 2, 1},
        };

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {}

    };

    /*
     * VECTOR
     */

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 1 && _domain.dim_ == 2 && _mapping == lolita::field::Mapping::Identity)
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{2, 1};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{0, 1, 1, 1},
        };

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {}

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 1 && _domain.dim_ == 3 && _mapping == lolita::field::Mapping::Identity)
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{3, 1};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{0, 1, 1, 1},
                lolita::core::field::MappingValues{0, 2, 2, 1},
        };

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {}

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 1 && _domain.dim_ == 2 && _mapping == lolita::field::Mapping::Gradient)
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{2, 2};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{0, 1, 1, 1},
                lolita::core::field::MappingValues{1, 0, 2, 1},
                lolita::core::field::MappingValues{1, 1, 3, 1},
        };

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {}

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 1 && _domain.dim_ == 3 && _mapping == lolita::field::Mapping::Gradient)
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{3, 3};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{0, 1, 1, 1},
                lolita::core::field::MappingValues{0, 2, 2, 1},
                lolita::core::field::MappingValues{1, 0, 3, 1},
                lolita::core::field::MappingValues{1, 1, 4, 1},
                lolita::core::field::MappingValues{1, 2, 5, 1},
                lolita::core::field::MappingValues{2, 0, 6, 1},
                lolita::core::field::MappingValues{2, 1, 7, 1},
                lolita::core::field::MappingValues{2, 2, 8, 1},
        };

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {}

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 1 && _mapping == lolita::field::Mapping::SmallStrain)
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{3, 3};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{0, 1, 1, 1},
                lolita::core::field::MappingValues{0, 2, 2, 1},
                lolita::core::field::MappingValues{1, 0, 3, 1},
                lolita::core::field::MappingValues{1, 1, 4, 1},
                lolita::core::field::MappingValues{1, 2, 5, 1},
                lolita::core::field::MappingValues{2, 0, 6, 1},
                lolita::core::field::MappingValues{2, 1, 7, 1},
                lolita::core::field::MappingValues{2, 2, 8, 1},
        };

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {}

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 1 && _domain.dim_ == 2 && _mapping == lolita::field::Mapping::SmallStrainPlane)
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{1, 4};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{1, 1, 1, 1},
                lolita::core::field::MappingValues{0, 0, 2, 0},
                lolita::core::field::MappingValues{0, 1, 3, lolita::numerics::sqrt_2},
        };

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {}

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 1 && _domain.dim_ == 3 && _mapping == lolita::field::Mapping::SmallStrainSolid)
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{1, 6};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{1, 1, 1, 1},
                lolita::core::field::MappingValues{2, 2, 2, 0},
                lolita::core::field::MappingValues{0, 1, 3, lolita::numerics::sqrt_2},
                lolita::core::field::MappingValues{0, 2, 4, lolita::numerics::sqrt_2},
                lolita::core::field::MappingValues{1, 2, 5, lolita::numerics::sqrt_2},
        };

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {}

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 1 && _domain.dim_ == 2 && _mapping == lolita::field::Mapping::LargeStrainPlane)
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{1, 5};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{1, 1, 1, 1},
                lolita::core::field::MappingValues{0, 0, 2, 0},
                lolita::core::field::MappingValues{0, 1, 3, 1},
                lolita::core::field::MappingValues{1, 0, 4, 1},
        };

        static
        void
        non_linear(
//                lolita::matrix::VectorConcept2<shape_.size_> auto & gradient
                auto gradient
        )
        {
            gradient(0) += 1.0;
            gradient(1) += 1.0;
            gradient(2) += 1.0;
        }

    };

    template<lolita::field::Tensor _tensor, lolita::geometry::Domain _domain, lolita::field::Mapping _mapping>
    requires(_tensor.ord_ == 1 && _domain.dim_ == 3 && _mapping == lolita::field::Mapping::LargeStrainSolid)
    struct MappingPolicy<_tensor, _domain, _mapping>
    {

        lolita::matrix::Shape const static constexpr shape_ = lolita::matrix::Shape{1, 9};

        std::array<lolita::core::field::MappingValues, shape_.size_> const static constexpr values_ = {
                lolita::core::field::MappingValues{0, 0, 0, 1},
                lolita::core::field::MappingValues{1, 1, 1, 1},
                lolita::core::field::MappingValues{2, 2, 2, 0},
                lolita::core::field::MappingValues{0, 1, 3, 1},
                lolita::core::field::MappingValues{0, 2, 4, 1},
                lolita::core::field::MappingValues{1, 2, 5, 1},
                lolita::core::field::MappingValues{1, 0, 6, 1},
                lolita::core::field::MappingValues{2, 0, 7, 1},
                lolita::core::field::MappingValues{2, 1, 8, 1},
        };

        static
        void
        non_linear(
                lolita::matrix::Vector<lolita::real, shape_.size_> & gradient
        )
        {
            gradient(0) += 1.0;
            gradient(1) += 1.0;
            gradient(2) += 1.0;
        }

    };

}

namespace lolita::core
{

    /**
     * @brief Basic structure to define an element
     */
    struct Element
    {

        /**
         * @brief Basic equality comparison operator
         * @param other
         * @return
         */
        constexpr
        lolita::boolean
        operator==(
                Element const & other
        )
        const = default;

        /**
         * @brief Basic equality comparison operator
         * @param other
         * @return
         */
        constexpr
        lolita::boolean
        operator!=(
                Element const & other
        )
        const = default;

        /**
         * @brief
         * @param os
         * @param quadrature
         * @return
         */
        friend
        std::ostream &
        operator<<(
                std::ostream & os,
                Element const & element
        )
        {
            os << lolita::utility::readLabel(element.tag_);
            return os;
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::core::Element
        Node()
        {
            return Element{"Node", 0, 0, 1};
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::core::Element
        LinearSegment()
        {
            return Element{"LinearSegment", 1, 1, 2};
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::core::Element
        LinearTriangle()
        {
            return Element{"LinearTriangle", 2, 1, 3};
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::core::Element
        LinearQuadrangle()
        {
            return Element{"LinearQuadrangle", 2, 1, 4};
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::core::Element
        LinearTetrahedron()
        {
            return Element{"LinearTetrahedron", 3, 1, 4};
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::boolean
        isPoint()
        const
        {
            return dim_ == 0;
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::boolean
        isCurve()
        const
        {
            return dim_ == 1;
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::boolean
        isFacet()
        const
        {
            return dim_ == 2;
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::boolean
        isSolid()
        const
        {
            return dim_ == 3;
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::boolean
        isSegment()
        const
        {
            return * this == LinearSegment();
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::boolean
        isTriangle()
        const
        {
            return * this == LinearTriangle();
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::boolean
        isQuadrangle()
        const
        {
            return * this == LinearQuadrangle();
        }

        /**
         * @brief
         * @return
         */
        constexpr
        lolita::boolean
        isTetrahedron()
        const
        {
            return * this == LinearTetrahedron();
        }

        /**
         * @brief The element tag, that fully defines it. Two elements cannot have the same tag
         */
        lolita::utility::Label tag_;

        /**
         * @brief The euclidean element dimension
         */
        lolita::index dim_;

        /**
         * @brief The element polynomial order
         */
        lolita::index ord_;

        /**
         * @brief The element number of nodes
         */
        lolita::index num_nodes_;

    };

    namespace element
    {

        /**
         * @brief
         * @tparam _element
         * @tparam _domain
         */
        template<lolita::core::Element _element, lolita::geometry::Domain _domain>
        struct ElementGeometry;

        /**
         * @brief
         * @tparam _element
         * @tparam _quadrature
         * @tparam _ord
         */
        template<lolita::core::Element _element, lolita::finite_element::Quadrature _quadrature, lolita::index _ord>
        struct ElementQuadrature;

    }

    namespace finite_element
    {

        /**
         * @brief Forward declaration
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
        struct FiniteElementGeometry;

        /**
         * @brief
         * @tparam _element
         * @tparam _domain
         * @tparam _finite_element
         */
        template<lolita::core::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
        struct FEObject;

        /**
         * @brief Forward declaration
         * @tparam _T
         * @tparam _element
         * @tparam _domain
         * @tparam _finite_element
         */
        template<
                lolita::core::Element _element,
                lolita::geometry::Domain _domain,
                auto... _finite_element
        >
        struct FiniteElementFinal;

        namespace basis
        {

            /**
             * @brief Implementation object for the definition of a polynomial basis function
             * @tparam _element the element on which the monomial is defined
             * @tparam _basis the basis enum
             * @tparam _ord the polynomial order of the basis
             */
            template<lolita::core::Element _element, lolita::finite_element::Basis _basis, lolita::index _ord>
            struct FiniteElementBasis;

        }

        namespace unknown
        {

//            /**
//             * @brief
//             */
//            enum struct UnknownType
//            {
//
//                /**
//                 * @brief Is involved in the linear system construction
//                 */
//                Structural,
//
//                /**
//                 * @brief Is needed for the computation of physical quantities, but does not intervene in the global linear system e.g. local unknowns
//                 */
//                Subsidiary,
//
//            };

            /**
             * @brief
             */
            struct UnknownType
            {

                /**
                 * @brief Basic equality comparison operator
                 * @param other
                 * @return
                 */
                constexpr
                lolita::boolean
                operator==(
                        UnknownType const & other
                )
                const = default;

                /**
                 * @brief Basic equality comparison operator
                 * @param other
                 * @return
                 */
                constexpr
                lolita::boolean
                operator!=(
                        UnknownType const & other
                )
                const = default;

                /**
                 * @brief
                 * @return
                 */
                static constexpr
                UnknownType
                Structural()
                {
                    return UnknownType{"Structural"};
                }

                /**
                 * @brief
                 * @return
                 */
                static constexpr
                UnknownType
                Subsidiary()
                {
                    return UnknownType{"Subsidiary"};
                }

                /**
                 * @brief
                 * @return
                 */
                static constexpr
                std::array<UnknownType, 2>
                UnknownTypes()
                {
                    return std::array<UnknownType, 2>{
                        UnknownType{"Structural"},
                        UnknownType{"Subsidiary"}
                    };
                }

                /**
                 * @brief
                 */
                lolita::utility::Label tag_;

            };

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
            struct Unknown;

            namespace detail
            {

                template<typename _T>
                struct _UnknownTrait : public std::false_type {};

                template<
                        lolita::geometry::Domain _domain,
                        lolita::field::Tensor _field,
                        lolita::integer _dim,
                        lolita::core::finite_element::unknown::UnknownType _unknown_type
                >
                struct _UnknownTrait<lolita::core::finite_element::unknown::Unknown<_domain, _field, _dim, _unknown_type>> : public std::true_type {};

            }

            /**
             * @brief
             * @tparam _T
             */
            template<typename _T>
            concept UnknownConcept = lolita::core::finite_element::unknown::detail::_UnknownTrait<_T>::value;

            /**
             * @brief
             * @tparam _T
             */
            template<lolita::core::finite_element::unknown::UnknownConcept... _T>
            struct Unknowns;

            namespace detail
            {

                template<typename _T>
                struct UnknownsTrait : public std::false_type {};

                template<lolita::core::finite_element::unknown::UnknownConcept... _T>
                struct UnknownsTrait<lolita::core::finite_element::unknown::Unknowns<_T...>> : public std::true_type {};

            }

            /**
             * @brief
             * @tparam _T
             */
            template<typename _T>
            concept UnknownsConcept = lolita::core::finite_element::unknown::detail::UnknownsTrait<_T>::value;

            template<lolita::core::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
            struct FiniteElementUnknowns;

        }

    }

    namespace mesh
    {

        /**
         * @brief
         * @tparam _element
         */
        template<lolita::core::Element _element>
        struct ElementInitializationData;

        /**
         * @brief
         * @tparam _domain
         * @tparam _finite_element
         */
        template<lolita::geometry::Domain _domain, auto... _finite_element>
        struct Mesh;

    }

    /**
     * @brief Holder for the position of an element in a mesh, or with respect to some other element
     */
    struct ElementCoordinates
    {

        /**
         * @brief The relative dimension of the element in the mesh, or with respect to some other element
         */
        lolita::integer dim_;

        /**
         * @brief The relative position of the element in the mesh, or with respect to some other element
         */
        lolita::integer tag_;

    };

    namespace detail
    {

        /**
         * @brief
         * @tparam _element
         * @tparam _domain
         * @tparam _args
         */
        template<lolita::core::Element _element, lolita::geometry::Domain _domain, auto... _args>
        struct _Span
        {

            lolita::core::Element const static constexpr element_ = _element;

        };

        /**
         * @brief
         */
        template<template<lolita::core::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _args>
        using _Points = std::tuple<
                _T<lolita::core::Element::Node(), _domain, _args...>
        >;

        /**
         * @brief
         */
        template<template<lolita::core::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _args>
        using _Curves = std::tuple<
                _T<lolita::core::Element::LinearSegment(), _domain, _args...>
        >;

        /**
         * @brief
         */
        template<template<lolita::core::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _args>
        using _Facets = std::tuple<
                _T<lolita::core::Element::LinearTriangle(), _domain, _args...>,
                _T<lolita::core::Element::LinearQuadrangle(), _domain, _args...>
        >;

        /**
         * @brief
         */
        template<template<lolita::core::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _args>
        using _Solids = std::tuple<
                _T<lolita::core::Element::LinearTetrahedron(), _domain, _args...>
        >;

        /**
         * @brief
         */
        template<template<lolita::core::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _args>
        using _Elements = std::tuple<
                _Points<_T, _domain, _args...>,
                _Curves<_T, _domain, _args...>,
                _Facets<_T, _domain, _args...>,
                _Solids<_T, _domain, _args...>
        >;

    }

    /**
     * @brief
     */
    template<template<lolita::core::Element, lolita::geometry::Domain, auto...> typename _T, lolita::geometry::Domain _domain, auto... _args>
    using Elements = decltype(lolita::utility::tupleSlice<0, _domain.dim_ + 1>(std::declval<lolita::core::detail::_Elements<_T, _domain, _args...>>()));

    /**
     * @brief
     * @tparam _element
     * @tparam _domain
     */
    template<lolita::core::Element _element, lolita::geometry::Domain _domain>
    struct ElementDescription
    {

        /**
         * @brief
         */
        lolita::core::Element const static constexpr element_ = _element;

        /**
         * @brief
         */
        lolita::geometry::Domain const static constexpr domain_= _domain;

    private:

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::core::ElementCoordinates
        _coordinates()
        {
            using _Elements = lolita::core::detail::_Elements<lolita::core::detail::_Span, _domain>;
            auto position = lolita::core::ElementCoordinates{-1, -1};
            auto f0 = [&] <lolita::index _i = 0u, lolita::index _j = 0u> (
                    auto & self
            )
                    constexpr mutable
            {
                using _Element = typename std::tuple_element_t<_j, std::tuple_element_t<_i, _Elements>>;
                if (_Element::element_ == _element) {
                    position.dim_ = _i;
                    position.tag_ = _j;
                }
                if constexpr (_j < std::tuple_size_v<std::tuple_element_t<_i, _Elements>> - 1) {
                    self.template operator()<_i, _j + 1u>(self);
                }
                else if constexpr (_i < std::tuple_size_v<_Elements> - 1) {
                    self.template operator()<_i + 1u, 0u>(self);
                }
            };
            f0(f0);
            return position;
        }

    public:

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::boolean
        isCell()
        {
            return _domain.dim_ - _element.dim_ == 0;
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::boolean
        isFace()
        {
            return _domain.dim_ - _element.dim_ == 1;
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::core::ElementCoordinates
        getCoordinates()
        {
            using _Elements = lolita::core::detail::_Elements<lolita::core::detail::_Span, _domain>;
            auto position = lolita::core::ElementCoordinates{-1, -1};
            auto f0 = [&] <lolita::index _i = 0u, lolita::index _j = 0u> (
                    auto & self
            )
                    constexpr mutable
            {
                using _Element = typename std::tuple_element_t<_j, std::tuple_element_t<_i, _Elements>>;
                if (_Element::element_ == _element) {
                    position.dim_ = _i;
                    position.tag_ = _j;
                }
                if constexpr (_j < std::tuple_size_v<std::tuple_element_t<_i, _Elements>> - 1) {
                    self.template operator()<_i, _j + 1u>(self);
                }
                else if constexpr (_i < std::tuple_size_v<_Elements> - 1) {
                    self.template operator()<_i + 1u, 0u>(self);
                }
            };
            f0(f0);
            return position;
        }

        /**
         * @brief
         * @tparam _component
         * @return
         */
        template<lolita::core::Element _component>
        static constexpr
        lolita::core::ElementCoordinates
        getComponentCoordinates()
        requires(_element != lolita::core::Element::Node())
        {
            using _Components = typename lolita::core::element::ElementGeometry<_element, _domain>::template Components<lolita::core::detail::_Span, _domain>;
            auto position = lolita::core::ElementCoordinates{-1, -1};
            auto f0 = [&] <lolita::index _i = 0u, lolita::index _j = 0u> (
                    auto & self
            )
                    constexpr mutable
            {
                using _Component = typename std::tuple_element_t<_j, std::tuple_element_t<_i, _Components>>::value_type;
                if (_Component::element_ == _component) {
                    position.dim_ = _i;
                    position.tag_ = _j;
                }
                if constexpr (_j < std::tuple_size_v<std::tuple_element_t<_i, _Components>> - 1) {
                    self.template operator()<_i, _j + 1u>(self);
                }
                else if constexpr (_i < std::tuple_size_v<_Components> - 1) {
                    self.template operator()<_i + 1u, 0u>(self);
                }
            };
            f0(f0);
            return position;
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        static constexpr
        lolita::core::Element
        getComponent()
        requires(_element != lolita::core::Element::Node())
        {
            using _Components = typename lolita::core::element::ElementGeometry<_element, _domain>::template Components<lolita::core::detail::_Span, _domain>;
            return std::tuple_element_t<_j, std::tuple_element_t<_i, _Components>>::value_type::element_;
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        static constexpr
        lolita::core::ElementDescription<getComponent<_i, _j>(), _domain>
        getComponentDescription()
        requires(_element != lolita::core::Element::Node())
        {
            return lolita::core::ElementDescription<getComponent<_i, _j>(), _domain>();
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        static constexpr
        lolita::index
        getNumComponents()
        requires(_element != lolita::core::Element::Node())
        {
            using _Components = typename lolita::core::element::ElementGeometry<_element, _domain>::template Components<lolita::core::detail::_Span, _domain>;
            return std::tuple_size_v<std::tuple_element_t<_j, std::tuple_element_t<_i, _Components>>>;
        }

        /**
         * @brief
         * @tparam _i
         * @return
         */
        template<lolita::index _i>
        static constexpr
        lolita::index
        getNumComponents()
        requires(_element != lolita::core::Element::Node())
        {
            using _Components = typename lolita::core::element::ElementGeometry<_element, _domain>::template Components<lolita::core::detail::_Span, _domain>;
            return std::tuple_size_v<std::tuple_element_t<_i, _Components>>;
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::index
        getNumComponents()
        requires(_element != lolita::core::Element::Node())
        {
            using _Components = typename lolita::core::element::ElementGeometry<_element, _domain>::template Components<lolita::core::detail::_Span, _domain>;
            return std::tuple_size_v<_Components>;
        }

        /**
         * @brief
         * @tparam _neighbour
         * @return
         */
        template<lolita::core::Element _neighbour>
        static constexpr
        lolita::core::ElementCoordinates
        getNeighbourCoordinates()
        {
            using _Neighbours = typename lolita::core::element::ElementGeometry<_element, _domain>::template Neighbours<lolita::core::detail::_Span, _domain>;
            auto position = lolita::core::ElementCoordinates{-1, -1};
            auto f0 = [&] <lolita::index _i = 0u, lolita::index _j = 0u> (
                    auto & self
            )
                    constexpr mutable
            {
                using _Neighbour = typename std::tuple_element_t<_j, std::tuple_element_t<_i, _Neighbours>>::value_type;
                if (_Neighbour::element_ == _neighbour) {
                    position.dim_ = _i;
                    position.tag_ = _j;
                }
                if constexpr (_j < std::tuple_size_v<std::tuple_element_t<_i, _Neighbours>> - 1) {
                    self.template operator()<_i, _j + 1u>(self);
                }
                else if constexpr (_i < std::tuple_size_v<_Neighbours> - 1) {
                    self.template operator()<_i + 1u, 0u>(self);
                }
            };
            f0(f0);
            return position;
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        static constexpr
        lolita::core::Element
        getNeighbour()
        {
            using _Neighbours = typename lolita::core::element::ElementGeometry<_element, _domain>::template Neighbours<lolita::core::detail::_Span, _domain>;
            return std::tuple_element_t<_j, std::tuple_element_t<_i, _Neighbours>>::value_type::element_;
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        static constexpr
        lolita::core::ElementDescription<getNeighbour<_i, _j>(), _domain>
        getNeighbourDescription()
        {
            return lolita::core::ElementDescription<getNeighbour<_i, _j>(), _domain>();
        }

        /**
         * @brief
         * @tparam _i
         * @return
         */
        template<lolita::index _i>
        static constexpr
        lolita::index
        getNumNeighbours()
        {
            using _Neighbours = typename lolita::core::element::ElementGeometry<_element, _domain>::template Neighbours<lolita::core::detail::_Span, _domain>;
            return std::tuple_size_v<std::tuple_element_t<_i, _Neighbours>>;
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::index
        getNumNeighbours()
        {
            using _Neighbours = typename lolita::core::element::ElementGeometry<_element, _domain>::template Neighbours<lolita::core::detail::_Span, _domain>;
            return std::tuple_size_v<_Neighbours>;
        }

    };

    namespace detail
    {

        template<typename _T>
        struct _ElementDescriptionConcept : public std::false_type {};

        template<lolita::core::Element _element, lolita::geometry::Domain _domain>
        struct _ElementDescriptionConcept<lolita::core::ElementDescription<_element, _domain>> : public std::true_type {};

    }

    /**
     * @brief
     * @tparam _T
     */
    template<typename _T>
    concept ElementDescriptionConcept = lolita::core::detail::_ElementDescriptionConcept<_T>::value;

    /**
     * @brief
     * @tparam _finite_element
     */
    template<lolita::core::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
    struct FiniteElementDescription
    {

    private:

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::integer
        _ord_quadrature()
        {
            auto ord_quadrature = _finite_element.ord_quadrature_;
            return _domain.frame_ == lolita::geometry::Frame::AxiSymmetric ? 2 * ord_quadrature + 1 : 2 * ord_quadrature;
        }

    public:

        /**
         * @brief
         */
        lolita::finite_element::FiniteElementConcept auto const static constexpr finite_element_ = _finite_element;

        /**
         * @brief
         */
        using ElementDescription = lolita::core::ElementDescription<_element, _domain>;

        /**
         * @brief
         */
        using Quadrature = lolita::core::element::ElementQuadrature<_element, _finite_element.quadrature_, _ord_quadrature()>;

        /**
         * @brief
         * @tparam _method
         * @return
         */
        template<lolita::finite_element::Method _method>
        static constexpr
        lolita::boolean
        hasMethod()
        {
            return finite_element_.discretization_.method_ == _method;
        }

        /**
         * @brief
         * @tparam _unknown_type
         * @return
         */
        template<lolita::core::finite_element::unknown::UnknownType _unknown_type>
        static constexpr
        lolita::integer
        getNumUnknowns()
        {
            auto get_num_element_unknowns = [&] <lolita::core::Element __element> () constexpr mutable {
                using __Unknowns = typename lolita::core::finite_element::unknown::FiniteElementUnknowns<__element, _domain, _finite_element>::Unknowns;
                return __Unknowns::template getNumUnknowns<_domain, _unknown_type>();
            };
            auto num_element_unknowns = get_num_element_unknowns.template operator()<_element>();
            auto get_num_components_unknowns = [&] <lolita::integer _i = 0, lolita::integer _j = 0> (auto & self) constexpr mutable {
                auto const constexpr _component = ElementDescription::template getComponent<_i, _j>();
                auto num_component_unknowns = get_num_element_unknowns.template operator()<_component>();
                auto num_components = ElementDescription::template getNumComponents<_i, _j>();
                num_element_unknowns += num_component_unknowns * num_components;
                if constexpr (_j < ElementDescription::template getNumComponents<_i>() - 1) {
                    self.template operator()<_i, _j + 1>(self);
                }
                else if constexpr (_i < ElementDescription::getNumComponents() - 1) {
                    self.template operator()<_i + 1, 0>(self);
                }
            };
            get_num_components_unknowns(get_num_components_unknowns);
            return num_element_unknowns;
        }

        static constexpr
        lolita::integer
        getNumUnknowns()
        {
            auto num_element_unknowns = lolita::integer(0);
            auto set_num_element_unknowns = [&] <lolita::integer _i = 0> (auto & self) constexpr mutable {
                num_element_unknowns += getNumUnknowns<lolita::core::finite_element::unknown::UnknownType::UnknownTypes()[_i]>();
                if constexpr (_i < lolita::core::finite_element::unknown::UnknownType::UnknownTypes().size() - 1) {
                    self.template operator()<_i + 1>(self);
                }
            };
            set_num_element_unknowns(set_num_element_unknowns);
            return num_element_unknowns;
        }

        /**
         * @brief
         * @tparam _unknown_type
         * @return
         */
        template<lolita::core::finite_element::unknown::UnknownType _unknown_type>
        static constexpr
        lolita::integer
        getNumOwnedUnknowns()
        {
            auto get_num_element_unknowns = [&] <lolita::core::Element __element> () constexpr mutable {
                using __Unknowns = typename lolita::core::finite_element::unknown::FiniteElementUnknowns<__element, _domain, _finite_element>::Unknowns;
                return __Unknowns::template getNumUnknowns<_domain, _unknown_type>();
            };
            auto num_element_unknowns = get_num_element_unknowns.template operator()<_element>();
            return num_element_unknowns;
        }

        static constexpr
        lolita::integer
        getNumOwnedUnknowns()
        {
            auto num_element_unknowns = lolita::integer(0);
            auto set_num_element_unknowns = [&] <lolita::integer _i = 0> (auto & self) constexpr mutable {
                num_element_unknowns += getNumOwnedUnknowns<lolita::core::finite_element::unknown::UnknownType::UnknownTypes()[_i]>();
                if constexpr (_i < lolita::core::finite_element::unknown::UnknownType::UnknownTypes().size() - 1) {
                    self.template operator()<_i + 1>(self);
                }
            };
            set_num_element_unknowns(set_num_element_unknowns);
            return num_element_unknowns;
        }

        /**
         * @brief
         * @tparam _mapping
         * @return
         */
        template<lolita::field::Mapping _mapping>
        static constexpr
        lolita::integer
        getMappingSize()
        {
            return lolita::core::field::MappingPolicy<_finite_element.unkown_.tensor_, _domain, _mapping>::shape_.size_;
        }

        /**
         * @brief
         * @tparam _mapping
         * @return
         */
        template<lolita::field::Mapping _mapping>
        static constexpr
        lolita::matrix::VectorBlock
        getMappingBlock()
        {

            auto mapping_row = lolita::integer(0);
            for (auto m : _finite_element.unknown_.mappings_) {
                if (m == _mapping) {
                    return lolita::matrix::VectorBlock(mapping_row, mapping_row + getMappingSize<_mapping>());
                }
                mapping_row += getMappingSize<_mapping>();
            }
            return lolita::matrix::VectorBlock();
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::integer
        getOperatorNumRows()
        {
            auto mapping_size = lolita::integer(0);
            auto set_dim_mapping = [&] <lolita::integer _i = 0> (auto & self) constexpr mutable {
                mapping_size += getMappingSize<_finite_element.unknown_.mappings_[_i]>();
                if constexpr (_i < _finite_element.unknown_.mappings_.size() - 1) {
                    self.template operator()<_i + 1>(self);
                }
            };
            set_dim_mapping(set_dim_mapping);
            return mapping_size;
        }

    };

    namespace detail
    {

        template<typename _T>
        struct _FiniteElementDescriptionConcept : public std::false_type {};

        template<lolita::core::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto _finite_element>
        struct _FiniteElementDescriptionConcept<lolita::core::FiniteElementDescription<_element, _domain, _finite_element>> : public std::true_type {};

    }

    /**
     * @brief
     * @tparam _T
     */
    template<typename _T>
    concept FiniteElementDescriptionConcept = lolita::core::detail::_FiniteElementDescriptionConcept<_T>::value;

    /**
     * @brief
     * @tparam _FiniteElement
     */
    template<lolita::core::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto... _finite_element>
    struct MixedElementDescription
    {

    private:

        /**
         * @brief A simple alias
         */
        using _FiniteElementDescriptions = std::tuple<lolita::core::FiniteElementDescription<_element, _domain, _finite_element>...>;

    public:

        /**
         * @brief
         */
        template<
                template<lolita::core::FiniteElementDescriptionConcept auto> typename _T
        >
        using MixedElement = std::tuple<_T<lolita::core::FiniteElementDescription<_element, _domain, _finite_element>{}>...>;

        /**
         * @brief
         */
        _FiniteElementDescriptions const static constexpr finite_elements_ = {
                lolita::core::FiniteElementDescription<_element, _domain, _finite_element>{}...
        };

        /**
         * @brief
         * @tparam _i
         * @return
         */
        template<lolita::integer _i>
        static constexpr
        std::tuple_element_t<_i, _FiniteElementDescriptions>
        getFiniteElementDescription()
        {
            return std::get<_i>(finite_elements_);
        }

        /**
         * @brief Fetch the finite element index within the _finite_element list
         * @tparam __finite_element the finite element object to find
         * @return the finite_element index if found, and the size of the _finite_element list otherwise
         */
        template<lolita::finite_element::FiniteElementConcept auto __finite_element>
        static constexpr
        lolita::index
        getFiniteElementIndex()
        {
            auto index = sizeof...(_finite_element);
            using _Elements = std::tuple<std::remove_cvref_t<decltype(_finite_element)>...>;
            auto const constexpr elements = _Elements{_finite_element...};
            auto set_index = [&] <lolita::index _i = 0u> (auto & self)
                    constexpr mutable
            {
                if constexpr (std::is_same_v<std::tuple_element_t<_i, _Elements>, std::remove_cvref_t<decltype(__finite_element)>>) {
                    if (__finite_element == std::get<_i>(elements)) {
                        index = _i;
                    }
                }
                if constexpr(_i < sizeof...(_finite_element) - 1) {
                    self.template operator()<_i + 1u>(self);
                }
            };
            set_index(set_index);
            return index;
        }

    };

    namespace detail
    {

        template<typename _T>
        struct _MixedElementDescriptionConcept : public std::false_type {};

        template<lolita::core::Element _element, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto... _finite_element>
        struct _MixedElementDescriptionConcept<lolita::core::MixedElementDescription<_element, _domain, _finite_element...>> : public std::true_type {};

    }

    /**
     * @brief
     * @tparam _T
     */
    template<typename _T>
    concept MixedElementDescriptionConcept = lolita::core::detail::_MixedElementDescriptionConcept<_T>::value;

    /**
     * @brief
     * @tparam _domain
     */
    template<lolita::geometry::Domain _domain>
    struct MeshDescription
    {

        /**
         * @brief
         */
        lolita::geometry::Domain const static constexpr domain_ = _domain;

        /**
         * @brief
         * @tparam _element
         * @return
         */
        template<lolita::core::Element _element>
        static constexpr
        lolita::core::ElementCoordinates
        getElementCoordinates()
        {
            using _Elements = lolita::core::Elements<lolita::core::detail::_Span, _domain>;
            auto position = lolita::core::ElementCoordinates{_element.dim_, -1};
            auto f0 = [&] <lolita::index _i = 0u> (auto & self) constexpr mutable {
                using _Element = std::tuple_element_t<_i, std::tuple_element_t<_element.dim_, _Elements>>;
                if (_Element::element_ == _element) {
                    position.tag_ = _i;
                }
                if constexpr (_i < std::tuple_size_v<std::tuple_element_t<_element.dim_, _Elements>> - 1) {
                    self.template operator()<_i + 1u>(self);
                }
            };
            f0(f0);
            return position;
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        static constexpr
        lolita::core::Element
        getElement()
        {
            return std::tuple_element_t<_j, std::tuple_element_t<_i, lolita::core::Elements<lolita::core::detail::_Span, _domain>>>::element_;
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        static constexpr
        lolita::core::ElementDescription<getElement<_i, _j>(), _domain>
        getElementDescription()
        {
            return lolita::core::ElementDescription<getElement<_i, _j>(), _domain>();
        }

        /**
         * @brief
         * @tparam _i
         * @return
         */
        template<lolita::index _i>
        static constexpr
        lolita::index
        getNumElements()
        {
            return std::tuple_size_v<std::tuple_element_t<_i, lolita::core::Elements<lolita::core::detail::_Span, _domain>>>;
        }

        /**
         * @brief
         * @return
         */
        static constexpr
        lolita::index
        getNumElements()
        {
            return std::tuple_size_v<lolita::core::Elements<lolita::core::detail::_Span, _domain>>;
        }

    };

    namespace detail
    {

        template<typename _T>
        struct MeshDescriptionTrait : public std::false_type {};

        template<lolita::geometry::Domain _domain>
        struct MeshDescriptionTrait<lolita::core::MeshDescription<_domain>> : public std::true_type {};

    }

    /**
     * @brief
     * @tparam _T
     */
    template<typename _T>
    concept MeshDescriptionConcept = lolita::core::detail::MeshDescriptionTrait<_T>::value;

}

#endif //LOLITA_NEW_HXX
