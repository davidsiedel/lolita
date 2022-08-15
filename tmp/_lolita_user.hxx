//
// Created by dsiedel on 06/05/22.
//

#ifndef LOLITA_LOLITA_USER_HXX
#define LOLITA_LOLITA_USER_HXX

#include <MGIS/Behaviour/BehaviourData.hxx>
#include <MGIS/Behaviour/BehaviourData.h>
#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/Behaviour.h>
#include <MGIS/Behaviour/MaterialStateManager.h>
#include <MGIS/Behaviour/MaterialStateManager.hxx>
#include <MGIS/Behaviour/Integrate.hxx>
#include <MGIS/Behaviour/Integrate.h>

#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"

namespace lolita
{

    namespace new_one
    {

        struct Field
        {

            Integer int_;

        };

    }

    namespace domain
    {
        /**
         * @brief 
         * 
         */
        using Point = lolita::matrix::Vector<Real, 3>;

        /**
         * @brief 
         * 
         * @tparam t_Matrix 
         */
        template<lolita::matrix::MatrixConcept t_Matrix>
        static inline
        lolita::matrix::Vector<typename t_Matrix::Scalar, lolita::matrix::rows<t_Matrix>()>
        getBarycenter(
                t_Matrix const & point_args
        )
        requires(lolita::matrix::rows<t_Matrix>() > 0)
        {
            auto barycenter = lolita::matrix::Vector<typename t_Matrix::Scalar, lolita::matrix::rows<t_Matrix>()>().setZero();
            for (auto i = 0; i < lolita::matrix::cols<t_Matrix>(); ++i) {
                barycenter += point_args.col(i);
            }
            barycenter /= Real(lolita::matrix::cols<t_Matrix>());
            return barycenter;
        }

        /**
         * @brief 
         * 
         */
        struct Frame : public lolita::utility::Enumeration<Frame>
        {

            constexpr
            Frame(
                    std::basic_string_view<Character> && tag
            )
            :
            lolita::utility::Enumeration<Frame>(std::forward<std::basic_string_view<Character>>(tag))
            {}

            static constexpr
            Frame
            AxiSymmetric()
            {
                return Frame("AxiSymmetric");
            }

            constexpr
            Boolean
            isAxiSymmetric()
            const
            {
                return * this == AxiSymmetric();
            }

            static constexpr
            Frame
            Cartesian()
            {
                return Frame("Cartesian");
            }

            constexpr
            Boolean
            isCartesian()
            const
            {
                return * this == Cartesian();
            }

        };

        struct Domain : public lolita::utility::Enumeration<Domain>
        {

            constexpr
            Domain(
                    std::basic_string_view<Character> && tag,
                    Integer dim,
                    lolita::domain::Frame && frame
            )
            :
            lolita::utility::Enumeration<Domain>(std::forward<std::basic_string_view<Character>>(tag)),
            dim_(dim),
            frame_(std::forward<lolita::domain::Frame>(frame))
            {}

            constexpr
            Boolean
            hasDim(
                    Integer dim
            )
            const
            {
                return dim_ == dim;
            }

            constexpr
            Boolean
            hasFrame(
                    lolita::domain::Frame const & frame
            )
            const
            {
                return frame_ == frame;
            }

            constexpr
            Boolean
            hasFrame(
                    lolita::domain::Frame && frame
            )
            const
            {
                return frame_ == std::forward<lolita::domain::Frame>(frame);
            }

            /**
             * @brief
             * @return
             */
            constexpr
            Integer
            dim()
            const
            {
                return dim_;
            }

            /**
             * @brief
             * @return
             */
            constexpr
            lolita::domain::Frame const &
            frame()
            const
            {
                return frame_;
            }

            Integer dim_;

            lolita::domain::Frame frame_;

        };

    }

    namespace field
    {

        struct Mapping : public lolita::utility::Enumeration<Mapping>
        {

            constexpr
            Mapping(
                    std::basic_string_view<Character> && tag
            )
            :
            lolita::utility::Enumeration<Mapping>(std::forward<std::basic_string_view<Character>>(tag))
            {}

            static constexpr
            Mapping
            Gradient()
            {
                return Mapping("Gradient");
            }

            constexpr
            Boolean
            isGradient()
            const
            {
                return * this == Gradient();
            }

            static constexpr
            Mapping
            Identity()
            {
                return Mapping("Identity");
            }

            constexpr
            Boolean
            isIdentity()
            const
            {
                return * this == Identity();
            }

            static constexpr
            Mapping
            Divergence()
            {
                return Mapping("Divergence");
            }

            constexpr
            Boolean
            isDivergence()
            const
            {
                return * this == Divergence();
            }

            static constexpr
            Mapping
            LargeStrain()
            {
                return Mapping("LargeStrain");
            }

            constexpr
            Boolean
            isLargeStrain()
            const
            {
                return * this == LargeStrain();
            }

            static constexpr
            Mapping
            SmallStrain()
            {
                return Mapping("SmallStrain");
            }

            constexpr
            Boolean
            isSmallStrain()
            const
            {
                return * this == SmallStrain();
            }

        };

        namespace detail
        {

            template<typename T>
            struct MappingConcept : public std::false_type {};

            template<>
            struct MappingConcept<Mapping> : public std::true_type {};

        }

        template<typename T>
        concept MappingConcept = detail::MappingConcept<T>::value;

        struct Field : public lolita::utility::Enumeration<Field>
        {

            constexpr
            Field(
                    std::basic_string_view<Character> && tag,
                    Integer dim
            )
            :
            lolita::utility::Enumeration<Field>(std::forward<std::basic_string_view<Character>>(tag)),
            ord_(dim)
            {}

            constexpr
            Field
            as(
                    std::basic_string_view<Character> && tag
            )
            const
            {
                return Field(std::forward<std::basic_string_view<Character>>(tag), ord_);
            }

            static constexpr
            Field
            Scalar()
            {
                return Field("Scalar", 0);
            }

            static constexpr
            Field
            Vector()
            {
                return Field("Vector", 1);
            }

            static constexpr
            Field
            Tensor()
            {
                return Field("Tensor", 2);
            }

            static constexpr
            Field
            Displacement()
            {
                return Field("Displacement", 1);
            }

            constexpr
            Boolean
            isDisplacement()
            const
            {
                return * this == Displacement();
            }

            static constexpr
            Field
            Damage()
            {
                return Field("Damage", 0);
            }

            constexpr
            Boolean
            isDamage()
            const
            {
                return * this == Damage();
            }

            /**
             * @brief
             * @return
             */
            constexpr
            Boolean
            isScalar()
            const
            {
                return ord_ == 0;
            }

            /**
             * @brief
             * @return
             */
            constexpr
            Boolean
            isVector()
            const
            {
                return ord_ == 1;
            }

            /**
             * @brief
             * @return
             */
            constexpr
            Boolean
            isTensor(
                    Integer ord
            )
            const
            {
                return ord_ == ord;
            }

            Integer ord_;

        };

        template<lolita::field::MappingConcept... t_Mapping>
        requires(sizeof...(t_Mapping) > 0)
        struct Unknown
        {

            constexpr
            Unknown(
                    lolita::field::Field const & tensor,
                    t_Mapping &&... mappings
            )
            :
            tensor_(tensor),
            mappings_({std::forward<lolita::field::Mapping>(mappings)...})
            {}

            constexpr
            Unknown(
                    lolita::field::Field && tensor,
                    t_Mapping &&... mappings
            )
            :
            tensor_(std::forward<lolita::field::Field>(tensor)),
            mappings_({std::forward<lolita::field::Mapping>(mappings)...})
            {}

            constexpr
            Unknown(
                    std::basic_string_view<Character> && tag,
                    Integer dim,
                    t_Mapping &&... mappings
            )
            :
            tensor_(std::forward<lolita::field::Field>(Field(std::forward<std::basic_string_view<Character>>(tag), dim))),
            mappings_({std::forward<lolita::field::Mapping>(mappings)...})
            {}

            constexpr
            Boolean
            operator==(
                    Unknown const &
                    other
            )
            const = default;

            constexpr
            Boolean
            operator!=(
                    Unknown const &
                    other
            )
            const = default;

            lolita::field::Field tensor_;

            std::array<lolita::field::Mapping, sizeof...(t_Mapping)> mappings_;

        };

        namespace detail
        {

            template<typename T>
            struct IsUnknownConcept : public std::false_type {};

            template<typename... MappingT>
            struct IsUnknownConcept<Unknown<MappingT...>> : public std::true_type {};

        }

        template<typename T>
        concept UnknownConcept = detail::IsUnknownConcept<T>::value;

    }

    namespace behaviour
    {

        using ParameterFunction = std::function<Real(lolita::domain::Point const &, Real const &)>;

        struct MgisParameter
        {

            MgisParameter(
                    std::basic_string<Character> &&
                    parameter_tag,
                    lolita::behaviour::ParameterFunction &&
                    function
            )
            :
            parameter_tag_(std::forward<std::basic_string<Character>>(parameter_tag)),
            function_(std::forward<lolita::behaviour::ParameterFunction>(function))
            {}

            constexpr
            Boolean
            operator==(
                    MgisParameter const & other
            )
            const
            {
                return other.parameter_tag_ == this->parameter_tag_;
            }

            constexpr
            Boolean
            operator!=(
                    MgisParameter const & other
            )
            const
            {
                return !(* this == other);
            }

            std::basic_string<Character> parameter_tag_;

            lolita::behaviour::ParameterFunction function_;

        };

        struct MgisBehaviourData
        {

            mgis::behaviour::Behaviour behaviour_;

            std::vector<lolita::behaviour::MgisParameter> parameters_;

        };

        struct MgisBehaviour
        {

            MgisBehaviour(
                    std::basic_string<Character> &&
                    unknown_tag,
                    std::basic_string<Character> &&
                    domain_tag,
                    std::string const &
                    path,
                    std::string const &
                    name,
                    mgis::behaviour::Hypothesis
                    hypothesis,
                    std::vector<lolita::behaviour::MgisParameter> &&
                    parameters
            )
            :
            unknown_tag_(std::forward<std::basic_string<Character>>(unknown_tag)),
            domain_tag_(std::forward<std::basic_string<Character>>(domain_tag)),
            behaviour_data_(std::make_shared<lolita::behaviour::MgisBehaviourData>(lolita::behaviour::MgisBehaviourData{
                mgis::behaviour::load(path, name, hypothesis),
                std::forward<std::vector<lolita::behaviour::MgisParameter>>(parameters)
            }))
            {}

            MgisBehaviour(
                    std::basic_string<Character> &&
                    unknown_tag,
                    std::basic_string<Character> &&
                    domain_tag,
                    std::string const &
                    path,
                    std::string const &
                    name,
                    mgis::behaviour::Hypothesis
                    hypothesis,
                    mgis::behaviour::FiniteStrainBehaviourOptions
                    finite_strain_behaviour_options,
                    std::vector<lolita::behaviour::MgisParameter> &&
                    parameters
            )
            :
            unknown_tag_(std::forward<std::basic_string<Character>>(unknown_tag)),
            domain_tag_(std::forward<std::basic_string<Character>>(domain_tag)),
            behaviour_data_(std::make_shared<lolita::behaviour::MgisBehaviourData>(lolita::behaviour::MgisBehaviourData{
                mgis::behaviour::load(finite_strain_behaviour_options, path, name, hypothesis),
                std::forward<std::vector<lolita::behaviour::MgisParameter>>(parameters)
            }))
            {}

            std::basic_string<Character> unknown_tag_;

            std::basic_string<Character> domain_tag_;

            std::shared_ptr<lolita::behaviour::MgisBehaviourData> behaviour_data_;

        };

    }

    namespace finite_element
    {

        /**
         * @brief
         */
        struct Quadrature : public lolita::utility::Enumeration<Quadrature>
        {

            /**
             * @brief
             * @param tag
             */
            constexpr
            Quadrature(
                    std::basic_string_view<Character> && tag
            )
            :
            lolita::utility::Enumeration<Quadrature>(std::forward<std::basic_string_view<Character>>(tag))
            {}

            /**
             * @brief
             * @return
             */
            static constexpr
            Quadrature
            Gauss()
            {
                return Quadrature("Gauss");
            }

            constexpr
            Boolean
            isGauss()
            const
            {
                return * this == Gauss();
            }

        };

        /**
         * @brief
         */
        struct Load : public lolita::utility::Enumeration<Load>
        {

            /**
             * @brief
             * @param tag
             */
            constexpr
            Load(
                    std::basic_string_view<Character> && tag
            )
            :
            lolita::utility::Enumeration<Load>(std::forward<std::basic_string_view<Character>>(tag))
            {}

            /**
             * @brief
             * @return
             */
            static constexpr
            Load
            Natural()
            {
                return Load("Natural");
            }

            constexpr
            Boolean
            isNatural()
            const
            {
                return * this == Natural();
            }

            /**
             * @brief
             * @return
             */
            static constexpr
            Load
            Constraint()
            {
                return Load("Constraint");
            }

            constexpr
            Boolean
            isConstraint()
            const
            {
                return * this == Constraint();
            }

        };

        /**
         * @brief
         */
        using Loading = std::function<Real(lolita::domain::Point const &, Real const &)>;

        /**
         * @brief
         */
        struct ScalarLoad
        {

            /**
             * @brief
             */
            auto const static constexpr zero = [] (auto const &, auto const &) constexpr { return Real(0); };

            /**
             * @brief
             */
            ScalarLoad()
            :
            loading_(zero),
            load_(lolita::finite_element::Load::Natural())
            {}

            /**
             * @brief
             * @param function
             * @param load
             */
            ScalarLoad(
                    lolita::finite_element::Loading &&
                    function,
                    lolita::finite_element::Load
                    load
            )
            :
            loading_(std::forward<lolita::finite_element::Loading>(function)),
            load_(load)
            {}

            /**
             * @brief
             * @param function
             * @param load
             */
            ScalarLoad(
                    lolita::finite_element::Loading const &
                    function,
                    lolita::finite_element::Load
                    load
            )
            :
            loading_(function),
            load_(load)
            {}

            /**
             * @brief
             * @param point
             * @param time
             * @return
             */
            Real
            getImposedValue(
                    lolita::domain::Point const &
                    point,
                    Real const &
                    time
            )
            const
            {
                return loading_(point, time);
            }

            /**
             * @brief
             * @return
             */
            Boolean
            isConstrained()
            const
            {
                return load_.isConstraint();
            }

            /**
             * @brief
             */
            lolita::finite_element::Load load_;

            /**
             * @brief
             */
            lolita::finite_element::Loading loading_;

        };

        /**
         * @brief
         */
        struct LoadEntry
        {

            /**
             * @brief
             * @param unknown_tag
             * @param domain_tag
             * @param row
             * @param col
             * @param function
             * @param loading
             */
            LoadEntry(
                    std::basic_string_view<Character> && unknown_tag,
                    std::basic_string_view<Character> && domain_tag,
                    Integer element_dim,
                    Integer row,
                    Integer col,
                    lolita::finite_element::Loading && function,
                    lolita::finite_element::Load && loading
            )
            :
            unknown_tag_(std::forward<std::basic_string_view<Character>>(unknown_tag)),
            domain_tag_(std::forward<std::basic_string_view<Character>>(domain_tag)),
            element_dim_(element_dim),
            components_(lolita::matrix::Coordinates{row, col}),
            scalar_load_(std::make_shared<lolita::finite_element::ScalarLoad>(lolita::finite_element::ScalarLoad(
                    std::forward<lolita::finite_element::Loading>(function),
                    std::forward<lolita::finite_element::Load>(loading))))
            {}

            /**
             * @brief
             * @param unknown_tag
             * @param domain_tag
             * @param row
             * @param col
             * @param function
             * @param loading
             */
            LoadEntry(
                    std::basic_string_view<Character> && unknown_tag,
                    std::basic_string_view<Character> && domain_tag,
                    Integer element_dim,
                    Integer row,
                    Integer col,
                    lolita::finite_element::Loading const & function,
                    lolita::finite_element::Load && loading
            )
            :
            unknown_tag_(std::forward<std::basic_string_view<Character>>(unknown_tag)),
            domain_tag_(std::forward<std::basic_string_view<Character>>(domain_tag)),
            element_dim_(element_dim),
            components_(lolita::matrix::Coordinates{row, col}),
            scalar_load_(std::make_shared<lolita::finite_element::ScalarLoad>(lolita::finite_element::ScalarLoad(
                    function,
                    std::forward<lolita::finite_element::Load>(loading))))
            {}

            /**
             * @brief
             */
            std::basic_string_view<Character> unknown_tag_;

            /**
             * @brief
             */
            std::basic_string_view<Character> domain_tag_;

            /**
             * @brief
             */
            Integer element_dim_;

            /**
             * @brief
             */
            lolita::matrix::Coordinates components_;

            /**
             * @brief
             */
            std::shared_ptr<lolita::finite_element::ScalarLoad> scalar_load_;

        };

        

        struct Basis
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
            Basis(
                lolita::finite_element::Basis::Type type,
                Integer ord
            )
            :
            type_(type),
            ord_(ord)
            {}
            
            /**
             * @brief 
             * 
             * @param other 
             * @return constexpr Boolean 
             */
            constexpr
            Boolean
            operator==(
                    Basis const & other
            )
            const = default;
            
            /**
             * @brief 
             * 
             * @param other 
             * @return constexpr Boolean 
             */
            constexpr
            Boolean
            operator!=(
                    Basis const & other
            )
            const = default;

        public:

            /**
             * @brief 
             * 
             * @param ord 
             * @return constexpr lolita::finite_element::Basis
             */
            static constexpr
            lolita::finite_element::Basis
            monomial(
                Integer ord
            )
            {
                return Basis(lolita::finite_element::Basis::Type::Monomial, ord);
            }

            /**
             * @brief 
             * 
             * @return constexpr Boolean 
             */
            constexpr
            Boolean
            isMonomial()
            const
            {
                return type_ == lolita::finite_element::Basis::Type::Monomial;
            }

            /**
             * @brief 
             * 
             */
            lolita::finite_element::Basis::Type type_;

            /**
             * @brief 
             * 
             */
            Integer ord_;

        };
        
        /**
         * @brief
         */
        struct FiniteElementMethod : public lolita::utility::Enumeration<FiniteElementMethod>
        {

            /**
             * @brief
             * @param tag
             */
            constexpr
            FiniteElementMethod(
                    std::basic_string_view<Character> && tag
            )
            :
            lolita::utility::Enumeration<FiniteElementMethod>(std::forward<std::basic_string_view<Character>>(tag))
            {}

            /**
             * @brief
             * @return
             */
            static constexpr
            FiniteElementMethod
            HHO()
            {
                return FiniteElementMethod("HHO");
            }

            constexpr
            Boolean
            isHHO()
            const
            {
                return * this == HHO();
            }

        };

        /**
         * @brief Implementation of the Hybrid Discontinuous Galerkin method
         */
        struct HybridHighOrder : public lolita::finite_element::FiniteElementMethod
        {

            /**
             * @brief
             */
            struct Stabilization : public lolita::utility::Enumeration<Stabilization>
            {

                /**
                 * @brief
                 * @param tag
                 */
                constexpr
                Stabilization(
                        std::basic_string_view<Character> && tag
                )
                :
                lolita::utility::Enumeration<Stabilization>(std::forward<std::basic_string_view<Character>>(tag))
                {}

                /**
                 * @brief
                 * @return
                 */
                static constexpr
                Stabilization
                HHO()
                {
                    return Stabilization("HHO");
                }

                constexpr
                Boolean
                isHHO()
                const
                {
                    return * this == HHO();
                }

                /**
                 * @brief
                 * @return
                 */
                static constexpr
                Stabilization
                HDG()
                {
                    return Stabilization("HDG");
                }

                constexpr
                Boolean
                isHDG()
                const
                {
                    return * this == HDG();
                }

            };

            /**
             * @brief
             */
            constexpr
            HybridHighOrder()
            :
            lolita::finite_element::FiniteElementMethod(lolita::finite_element::FiniteElementMethod::HHO()),
            ord_cell_(-1),
            ord_face_(-1)
            {}

            /**
             * @brief
             * @param ord_cell
             * @param ord_face
             */
            constexpr
            HybridHighOrder(
                    Integer ord_cell,
                    Integer ord_face
            )
            :
            lolita::finite_element::FiniteElementMethod(lolita::finite_element::FiniteElementMethod::HHO()),
            ord_cell_(ord_cell),
            ord_face_(ord_face)
            {}

            /**
             * @brief
             * @param mapping
             * @return
             */
            constexpr
            Integer
            ordMapping(
                    lolita::field::Mapping mapping
            )
            const
            {
                return mapping == lolita::field::Mapping::Identity() ? ord_cell_ : ord_face_;
            }

            /**
             * @brief
             * @return
             */
            constexpr
            Integer
            ordCell()
            const
            {
                return ord_cell_;
            }

            /**
             * @brief
             * @return
             */
            constexpr
            Integer
            ordFace()
            const
            {
                return ord_face_;
            }

            /**
             * @brief
             */
            Integer ord_cell_;

            /**
             * @brief
             */
            Integer ord_face_;

        };

        template<typename _T>
        concept FiniteElementMethodConcept = std::derived_from<_T, lolita::finite_element::FiniteElementMethod> && requires(
                std::remove_reference_t<_T> const arg,
                lolita::field::Mapping mapping
        )
        {
            { arg.ord_cell_ } -> std::convertible_to<Integer>;
            { arg.ordMapping(mapping) } -> std::convertible_to<Integer>;
        };

        template<lolita::field::UnknownConcept auto _unknown, auto _behaviour, lolita::finite_element::FiniteElementMethodConcept auto _finite_element_method>
        struct FiniteElement
        {

            using UnknownType = std::remove_cvref_t<decltype(_unknown)>;

            using BehaviourType = std::remove_cvref_t<decltype(_behaviour)>;

            using FiniteElementMethodType = std::remove_cvref_t<decltype(_finite_element_method)>;

            constexpr
            FiniteElement(
                    lolita::finite_element::Quadrature quadrature,
                    Integer ord_quadrature
            )
            :
            unknown_(_unknown),
            behaviour_(_behaviour),
            discretization_(_finite_element_method),
            quadrature_(quadrature),
            ord_quadrature_(ord_quadrature)
            {}

            constexpr
            FiniteElement(
                    lolita::finite_element::Quadrature quadrature
            )
            :
            unknown_(_unknown),
            behaviour_(_behaviour),
            discretization_(_finite_element_method),
            quadrature_(quadrature),
            ord_quadrature_()
            {}

            constexpr
            Boolean
            operator==(
                    FiniteElement const & other
            )
            const = default;

            constexpr
            Boolean
            operator!=(
                    FiniteElement const & other
            )
            const = default;

            constexpr
            UnknownType const &
            unknown()
            const
            {
                return unknown_;
            }

            constexpr
            FiniteElementMethodType const &
            discretization()
            const
            {
                return discretization_;
            }

            constexpr
            lolita::finite_element::Quadrature const &
            quadrature()
            const
            {
                return quadrature_;
            }

            constexpr
            Integer
            ordQuadrature()
            const
            {
                return ord_quadrature_;
            }

            constexpr
            lolita::field::Field const &
            field()
            const
            {
                return unknown_.field();
            }

            constexpr
            Boolean
            hasMethod(
                    lolita::finite_element::FiniteElementMethod method
            )
            const
            {
                return discretization_.tag_ == method.tag_;
            }

            UnknownType unknown_;

            BehaviourType behaviour_;

            FiniteElementMethodType discretization_;

            lolita::finite_element::Quadrature quadrature_;

            Integer ord_quadrature_;

        };

        namespace detail
        {

            template<typename T>
            struct IsFiniteElementConcept : std::false_type {};

            template<auto _unknown, auto _behaviour, auto _finite_element_method>
            struct IsFiniteElementConcept<FiniteElement<_unknown, _behaviour, _finite_element_method>> : std::true_type {};

        }

        template<typename _T>
        concept FiniteElementConcept = detail::IsFiniteElementConcept<_T>::value;

        template<auto _arg>
        concept HybridHighOrderFiniteElementConcept = FiniteElementConcept<std::remove_cvref_t<decltype(_arg)>> && requires(
                std::remove_reference_t<decltype(_arg)> const arg
        )
        {
            { arg.discretization_.tag_ == lolita::finite_element::FiniteElementMethod::HHO().tag_ };
        };

        template<auto... _finite_elements>
        struct ElementGroup;

        namespace detail
        {

            template<typename T>
            struct IsMixedElementConcept : std::false_type {};

            template<auto... _finite_elements>
            struct IsMixedElementConcept<ElementGroup<_finite_elements...>> : std::true_type {};

        }

        template<typename _T>
        concept ElementGroupConcept = lolita::finite_element::detail::IsMixedElementConcept<_T>::value;

        template<auto... _finite_elements>
        struct ElementGroup
        {

        private:

            /**
             * @brief A simple alias
             */
            using t_FiniteElements = std::tuple<std::remove_cvref_t<decltype(_finite_elements)>...>;

            /**
             * @brief A simple alias
             */
            using _Self = lolita::finite_element::ElementGroup<_finite_elements...>;

        public:

            /**
             * @brief
             */
            template<template<auto, auto, auto> typename _T, auto _element, auto _domain>
            using ElementPointers = std::tuple<std::shared_ptr<_T<_element, _domain, _finite_elements>>...>;

            /**
             * @brief
             */
            template<template<auto, auto, auto> typename _T, auto _element, auto _domain>
            using Elements = std::tuple<_T<_element, _domain, _finite_elements>...>;

            /**
             * @brief
             */
            t_FiniteElements const static constexpr finite_elements_ = {_finite_elements...};

            static constexpr
            Boolean
            has()
            {
                return (lolita::finite_element::FiniteElementConcept<std::remove_cvref_t<decltype(_finite_elements)>> && ...);
            }

            static constexpr
            Integer
            size()
            {
                return sizeof...(_finite_elements);
            }

            /**
             * @brief
             * @tparam _i
             * @return
             */
            template<Integer _i>
            static constexpr
            std::tuple_element_t<_i, t_FiniteElements>
            getFiniteElement()
            {
                return std::get<_i>(finite_elements_);
            }

            /**
             * @brief
             * @return
             */
            static constexpr
            Integer
            count()
            {
                if constexpr (has()) {
                    return size();
                }
                else {
                    return lolita::numerics::sum(_finite_elements.count()...);
                }
            }

            template<auto __finite_element>
            static constexpr
            Boolean
            hasFiniteElement()
            {
                if constexpr(has()) {
                    auto index = false;
                    using _Elements = std::tuple<std::remove_cvref_t<decltype(_finite_elements)>...>;
                    auto const constexpr elements = _Elements{_finite_elements...};
                    auto set_index = [&] <Integer _i = 0u> (auto & self) constexpr mutable {
                        if constexpr (std::is_same_v<std::tuple_element_t<_i, _Elements>, std::remove_cvref_t<decltype(__finite_element)>>) {
                            if (__finite_element == std::get<_i>(elements)) {
                                index = true;
                            }
                        }
                        if constexpr(_i < sizeof...(_finite_elements) - 1) {
                            self.template operator()<_i + 1u>(self);
                        }
                    };
                    set_index(set_index);
                    return index;
                }
                else {
                    return (_finite_elements.template hasFiniteElement<__finite_element>() || ...);
                }
            }

            template<auto __finite_element>
            static constexpr
            void
            getFiniteElementIndex2(
                    Integer & index,
                    Boolean & found
            )
            {
                if constexpr(has()) {
                    if constexpr (hasFiniteElement<__finite_element>()) {
                        auto index2 = 0;
                        using _Elements = std::tuple<std::remove_cvref_t<decltype(_finite_elements)>...>;
                        auto const constexpr elements = _Elements{_finite_elements...};
                        auto set_index = [&] <Integer _i = 0u> (auto & self) constexpr mutable {
                            if constexpr (std::is_same_v<std::tuple_element_t<_i, _Elements>, std::remove_cvref_t<decltype(__finite_element)>>) {
                                if (__finite_element == std::get<_i>(elements)) {
                                    index2 = _i;
                                }
                            }
                            if constexpr(_i < sizeof...(_finite_elements) - 1) {
                                self.template operator()<_i + 1u>(self);
                            }
                        };
                        set_index(set_index);
                        index += index2;
                        found = true;
                    }
                    else {
                        if (!found) {
                            index += size();
                        }
                    };
                }
                else {
                    (_finite_elements.template getFiniteElementIndex2<__finite_element>(index, found), ...);
                }
            }

            template<auto __finite_element>
            static constexpr
            Integer
            getFiniteElementIndex()
            {
                auto index = 0;
                auto found = false;
//                auto mli = [&] <typename _T, Integer __i = 0> (auto & self1) constexpr mutable {
//                    if constexpr(_T::has()) {
//                        if constexpr (_T::template hasFiniteElement<__finite_element>()) {
//                            auto index2 = 0;
//                            using _Elements = std::tuple<std::remove_cvref_t<decltype(_finite_elements)>...>;
//                            auto const constexpr elements = _Elements{_finite_elements...};
//                            auto set_index = [&] <Integer _i = 0u> (auto & self) constexpr mutable {
//                                if constexpr (std::is_same_v<std::tuple_element_t<_i, _Elements>, std::remove_cvref_t<decltype(__finite_element)>>) {
//                                    if (__finite_element == std::get<_i>(elements)) {
//                                        index2 = _i;
//                                    }
//                                }
//                                if constexpr(_i < sizeof...(_finite_elements) - 1) {
//                                    self.template operator()<_i + 1u>(self);
//                                }
//                            };
//                            set_index(set_index);
//                            index += index2;
//                            found = true;
//                        }
//                        else {
//                            if (!found) {
//                                index += size();
//                            }
//                        };
//                    }
//                    else {
//                        if constexpr (__i < std::tuple_size_v<_T> - 1) {
//                            self1.template operator ()<std::tuple_element_t<__i + 1, typename _T::_FiniteElements>, __i + 1>(self1);
//                        }
//                    }
//                };
                getFiniteElementIndex2<__finite_element>(index, found);
                return index;
            }

//            template<auto __finite_element>
//            static constexpr
//            Integer
//            getFiniteElementIndex2()
//            {
//                auto index = -1;
//                auto mli = [&] <auto fe = lolita::finite_element::ElementGroup<_finite_elements...>{}> (auto & self) constexpr mutable {
//                    if constexpr (fe.has()) {
//                        if (index < 0) {
//                            if (fe.template hasFiniteElement<__finite_element>()) {
//                                index = lolita::numerics::abs(index) + 1;
//                            }
//                            else {
//                                index -= size();
//                            }
//                        }
//                    }
//                    else {
//                        (self.template operator()<_finite_elements>(self), ...);
//                    }
//                };
//                mli(mli);
//                return index;
//            }

            /**
             * @brief Fetch the finite element index within the _finite_element list
             * @tparam __finite_element the finite element object to find
             * @return the finite_element index if found, and the size of the _finite_element list otherwise
             */
//            template<auto __finite_element>
//            static constexpr
//            Integer
//            getFiniteElementIndex()
//            requires((lolita::finite_element::FiniteElementConcept<std::remove_cvref_t<decltype(_finite_elements)>>) && ...)
//            {
//                auto index = sizeof...(_finite_elements);
//                using _Elements = std::tuple<std::remove_cvref_t<decltype(_finite_elements)>...>;
//                auto const constexpr elements = _Elements{_finite_elements...};
//                auto set_index = [&] <Integer _i = 0u> (auto & self)
//                        constexpr mutable
//                {
//                    if constexpr (std::is_same_v<std::tuple_element_t<_i, _Elements>, std::remove_cvref_t<decltype(__finite_element)>>) {
//                        if (__finite_element == std::get<_i>(elements)) {
//                            index = _i;
//                        }
//                    }
//                    if constexpr(_i < sizeof...(_finite_elements) - 1) {
//                        self.template operator()<_i + 1u>(self);
//                    }
//                };
//                set_index(set_index);
//                return index;
//            }

        };

    }

    namespace mesh
    {

        enum struct Format
        {

            Gmsh,

        };

        template<lolita::domain::Domain domain, lolita::field::UnknownConcept auto... unknowns>
        struct Mesh
        {

        };

    }

}

#endif //LOLITA_LOLITA_USER_HXX
