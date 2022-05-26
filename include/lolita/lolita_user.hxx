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

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"

namespace lolita
{

    namespace geometry
    {

        using Point = lolita::matrix::Vector<lolita::real, 3>;

        /*
         *
         */

        template<lolita::matrix::MatrixConcept _Matrix>
        static inline
        lolita::matrix::Vector<typename _Matrix::Scalar, lolita::matrix::rows<_Matrix>()>
        getBarycenter(
                _Matrix const &
                point_args
        )
        requires(lolita::matrix::rows<_Matrix>() > 0)
        {
            auto barycenter = lolita::matrix::Vector<typename _Matrix::Scalar, lolita::matrix::rows<_Matrix>()>().setZero();
            for (auto i = 0; i < lolita::matrix::cols<_Matrix>(); ++i) {
                barycenter += point_args.col(i);
            }
            barycenter /= lolita::real(lolita::matrix::cols<_Matrix>());
            return barycenter;
        }

        /*
         *
         */

        template<lolita::matrix::MatrixConcept _Matrix>
        static inline
        auto
        getNormalVector2(
                _Matrix const &
                vector_args
        );

        template<lolita::matrix::MatrixConcept _Matrix>
        static inline
        auto
        getNormalVector2(
                _Matrix const &
                vector_args
        )
        requires(lolita::matrix::rows<_Matrix>() == 1 && lolita::matrix::cols<_Matrix>() == 3)
        {
            auto norm = vector_args.norm();
            return lolita::matrix::Vector<typename _Matrix::Scalar, 3>{vector_args(1)/norm, -vector_args(0)/norm, 0};
        }

        template<lolita::matrix::MatrixConcept _Matrix>
        static inline
        auto
        getNormalVector2(
                _Matrix const &
                vector_args
        )
        requires(lolita::matrix::rows<_Matrix>() == 3 && lolita::matrix::cols<_Matrix>() == 3)
        {
            auto v0 = vector_args.col(0) / vector_args.col(0).norm();
            auto v1 = vector_args.col(1) / vector_args.col(1).norm();
            return v0.cross(v1);
        }

        /*
         *
         */

        template<lolita::matrix::MatrixConcept _Matrix>
        static inline
        auto
        getNormalVector(
                _Matrix const &
                vector_args
        );

        template<lolita::matrix::MatrixConcept _Matrix>
        static inline
        auto
        getNormalVector(
                _Matrix const &
                vector_args
        )
        requires(lolita::matrix::rows<_Matrix>() == 1 && lolita::matrix::cols<_Matrix>() == 1)
        {
            return vector_args;
        }

        template<lolita::matrix::MatrixConcept _Matrix>
        static inline
        auto
        getNormalVector(
                _Matrix const &
                vector_args
        )
        requires(lolita::matrix::rows<_Matrix>() == 2 && lolita::matrix::cols<_Matrix>() == 1)
        {
            auto norm = vector_args.norm();
            return lolita::matrix::Vector<typename _Matrix::Scalar, 2>{vector_args(1)/norm, -vector_args(0)/norm};
        }

        template<lolita::matrix::MatrixConcept _Matrix>
        static inline
        auto
        getNormalVector(
                _Matrix const &
                vector_args
        )
        requires(lolita::matrix::rows<_Matrix>() == 3 && lolita::matrix::cols<_Matrix>() == 2)
        {
            auto v0 = vector_args.col(0);
            auto v1 = vector_args.col(1);
            v0 /= v0.norm();
            v1 /= v1.norm();
            return v0.cross(v1);
        }

        /*
         *
         */

        template<lolita::matrix::MatrixConcept _Matrix>
        static inline
        auto
        getRotationMatrix(
                _Matrix const &
                point_args
        );

        template<lolita::matrix::MatrixConcept _Matrix>
        static inline
        auto
        getRotationMatrix(
                _Matrix const &
                point_args
        )
        requires(lolita::matrix::rows<_Matrix>() == 1 && lolita::matrix::cols<_Matrix>() == 1)
        {
            return point_args;
        }

        template<lolita::matrix::MatrixConcept _Matrix>
        static inline
        auto
        getRotationMatrix(
                _Matrix const &
                point_args
        )
        requires(lolita::matrix::rows<_Matrix>() == 2 && lolita::matrix::cols<_Matrix>() == 2)
        {
            auto rotation_matrix = lolita::matrix::Matrix<typename _Matrix::Scalar, 2, 2>();
            auto edge = point_args.col(1) - point_args.col(0);
            edge /= edge.norm();
            rotation_matrix(0, 0) = edge(0);
            rotation_matrix(0, 1) = edge(1);
            rotation_matrix(1, 0) = edge(1);
            rotation_matrix(1, 1) = -edge(0);
            return rotation_matrix;
        }

        template<lolita::matrix::MatrixConcept _Matrix>
        static inline
        auto
        getRotationMatrix(
                _Matrix const &
                point_args
        )
        requires(lolita::matrix::rows<_Matrix>() == 3 && lolita::matrix::cols<_Matrix>() == 3)
        {
            auto rotation_matrix = lolita::matrix::Matrix<typename _Matrix::Scalar, 3, 3>();
            auto x_axis = point_args.col(2) - point_args.col(0);
            auto vector = point_args.col(1) - point_args.col(0);
            auto z_axis = x_axis.cross(vector);
            auto y_axis = z_axis.cross(x_axis);
            x_axis /= x_axis.norm();
            y_axis /= y_axis.norm();
            z_axis /= z_axis.norm();
            rotation_matrix.row(0) = x_axis;
            rotation_matrix.row(1) = y_axis;
            rotation_matrix.row(2) = z_axis;
            return rotation_matrix;
        }

        enum struct Frame
        {

            AxiSymmetric,
            Cartesian,

        };

        static constexpr
        std::basic_string_view<lolita::character>
        readFrame(
                lolita::geometry::Frame
                frame
        )
        {
            if (frame == lolita::geometry::Frame::AxiSymmetric) {
                return std::basic_string_view<lolita::character>("AxiSymmetric");
            }
            else {
                return std::basic_string_view<lolita::character>("Cartesian");
            }
        }

        struct Domain
        {

            constexpr
            Domain(
                    std::basic_string_view<lolita::character> &&
                    tag,
                    lolita::index
                    dim,
                    lolita::geometry::Frame
                    frame
            )
            :
            tag_(lolita::utility::label(std::forward<std::basic_string_view<lolita::character>>(tag))),
            dim_(dim),
            frame_(frame)
            {}

            constexpr
            lolita::boolean
            operator==(
                    Domain const &
                    other
            )
            const = default;

            constexpr
            lolita::boolean
            operator!=(
                    Domain const &
                    other
            )
            const = default;

            constexpr
            lolita::index
            ordIntegration(
                    lolita::index
                    ord
            )
            const
            {
                return frame_ == lolita::geometry::Frame::AxiSymmetric ? 2 * ord + 1 : 2 * ord;
            }

            lolita::utility::Label tag_;

            lolita::index dim_;

            lolita::geometry::Frame frame_;

        };

    }

    namespace field
    {

        enum struct Mapping
        {

            Gradient,
            Identity,
            Divergence,
            LargeStrain,
            SmallStrain,
            LargeStrainPlane,
            SmallStrainPlane,
            LargeStrainSolid,
            SmallStrainSolid

        };

//        static constexpr
//        std::basic_string_view<lolita::character>
//        readMapping(
//                lolita::field::Mapping
//                mapping
//        )
//        {
//            if (mapping == lolita::field::Mapping::Gradient) {
//                return std::basic_string_view<lolita::character>("Gradient");
//            }
//            else if (mapping == lolita::field::Mapping::Identity) {
//                return std::basic_string_view<lolita::character>("Identity");
//            }
//            else if (mapping == lolita::field::Mapping::Divergence) {
//                return std::basic_string_view<lolita::character>("Divergence");
//            }
//            else if (mapping == lolita::field::Mapping::LargeStrain) {
//                return std::basic_string_view<lolita::character>("LargeStrain");
//            }
//            else {
//                return std::basic_string_view<lolita::character>("SmallStrain");
//            }
//        }

        struct Tensor
        {

            constexpr
            Tensor(
                    std::basic_string_view<lolita::character> &&
                    tag,
                    lolita::index
                    dim
            )
            :
            tag_(lolita::utility::label(std::forward<std::basic_string_view<lolita::character>>(tag))),
            ord_(dim)
            {}

            constexpr
            lolita::boolean
            operator==(
                    Tensor const &
                    other
            )
            const = default;

            constexpr
            lolita::boolean
            operator!=(
                    Tensor const &
                    other
            )
            const = default;

            lolita::utility::Label tag_;

            lolita::index ord_;

        };

//        namespace detail {
//
//            static constexpr
//            lolita::utility::Label
//            getOutputFieldLabel(
//                    lolita::field::Tensor const &
//                    tensor,
//                    lolita::field::Mapping
//                    mapping
//            )
//            {
//                return lolita::utility::label(lolita::utility::readLabel(tensor.tag_), lolita::field::readMapping(mapping));
//            }
//
//        }

        template<typename... _Mapping>
        requires((std::same_as<_Mapping, lolita::field::Mapping> && ...) && sizeof...(_Mapping) > 0)
        struct Unknown
        {

            constexpr
            Unknown(
                    lolita::field::Tensor const &
                    tensor,
                    _Mapping &&...
                    mappings
            )
            :
            tensor_(tensor),
            mappings_({std::forward<lolita::field::Mapping>(mappings)...})
            {}

            constexpr
            Unknown(
                    lolita::field::Tensor &&
                    tensor,
                    _Mapping &&...
                    mappings
            )
            :
            tensor_(std::forward<lolita::field::Tensor>(tensor)),
            mappings_({std::forward<lolita::field::Mapping>(mappings)...})
            {}

            constexpr
            Unknown(
                    std::basic_string_view<lolita::character> &&
                    tag,
                    lolita::index
                    dim,
                    _Mapping &&...
                    mappings
            )
            :
            tensor_(std::forward<lolita::field::Tensor>(Tensor(std::forward<std::basic_string_view<lolita::character>>(tag), dim))),
            mappings_({std::forward<lolita::field::Mapping>(mappings)...})
            {}

            constexpr
            lolita::boolean
            operator==(
                    Unknown const &
                    other
            )
            const = default;

            constexpr
            lolita::boolean
            operator!=(
                    Unknown const &
                    other
            )
            const = default;

            lolita::field::Tensor tensor_;

            std::array<lolita::field::Mapping, sizeof...(_Mapping)> mappings_;

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

//        struct MaterialProperty : public Tensor
//        {
//
//            constexpr
//            MaterialProperty(
//                    std::basic_string_view<lolita::character> &&
//                    tag
//            )
//            :
//            Tensor(std::forward<std::basic_string_view<lolita::character>>(tag), 0)
//            {}
//
//        };
//
//        struct InternalVariable : public Tensor
//        {
//
//            constexpr
//            InternalVariable(
//                    std::basic_string_view<lolita::character> &&
//                    tag,
//                    lolita::index
//                    ord
//            )
//            :
//            Tensor(std::forward<std::basic_string_view<lolita::character>>(tag), ord)
//            {}
//
//        };
//
//        struct ExternalVariable : public Tensor
//        {
//
//            constexpr
//            ExternalVariable(
//                    std::basic_string_view<lolita::character> &&
//                    tag,
//                    lolita::index
//                    ord
//            )
//            :
//            Tensor(std::forward<std::basic_string_view<lolita::character>>(tag), ord)
//            {}
//
//        };
//
//        struct SomeField : public Tensor
//        {
//
//            constexpr
//            SomeField(
//                    std::basic_string_view<lolita::character> &&
//                    tag,
//                    lolita::index
//                    ord
//            )
//            :
//            Tensor(std::forward<std::basic_string_view<lolita::character>>(tag), ord)
//            {}
//
//        };

    }

    namespace behaviour
    {

        template<lolita::field::UnknownConcept auto... _unknown>
        struct MgisBehaviour3
        {

//            using MaterialPoint = mgis::behaviour::BehaviourData;

            constexpr
            MgisBehaviour3()
            {}

            constexpr
            lolita::boolean
            operator==(
                    MgisBehaviour3 const & other
            )
            const = default;

            constexpr
            lolita::boolean
            operator!=(
                    MgisBehaviour3 const & other
            )
            const = default;

            template<lolita::field::UnknownConcept auto __unknown>
            constexpr
            lolita::index
            getUnknownIndex()
            const
            {
                auto idx = lolita::index(sizeof...(_unknown));
                auto set_index = [&] <lolita::index _i = 0> (
                        auto & self
                )
                        constexpr mutable
                {
                    using Unknowns = std::tuple<std::remove_cvref_t<decltype(_unknown)>...>;
                    auto constexpr unknowns = Unknowns{_unknown...};
                    if constexpr (std::is_same_v<std::tuple_element_t<_i, Unknowns>, std::remove_cvref_t<decltype(__unknown)>>) {
                        if (lolita::utility::template get<_i>(unknowns) == __unknown) {
                            idx = _i;
                        }
                    }
                    if constexpr (_i < sizeof...(_unknown) - 1) {
                        self.template operator()<_i + 1>(self);
                    }
                };
                set_index(set_index);
                return idx;
            }

        };

        namespace detail
        {

            template<typename _T>
            struct IsMgisBehaviour3 : std::false_type {};


            template<lolita::field::UnknownConcept auto... _unknown>
            struct IsMgisBehaviour3<MgisBehaviour3<_unknown...>> : std::true_type {};

        }

        template<typename _T>
        concept MgisBehaviour3Concept = detail::IsMgisBehaviour3<_T>::value;

        template<lolita::behaviour::MgisBehaviour3Concept auto _behaviour>
        struct Behaviour
        {

            lolita::behaviour::MgisBehaviour3Concept auto const static constexpr b = _behaviour;

            Behaviour()
            {}

            mgis::behaviour::Behaviour behaviour_;

        };

        template<lolita::behaviour::MgisBehaviour3Concept auto... _behaviour>
        struct BehaviourHolder
        {

            using _Behaviours = std::tuple<lolita::behaviour::Behaviour<_behaviour>...>;

            BehaviourHolder(
                    lolita::behaviour::Behaviour<_behaviour> const &...
                    behaviours
            )
            :
            behaviours_(behaviours...)
            {}

            constexpr
            lolita::boolean
            operator==(
                    BehaviourHolder const & other
            )
            const = default;

            constexpr
            lolita::boolean
            operator!=(
                    BehaviourHolder const & other
            )
            const = default;

        private:

            template<lolita::behaviour::MgisBehaviour3Concept auto __behaviour>
            static constexpr
            lolita::index
            getBehaviourIndex()
            {
                auto index = std::tuple_size_v<_Behaviours>;
                auto set_index = [&] <lolita::index _i = 0> (auto & self)
                constexpr mutable
                {
                    if constexpr (std::is_same_v<lolita::behaviour::Behaviour<__behaviour>, std::tuple_element_t<_i, _Behaviours>>) {
                        index = _i;
                    }
                    if constexpr (_i < std::tuple_size_v<_Behaviours> - 1) {
                        self.template operator()<_i + 1>(self);
                    }
                };
                set_index(set_index);
                return index;
            }

        public:

            template<lolita::behaviour::MgisBehaviour3Concept auto __behaviour>
            std::tuple_element_t<getBehaviourIndex<__behaviour>(), _Behaviours> const &
            getBehaviour()
            const
            {
                return std::get<getBehaviourIndex<__behaviour>()>(behaviours_);
            }

            template<lolita::behaviour::MgisBehaviour3Concept auto __behaviour>
            std::tuple_element_t<getBehaviourIndex<__behaviour>(), _Behaviours> &
            getBehaviour()
            {
                return std::get<getBehaviourIndex<__behaviour>()>(behaviours_);
            }

            _Behaviours behaviours_;

        };

        using ParameterFunction = std::function<lolita::real(lolita::geometry::Point const &, lolita::real const &)>;

        struct MgisParameter
        {

            MgisParameter(
                    std::basic_string<lolita::character> &&
                    parameter_tag,
                    lolita::behaviour::ParameterFunction &&
                    function
            )
            :
            parameter_tag_(std::forward<std::basic_string<lolita::character>>(parameter_tag)),
            function_(std::forward<lolita::behaviour::ParameterFunction>(function))
            {}

            constexpr
            lolita::boolean
            operator==(
                    MgisParameter const & other
            )
            const
            {
                return other.parameter_tag_ == this->parameter_tag_;
            }

            constexpr
            lolita::boolean
            operator!=(
                    MgisParameter const & other
            )
            const
            {
                return !(* this == other);
            }

            std::basic_string<lolita::character> parameter_tag_;

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
                    std::basic_string<lolita::character> &&
                    unknown_tag,
                    std::basic_string<lolita::character> &&
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
            unknown_tag_(std::forward<std::basic_string<lolita::character>>(unknown_tag)),
            domain_tag_(std::forward<std::basic_string<lolita::character>>(domain_tag)),
            behaviour_data_(std::make_shared<lolita::behaviour::MgisBehaviourData>(lolita::behaviour::MgisBehaviourData{
                mgis::behaviour::load(path, name, hypothesis),
                std::forward<std::vector<lolita::behaviour::MgisParameter>>(parameters)
            }))
            {}

            MgisBehaviour(
                    std::basic_string<lolita::character> &&
                    unknown_tag,
                    std::basic_string<lolita::character> &&
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
            unknown_tag_(std::forward<std::basic_string<lolita::character>>(unknown_tag)),
            domain_tag_(std::forward<std::basic_string<lolita::character>>(domain_tag)),
            behaviour_data_(std::make_shared<lolita::behaviour::MgisBehaviourData>(lolita::behaviour::MgisBehaviourData{
                mgis::behaviour::load(finite_strain_behaviour_options, path, name, hypothesis),
                std::forward<std::vector<lolita::behaviour::MgisParameter>>(parameters)
            }))
            {}

            std::basic_string<lolita::character> unknown_tag_;

            std::basic_string<lolita::character> domain_tag_;

            std::shared_ptr<lolita::behaviour::MgisBehaviourData> behaviour_data_;

        };

    }

    namespace finite_element
    {

        enum struct Method
        {

            HHO,
            Lagrange

        };

        enum struct Loading
        {

            Natural,
            Constraint,

        };

        enum struct Basis
        {

            Monomial,
            Lagrange,

        };

        enum struct Quadrature
        {

            Gauss,

        };

        using LoadFunction = std::function<lolita::real(lolita::geometry::Point const &, lolita::real const &)>;

        struct LoadComponent
        {

            auto const static constexpr zero = [] (auto const &, auto const &) constexpr { return lolita::real(0); };

            LoadComponent()
            :
            function_(zero),
            loading_(lolita::finite_element::Loading::Natural)
            {}

            LoadComponent(
                    lolita::finite_element::LoadFunction &&
                    function,
                    lolita::finite_element::Loading
                    loading
            )
            :
            function_(std::forward<lolita::finite_element::LoadFunction>(function)),
            loading_(loading)
            {}

            LoadComponent(
                    lolita::finite_element::LoadFunction const &
                    function,
                    lolita::finite_element::Loading
                    loading
            )
            :
            function_(function),
            loading_(loading)
            {}

            lolita::real
            getImposedValue(
                    lolita::geometry::Point const &
                    point,
                    lolita::real const &
                    time
            )
            const
            {
                return function_(point, time);
            }

            lolita::finite_element::Loading loading_;

            lolita::finite_element::LoadFunction function_;

        };

        struct Load
        {

            Load(
                    std::basic_string_view<lolita::character> &&
                    unknown_tag,
                    std::basic_string_view<lolita::character> &&
                    domain_tag,
                    lolita::index
                    row,
                    lolita::index
                    col,
                    lolita::finite_element::LoadFunction &&
                    function,
                    lolita::finite_element::Loading
                    loading
            )
            :
            unknown_tag_(std::forward<std::basic_string_view<lolita::character>>(unknown_tag)),
            domain_tag_(std::forward<std::basic_string_view<lolita::character>>(domain_tag)),
            components_(lolita::matrix::Coordinates{row, col}),
            load_(std::make_shared<lolita::finite_element::LoadComponent>(lolita::finite_element::LoadComponent(std::forward<lolita::finite_element::LoadFunction>(function), loading)))
            {}

            Load(
                    std::basic_string_view<lolita::character> && unknown_tag,
                    std::basic_string_view<lolita::character> && domain_tag,
                    lolita::index row,
                    lolita::index col,
                    lolita::finite_element::LoadFunction const & function,
                    lolita::finite_element::Loading loading
            )
            :
            unknown_tag_(std::forward<std::basic_string_view<lolita::character>>(unknown_tag)),
            domain_tag_(std::forward<std::basic_string_view<lolita::character>>(domain_tag)),
            components_(lolita::matrix::Coordinates{row, col}),
            load_(std::make_shared<lolita::finite_element::LoadComponent>(lolita::finite_element::LoadComponent(function, loading)))
            {}

            std::basic_string_view<lolita::character> unknown_tag_;

            std::basic_string_view<lolita::character> domain_tag_;

            lolita::matrix::Coordinates components_;

            std::shared_ptr<lolita::finite_element::LoadComponent> load_;

        };

        struct FiniteElementMethod
        {

            constexpr
            FiniteElementMethod(
                    lolita::finite_element::Method method
            )
            :
            method_(method)
            {}

            constexpr
            lolita::boolean
            operator==(
                    FiniteElementMethod const & other
            )
            const = default;

            constexpr
            lolita::boolean
            operator!=(
                    FiniteElementMethod const & other
            )
            const = default;

            lolita::finite_element::Method method_;

        };

        struct HybridHighOrder : public lolita::finite_element::FiniteElementMethod
        {

            constexpr
            HybridHighOrder()
            :
            FiniteElementMethod(lolita::finite_element::Method::HHO),
            ord_cell_(-1),
            ord_face_(-1)
            {}

            constexpr
            HybridHighOrder(
                    lolita::integer ord_cell,
                    lolita::integer ord_face
            )
            :
            FiniteElementMethod(lolita::finite_element::Method::HHO),
            ord_cell_(ord_cell),
            ord_face_(ord_face)
            {}

            constexpr
            lolita::integer
            ordMapping(
                    lolita::field::Mapping mapping
            )
            const
            {
                return mapping == lolita::field::Mapping::Identity ? ord_cell_ : ord_face_;
            }

            lolita::integer ord_cell_;

            lolita::integer ord_face_;

        };

        template<typename _T>
        concept FiniteElementMethodConcept = std::derived_from<_T, lolita::finite_element::FiniteElementMethod> && requires(
                std::remove_reference_t<_T> const arg,
                lolita::field::Mapping mapping
        )
        {
            { arg.ord_cell_ } -> std::convertible_to<lolita::index>;
            { arg.ordMapping(mapping) } -> std::convertible_to<lolita::index>;
        };

        template<lolita::field::UnknownConcept auto _unknown, lolita::behaviour::MgisBehaviour3Concept auto _behaviour, lolita::finite_element::FiniteElementMethodConcept auto _finite_element_method>
        struct FiniteElement
        {

            using UnknownType = std::remove_cvref_t<decltype(_unknown)>;

            using BehaviourType = std::remove_cvref_t<decltype(_behaviour)>;

            using FiniteElementMethodType = std::remove_cvref_t<decltype(_finite_element_method)>;

            constexpr
            FiniteElement(
                    lolita::finite_element::Quadrature quadrature,
                    lolita::index ord_quadrature
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
            lolita::boolean
            operator==(
                    FiniteElement const & other
            )
            const = default;

            constexpr
            lolita::boolean
            operator!=(
                    FiniteElement const & other
            )
            const = default;

            UnknownType unknown_;

            BehaviourType behaviour_;

            FiniteElementMethodType discretization_;

            lolita::finite_element::Quadrature quadrature_;

            lolita::index ord_quadrature_;

        };

        namespace detail
        {

            template<typename T>
            struct IsFiniteElementConcept : std::false_type {};

            template<lolita::field::UnknownConcept auto _unknown, auto _behaviour, lolita::finite_element::FiniteElementMethodConcept auto _finite_element_method>
            struct IsFiniteElementConcept<FiniteElement<_unknown, _behaviour, _finite_element_method>> : std::true_type {};

        }

        template<typename _T>
        concept FiniteElementConcept = detail::IsFiniteElementConcept<_T>::value;

        template<lolita::finite_element::FiniteElementConcept... FiniteElementT>
        struct MixedElement
        {



        };

    }

    namespace mesh
    {

        enum struct Format
        {

            Gmsh,

        };

        template<lolita::geometry::Domain domain, lolita::field::UnknownConcept auto... unknowns>
        struct Mesh
        {

        };

    }

}

#endif //LOLITA_LOLITA_USER_HXX
