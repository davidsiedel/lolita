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

        template<lolita::matrix::MatrixConcept _Matrix>
        static inline
        lolita::matrix::Vector<typename _Matrix::Scalar, lolita::matrix::rows<_Matrix>()>
        getBarycenter(
                _Matrix const &
                point_args
        )
        requires(lolita::matrix::rows<_Matrix>() > 0)
        {
            auto barycenter = lolita::matrix::Vector<typename _Matrix::Scalar, lolita::matrix::rows<_Matrix>()>();
            for (auto i = 0; i < lolita::matrix::cols<_Matrix>(); ++i) {
                barycenter += point_args.col(i);
            }
            barycenter /= lolita::real(lolita::matrix::cols<_Matrix>());
            return barycenter;
        }

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
            SmallStrain

        };

        static constexpr
        std::basic_string_view<lolita::character>
        readMapping(
                lolita::field::Mapping
                mapping
        )
        {
            if (mapping == lolita::field::Mapping::Gradient) {
                return std::basic_string_view<lolita::character>("Gradient");
            }
            else if (mapping == lolita::field::Mapping::Identity) {
                return std::basic_string_view<lolita::character>("Identity");
            }
            else if (mapping == lolita::field::Mapping::Divergence) {
                return std::basic_string_view<lolita::character>("Divergence");
            }
            else if (mapping == lolita::field::Mapping::LargeStrain) {
                return std::basic_string_view<lolita::character>("LargeStrain");
            }
            else {
                return std::basic_string_view<lolita::character>("SmallStrain");
            }
        }

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

            constexpr
            lolita::matrix::Cardinality
            cardinality(
                    lolita::geometry::Domain const &
                    domain
            )
            const
            {
                if (ord_ == 0) {
                    return lolita::matrix::Cardinality(lolita::numerics::pow(domain.dim_, 0), lolita::numerics::pow(domain.dim_, 0));
                }
                else if (ord_ == 1) {
                    return lolita::matrix::Cardinality(lolita::numerics::pow(domain.dim_, 1), lolita::numerics::pow(domain.dim_, 0));
                }
                else if (ord_ == 2) {
                    return lolita::matrix::Cardinality(lolita::numerics::pow(domain.dim_, 1), lolita::numerics::pow(domain.dim_, 1));
                }
                else if (ord_ == 3) {
                    return lolita::matrix::Cardinality(lolita::numerics::pow(domain.dim_, 2), lolita::numerics::pow(domain.dim_, 1));
                }
                else {
                    return lolita::matrix::Cardinality(lolita::numerics::pow(domain.dim_, 2), lolita::numerics::pow(domain.dim_, 2));
                }
            }

            lolita::utility::Label tag_;

            lolita::index ord_;

        };

        namespace detail
        {

            static constexpr
            lolita::utility::Label
            getOutputFieldLabel(
                    lolita::field::Tensor const &
                    tensor,
                    lolita::field::Mapping
                    mapping
            )
            {
                return lolita::utility::label(lolita::utility::readLabel(tensor.tag_), lolita::field::readMapping(mapping));
            }

            template<lolita::field::Tensor tensor, lolita::geometry::Domain domain, lolita::field::Mapping mapping>
            struct MappingPolicy;

            template<lolita::field::Tensor tensor, lolita::geometry::Domain domain, lolita::field::Mapping mapping>
            requires(mapping == lolita::field::Mapping::Identity)
            struct MappingPolicy<tensor, domain, mapping>
            {

                lolita::field::Tensor const static constexpr tensor_output_ = lolita::field::Tensor(getOutputFieldLabel(tensor, mapping), tensor.ord_);

            };

            template<lolita::field::Tensor tensor, lolita::geometry::Domain domain, lolita::field::Mapping mapping>
            requires(mapping == lolita::field::Mapping::Gradient)
            struct MappingPolicy<tensor, domain, mapping>
            {

                lolita::field::Tensor const static constexpr tensor_output_ = lolita::field::Tensor(getOutputFieldLabel(tensor, mapping), tensor.ord_ + 1);

            };

        }

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

        template<typename... MappingT>
        requires((std::same_as<MappingT, lolita::field::Mapping> && ...) && sizeof...(MappingT) > 0)
        struct DegreeOfFreedom : public Tensor
        {

            constexpr
            DegreeOfFreedom(
                    std::basic_string_view<lolita::character> &&
                    tag,
                    lolita::index
                    dim,
                    MappingT &&...
                    mappings
            )
            :
            Tensor(std::forward<std::basic_string_view<lolita::character>>(tag), dim),
            mappings_({std::forward<lolita::field::Mapping>(mappings)...})
            {}

            std::array<lolita::field::Mapping, sizeof...(MappingT)> mappings_;

        };

        namespace detail
        {

            template<typename T>
            struct IsDegreeOfFreedomConcept : public std::false_type {};

            template<typename... MappingT>
            struct IsDegreeOfFreedomConcept<DegreeOfFreedom<MappingT...>> : public std::true_type {};

            template<typename T>
            struct IsUnknownConcept : public std::false_type {};

            template<typename... MappingT>
            struct IsUnknownConcept<Unknown<MappingT...>> : public std::true_type {};

        }

        template<typename T>
        concept DegreeOfFreedomConcept = detail::IsDegreeOfFreedomConcept<T>::value;

        template<typename T>
        concept UnknownConcept = detail::IsUnknownConcept<T>::value;

        struct MaterialProperty : public Tensor
        {

            constexpr
            MaterialProperty(
                    std::basic_string_view<lolita::character> &&
                    tag
            )
            :
            Tensor(std::forward<std::basic_string_view<lolita::character>>(tag), 0)
            {}

        };

        struct InternalVariable : public Tensor
        {

            constexpr
            InternalVariable(
                    std::basic_string_view<lolita::character> &&
                    tag,
                    lolita::index
                    ord
            )
            :
            Tensor(std::forward<std::basic_string_view<lolita::character>>(tag), ord)
            {}

        };

        struct ExternalVariable : public Tensor
        {

            constexpr
            ExternalVariable(
                    std::basic_string_view<lolita::character> &&
                    tag,
                    lolita::index
                    ord
            )
            :
            Tensor(std::forward<std::basic_string_view<lolita::character>>(tag), ord)
            {}

        };

        struct SomeField : public Tensor
        {

            constexpr
            SomeField(
                    std::basic_string_view<lolita::character> &&
                    tag,
                    lolita::index
                    ord
            )
            :
            Tensor(std::forward<std::basic_string_view<lolita::character>>(tag), ord)
            {}

        };

    }

    namespace behaviour
    {

        template<lolita::field::UnknownConcept... UnknownT>
        struct BHVData
        {

            BHVData(
                    UnknownT const &...
                    unknowns
            )
            :
            unknowns_({unknowns...})
            {}

            std::array<lolita::field::Tensor, 2> auxiliary_fields_;

            lolita::utility::Aggregate<UnknownT...> unknowns_;

        };

        struct MgisBehaviour
        {



        };

        template<lolita::field::UnknownConcept... UnknownT>
        struct MgisBehaviour2
        {

            using MaterialPoint = mgis::behaviour::BehaviourData;

            constexpr
            MgisBehaviour2(
                    UnknownT const &...
                    unknowns
            )
            :
            unknowns_({unknowns...})
            {}

            template<lolita::field::UnknownConcept auto unknown>
            constexpr
            lolita::index
            getUnknownIndex()
            const
            {
                auto idx = lolita::index(sizeof...(UnknownT));
                auto set_index = [&] <lolita::index I = 0> (
                        auto & self
                )
                constexpr mutable
                {
                    if constexpr (std::is_same_v<std::tuple_element_t<I, std::tuple<UnknownT...>>, std::remove_cvref_t<decltype(unknown)>>) {
                        if (lolita::utility::template get<I>(unknowns_) == unknown) {
                            idx = I;
                        }
//                        if (std::template get<I>(unknowns_) == unknown) {
//                            idx = I;
//                        }
                    }
                    if constexpr (I < sizeof...(UnknownT) - 1) {
                        self.template operator()<I + 1>(self);
                    }
                };
                set_index(set_index);
                return idx;
            }

            lolita::utility::Aggregate<UnknownT...> unknowns_;
//            std::tuple<UnknownT...> unknowns_;

        };

        template<lolita::field::UnknownConcept auto... _unknown>
        struct MgisBehaviour3
        {

            using MaterialPoint = mgis::behaviour::BehaviourData;

            constexpr
            MgisBehaviour3()
            {}

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
                    auto constexpr unknowns = Unknowns();
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

    }

    namespace finite_element
    {

        enum struct Quadrature
        {

            Gauss,

        };

        static constexpr
        std::basic_string_view<lolita::character>
        readQuadrature(
                lolita::finite_element::Quadrature
                quadrature
        )
        {
            if (quadrature == lolita::finite_element::Quadrature::Gauss) {
                return std::basic_string_view<lolita::character>("Gauss");
            }
            else {
                return std::basic_string_view<lolita::character>("SmallStrain");
            }
        }

        enum struct Method
        {

            HHO,
            Lagrange

        };

        struct FiniteElementMethod
        {

            constexpr
            FiniteElementMethod(
                    lolita::finite_element::Method
                    method
            )
            :
            method_(method)
            {}

            lolita::finite_element::Method method_;

        };

        struct HybridHighOrder : public lolita::finite_element::FiniteElementMethod
        {

            constexpr
            HybridHighOrder(
                    lolita::index
                    ord_cell,
                    lolita::index
                    ord_face
            )
            :
            FiniteElementMethod(lolita::finite_element::Method::HHO),
            ord_cell_(ord_cell),
            ord_face_(ord_face)
            {}

            constexpr
            lolita::index
            ordMapping(
                    lolita::field::Mapping
                    mapping
            )
            const
            {
                return mapping == lolita::field::Mapping::Identity ? ord_cell_ : ord_face_;
            }

            lolita::index ord_cell_;

            lolita::index ord_face_;

        };

        template<typename _T>
        concept FiniteElementMethodConcept = std::derived_from<_T, lolita::finite_element::FiniteElementMethod> && requires(
                std::remove_reference_t<_T> const
                arg,
                lolita::field::Mapping
                mapping
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
                    lolita::finite_element::Quadrature
                    quadrature,
                    lolita::index
                    ord_quadrature
            )
            :
            unknown_(_unknown),
            behaviour_(_behaviour),
            finite_element_method_(_finite_element_method),
            quadrature_(quadrature),
            ord_quadrature_(ord_quadrature)
            {}

            constexpr
            FiniteElement(
                    lolita::finite_element::Quadrature
                    quadrature
            )
            :
            unknown_(_unknown),
            behaviour_(_behaviour),
            finite_element_method_(_finite_element_method),
            quadrature_(quadrature),
            ord_quadrature_()
            {}

            UnknownType unknown_;

            BehaviourType behaviour_;

            FiniteElementMethodType finite_element_method_;

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

        enum struct Loading
        {

            Natural,
            Constraint,

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
                    i,
                    lolita::index
                    j,
                    lolita::finite_element::LoadComponent &&
                    fun
            )
            :
                    unknown_tag_(std::forward<std::basic_string_view<lolita::character>>(unknown_tag)),
                    domain_tag_(std::forward<std::basic_string_view<lolita::character>>(domain_tag)),
                    components_(lolita::matrix::Coordinates(i, j)),
                    load_(std::make_shared<lolita::finite_element::LoadComponent>(fun))
            {}

            std::basic_string_view<lolita::character> unknown_tag_;

            std::basic_string_view<lolita::character> domain_tag_;

            lolita::matrix::Coordinates components_;

            std::shared_ptr<lolita::finite_element::LoadComponent> load_;

        };

        template<lolita::finite_element::FiniteElementConcept... FiniteElementT>
        struct MixedElement
        {



        };

    }

    namespace mesh
    {

        enum struct MeshFormatType
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
