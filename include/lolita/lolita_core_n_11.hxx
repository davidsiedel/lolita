#ifndef CC79BDC5_49DB_4A81_8A93_18ABD6551AF1
#define CC79BDC5_49DB_4A81_8A93_18ABD6551AF1

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_defs.hxx"
#include "lolita/lolita_core_n_00.hxx"
#include "lolita/lolita_core_n_0.hxx"

namespace lolita2::geometry
{

    template<template<Element, Domain, auto...> typename t_FiniteElement, Element t_element, Domain t_domain, auto... t_args>
    struct FiniteElementGeometry;

    namespace quadrature
    {
    
        template<Element t_element, Quadrature t_quadrature>
        struct ElementQuadratureRuleTraits;
        
        template<Element t_element, Quadrature t_quadrature>
        requires(t_element.isNode() || !t_element.isNode())
        struct ElementQuadratureRuleTraits<t_element, t_quadrature>
        {
            
            lolita::integer static constexpr dim_ = 1;

            std::array<std::array<lolita::real, 3>, dim_> static constexpr reference_points_ = {
                    +0.0000000000000000, +0.0000000000000000, +0.0000000000000000
            };
            
            std::array<lolita::real, dim_> static constexpr reference_weights_ = {
                    +1.0000000000000000
            };

        };

    }

    namespace basis
    {
        
        template<Basis t_basis, auto t_arg>
        struct FiniteElementBasisTraits;

        template<Basis t_basis, lolita::integer t_dim>
        requires(t_basis.isMonomial())
        struct FiniteElementBasisTraits<t_basis, t_dim>
        {
            
            lolita::integer static constexpr dim_ = lolita::numerics::binomial(t_dim + t_basis.ord_, t_dim);
            
            lolita::integer static constexpr num_components_ = lolita::numerics::binomial(t_dim + t_basis.ord_, t_dim);
            
            lolita::integer static constexpr size_ = lolita::numerics::binomial(t_dim + t_basis.ord_, t_dim);
            
            template<template<Element, Domain, auto...> typename t_FiniteElement, Element t_element, Domain t_domain, auto... t_args>
            struct Implementation : FiniteElementGeometry<t_FiniteElement, t_element, t_domain, t_args...>
            {

            private:
            
                static constexpr
                std::array<std::array<lolita::integer, 3>, dim_>
                getExponents()
                {
                    auto exponents = std::array<std::array<lolita::integer, 3>, dim_>();
                    auto row = lolita::integer(0);
                    if constexpr (t_dim == 0)
                    {
                        exponents[row][0] = 0;
                        exponents[row][1] = 0;
                        exponents[row][2] = 0;
                    }
                    else if constexpr (t_dim == 1)
                    {
                        for (lolita::integer i = 0; i < t_basis.ord_ + 1; ++i)
                        {
                            exponents[row][0] = i;
                            exponents[row][1] = 0;
                            exponents[row][2] = 0;
                            row += 1;
                        }
                    }
                    else if constexpr (t_dim == 2)
                    {
                        for (lolita::integer i = 0; i < t_basis.ord_ + 1; ++i)
                        {
                            for (lolita::integer j = 0; j < i + 1; ++j)
                            {
                                exponents[row][0] = i - j;
                                exponents[row][1] = j;
                                exponents[row][2] = 0;
                                row += 1;
                            }
                        }
                    }
                    else if constexpr (t_dim == 3)
                    {
                        for (lolita::integer i = 0; i < t_basis.ord_ + 1; ++i)
                        {
                            for (lolita::integer j = 0; j < i + 1; ++j)
                            {
                                for (lolita::integer k = 0; k < i + 1; ++k)
                                {
                                    if (j + k < i + 1)
                                    {
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
                
                std::array<std::array<lolita::integer, 3>, dim_> static constexpr exponents_ = getExponents();

            public:
            
                lolita::matrix::Vector<lolita::real, dim_>
                getBasisEvaluation(
                    Point const & point
                )
                const
                {
                    auto basis_vector_values = lolita::matrix::Vector<lolita::real, dim_>();
                    auto const centroid = this->getReferenceCentroid();
                    auto const diameters = this->getCurrentDiameters();
                    for (lolita::integer i = 0; i < dim_; ++i)
                    {
                        auto value = lolita::real(1);
                        for (lolita::integer j = 0; j < t_dim; ++j)
                        {
                            auto dist = this->getRiemannianDistance(centroid, point, j);
                            value *= std::pow(2.0 * dist / diameters(j), exponents_[i][j]);
                        }
                        basis_vector_values(i) = value;
                    }
                    return basis_vector_values;
                }
                
                lolita::matrix::Vector<lolita::real, dim_>
                getBasisDerivative(
                    Point const & point,
                    lolita::integer derivative_direction
                )
                const
                {
                    auto basis_vector_values = lolita::matrix::Vector<lolita::real, dim_>();
                    auto const centroid = this->getReferenceCentroid();
                    auto const diameters = this->getCurrentDiameters();
                    for (lolita::integer i = 0; i < dim_; ++i)
                    {
                        auto value = lolita::real(1);
                        for (lolita::integer j = 0; j < t_dim; ++j)
                        {
                            if (j != derivative_direction)
                            {
                                auto dist = this->getRiemannianDistance(centroid, point, j);
                                value *= std::pow(2.0 * (dist) / diameters(j), exponents_[i][j]);
                            }
                            else
                            {
                                if (exponents_[i][j] > 0)
                                {
                                    auto c = 2.0 * exponents_[i][j] / diameters(j);
                                    auto dist = this->getRiemannianDistance(centroid, point, j);
                                    value *= c * std::pow(2.0 * (dist) / diameters(j), exponents_[i][j] - 1);
                                }
                                else
                                {
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

        template<Basis t_basis, Element t_element>
        requires(t_basis.isMonomial())
        struct FiniteElementBasisTraits<t_basis, t_element> : FiniteElementBasisTraits<t_basis, t_element.dim_>
        {};

    }

    namespace discretization
    {

        namespace hybrid_discontinuous_galerkin
        {

            struct Stabilization
            {

                enum Type{

                    Hdg,
                    Hho,

                };

                constexpr
                Stabilization(
                    Type type
                )
                :
                type_(type)
                {}

                constexpr
                lolita::boolean
                operator==(
                    Stabilization const & other
                )
                const = default;

                constexpr
                lolita::boolean
                operator!=(
                    Stabilization const & other
                )
                const = default;

                constexpr
                lolita::boolean
                isHdg()
                const
                {
                    return type_ == Type::Hdg;
                }

                constexpr
                lolita::boolean
                isHho()
                const
                {
                    return type_ == Type::Hho;
                }

                Type type_;

            };

            template<Basis t_cell_basis, Basis t_face_basis>
            struct HybridDiscontinuousGalerkinElementTraits
            {

                template<Element t_element, Domain t_domain, Field t_field>
                static constexpr
                lolita::integer
                getNumUnknowns()
                {
                    // FieldTraits<t_domain, t_field>::size_;
                    return 0;
                }

            };

        }

    }
    
    template<template<Element, Domain, auto...> typename t_T, Element t_element, Domain t_domain, auto... t_args>
    struct FiniteElementGeometry
    {

    private:
    
        using t_ElementTraits = ElementTraits<t_element, t_domain>;
        
        template<Element t__element, Domain t__domain>
        using t_ElementPointer = std::shared_ptr<t_T<t__element, t__domain, t_args...>>;

    public:
    
        using t_InnerNeighbors = typename t_ElementTraits::template InnerConnectivity<t_ElementPointer>;
        
        using t_OuterNeighbors = typename t_ElementTraits::template OuterConnectivity<t_ElementPointer>;
        
        lolita::boolean
        operator==(
            FiniteElementGeometry const & other
        )
        const = default;
        
        lolita::boolean
        operator!=(
            FiniteElementGeometry const & other
        )
        const = default;
        
        std::basic_string<lolita::character>
        getHash()
        const
        requires(t_element.isNode())
        {
            return std::to_string(this->tag_);
        }
        
        std::basic_string<lolita::character>
        getHash()
        const
        {
            std::basic_stringstream<lolita::character> hash;
            auto const & nodes = getInnerNeighbors<t_element.dim_ - 1, 0>();
            for (auto const & node : nodes)
            {
                hash << node->getHash();
            }
            return hash.str();
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        static constexpr
        lolita::integer
        getInnerNeighborNodeConnection(
            lolita::integer i,
            lolita::integer j
        )
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(ElementTraits<t_element, t_domain>::node_connectivity_))[i][j];
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_OuterNeighbors>> &
        getOuterNeighbors()
        {
            return std::get<t_j>(std::get<t_i>(outer_neighbors_));
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_OuterNeighbors>> const &
        getOuterNeighbors()
        const
        {
            return std::get<t_j>(std::get<t_i>(outer_neighbors_));
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_InnerNeighbors>> &
        getInnerNeighbors()
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(inner_neighbors_));
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_InnerNeighbors>> const &
        getInnerNeighbors()
        const
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(inner_neighbors_));
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        lolita::integer
        getInnerNeighborIndex(
            lolita::integer i
        )
        const
        requires(!t_element.isNode())
        {
            auto constexpr t_component = t_ElementTraits::template getInnerNeighbor<t_i, t_j>();
            using t_NeighbourTraits = ElementTraits<t_component, t_domain>;
            auto constexpr t_coordinates = t_NeighbourTraits::template getOuterNeighborCoordinates<t_element>();
            auto const & items = getInnerNeighbors<t_i, t_j>()[i]->template getOuterNeighbors<t_coordinates.dim_, t_coordinates.tag_>();
            auto is_equal = [&] (t_ElementPointer<t_element, t_domain> const & ptr_element)
            {
                return * ptr_element == * this;
            };
            return std::distance(items.begin(), std::find_if(items.begin(), items.end(), is_equal));
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        lolita::integer
        getInnerNeighborOrientation(
            lolita::integer i
        )
        const
        requires(!t_element.isNode())
        {
            return getInnerNeighborIndex<t_i, t_j>(i) == 0 ? 1 : -1;
        }
    
        Point &
        getCurrentCoordinates()
        requires(t_element.isNode())
        {
            return * this->coordinates_;
        }
        
        Point const &
        getCurrentCoordinates()
        const
        requires(t_element.isNode())
        {
            return * this->coordinates_;
        }
        
        lolita::matrix::Matrix<lolita::real, 3, t_element.num_nodes_>
        getCurrentCoordinates()
        const
        requires(!t_element.isNode())
        {
            auto current_nodes_coordinates = lolita::matrix::Matrix<lolita::real, 3, t_element.num_nodes_>();
            auto count = lolita::integer(0);
            for (auto const & node : this->template getInnerNeighbors<t_element.dim_ - 1, 0>())
            {
                current_nodes_coordinates.col(count) = node->getCurrentCoordinates();
                count ++;
            }
            return current_nodes_coordinates;
        }
        
        static
        lolita::matrix::Span<lolita::matrix::Matrix<lolita::real, 3, t_element.num_nodes_, lolita::matrix::col_major> const>
        getReferenceCoordinates()
        {
            using t_ReferenceCoordinates = lolita::matrix::Span<lolita::matrix::Matrix<lolita::real, 3, t_element.num_nodes_, lolita::matrix::col_major> const>;
            return t_ReferenceCoordinates(ElementTraits<t_element, t_domain>::reference_nodes_.begin()->begin());
        }
        
        static
        lolita::real
        getShapeMappingEvaluation(
            lolita::matrix::Vector<lolita::real, t_element.num_nodes_> const & nodal_field_values,
            Point const & reference_point
        )
        {
            return t_ElementTraits::getShapeMappingEvaluation(nodal_field_values, reference_point);
        }
        
        static
        lolita::real
        getShapeMappingDerivative(
            lolita::matrix::Vector<lolita::real, t_element.num_nodes_> const & nodal_field_values,
            Point const & reference_point,
            lolita::integer derivative_direction
        )
        {
            return t_ElementTraits::getShapeMappingDerivative(nodal_field_values, reference_point, derivative_direction);
        }
        
        lolita::real
        getShapeMappingDifferential(
            Point const & point
        )
        const
        {
            auto const current_coordinates = this->getCurrentCoordinates();
            auto ru = lolita::matrix::Matrix<lolita::real, 3, 3>();
            auto du = lolita::real(0);
            ru.setZero();
            for (lolita::integer i = 0; i < t_domain.dim_; ++i)
            {
                for (lolita::integer j = 0; j < t_element.dim_; ++j)
                {
                    ru(i, j) = FiniteElementGeometry::getShapeMappingDerivative(current_coordinates.row(i), point, j);
                }
            }
            if constexpr (t_element.dim_ == 3)
            {
                du = lolita::numerics::abs((ru.col(0).template cross(ru.col(1))).template dot(ru.col(2)));
            }
            else if constexpr (t_element.dim_ == 2)
            {
                du = lolita::numerics::abs((ru.col(0).template cross(ru.col(1))).norm());
            }
            else
            {
                du = lolita::numerics::abs(ru.col(0).norm());
            }
            if constexpr (t_domain.frame_ == Domain::Frame::AxiSymmetric)
            {
                lolita::real r0 = getShapeMappingEvaluation(current_coordinates.row(0), point);
                if (r0 < 1.e-10)
                {
                    r0 = 1.e-10;
                }
                du *= 2.0 * lolita::numerics::pi * r0;
            }
            return du;
        }
        
        lolita::real
        getRiemannianDistance(
            Point const & first_point,
            Point const & second_point,
            lolita::integer direction = -1
        )
        const
        {
            if constexpr (t_ElementTraits::hasDim(0))
            {
                auto const & current_coordinates = this->getCurrentCoordinates();
                auto distance = lolita::real();
                auto mp0 = Point();
                auto mp1 = Point();
                for (lolita::integer i = 0; i < t_element.dim_; ++i)
                {
                    mp0(i) = FiniteElementGeometry::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
                    mp1(i) = FiniteElementGeometry::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
                }
                direction == -1 ? distance = (mp1 - mp0).norm() : distance = (mp1 - mp0)(direction);
                return distance;
            }
            else
            {
                using SegmentQuadrature = quadrature::ElementQuadratureRuleTraits<Element::segment(1), Quadrature::gauss(4)>;
                auto distance = lolita::real(0);
                auto dt = lolita::real();
                auto const current_nodes_coordinates = this->getCurrentCoordinates();
                for (lolita::integer q = 0; q < SegmentQuadrature::dim_; ++q)
                {
                    auto pq = SegmentQuadrature::reference_points_[q][0];
                    auto wq = SegmentQuadrature::reference_weights_[q];
                    auto ru = lolita::matrix::Matrix<lolita::real, 3, 3>().setZero();
                    auto difference = second_point - first_point;
                    auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
                    for (lolita::integer i = 0; i < t_domain.dim_; ++i)
                    {
                        for (lolita::integer j = 0; j < t_element.dim_; ++j)
                        {
                            if (direction == -1 || i == static_cast<lolita::integer>(direction))
                            {
                                auto du = (1.0 / 2.0) * (second_point(j) - first_point(j));
                                auto dx = FiniteElementGeometry::getShapeMappingDerivative(current_nodes_coordinates.row(i), uq, j);
                                ru(i, j) = dx * du;
                            }
                        }
                    }
                    if constexpr (t_element.isFacet())
                    {
                        auto Eff = ru.col(0).template dot(ru.col(0));
                        auto Fff = ru.col(0).template dot(ru.col(1));
                        auto Gff = ru.col(1).template dot(ru.col(1));
                        dt = std::sqrt(Eff + 2.0 * Fff + Gff);
                    }
                    else if constexpr (t_element.isCurve())
                    {
                        auto Eff = ru.col(0).template dot(ru.col(0));
                        dt = std::sqrt(Eff);
                    }
                    else
                    {
                        dt = 0;
                    }
                    distance += wq * dt;
                }
                return distance;
            }
        }
        
        Point
        getCurrentDiameters()
        const
        {
            auto reference_coordinates = FiniteElementGeometry::getReferenceCoordinates();
            auto current_diameters = Point().setZero();
            for (lolita::integer i = 0; i < t_element.num_nodes_; ++i)
            {
                for (lolita::integer j = i + 1; j < t_element.num_nodes_; ++j)
                {
                    auto const & pt0 = reference_coordinates.col(i);
                    auto const & pt1 = reference_coordinates.col(j);
                    for (lolita::integer k = 0; k < 3; ++k)
                    {
                        auto new_value = lolita::numerics::abs(getRiemannianDistance(pt0, pt1, k));
                        auto & current_value = current_diameters(k);
                        if (new_value > current_value)
                        {
                            current_value = new_value;
                        }
                    }
                }
            }
            return current_diameters;
        }
        
        Point
        getCurrentCentroid()
        const
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto barycenter = Point();
            barycenter.setZero();
            for (auto i = 0; i < t_element.num_nodes_; ++i) {
                barycenter += current_nodes_coordinates.col(i);
            }
            barycenter /= lolita::real(t_element.num_nodes_);
            return barycenter;
        }
        
        static
        Point
        getReferenceCentroid()
        {
            auto reference_nodes_coordinates = FiniteElementGeometry::getReferenceCoordinates();
            auto barycenter = Point();
            barycenter.setZero();
            for (auto i = 0; i < t_element.num_nodes_; ++i) {
                barycenter += reference_nodes_coordinates.col(i);
            }
            barycenter /= lolita::real(t_element.num_nodes_);
            return barycenter;
        }
        
        static
        Point
        getReferenceDiameters()
        {
            auto dts = Point();
            auto nds = FiniteElementGeometry::getReferenceCoordinates();
            dts.setZero();
            for (lolita::integer i = 0; i < t_element.num_nodes_; ++i)
            {
                for (lolita::integer j = i + 1; j < t_element.num_nodes_; ++j)
                {
                    for (lolita::integer k = 0; k < 3; ++k)
                    {
                        auto new_value = lolita::numerics::abs(nds(k, i) - nds(k, j));
                        auto & current_value = dts(k);
                        if (new_value > current_value)
                        {
                            current_value = new_value;
                        }
                    }
                }
            }
            return dts;
        }
        
        Point
        getNormalVector(
            Point const & point
        )
        const
        requires(t_ElementTraits::isFace())
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto ru = lolita::matrix::Matrix<lolita::real, 3, t_element.dim_>();
            ru.setZero();
            for (lolita::integer i = 0; i < 3; ++i)
            {
                for (lolita::integer j = 0; j < t_element.dim_; ++j)
                {
                    ru(i, j) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, j);
                }
            }
            if constexpr (t_element.isNode()) {
                return Point{0, 0, 0};
            }
            else if constexpr (t_element.isCurve())
            {
                return Point{ru(1)/ru.norm(), -ru(0)/ru.norm(), 0};
            }
            else
            {
                return (ru.col(0) / ru.col(0).norm()).cross((ru.col(1) / ru.col(1).norm()));
            }
        }
        
        Point
        getTangentVector(
            Point const & point,
            lolita::integer direction
        )
        const
        requires(t_element.isCurve())
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto tangent_vector = Point();
            for (lolita::integer i = 0; i < 3; ++i)
            {
                tangent_vector(i) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, direction);
            }
            return tangent_vector;
        }
        
        template<Quadrature t_quadrature>
        static
        lolita::real
        getReferenceQuadratureWeight(
            lolita::integer index
        )
        {
            return quadrature::ElementQuadratureRuleTraits<t_element, t_quadrature>::reference_weights_[index];
        }
        
        template<Quadrature t_quadrature>
        static
        lolita::matrix::Span<Point const>
        getReferenceQuadraturePoint(
            lolita::integer index
        )
        {
            return lolita::matrix::Span<Point const>(
                    quadrature::ElementQuadratureRuleTraits<t_element, t_quadrature>::reference_points_[index].begin()
            );
        }
        
        template<Quadrature t_quadrature>
        lolita::real
        getCurrentQuadratureWeight(
            lolita::integer index
        )
        const
        {
            auto w = getReferenceQuadratureWeight<t_quadrature>(index);
            return w * getShapeMappingDifferential(getReferenceQuadraturePoint<t_quadrature>(index));
        }
        
        template<Quadrature t_quadrature>
        Point
        getCurrentQuadraturePoint(
            lolita::integer index
        )
        const
        {
            auto p = Point();
            auto const nds = this->getCurrentCoordinates();
            for (lolita::integer j = 0; j < 3; ++j)
            {
                p(j) = FiniteElementGeometry::getShapeMappingEvaluation(nds.row(j), getReferenceQuadraturePoint<t_quadrature>(index));
            }
            return p;
        }
        
        template<Quadrature t_quadrature, lolita::integer _i, lolita::integer _j>
        static
        lolita::real
        getInnerNeighborReferenceQuadratureWeight(
            lolita::integer index
        )
        requires(!t_element.isNode())
        {
            auto const constexpr _component = t_ElementTraits::template getInnerNeighbor<_i, _j>();
            using ComponentGeometry = FiniteElementGeometry<t_T, _component, t_domain, t_args...>;
            return ComponentGeometry::template getReferenceQuadratureWeight<t_quadrature>(index);
        }
        
        template<Quadrature t_quadrature, lolita::integer _i, lolita::integer _j>
        static
        Point
        getInnerNeighborReferenceQuadraturePoint(
            lolita::integer component_index,
            lolita::integer index
        )
        requires(!t_element.isNode())
        {
            auto const constexpr _component = t_ElementTraits ::template getInnerNeighbor<_i, _j>();
            auto p = Point();
            using ComponentGeometry = FiniteElementGeometry<t_T, _component, t_domain, t_args...>;
            auto const & elt_reference_nodes = ElementTraits<t_element, t_domain>::reference_nodes_;
            for (lolita::integer i = 0; i < 3; ++i)
            {
                auto cpt_coordinates = lolita::matrix::Vector<lolita::real, _component.num_nodes_>();
                for (lolita::integer j = 0; j < _component.num_nodes_; ++j)
                {
                    auto const node_tag = getInnerNeighborNodeConnection<_i, _j>(component_index, j);//.get(component_index).get(j);
                    cpt_coordinates(j) = elt_reference_nodes[node_tag][i];
                }
                auto cpt_reference_point = ComponentGeometry::template getReferenceQuadraturePoint<t_quadrature>(index);
                p(i) = ComponentGeometry::getShapeMappingEvaluation(cpt_coordinates, cpt_reference_point);
            }
            return p;
        }
        
        template<Quadrature t_quadrature, lolita::integer t_i, lolita::integer t_j>
        lolita::real
        getInnerNeighborCurrentQuadratureWeight(
            lolita::integer component_index,
            lolita::integer index
        )
        const
        requires(!t_element.isNode())
        {
            auto const & cmp =  this->template getInnerNeighbors<t_i, t_j>()[component_index];//.template get<I>().template get<J>().get(component_index).get();
            return cmp->template getCurrentQuadratureWeight<t_quadrature>(index);
        }
        
        template<Quadrature t_quadrature, lolita::integer t_i, lolita::integer t_j>
        Point
        getInnerNeighborCurrentQuadraturePoint(
            lolita::integer component_index,
            lolita::integer index
        )
        const
        requires(!t_element.isNode())
        {
            auto p = Point();
            auto const cpt_ref_pnt = getInnerNeighborReferenceQuadraturePoint<t_quadrature, t_i, t_j>(component_index, index);
            auto const nds = this->getCurrentCoordinates();
            for (lolita::integer j = 0; j < 3; ++j)
            {
                p(j) = FiniteElementGeometry::getShapeMappingEvaluation(nds.row(j), cpt_ref_pnt);
            }
            return p;
        }
        
        template<Basis t_basis, auto t_arg>
        lolita::matrix::Vector<lolita::real, basis::FiniteElementBasisTraits<t_basis, t_arg>::dim_>
        getBasisEvaluation(
            Point const & point
        )
        const
        {
            using t_FiniteElementBasisTraits = typename basis::FiniteElementBasisTraits<t_basis, t_arg>;
            using t_Implementation = t_FiniteElementBasisTraits::template Implementation<lolita2::geometry::FiniteElementGeometry, t_element, t_domain>;
            return static_cast<t_Implementation const *>(this)->getBasisEvaluation(point);
        }
        
        template<Basis t_basis, auto t_arg>
        lolita::matrix::Vector<lolita::real, basis::FiniteElementBasisTraits<t_basis, t_arg>::dim_>
        getBasisDerivative(
            Point const & point,
            lolita::integer derivative_direction
        )
        const
        {
            using t_FiniteElementBasisTraits = typename basis::FiniteElementBasisTraits<t_basis, t_arg>;
            using t_Implementation = t_FiniteElementBasisTraits::template Implementation<lolita2::geometry::FiniteElementGeometry, t_element, t_domain>;
            return static_cast<t_Implementation const *>(this)->getBasisDerivative(point, derivative_direction);
        }
        
        template<Basis t_basis>
        lolita::matrix::Vector<lolita::real, basis::FiniteElementBasisTraits<t_basis, t_element>::dim_>
        getBasisEvaluation(
            Point const & point
        )
        const
        {
            return getBasisEvaluation<t_basis, t_element>(point);
        }
        
        template<Basis t_basis>
        lolita::matrix::Vector<lolita::real, basis::FiniteElementBasisTraits<t_basis, t_element>::dim_>
        getBasisDerivative(
            Point const & point,
            lolita::integer derivative_direction
        )
        const
        {
            return getBasisDerivative<t_basis, t_element>(point, derivative_direction);
        }
        
        t_OuterNeighbors outer_neighbors_;
        
        t_InnerNeighbors inner_neighbors_;
        
        lolita::natural tag_;

        std::vector<std::shared_ptr<MeshDomain>> domains_;
        
        std::shared_ptr<Point> coordinates_;

    };

    template<Element t_element, Domain t_domain, Field t_field, Basis t_basis>
    struct FiniteElementUnknown
    {

        lolita::integer static constexpr num_basis_components_ = basis::FiniteElementBasisTraits<t_basis, t_element>::dim_;

        lolita::integer static constexpr num_field_components_ = FieldTraits<t_field, t_domain>::size__;

        lolita::integer static constexpr num_unknown_components_ = num_basis_components_ * num_field_components_;

        void
        setUnknown(
            std::shared_ptr<lolita::matrix::Vector<lolita::real>> & unknown
        )
        {
            unknown_ = unknown;
            index_ = unknown->size();
            unknown->resize(unknown->size() + num_unknown_components_);
        }

        lolita::matrix::Span<lolita::matrix::Vector<lolita::real, num_unknown_components_> const>
        getCoefficients()
        const
        {
            return lolita::matrix::Span<lolita::matrix::Vector<lolita::real, num_unknown_components_> const>(unknown_->data() + 1);
        }

        std::shared_ptr<lolita::matrix::Vector<lolita::real>> unknown_;

        lolita::integer index_;

    };

    template<Element t_element, Domain t_domain, auto t_arg>
    struct FiniteElement : FiniteElementGeometry<FiniteElement, t_element, t_domain, t_arg>
    {

        void
        make()
        {

        }

    };

    template<Element t_element, Domain t_domain, auto... t_args>
    struct FiniteElementHolder : FiniteElementGeometry<FiniteElementHolder, t_element, t_domain, t_args...>
    {

        using t_FiniteElements = std::tuple<std::shared_ptr<FiniteElement<t_element, t_domain, t_args>>...>;

        template<lolita::integer t_i>
        using t_FiniteElement = typename std::tuple_element_t<t_i, t_FiniteElements>::element_type;
          
        template<lolita::integer t_i>
        std::tuple_element_t<t_i, t_FiniteElements> const &
        getFiniteElement()
        const
        {
            return std::get<t_i>(finite_elements_);
        }
        
        template<lolita::integer t_i>
        std::tuple_element_t<t_i, t_FiniteElements> &
        getFiniteElement()
        {
            return std::get<t_i>(finite_elements_);
        }

        t_FiniteElements finite_elements_;

    };

}

#endif /* CC79BDC5_49DB_4A81_8A93_18ABD6551AF1 */
