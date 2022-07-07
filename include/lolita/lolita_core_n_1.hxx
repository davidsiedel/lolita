#ifndef A0078D95_FD87_4075_93DB_523A0055E16A
#define A0078D95_FD87_4075_93DB_523A0055E16A

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_defs.hxx"
#include "lolita/lolita_core_n_00.hxx"
#include "lolita/lolita_core_n_0.hxx"

namespace lolita2::geometry
{

    template<template<Element, Domain, auto...> typename t_FiniteElement, Element t_element, Domain t_domain, auto... t_args>
    struct FiniteElementConnectivity
    {

    private:
    
        using t_ElementTraits = ElementTraits<t_element, t_domain>;
        
        template<Element t__element, Domain t__domain>
        using t_ElementPointer = std::shared_ptr<FiniteElementConnectivity<t_FiniteElement, t__element, t__domain, t_args...>>;

    public:
    
        using Components = typename t_ElementTraits::template InnerConnectivity<t_ElementPointer>;
        
        using Neighbours = typename t_ElementTraits::template OuterConnectivity<t_ElementPointer>;
        
        lolita::natural tag_;

        // std::basic_string_view<lolita::character> hash_;
        
        Neighbours neighbours_;
        
        Components components_;
        
        std::vector<std::shared_ptr<std::basic_string<lolita::character>>> domains_;
        
        std::shared_ptr<Point> coordinates_;
        
        lolita::boolean
        operator==(
            FiniteElementConnectivity const & other
        )
        const = default;
        
        lolita::boolean
        operator!=(
            FiniteElementConnectivity const & other
        )
        const = default;
        
        std::basic_string<lolita::character>
        hash()
        const
        requires(t_element.isNode())
        {
            return std::to_string(this->tag_);
        }
        
        std::basic_string<lolita::character>
        hash()
        const
        {
            std::basic_stringstream<lolita::character> hash;
            auto const & nodes = getComponents<t_element.dim_ - 1, 0>();
            for (auto const & node : nodes)
            {
                hash << node->hash();
            }
            return hash.str();
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        static constexpr
        lolita::integer
        getComponentNodeConnection(
            lolita::integer i,
            lolita::integer j
        )
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(ElementTraits<t_element, t_domain>::node_connectivity_))[i][j];
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, Neighbours>> &
        getNeighbours()
        {
            return std::get<t_j>(std::get<t_i>(neighbours_));
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, Neighbours>> const &
        getNeighbours()
        const
        {
            return std::get<t_j>(std::get<t_i>(neighbours_));
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, Components>> &
        getComponents()
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(components_));
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, Components>> const &
        getComponents()
        const
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(components_));
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        lolita::integer
        getComponentIndex(
            lolita::integer i
        )
        const
        requires(!t_element.isNode())
        {
            auto constexpr t_component = t_ElementTraits::template getComponent<t_i, t_j>();
            using t_NeighbourTraits = ElementTraits<t_component, t_domain>;
            auto constexpr t_coordinates = t_NeighbourTraits::template getNeighbourCoordinates<t_element>();
            auto const & items = getComponents<t_i, t_j>()[i]->template getNeighbours<t_coordinates.dim_, t_coordinates.tag_>();
            auto is_equal = [&] (t_ElementPointer<t_element, t_domain> const & ptr_element)
            {
                return * ptr_element == * this;
            };
            return std::distance(items.begin(), std::find_if(items.begin(), items.end(), is_equal));
        }
        
        template<lolita::integer t_i, lolita::integer t_j>
        lolita::integer
        getComponentOrientation(
            lolita::integer i
        )
        const
        requires(!t_element.isNode())
        {
            return getComponentIndex<t_i, t_j>(i) == 0 ? 1 : -1;
        }

    };

    template<Element t_element, Domain t_domain>
    struct FiniteElementInput : FiniteElementConnectivity<FiniteElementInput, t_element, t_domain>
    {

        static
        std::basic_string<lolita::character>
        getElementHash(
            std::array<lolita::integer, t_element.num_nodes_> node_tags
        )
        requires(!t_element.isNode())
        {
            auto element_hash = std::basic_stringstream<lolita::character>();
            std::sort(std::execution::par_unseq, node_tags.begin(), node_tags.end());
            for (auto node_tag : node_tags)
            {
                element_hash << node_tag;
            }
            return element_hash.str();
        }

        static
        std::basic_string<lolita::character>
        getElementHash(
            lolita::natural node_tag
        )
        requires(t_element.isNode())
        {
            return std::to_string(node_tag);
        }
        
        template<lolita::integer t_i = 0, lolita::integer t_j = 0>
        static
        void
        makeElement(
            auto & list,
            std::shared_ptr<FiniteElementInput> & ptr_element,
            std::array<lolita::integer, t_element.num_nodes_> node_tags
        )
        requires(!t_element.isNode())
        {
            auto const constexpr _component = ElementTraits<t_element, t_domain>::template getComponent<t_i, t_j>();
            auto const constexpr _is_initialized = t_i == 0 && t_j == 0;
            auto const constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            auto const constexpr t_node_coordinates = ElementTraits<t_element, t_domain>::template getComponentCoordinates<Element::node()>();
            auto const constexpr _component_coordinates = DomainTraits<t_domain>::template getElementCoordinates<_component>();
            auto const constexpr _neighbour_coordinates = ElementTraits<_component, t_domain>::template getNeighbourCoordinates<t_element>();
            auto & components = list.template getElements<_component_coordinates.dim_, _component_coordinates.tag_>();
            auto & element_component_array = ptr_element->template getComponents<t_i, t_j>();
            for (auto i = 0; i < element_component_array.size(); ++i)
            {
                auto component_hash = std::basic_string<lolita::character>();
                if constexpr(!_component.isNode())
                {
                    auto component_node_tags = std::array<lolita::integer, _component.num_nodes_>();
                    for (lolita::integer j = 0; j < _component.num_nodes_; ++j)
                    {
                        auto const k = FiniteElementInput::template getComponentNodeConnection<t_i, t_j>(i, j);
                        component_node_tags[j] = node_tags[k];
                    }
                    component_hash = lolita2::geometry::FiniteElementInput<_component, t_domain>::getElementHash(component_node_tags);
                    if (!components.contains(component_hash))
                    {
                        using t_Component = lolita2::geometry::FiniteElementInput<_component, t_domain>;
                        auto ptr_component = std::make_shared<t_Component>(t_Component());
                        lolita2::geometry::FiniteElementInput<_component, t_domain>::template makeElement<0, 0>(ptr_component, component_node_tags);
                    }
                }
                else
                {
                    component_hash = std::to_string(node_tags[FiniteElementInput::template getComponentNodeConnection<t_i, t_j>(i, 0)]);
                }
                element_component_array[i] = components[component_hash];
                components[component_hash]->template getNeighbours<_neighbour_coordinates.dim_, _neighbour_coordinates.tag_>().push_back(ptr_element);
            }
            if constexpr (t_j < ElementTraits<t_element, t_domain>::template getNumComponents<t_i>() - 1)
            {
                makeElement<t_i, t_j + 1u>(ptr_element, node_tags);
            }
            else if constexpr (t_i < ElementTraits<t_element, t_domain>::getNumComponents() - 1)
            {
                makeElement<t_i + 1u, 0u>(ptr_element, node_tags);
            }
            if constexpr (_is_initialized)
            {
                auto const & nodes = ptr_element->template getComponents<t_node_coordinates.dim_, t_node_coordinates.tag_>();
                auto domains = std::set<std::shared_ptr<std::basic_string<lolita::character>>>();
                for (auto const & domain : nodes[0]->domains_)
                {
                    auto has_domain = true;
                    for (lolita::integer j = 1; j < t_element.num_nodes_; ++j)
                    {
                        has_domain = std::find(nodes[j]->domains_.begin(), nodes[j]->domains_.end(), domain) != nodes[j]->domains_.end();
                    }
                    if (has_domain)
                    {
                        domains.insert(domain);
                    }
                }
                ptr_element->tag_ = list.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>().size();
                ptr_element->domains_.assign(domains.begin(), domains.end());
                list.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>()[getElementHash(node_tags)] = ptr_element;
            }
        }

        static
        void
        makeElement(
            auto & list,
            std::shared_ptr<FiniteElementInput> & ptr_element,
            lolita::natural && tag,
            Point && coordinates,
            std::vector<std::shared_ptr<std::basic_string<lolita::character>>> && domains
        )
        requires(t_element.isNode())
        {
            auto const constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            ptr_element->tag_ = std::forward<lolita::natural>(tag);
            ptr_element->domains_ = std::forward<std::vector<std::shared_ptr<std::basic_string<lolita::character>>>>(domains);
            ptr_element->coordinates_ = std::make_shared<Point>(std::forward<Point>(coordinates));
            list.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>()[std::to_string(tag)] = ptr_element;
        }

    };

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
        
        template<Element t_element, Basis t_basis>
        struct FiniteElementBasisTraits;

        template<Element t_element, Basis t_basis>
        requires(t_basis.isMonomial())
        struct FiniteElementBasisTraits<t_element, t_basis>
        {
            
            lolita::integer static constexpr dim_ = lolita::numerics::binomial(t_element.dim_ + t_basis.ord_, t_element.dim_);
            
            template<template<Element, Domain, auto...> typename t_FiniteElement, Domain t_domain, auto... t_args>
            struct Implementation : FiniteElementConnectivity<t_FiniteElement, t_element, t_domain, t_args...>
            {

            private:
            
                static constexpr
                std::array<std::array<lolita::integer, 3>, dim_>
                getExponents()
                {
                    auto exponents = std::array<std::array<lolita::integer, 3>, dim_>();
                    auto row = lolita::integer(0);
                    if constexpr (t_element.dim_ == 0)
                    {
                        exponents[row][0] = 0;
                        exponents[row][1] = 0;
                        exponents[row][2] = 0;
                    }
                    else if constexpr (t_element.dim_ == 1)
                    {
                        for (lolita::integer i = 0; i < t_basis.ord_ + 1; ++i)
                        {
                            exponents[row][0] = i;
                            exponents[row][1] = 0;
                            exponents[row][2] = 0;
                            row += 1;
                        }
                    }
                    else if constexpr (t_element.dim_ == 2)
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
                    else if constexpr (t_element.dim_ == 3)
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
                        for (lolita::integer j = 0; j < t_element.dim_; ++j)
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
                        for (lolita::integer j = 0; j < t_element.dim_; ++j)
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

    }

    template<Element t_element, Domain t_domain>
    struct FiniteElementGeometry : FiniteElementConnectivity<FiniteElementGeometry, t_element, t_domain>
    {

    private:
    
        using t_ElementTraits = ElementTraits<t_element, t_domain>;

    public:
    
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
            for (auto const & node : this->template getComponents<t_element.dim_ - 1, 0>())
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
            using _ReferenceCoordinates = lolita::matrix::Span<lolita::matrix::Matrix<lolita::real, 3, t_element.num_nodes_, lolita::matrix::col_major> const>;
            return _ReferenceCoordinates(ElementTraits<t_element, t_domain>::reference_nodes_.begin()->begin());
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
                // auto constexpr _quadrature = Quadrature::gauss(4);
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
        getComponentReferenceQuadratureWeight(
            lolita::integer index
        )
        requires(!t_element.isNode())
        {
            auto const constexpr _component = t_ElementTraits::template getComponent<_i, _j>();
            using ComponentGeometry = FiniteElementGeometry<_component, t_domain>;
            return ComponentGeometry::template getReferenceQuadratureWeight<t_quadrature>(index);
        }
        
        template<Quadrature t_quadrature, lolita::integer _i, lolita::integer _j>
        static
        Point
        getComponentReferenceQuadraturePoint(
            lolita::integer component_index,
            lolita::integer index
        )
        requires(!t_element.isNode())
        {
            auto const constexpr _component = t_ElementTraits ::template getComponent<_i, _j>();
            auto p = Point();
            using ComponentGeometry = FiniteElementGeometry<_component, t_domain>;
            auto const & elt_reference_nodes = ElementTraits<t_element, t_domain>::reference_nodes_;
            for (lolita::integer i = 0; i < 3; ++i)
            {
                auto cpt_coordinates = lolita::matrix::Vector<lolita::real, _component.num_nodes_>();
                for (lolita::integer j = 0; j < _component.num_nodes_; ++j)
                {
                    auto const node_tag = getComponentNodeConnection<_i, _j>(component_index, j);//.get(component_index).get(j);
                    cpt_coordinates(j) = elt_reference_nodes[node_tag][i];
                }
                auto cpt_reference_point = ComponentGeometry::template getReferenceQuadraturePoint<t_quadrature>(index);
                p(i) = ComponentGeometry::getShapeMappingEvaluation(cpt_coordinates, cpt_reference_point);
            }
            return p;
        }
        
        template<Quadrature t_quadrature, lolita::integer _i, lolita::integer _j>
        lolita::real
        getComponentCurrentQuadratureWeight(
            lolita::integer component_index,
            lolita::integer index
        )
        const
        requires(!t_element.isNode())
        {
            auto const & cmp =  this->template getComponents<_i, _j>()[component_index];//.template get<I>().template get<J>().get(component_index).get();
            return cmp->template getCurrentQuadratureWeight<t_quadrature>(index);
        }
        
        template<Quadrature t_quadrature, lolita::integer _i, lolita::integer _j>
        Point
        getComponentCurrentQuadraturePoint(
            lolita::integer component_index,
            lolita::integer index
        )
        const
        requires(!t_element.isNode())
        {
            auto p = Point();
            auto const cpt_ref_pnt = getComponentReferenceQuadraturePoint<t_quadrature, _i, _j>(component_index, index);
            auto const nds = this->getCurrentCoordinates();
            for (lolita::integer j = 0; j < 3; ++j)
            {
                p(j) = FiniteElementGeometry::getShapeMappingEvaluation(nds.row(j), cpt_ref_pnt);
            }
            return p;
        }
        
        template<Basis t_basis>
        lolita::matrix::Vector<lolita::real, basis::FiniteElementBasisTraits<t_element, t_basis>::dim_>
        getBasisEvaluation(
            Point const & point
        )
        const
        {
            using t_FiniteElementBasisTraits = typename basis::FiniteElementBasisTraits<t_element, t_basis>;
            using t_Implementation = t_FiniteElementBasisTraits::template Implementation<lolita2::geometry::FiniteElementGeometry, t_domain>;
            return static_cast<t_Implementation const *>(this)->getBasisEvaluation(point);
        }
        
        template<Basis t_basis>
        lolita::matrix::Vector<lolita::real, basis::FiniteElementBasisTraits<t_element, t_basis>::dim_>
        getBasisDerivative(
            Point const & point,
            lolita::integer derivative_direction
        )
        const
        {
            using t_FiniteElementBasisTraits = typename basis::FiniteElementBasisTraits<t_element, t_basis>;
            using t_Implementation = t_FiniteElementBasisTraits::template Implementation<lolita2::geometry::FiniteElementGeometry, t_domain>;
            return static_cast<t_Implementation const *>(this)->getBasisDerivative(point, derivative_direction);
        }

    };

}

#endif /* A0078D95_FD87_4075_93DB_523A0055E16A */
