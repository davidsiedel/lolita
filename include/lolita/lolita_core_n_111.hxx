#ifndef F548FA01_2D9F_4847_B854_824601F31371
#define F548FA01_2D9F_4847_B854_824601F31371

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_defs.hxx"
#include "lolita/lolita_core_n_00.hxx"
#include "lolita/lolita_core_n_0.hxx"

namespace lolita2::geometry
{

    namespace detail
    {

        template<template<Element, Domain, auto...> typename t_T, auto... t_args>
        struct ElementMapCollectionTraits
        {

        private:

            template<Element t_element, Domain t_domain>
            using t_Elements = std::map<std::basic_string<lolita::character>, std::shared_ptr<t_T<t_element, t_domain, t_args...>>>;

        public:

            template<Domain t_domain>
            using ElementMapCollection = ElementCollection<t_Elements, t_domain>;

        };

    }

    template<template<Element, Domain, auto...> typename t_T, Domain t_domain, auto... t_args>
    using ElementMapCollection = typename detail::ElementMapCollectionTraits<t_T, t_args...>::template ElementMapCollection<t_domain>;

    struct MeshDomain
    {

        MeshDomain(
            lolita::integer dim,
            std::basic_string<lolita::character> const & tag
        )
        :
        dim_(dim),
        tag_(tag)
        {}

        lolita::integer dim_;

        std::basic_string<lolita::character> tag_;

    };

    using MeshDomains = std::vector<std::shared_ptr<MeshDomain>>;

    template<template<Element, Domain, auto...> typename t_FiniteElement, Domain t_domain, auto... t_args>
    struct ElementSet
    {

        ElementSet()
        {}
        
        friend
        std::ostream &
        operator<<(
            std::ostream & os,
            ElementSet const & mesh
        )
        {
            auto print_element_components = [&] <Element t_element, lolita::integer t_i = 0, lolita::integer t_j = 0> (
                    auto const & elem_arg,
                    auto & self
            )
            mutable
            {
                if constexpr (!t_element.isNode())
                {
                    auto const constexpr _component = ElementTraits<t_element, t_domain>::template getComponent<t_i, t_j>();
                    for (auto const & c_ : elem_arg->template getComponents<t_i, t_j>())
                    {
                        os << "layer : " << t_i << " type : " << t_j << " <-- " << _component << " " << c_->hash() << std::endl;
                    }
                    if constexpr (t_j < ElementTraits<t_element, t_domain>::template getNumComponents<t_i>() - 1)
                    {
                        self.template operator()<t_element, t_i, t_j + 1>(elem_arg, self);
                    }
                    else if constexpr (t_i < ElementTraits<t_element, t_domain>::getNumComponents() - 1)
                    {
                        self.template operator()<t_element, t_i + 1, 0>(elem_arg, self);
                    }
                }
            };
            auto print_element_neighbours = [&] <Element t_element, lolita::integer t_i = 0, lolita::integer t_j = 0> (
                    auto const & elem_arg,
                    auto & self
            )
            mutable
            {
                auto const constexpr _neighbour = ElementTraits<t_element, t_domain>::template getNeighbour<t_i, t_j>();
                for (auto const & c_ : elem_arg->template getNeighbours<t_i, t_j>())
                {
                    if constexpr (!t_element.isNode() && t_i == 0)
                    {
                        os << "layer : " << t_i << " type : " << t_j << " <-> " << _neighbour << " " << c_->hash() << std::endl;
                    }
                    else
                    {
                        os << "layer : " << t_i << " type : " << t_j << " --> " << _neighbour << " " << c_->hash() << std::endl;
                    }
                }
                if constexpr (t_j < ElementTraits<t_element, t_domain>::template getNumNeighbours<t_i>() - 1)
                {
                    self.template operator()<t_element, t_i, t_j + 1>(elem_arg, self);
                }
                else if constexpr (t_i < ElementTraits<t_element, t_domain>::getNumNeighbours() - 1)
                {
                    self.template operator()<t_element, t_i + 1, 0>(elem_arg, self);
                }
            };
            auto print_elements2 = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                if constexpr (t_i == 0 && t_j == 0)
                {
                    os << "*** Elements : " << std::endl;
                }
                for (auto const & element : mesh.elements_.template getElements<t_i, t_j>())
                {
                    os << "* Element : " << t_element << " " << element.second->hash() << std::endl;
                    os << "domains : ";
                    for (auto const & domain : element.second->domains_)
                    {
                        os << domain->tag_ << " ";
                    }
                    os << std::endl;
                    print_element_components.template operator()<t_element>(element.second, print_element_components);
                    print_element_neighbours.template operator()<t_element>(element.second, print_element_neighbours);
                }                
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < DomainTraits<t_domain>::getNumElements() - 1)
                {
                    self.template operator()<t_i + 1, 0>(self);
                }
            };
            print_elements2(print_elements2);
            return os;
        }

        ElementMapCollection<t_FiniteElement, t_domain, t_args...> elements_;

        MeshDomains domains_;

    };

    template<Element t_element, Domain t_domain>
    struct MeshElement;
    
    template<Element t_element, Domain t_domain>
    requires(!t_element.isNode())
    struct MeshElement<t_element, t_domain>
    {
        
        template<lolita::integer t_i, lolita::integer t_j>
        static constexpr
        lolita::integer
        getComponentNodeConnection(
            lolita::integer i,
            lolita::integer j
        )
        {
            return std::get<t_j>(std::get<t_i>(ElementTraits<t_element, t_domain>::node_connectivity_))[i][j];
        }

        static
        std::basic_string<lolita::character>
        getHash(
            std::array<lolita::natural, t_element.num_nodes_> node_tags
        )
        {
            auto hash = std::basic_stringstream<lolita::character>();
            std::sort(std::execution::par_unseq, node_tags.begin(), node_tags.end());
            for (auto node_tag : node_tags)
            {
                hash << node_tag;
            }
            return hash.str();
        }

        MeshElement()
        :
        node_tags_()
        {}

        template<lolita::integer t_i, lolita::integer t_j, template<Element, Domain, auto...> typename t_FiniteElement, auto... t_args>
        void
        setFiniteElement(
            ElementSet<t_FiniteElement, t_domain, t_args...> & list,
            std::shared_ptr<t_FiniteElement<t_element, t_domain, t_args...>> & ptr_element
        )
        const
        {
            auto const constexpr t_component = ElementTraits<t_element, t_domain>::template getComponent<t_i, t_j>();
            auto const constexpr t_is_initialized = t_i == 0 && t_j == 0;
            auto const constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            auto const constexpr t_node_coordinates = ElementTraits<t_element, t_domain>::template getComponentCoordinates<Element::node()>();
            auto const constexpr t_component_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_component>();
            auto const constexpr t_neighbour_coordinates = ElementTraits<t_component, t_domain>::template getNeighbourCoordinates<t_element>();
            auto & components = list.elements_.template getElements<t_component_coordinates.dim_, t_component_coordinates.tag_>();
            for (lolita::integer i = 0; i < ElementTraits<t_element, t_domain>::template getNumComponents<t_i, t_j>(); ++i)
            {
                auto component_hash = std::basic_string<lolita::character>();
                if constexpr(!t_component.isNode())
                {
                    auto cpt = MeshElement<t_component, t_domain>();
                    for (lolita::integer j = 0; j < t_component.num_nodes_; ++j)
                    {
                        cpt.node_tags_[j] = node_tags_[getComponentNodeConnection<t_i, t_j>(i, j)];
                    }
                    component_hash = MeshElement<t_component, t_domain>::getHash(cpt.node_tags_);
                    if (!components.contains(component_hash))
                    {
                        using t_Component = t_FiniteElement<t_component, t_domain, t_args...>;
                        auto ptr_component = std::make_shared<t_Component>(t_Component());
                        cpt.template setFiniteElement<0, 0, t_FiniteElement, t_args...>(list, ptr_component);
                    }
                }
                else
                {
                    component_hash = MeshElement<t_component, t_domain>::getHash(node_tags_[getComponentNodeConnection<t_i, t_j>(i, 0)]);
                }
                ptr_element->template getComponents<t_i, t_j>()[i] = components[component_hash];
                components[component_hash]->template getNeighbours<t_neighbour_coordinates.dim_, t_neighbour_coordinates.tag_>().push_back(ptr_element);
            }
            if constexpr (t_j < ElementTraits<t_element, t_domain>::template getNumComponents<t_i>() - 1)
            {
                setFiniteElement<t_i, t_j + 1, t_FiniteElement, t_args...>(list, ptr_element);
            }
            else if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumComponents<>() - 1)
            {
                setFiniteElement<t_i + 1, 0, t_FiniteElement, t_args...>(list, ptr_element);
            }
            if constexpr (t_is_initialized)
            {
                auto const & nodes = ptr_element->template getComponents<t_node_coordinates.dim_, t_node_coordinates.tag_>();
                auto domains = std::set<std::shared_ptr<MeshDomain>>();
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
                ptr_element->tag_ = list.elements_.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>().size();
                ptr_element->domains_.assign(domains.begin(), domains.end());
                list.elements_.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>()[getHash(node_tags_)] = ptr_element;
            }
        }

        template<template<Element, Domain, auto...> typename t_FiniteElement, auto... t_args>
        void
        makeElement(
            ElementSet<t_FiniteElement, t_domain, t_args...> & list
        )
        const
        {
            auto ptr_element = std::make_shared<t_FiniteElement<t_element, t_domain, t_args...>>(t_FiniteElement<t_element, t_domain, t_args...>());
            setFiniteElement<0, 0, t_FiniteElement, t_args...>(list, ptr_element);
        }

        template<lolita::integer _k, template<Element, Domain, auto...> typename t_FiniteElement, auto... t_args>
        static
        void
        buildElementNN(
            std::shared_ptr<t_FiniteElement<t_element, t_domain, t_args...>> & ptr_element
        )
        {
            auto const constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            auto const constexpr t_node_coordinates = ElementTraits<t_element, t_domain>::template getComponentCoordinates<Element::node()>();
            auto const & element_nds = ptr_element->template getComponents<t_node_coordinates.dim_, t_node_coordinates.tag_>();
            for (auto i = 0; i < ElementTraits<t_element, t_domain>::template getNumComponents<t_node_coordinates.dim_, t_node_coordinates.tag_>(); ++i)
            {
                auto const & nde = element_nds[i];
                auto const & ngs = nde->template getNeighbours<t_element_coordinates.dim_ - 1, _k>();
                for (auto const & neighbour : ngs)
                {
                    if (((neighbour->tag_ != ptr_element->tag_) && _k == t_element_coordinates.tag_) || (_k != t_element_coordinates.tag_))
                    {
                        auto & element_ngs = ptr_element->template getNeighbours<0, _k>();
                        auto found = false;
                        for (auto const & ngb: element_ngs)
                        {
                            if (ngb->tag_ == neighbour->tag_)
                            {
                                found = true;
                                break;
                            }
                        }
                        if (!found)
                        {
                            element_ngs.push_back(neighbour);
                        }
                    }
                }
            }
            if constexpr (_k < DomainTraits<t_domain>::template getNumElements<t_element_coordinates.dim_>() - 1)
            {
                buildElementNN<_k + 1, t_FiniteElement, t_args...>(ptr_element);
            }
        }

        template<template<Element, Domain, auto...> typename t_FiniteElement, auto... t_args>
        static
        void
        buildElementN(
            std::shared_ptr<t_FiniteElement<t_element, t_domain, t_args...>> & ptr_element
        )
        {
            buildElementNN<0, t_FiniteElement, t_args...>(ptr_element);
        }

        std::array<lolita::natural, t_element.num_nodes_> node_tags_;

    };

    template<Element t_element, Domain t_domain>
    requires(t_element.isNode())
    struct MeshElement<t_element, t_domain>
    {

        static
        std::basic_string<lolita::character>
        getHash(
            lolita::natural tag
        )
        {
            return std::to_string(tag);
        }

        MeshElement()
        :
        tag_(),
        coordinates_(std::make_shared<Point>(Point::Zero())),
        domains_()
        {}

        template<template<Element, Domain, auto...> typename t_FiniteElement, auto... t_args>
        void
        setFiniteElement(
            ElementSet<t_FiniteElement, t_domain, t_args...> & list,
            std::shared_ptr<t_FiniteElement<t_element, t_domain, t_args...>> & ptr_element
        )
        const
        {
            auto constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            ptr_element->tag_ = tag_;
            ptr_element->domains_ = domains_;
            ptr_element->coordinates_ = coordinates_;
            list.elements_.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>()[getHash(tag_)] = ptr_element;
        }

        template<template<Element, Domain, auto...> typename t_FiniteElement, auto... t_args>
        void
        makeElement(
            ElementSet<t_FiniteElement, t_domain, t_args...> & list
        )
        const
        {
            auto ptr_element = std::make_shared<t_FiniteElement<t_element, t_domain, t_args...>>(t_FiniteElement<t_element, t_domain, t_args...>());
            setFiniteElement<t_FiniteElement, t_args...>(list, ptr_element);
        }

        template<template<Element, Domain, auto...> typename t_FiniteElement, auto... t_args>
        static
        void
        buildElementN(
            std::shared_ptr<t_FiniteElement<t_element, t_domain, t_args...>> & ptr_element
        )
        {}
        
        lolita::natural tag_;
        
        std::shared_ptr<Point> coordinates_;

        std::vector<std::shared_ptr<MeshDomain>> domains_;

    };
    
    template<Domain t_domain>
    struct MeshData
    {

        template<template<Element, Domain, auto...> typename t_FiniteElement, auto... t_args>
        ElementSet<t_FiniteElement, t_domain, t_args...>
        make()
        const
        {
            auto mmm = ElementSet<t_FiniteElement, t_domain, t_args...>();
            auto mkl = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & element : elements_.template getElements<t_i, t_j>())
                {
                    element.second->template makeElement<t_FiniteElement, t_args...>(mmm);
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_i, t_j + 1u>(self);
                }
                else if constexpr (t_i < DomainTraits<t_domain>::getNumElements() - 1)
                {
                    self.template operator()<t_i + 1u, 0u>(self);
                }
            };
            mkl(mkl);
            auto make_elements = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto const constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                for (auto & element : mmm.elements_.template getElements<t_i, t_j>())
                {
                    MeshElement<t_element, t_domain>::template buildElementN<t_FiniteElement, t_args...>(element.second);
                }
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < DomainTraits<t_domain>::getNumElements() - 1)
                {
                    self.template operator()<t_i + 1, 0>(self);
                }
            };
            make_elements(make_elements);
            mmm.domains_ = domains_;
            return mmm;
        }

        MeshData()
        {}

        ElementMapCollection<MeshElement, t_domain> elements_;

        MeshDomains domains_;

    };    
    
    struct MeshFileParserBase
    {

        MeshFileParserBase(
            std::basic_string_view<lolita::character> str
        )
        :
        file_(str)
        {}

        lolita::utility::File file_;
        
    };

    struct GmshFileParser : MeshFileParserBase
    {

    private:

        static constexpr
        Element
        getElemType(
                lolita::integer tag
        )
        {
            if (tag == 15) return Element::node();
            else if (tag == 01) return Element::segment(1);
            else if (tag == 02) return Element::triangle(1);
            else if (tag == 03) return Element::quadrangle(1);
            else return Element::node();
        }

        struct PhysicalEntity
        {

            static
            std::basic_string<lolita::character>
            getHash(
                lolita::integer tag
            )
            {
                return std::to_string(tag);
            }

            PhysicalEntity(
                std::basic_string<lolita::character> const & name,
                lolita::integer dim
            )
            :
            mesh_domain_(std::make_shared<MeshDomain>(dim, name))
            {}

            std::shared_ptr<MeshDomain> mesh_domain_;

        };

        struct GeometricalEntity
        {

            static
            std::basic_string<lolita::character>
            getHash(
                lolita::integer dim,
                lolita::integer tag
            )
            {
                auto hash = std::basic_stringstream<lolita::character>();
                hash << dim;
                hash << tag;
                return hash.str();
            }

            GeometricalEntity()
            :
            physical_entities_(),
            bounding_entities_()
            {}

            void
            addPhysicalEntity(
                std::shared_ptr<PhysicalEntity> const & physical_entity
            )
            {
                physical_entities_.insert(physical_entity);
                for (auto & bounding_entity : bounding_entities_)
                {
                    bounding_entity->addPhysicalEntity(physical_entity);
                }
            }

            void
            addBoundingEntity(
                std::shared_ptr<GeometricalEntity> & bounding_entity
            )
            {
                bounding_entities_.insert(bounding_entity);
                for (auto const & physical_entity : physical_entities_)
                {
                    bounding_entity->addPhysicalEntity(physical_entity);
                }
            }

            std::set<std::shared_ptr<PhysicalEntity>> physical_entities_;

            std::set<std::shared_ptr<GeometricalEntity>> bounding_entities_;

        };
        
        void
        setPhysicalEntities()
        {
            auto const & file_lines = this->file_.lines_;
            auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$PhysicalNames"));
            auto offset = 1;
            auto line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
            auto num_physical_names = lolita::integer();
            line_stream >> num_physical_names;
            offset += 1;
            for (lolita::integer i = 0; i < num_physical_names; ++i)
            {
                line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                auto dim = lolita::integer();
                auto name = std::basic_string<lolita::character>();
                auto physical_entity_tag = lolita::integer();
                line_stream >> dim >> physical_entity_tag >> name;
                lolita::utility::removeCharacter(name, '"');
                auto hash = PhysicalEntity::getHash(physical_entity_tag);
                physical_entities_[hash] = std::make_shared<PhysicalEntity>(PhysicalEntity(name, dim));
                offset += 1;
            }
        }

        void
        setGeometricalEntities()
        {
            auto const & file_lines = this->file_.lines_;
            auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Entities"));
            auto offset = 1;
            auto line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
            auto num_points = lolita::integer();
            auto num_curves = lolita::integer();
            auto num_surfaces = lolita::integer();
            auto num_volumes = lolita::integer();
            line_stream >> num_points >> num_curves >> num_surfaces >> num_volumes;
            auto num_domains = std::array<lolita::integer, 4>{num_points, num_curves, num_surfaces, num_volumes};
            offset += 1;
            for (lolita::integer i = 0; i < 4; ++i)
            {
                for (lolita::integer j = 0; j < num_domains[i]; ++j)
                {
                    line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                    auto tag = lolita::integer();
                    line_stream >> tag;
                    auto geometrical_entity_hash = GeometricalEntity::getHash(i, tag);
                    geometrical_entities_[geometrical_entity_hash] = std::make_shared<GeometricalEntity>(GeometricalEntity());
                    if (i == 0)
                    {
                        for (lolita::integer k = 0; k < 3; ++k)
                        {
                            auto a = lolita::real();
                            line_stream >> a;
                        }
                    }
                    else
                    {
                        for (lolita::integer k = 0; k < 6; ++k)
                        {
                            auto a = lolita::real();
                            line_stream >> a;
                        }
                    }
                    auto num_physical_entities = lolita::integer();
                    line_stream >> num_physical_entities;
                    for (lolita::integer k = 0; k < num_physical_entities; ++k)
                    {
                        auto physical_entity_tag = lolita::integer();
                        line_stream >> physical_entity_tag;
                        auto physical_entity_hash = PhysicalEntity::getHash(physical_entity_tag);
                        geometrical_entities_[geometrical_entity_hash]->addPhysicalEntity(physical_entities_[physical_entity_hash]);
                    }
                    if (i > 0)
                    {
                        auto num_bounding_entities = lolita::integer();
                        line_stream >> num_bounding_entities;
                        for (lolita::integer k = 0; k < num_bounding_entities; ++k)
                        {
                            auto bounding_entity_tag = lolita::integer();
                            line_stream >> bounding_entity_tag;
                            auto bounding_entity_hash = GeometricalEntity::getHash(i - 1, std::abs(bounding_entity_tag));
                            geometrical_entities_[geometrical_entity_hash]->addBoundingEntity(geometrical_entities_[bounding_entity_hash]);
                        }
                    }
                    offset += 1;
                }
            }
        }

    public:

        GmshFileParser(
            std::basic_string_view<lolita::character> str
        )
        :
        MeshFileParserBase(str)
        {
            setPhysicalEntities();
            setGeometricalEntities();
        }

        template<Domain t_domain>
        void
        setMeshDomains(
            MeshData<t_domain> & mesh_data
        )
        const
        {
            for (auto const & physical_entity : physical_entities_)
            {
                mesh_data.domains_.push_back(physical_entity.second->mesh_domain_);
            }
        }

        template<Domain t_domain, Element t_element>
        void
        setMeshNodes(
            MeshData<t_domain> & mesh_data
        )
        const
        requires(t_element.isNode())
        {
            auto constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            auto & elements = mesh_data.elements_.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>();
            auto const & file_lines = this->file_.lines_;
            auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Nodes"));
            auto offset = 1;
            auto line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
            auto num_entity_block = lolita::integer();
            auto num_nodes = lolita::integer();
            auto min_node_tag = lolita::integer();
            auto max_node_tag = lolita::integer();
            line_stream >> num_entity_block >> num_nodes >> min_node_tag >> max_node_tag;
            offset += 1;
            for (lolita::integer i = 0; i < num_entity_block; ++i)
            {
                line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                auto entity_dim = lolita::integer();
                auto entity_tag = lolita::integer();
                auto parametric = lolita::integer();
                auto num_nodes_in_block = lolita::integer();
                line_stream >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
                auto const entity_hash = GeometricalEntity::getHash(entity_dim, entity_tag);
                offset += 1;
                for (lolita::integer j = 0; j < num_nodes_in_block; ++j)
                {
                    auto mesh_element = std::make_shared<MeshElement<t_element, t_domain>>(MeshElement<t_element, t_domain>());
                    auto tag = lolita::natural();
                    auto coordinates_ = Point();
                    auto domains_ = std::vector<std::shared_ptr<std::basic_string<lolita::character>>>();
                    line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                    line_stream >> tag;
                    mesh_element->tag_ = tag - 1;
                    line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset + num_nodes_in_block]);
                    for (lolita::integer k = 0; k < 3; ++k)
                    {
                        line_stream >> mesh_element->coordinates_->operator()(k);
                    }
                    for (auto const & physical_entity : geometrical_entities_.at(entity_hash)->physical_entities_)
                    {
                        mesh_element->domains_.push_back(physical_entity->mesh_domain_);
                    }
                    auto const element_hash = MeshElement<t_element, t_domain>::getHash(mesh_element->tag_);
                    elements[element_hash] = mesh_element;
                    offset += 1;
                }
                offset += num_nodes_in_block;
            }
        }
        
        template<Domain t_domain, Element t_element>
        void
        setMeshNodes(
            MeshData<t_domain> & mesh_data
        )
        const
        requires(!t_element.isNode())
        {
            auto constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            auto & elements = mesh_data.elements_.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>();
            auto const & file_lines = this->file_.lines_;
            auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Elements"));
            auto offset = 1;
            auto line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
            auto num_entity_blocks = 0;
            auto num_elements = 0;
            auto min_element_tag = 0;
            auto max_element_tag = 0;
            line_stream >> num_entity_blocks >> num_elements >> min_element_tag >> max_element_tag;
            offset += 1;
            for (lolita::integer i = 0; i < num_entity_blocks; ++i)
            {
                line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                auto entity_dim = 0;
                auto entity_tag = 0;
                auto element_type_tag = 0;
                auto num_elements_in_block = 0;
                line_stream >> entity_dim >> entity_tag >> element_type_tag >> num_elements_in_block;
                offset += 1;
                auto const element = getElemType(element_type_tag);
                if (element == t_element)
                {
                    for (lolita::integer j = 0; j < num_elements_in_block; ++j)
                    {
                        auto mesh_element = std::make_shared<MeshElement<t_element, t_domain>>(MeshElement<t_element, t_domain>());
                        line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                        auto tag = lolita::integer();
                        line_stream >> tag;
                        for (lolita::integer k = 0; k < t_element.num_nodes_; ++k)
                        {
                            auto node_tag = lolita::integer();
                            line_stream >> node_tag;
                            mesh_element->node_tags_[k] = node_tag - 1;
                        }
                        auto const element_hash = MeshElement<t_element, t_domain>::getHash(mesh_element->node_tags_);
                        elements[element_hash] = mesh_element;
                        offset += 1;
                    }
                }
                else
                {
                    offset += num_elements_in_block;
                }
            }
        }
        
    private:
        
        std::map<std::basic_string<lolita::character>, std::shared_ptr<PhysicalEntity>> physical_entities_;

        std::map<std::basic_string<lolita::character>, std::shared_ptr<GeometricalEntity>> geometrical_entities_;

    };
    
    template<template<Element, Domain, auto...> typename t_FiniteElement, Domain t_domain, auto... t_args>
    struct MeshFileParser
    {

        static
        ElementSet<t_FiniteElement, t_domain, t_args...>
        makeIt(
            std::basic_string_view<lolita::character> str
        )
        {
            auto gmsh = GmshFileParser(str);
            auto mshh = MeshData<t_domain>();
            gmsh.setMeshDomains(mshh);
            auto mkl = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (
                auto & self
            )
            mutable
            {
                auto const constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                gmsh.template setMeshNodes<t_domain, t_element>(mshh);
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1)
                {
                    self.template operator()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < DomainTraits<t_domain>::getNumElements() - 1)
                {
                    self.template operator()<t_i + 1, 0>(self);
                }
            };
            mkl(mkl);
            return mshh.template make<t_FiniteElement, t_args...>();
        }

    };
    
}


#endif /* F548FA01_2D9F_4847_B854_824601F31371 */
