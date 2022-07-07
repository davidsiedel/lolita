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

    template<template<Element, Domain, auto...> typename t_FiniteElement, Domain t_domain, auto... t_args>
    struct MeshSetBase2
    {

    private:

        template<template<Element, Domain, auto...> typename t__FiniteElement, auto... t__args>
        struct FiniteElementCollectionTraits
        {

            template<Element t_element, Domain t__domain>
            using t_FiniteElements = std::map<std::basic_string<lolita::character>, std::shared_ptr<t__FiniteElement<t_element, t__domain, t__args...>>>;

            using t_FiniteElementCollection = ElementCollection<t_FiniteElements, t_domain>;

        };

        template<template<Element, Domain, auto...> typename t__FiniteElement, auto... t__args>
        using t_FiniteElements = typename FiniteElementCollectionTraits<t__FiniteElement, t__args...>::t_FiniteElementCollection;

    public:

        MeshSetBase2()
        {}

        using FiniteElements = typename FiniteElementCollectionTraits<t_FiniteElement, t_args...>::t_FiniteElementCollection;

        using Domains = std::vector<std::shared_ptr<MeshDomain>>;
        
        friend
        std::ostream &
        operator<<(
            std::ostream & os,
            MeshSetBase2 const & mesh
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

        FiniteElements elements_;

        Domains domains_;

    };

    template<Element t_element>
    struct MeshElement;
    
    template<Element t_element>
    requires(!t_element.isNode())
    struct MeshElement<t_element>
    {
        
        template<Domain t_domain, lolita::integer t_i, lolita::integer t_j>
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

        template<lolita::integer t_i, lolita::integer t_j, template<Element, Domain, auto...> typename t_FiniteElement, Domain t_domain, auto... t_args>
        void
        setFiniteElement(
            MeshSetBase2<t_FiniteElement, t_domain, t_args...> & list,
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
                    auto cpt = MeshElement<t_component>();
                    for (lolita::integer j = 0; j < t_component.num_nodes_; ++j)
                    {
                        cpt.node_tags_[j] = node_tags_[getComponentNodeConnection<t_domain, t_i, t_j>(i, j)];
                    }
                    component_hash = MeshElement<t_component>::getHash(cpt.node_tags_);
                    if (!components.contains(component_hash))
                    {
                        using t_Component = t_FiniteElement<t_component, t_domain, t_args...>;
                        auto ptr_component = std::make_shared<t_Component>(t_Component());
                        cpt.template setFiniteElement<0, 0, t_FiniteElement, t_domain, t_args...>(list, ptr_component);
                    }
                }
                else
                {
                    component_hash = MeshElement<t_component>::getHash(node_tags_[getComponentNodeConnection<t_domain, t_i, t_j>(i, 0)]);
                }
                ptr_element->template getComponents<t_i, t_j>()[i] = components[component_hash];
                components[component_hash]->template getNeighbours<t_neighbour_coordinates.dim_, t_neighbour_coordinates.tag_>().push_back(ptr_element);
            }
            if constexpr (t_j < ElementTraits<t_element, t_domain>::template getNumComponents<t_i>() - 1)
            {
                setFiniteElement<t_i, t_j + 1, t_FiniteElement, t_domain, t_args...>(list, ptr_element);
            }
            else if constexpr (t_i < ElementTraits<t_element, t_domain>::template getNumComponents<>() - 1)
            {
                setFiniteElement<t_i + 1, 0, t_FiniteElement, t_domain, t_args...>(list, ptr_element);
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

        template<template<Element, Domain, auto...> typename t_FiniteElement, Domain t_domain, auto... t_args>
        void
        makeElement(
            MeshSetBase2<t_FiniteElement, t_domain, t_args...> & list
        )
        const
        {
            auto ptr_element = std::make_shared<t_FiniteElement<t_element, t_domain, t_args...>>(t_FiniteElement<t_element, t_domain, t_args...>());
            setFiniteElement<0, 0, t_FiniteElement, t_domain, t_args...>(list, ptr_element);
        }

        std::array<lolita::natural, t_element.num_nodes_> node_tags_;

    };

    template<Element t_element>
    requires(t_element.isNode())
    struct MeshElement<t_element>
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

        template<template<Element, Domain, auto...> typename t_FiniteElement, Domain t_domain, auto... t_args>
        void
        setFiniteElement(
            MeshSetBase2<t_FiniteElement, t_domain, t_args...> & list,
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

        template<template<Element, Domain, auto...> typename t_FiniteElement, Domain t_domain, auto... t_args>
        void
        makeElement(
            MeshSetBase2<t_FiniteElement, t_domain, t_args...> & list
        )
        const
        {
            auto ptr_element = std::make_shared<t_FiniteElement<t_element, t_domain, t_args...>>(t_FiniteElement<t_element, t_domain, t_args...>());
            setFiniteElement<t_FiniteElement, t_domain, t_args...>(list, ptr_element);
        }        
        
        lolita::natural tag_;
        
        std::shared_ptr<Point> coordinates_;

        std::vector<std::shared_ptr<MeshDomain>> domains_;

    };
    
    template<Domain t_domain>
    struct MeshData
    {

    private:

        template<Element t_element, Domain>
        using t_MeshElements = std::vector<MeshElement<t_element>>;

    public:

        template<template<Element, Domain, auto...> typename t_FiniteElement, auto... t_args>
        MeshSetBase2<t_FiniteElement, t_domain, t_args...>
        make()
        const
        {
            auto mmm = MeshSetBase2<t_FiniteElement, t_domain, t_args...>();
            auto mkl = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (
                auto & self
            )
            mutable
            {
                for (auto const & element : elements_.template getElements<t_i, t_j>())
                {
                    auto k = 9;
                    element.template makeElement<t_FiniteElement, t_domain, t_args...>(mmm);
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
            mmm.domains_ = domains_;
            return mmm;
        }

        MeshData()
        {}

        using MeshElements = ElementCollection<t_MeshElements, t_domain>;

        using MeshDomains = std::vector<std::shared_ptr<MeshDomain>>;

        MeshElements elements_;

        MeshDomains domains_;

    };    
    
    template<Domain t_domain>
    struct MeshFileBase
    {

        MeshFileBase(
            std::basic_string_view<lolita::character> str
        )
        :
        file_(str)
        {}

        lolita::utility::File file_;
        
    };

    struct NewGmsh
    {

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

            PhysicalEntity(
                std::basic_string<lolita::character> const & name
            )
            :
            name_(name)
            {}

            std::basic_string<lolita::character> name_;

        };

        struct GeometricalEntity
        {

            GeometricalEntity()
            :
            physical_entities2_(),
            bounding_entities_()
            {}

            void
            addPhysicalEntity(
                std::shared_ptr<PhysicalEntity> const & pe
            )
            {
                physical_entities2_.insert(pe);
                for (auto & be : bounding_entities_)
                {
                    be->addPhysicalEntity(pe);
                }
            }

            void
            addBoundingEntity(
                std::shared_ptr<GeometricalEntity> & be
            )
            {
                bounding_entities_.insert(be);
                for (auto const & pe : physical_entities2_)
                {
                    be->addPhysicalEntity(pe);
                }
            }

            std::set<std::shared_ptr<PhysicalEntity>> physical_entities2_;

            std::set<std::shared_ptr<GeometricalEntity>> bounding_entities_;

        };

        NewGmsh(
            std::basic_string_view<lolita::character> str
        )
        :
        file_(str)
        {
            setPhysicalEntities();
            setGeometricalEntities();
        }

        lolita::utility::File file_;

        void
        setPhysicalEntities()
        {
            auto const & file_lines = file_.lines_;
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
                auto physical_entity_tag = std::array<lolita::integer, 1>();
                line_stream >> dim >> physical_entity_tag[0] >> name;
                lolita::utility::removeCharacter(name, '"');
                physical_entities_[physical_entity_tag] = std::make_shared<PhysicalEntity>(PhysicalEntity(name));
                offset += 1;
            }
        }

        void
        setGeometricalEntities()
        {
            auto const & file_lines = file_.lines_;
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
                    auto geometrical_entity_tag = std::array<lolita::integer, 2>{i, lolita::integer()};
                    line_stream >> geometrical_entity_tag[1];
                    std::cout << "making entity " << geometrical_entity_tag[0] << geometrical_entity_tag[1] << std::endl;
                    geometrical_entities_[geometrical_entity_tag] = std::make_shared<GeometricalEntity>(GeometricalEntity());
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
                        auto physical_entity_tag = std::array<lolita::integer, 1>();
                        line_stream >> physical_entity_tag[0];
                        geometrical_entities_[geometrical_entity_tag]->addPhysicalEntity(physical_entities_[physical_entity_tag]);
                    }
                    if (i > 0)
                    {
                        auto num_bounding_entities = lolita::integer();
                        line_stream >> num_bounding_entities;
                        for (lolita::integer k = 0; k < num_bounding_entities; ++k)
                        {
                            auto bounding_entity_tag = std::array<lolita::integer, 2>{i - 1, lolita::integer()};
                            line_stream >> bounding_entity_tag[1];
                            bounding_entity_tag[1] = std::abs(bounding_entity_tag[1]);
                            geometrical_entities_[geometrical_entity_tag]->addBoundingEntity(geometrical_entities_[bounding_entity_tag]);
                        }
                    }
                    offset += 1;
                }
            }
        }

        std::map<std::array<lolita::integer, 1>, std::shared_ptr<PhysicalEntity>> physical_entities_;

        std::map<std::array<lolita::integer, 2>, std::shared_ptr<GeometricalEntity>> geometrical_entities_;

    };
    
    template<MeshFileFormat t_mesh_data, Domain t_domain>
    struct MeshFileDataStructure;
    
    template<MeshFileFormat t_mesh_data, Domain t_domain>
    struct MeshFileParser
    {

        using t_Implementation = typename MeshFileDataStructure<t_mesh_data, t_domain>::Implementation;

        using t_Module = typename MeshFileDataStructure<t_mesh_data, t_domain>::Module;

        MeshFileParser(
            std::basic_string_view<lolita::character> mesh_file_path
        )
        :
        file_(mesh_file_path)
        {
            auto const module = t_Module(file_);
            makeDomains(module);
            auto set_elements = [&] <lolita::integer t_i = 0u, lolita::integer t_j = 0u> (
                auto & set_elements_imp
            )
            mutable
            {
                auto const constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                makeElement<t_element>(module);
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1) {
                    set_elements_imp.template operator()<t_i, t_j + 1u>(set_elements_imp);
                }
                else if constexpr (t_i < DomainTraits<t_domain>::getNumElements() - 1) {
                    set_elements_imp.template operator()<t_i + 1u, 0u>(set_elements_imp);
                }
            };
            set_elements(set_elements);
        }

        static
        MeshData<t_domain>
        readMesh(
            std::basic_string_view<lolita::character> mesh_file_path
        )
        {
            return MeshFileParser(mesh_file_path).mesh_data_;
        }
    
        void
        makeDomains(
            t_Module const & module
        )
        {
            static_cast<t_Implementation *>(this)->makeDomains(module);
        }
        
        template<Element t_element>
        void
        makeElement(
            t_Module const & module
        )
        {
            static_cast<t_Implementation *>(this)->template makeElement<t_element>(module);
        }

        lolita::utility::File file_;

        MeshData<t_domain> mesh_data_;

    };
    
    template<MeshFileFormat t_mesh_data, Domain t_domain>
    requires(t_mesh_data.isGmsh())
    struct MeshFileDataStructure<t_mesh_data, t_domain>
    {

        lolita::integer static constexpr node_index_offset = 1;

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

        struct Module
        {

            lolita::integer static constexpr size_ = 4;

            struct PhysicalEntity
            {

                lolita::boolean
                operator==(
                    PhysicalEntity const &
                    other
                )
                const = default;

                lolita::boolean
                operator!=(
                    PhysicalEntity const &
                    other
                )
                const = default;

                lolita::integer tag_;

                lolita::integer dim_;

                std::basic_string<lolita::character> name_;

            };

            struct GeometricalEntity
            {

                lolita::boolean
                operator==(
                    GeometricalEntity const &
                    other
                )
                const = default;

                lolita::boolean
                operator!=(
                    GeometricalEntity const &
                    other
                )
                const = default;

                lolita::integer tag_;

                lolita::integer dim_;

                std::vector<lolita::integer> physical_entities_tags_;

                std::vector<lolita::integer> bounding_entities_tags_;

            };

            struct PhysicalGroup
            {

                lolita::boolean
                operator==(
                    PhysicalGroup const &
                    other
                )
                const = default;

                lolita::boolean
                operator!=(
                    PhysicalGroup const &
                    other
                )
                const = default;

                std::basic_string<lolita::character> name_;

                lolita::integer dim_;

                std::array<std::vector<lolita::integer>, size_> geometrical_entities_tags_;

            };

            explicit
            Module(
                lolita::utility::File const & file
            )
            {
                setGeometricalEntities(file);
                setPhysicalEntities(file);
                setPhysicalGroups();
            }

            std::array<std::set<lolita::integer>, size_>
            getSubGeometricalEntities(
                lolita::integer d,
                lolita::integer t,
                std::array<std::set<lolita::integer>, size_> & a
            )
            {
                a[d].insert(t);
                for (lolita::integer i = 0; i < geometrical_entities_[d][t - 1].bounding_entities_tags_.size(); ++i) {
                    lolita::integer const d2 = d - 1;
                    lolita::integer const t2 = geometrical_entities_[d][t - 1].bounding_entities_tags_[i];
                    a = getSubGeometricalEntities(d2, t2, a);
                }
                return a;
            }

            void
            setGeometricalEntities(
                lolita::utility::File const & file_
            )
            {
                auto const & file_lines = file_.lines_;
                lolita::integer line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Entities"));
                lolita::integer offset = 1;
                std::basic_stringstream<lolita::character> line_stream(file_lines[line_start + offset]);
                lolita::integer num_points;
                lolita::integer num_curves;
                lolita::integer num_surfaces;
                lolita::integer num_volumes;
                line_stream >> num_points >> num_curves >> num_surfaces >> num_volumes;
                auto num_domains = std::array<lolita::integer, size_>{num_points, num_curves, num_surfaces, num_volumes};
                offset += 1;
                for (lolita::integer i = 0; i < size_; ++i) {
                    for (lolita::integer j = 0; j < num_domains[i]; ++j) {
                        line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                        lolita::integer tag;
                        line_stream >> tag;
                        if (i == 0) {
                            for (lolita::integer k = 0; k < 3; ++k) {
                                lolita::real a;
                                line_stream >> a;
                            }
                        } else {
                            for (lolita::integer k = 0; k < 6; ++k) {
                                lolita::real a;
                                line_stream >> a;
                            }
                        }
                        lolita::integer num_physical_entities;
                        line_stream >> num_physical_entities;
                        std::vector<lolita::integer> physical_entities_tags;
                        for (lolita::integer k = 0; k < num_physical_entities; ++k) {
                            lolita::integer physical_entity_tag;
                            line_stream >> physical_entity_tag;
                            physical_entities_tags.push_back(physical_entity_tag);
                        }
                        std::vector<lolita::integer> bounding_entities_tags = {};
                        if (i > 0) {
                            lolita::integer num_bounding_entities;
                            line_stream >> num_bounding_entities;
                            for (lolita::integer k = 0; k < num_bounding_entities; ++k) {
                                lolita::integer bounding_entity_tag;
                                line_stream >> bounding_entity_tag;
                                bounding_entities_tags.push_back(std::abs(bounding_entity_tag));
                            }
                        }
                        geometrical_entities_[i].push_back(GeometricalEntity{
                                tag,
                                i,
                                physical_entities_tags,
                                bounding_entities_tags
                        });
                        offset += 1;
                    }
                }
            }

            void
            setPhysicalEntities(
                lolita::utility::File const & file_
            )
            {
                auto const & file_lines = file_.lines_;
                lolita::integer line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$PhysicalNames"));
                lolita::integer offset = 1;
                std::basic_stringstream<lolita::character> line_stream(file_lines[line_start + offset]);
                lolita::integer num_physical_names;
                line_stream >> num_physical_names;
                offset += 1;
                for (lolita::integer i = 0; i < num_physical_names; ++i) {
                    line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                    lolita::integer dim;
                    lolita::integer tag;
                    std::basic_string<lolita::character> name;
                    line_stream >> dim >> tag >> name;
                    lolita::utility::removeCharacter(name, '"');
                    physical_entities_.push_back(PhysicalEntity{tag, dim, name});
                    offset += 1;
                }
            }

            void
            setPhysicalGroups()
            {
                for (lolita::integer i = 0; i < physical_entities_.size(); ++i) {
                    std::array<std::set<lolita::integer>, size_> tags;
                    for (lolita::integer j = 0; j < size_; ++j) {
                        for (lolita::integer k = 0; k < geometrical_entities_[j].size(); ++k) {
                            for (lolita::integer l = 0; l < geometrical_entities_[j][k].physical_entities_tags_.size(); ++l) {
                                lolita::integer const & t0 = geometrical_entities_[j][k].physical_entities_tags_[l];
                                lolita::integer const & t1 = physical_entities_[i].tag_;
                                if (t0 == t1) {
                                    tags = getSubGeometricalEntities(j, k + 1, tags);
                                }
                            }
                        }
                    }
                    std::array<std::vector<lolita::integer>, size_> group_tags;
                    for (lolita::integer j = 0; j < size_; ++j) {
                        group_tags[j].assign(tags[j].begin(), tags[j].end());
                    }
                    physical_groups_.push_back(PhysicalGroup{physical_entities_[i].name_, physical_entities_[i].dim_, group_tags});
                }
            }
            
            std::array<std::vector<GeometricalEntity>, size_> geometrical_entities_;

            std::vector<PhysicalEntity> physical_entities_;

            std::vector<PhysicalGroup> physical_groups_;

        };
        
        struct Implementation : MeshFileParser<t_mesh_data, t_domain>
        {

        private:

            using t_MeshParser = MeshFileParser<t_mesh_data, t_domain>;

        public:

            void
            makeDomains(
                Module const & module
            )
            {
                auto const & physical_groups = module.physical_groups_;
                for (lolita::integer k = 0; k < physical_groups.size(); ++k) {
                    this->mesh_data_.domains_.push_back(std::make_shared<MeshDomain>(
                        MeshDomain(physical_groups[k].dim_, physical_groups[k].name_)
                    ));
                }
            }

            template<Element t_element>
            void
            makeElement(
                Module const & module
            )
            requires(t_element.isNode())
            {
                auto constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
                auto const & file_lines = this->file_.lines_;
                auto const & physical_groups = module.physical_groups_;
                auto & elements = this->mesh_data_.elements_.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>();
                auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Nodes"));
                auto offset = 1;
                auto line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                auto num_entity_block = lolita::integer();
                auto num_nodes = lolita::integer();
                auto min_node_tag = lolita::integer();
                auto max_node_tag = lolita::integer();
                line_stream >> num_entity_block >> num_nodes >> min_node_tag >> max_node_tag;
                offset += 1;
                for (lolita::integer i = 0; i < num_entity_block; ++i) {
                    line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                    auto entity_dim = lolita::integer();
                    auto entity_tag = lolita::integer();
                    auto parametric = lolita::integer();
                    auto num_nodes_in_block = lolita::integer();
                    line_stream >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
                    offset += 1;
                    for (lolita::integer j = 0; j < num_nodes_in_block; ++j) {
                        auto mesh_element = MeshElement<t_element>();
                        line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                        line_stream >> mesh_element.tag_;
                        mesh_element.tag_ --;
                        line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset + num_nodes_in_block]);
                        for (lolita::integer k = 0; k < 3; ++k) {
                            line_stream >> mesh_element.coordinates_->operator ()(k);
                        }
                        for (lolita::integer k = 0; k < physical_groups.size(); ++k) {
                            auto const & node_set_name = physical_groups[k].name_;
                            for (lolita::integer l = 0; l < physical_groups[k].geometrical_entities_tags_[entity_dim].size(); ++l) {
                                lolita::integer tag = physical_groups[k].geometrical_entities_tags_[entity_dim][l];
                                if (entity_tag == tag) {
                                    auto is_equal = [&] (std::shared_ptr<MeshDomain> const & ptr_domain) {
                                        return ptr_domain->tag_ == node_set_name;
                                    };
                                    auto const & domains = this->mesh_data_.domains_;
                                    auto domain_index = std::distance(domains.begin(), std::find_if(domains.begin(), domains.end(), is_equal));
                                    mesh_element.domains_.push_back(domains[domain_index]);
                                }
                            }
                        }
                        elements.push_back(mesh_element);
                        offset += 1;
                    }
                    offset += num_nodes_in_block;
                }
            }

            template<Element t_element>
            void
            makeElement(
                Module const & module
            )
            requires(!t_element.isNode())
            {
                auto const constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
                auto const & file_lines = this->file_.lines_;
                auto & elements = this->mesh_data_.elements_.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>();
                auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Elements"));
                auto offset = 1;
                auto line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                auto num_entity_blocks = 0;
                auto num_elements = 0;
                auto min_element_tag = 0;
                auto max_element_tag = 0;
                line_stream >> num_entity_blocks >> num_elements >> min_element_tag >> max_element_tag;
                offset += 1;
                for (lolita::integer i = 0; i < num_entity_blocks; ++i) {
                    line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                    auto entity_dim = 0;
                    auto entity_tag = 0;
                    auto element_type_tag = 0;
                    auto num_elements_in_block = 0;
                    line_stream >> entity_dim >> entity_tag >> element_type_tag >> num_elements_in_block;
                    offset += 1;
                    auto const element = getElemType(element_type_tag);
                    if (entity_dim == t_domain.dim_ && element == t_element) {
                        for (lolita::integer j = 0; j < num_elements_in_block; ++j) {
                            auto mesh_element = MeshElement<t_element>();
                            line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                            auto tag = lolita::integer();
                            line_stream >> tag;
                            for (lolita::integer k = 0; k < element.num_nodes_; ++k) {
                                auto node_tag = lolita::natural();
                                line_stream >> node_tag;
                                mesh_element.node_tags_[k] = node_tag - 1;
                            }
                            elements.push_back(mesh_element);
                            offset += 1;
                        }
                    }
                    else {
                        offset += num_elements_in_block;
                    }
                }
            }

        };

    };



    template<template<Element, Domain, auto...> typename t_T, Element t_element, Domain t_domain, auto... t_args>
    struct FiniteElementConnectivity
    {

    private:
    
        using t_ElementTraits = ElementTraits<t_element, t_domain>;
        
        template<Element t__element, Domain t__domain>
        using t_ElementPointer = std::shared_ptr<t_T<t__element, t__domain, t_args...>>;

    public:
    
        using Components = typename t_ElementTraits::template InnerConnectivity<t_ElementPointer>;
        
        using Neighbours = typename t_ElementTraits::template OuterConnectivity<t_ElementPointer>;
        
        lolita::natural tag_;

        // std::basic_string_view<lolita::character> hash_;
        
        Neighbours neighbours_;
        
        Components components_;
        
        // std::vector<std::shared_ptr<std::basic_string<lolita::character>>> domains_;

        std::vector<std::shared_ptr<MeshDomain>> domains_;
        
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

    template<Element t_element, Domain t_domain, auto t_arg>
    struct Fefinal : FiniteElementConnectivity<Fefinal, t_element, t_domain, t_arg>
    {

    };

    template<Element t_element, Domain t_domain, ElementGroupConcept auto... t_element_group>
    struct FiniteElementHolder : FiniteElementConnectivity<FiniteElementHolder, t_element, t_domain, t_element_group...>
    {

        using Elements = std::tuple<std::shared_ptr<Fefinal<t_element, t_domain, t_element_group>>...>;

        Elements elements_;
    
        template<lolita::integer t_i>
        std::tuple_element_t<t_i, Elements> const &
        getElement()
        const
        {
            return std::get<t_i>(elements_);
        }
        
        template<lolita::integer t_i>
        std::tuple_element_t<t_i, Elements> &
        getElement()
        {
            return std::get<t_i>(elements_);
        }

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
            std::shared_ptr<FiniteElementHolder> & ptr_element,
            std::array<lolita::integer, t_element.num_nodes_> const & node_tags
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
                        auto const k = FiniteElementHolder::template getComponentNodeConnection<t_i, t_j>(i, j);
                        component_node_tags[j] = node_tags[k];
                    }
                    component_hash = lolita2::geometry::FiniteElementHolder<_component, t_domain, t_element_group...>::getElementHash(component_node_tags);
                    if (!components.contains(component_hash))
                    {
                        using t_Component = lolita2::geometry::FiniteElementHolder<_component, t_domain, t_element_group...>;
                        auto ptr_component = std::make_shared<t_Component>(t_Component());
                        lolita2::geometry::FiniteElementHolder<_component, t_domain, t_element_group...>::template makeElement<0, 0>(ptr_component, component_node_tags);
                    }
                }
                else
                {
                    component_hash = std::to_string(node_tags[FiniteElementHolder::template getComponentNodeConnection<t_i, t_j>(i, 0)]);
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

        // static
        // void
        // makeElement(
        //     auto & list,
        //     std::shared_ptr<FiniteElementHolder> & ptr_element,
        //     lolita::natural const & tag,
        //     Point const & coordinates,
        //     std::vector<std::shared_ptr<std::basic_string<lolita::character>>> const & domains
        // )
        // requires(t_element.isNode())
        // {
        //     auto const constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
        //     auto hklm = std::forward<std::shared_ptr<FiniteElementHolder>>(ptr_element);
        //     hklm->tag_ = std::forward<lolita::natural>(tag);
        //     hklm->domains_ = std::forward<std::vector<std::shared_ptr<std::basic_string<lolita::character>>>>(domains);
        //     hklm->coordinates_ = std::make_shared<Point>(std::forward<Point>(coordinates));
        //     list.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>()[std::to_string(tag)] = hklm;
        // }

        template<lolita::integer t_i = 1>
        static
        void
        hello(
            auto & list,
            std::shared_ptr<FiniteElementHolder> ptr_element,
            lolita::natural const & tag,
            Point const & coordinates,
            std::vector<std::shared_ptr<std::basic_string<lolita::character>>> const & domains
        )
        {
            auto const constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            list.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>();
            std::cout << "yo !" << std::endl;
            // auto hklm = std::forward<std::shared_ptr<FiniteElementHolder>>(ptr_element);
            if constexpr (t_i > 0)
            {
                hello<t_i - 1>(list, ptr_element, tag, coordinates, domains);
            }
            // auto tagg = std::forward<lolita::natural>(tag);
            // auto coordinatess = std::forward<Point>(coordinates);
            // auto domainss = std::forward<std::vector<std::shared_ptr<std::basic_string<lolita::character>>>>(domains);
        }

    };

    template<template<Element, Domain, auto...> typename t_FiniteElement, Domain t_domain, auto... t_args>
    struct MeshSetBase
    {

    private:

        template<template<Element, Domain, auto...> typename t__FiniteElement, auto... t__args>
        struct FiniteElementCollectionTraits
        {

            template<Element t_element, Domain t__domain>
            using t_FiniteElements = std::map<std::basic_string<lolita::character>, std::shared_ptr<t__FiniteElement<t_element, t__domain, t__args...>>>;

            using t_FiniteElementCollection = ElementCollection<t_FiniteElements, t_domain>;

        };

        template<template<Element, Domain, auto...> typename t__FiniteElement, auto... t__args>
        using t_FiniteElements = typename FiniteElementCollectionTraits<t__FiniteElement, t__args...>::t_FiniteElementCollection;

    public:

        using FiniteElements = typename FiniteElementCollectionTraits<t_FiniteElement, t_args...>::t_FiniteElementCollection;

        using Domains = std::map<std::basic_string<lolita::character>, std::shared_ptr<std::basic_string<lolita::character>>>;
        
        friend
        std::ostream &
        operator<<(
            std::ostream & os,
            MeshSetBase const & mesh
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

        FiniteElements elements_;

        Domains domains_;

    };

    // template<Domain t_domain, auto... t_args>
    // using Mesh = MeshSetBase<FiniteElementHolder, t_domain, t_args...>;
    
    template<template<Element, Domain, auto...> typename t_FiniteElement, MeshFileFormat t_mesh, Domain t_domain, auto... t_args>
    struct MeshParserModule;
    
    template<template<Element, Domain, auto...> typename t_FiniteElement, MeshFileFormat t_mesh, Domain t_domain, auto... t_args>
    struct MeshParser
    {

    private:

        using t_Module = typename MeshParserModule<t_FiniteElement, t_mesh, t_domain, t_args...>::Module;

        using t_Implementation = typename MeshParserModule<t_FiniteElement, t_mesh, t_domain, t_args...>::Implementation;

    public:
    
        MeshParser(
            std::basic_string_view<lolita::character> mesh_file_path
        )
        {
            
            auto const module = t_Module(mesh_file_path);
            makeDomains(module);
            auto set_elements = [&] <lolita::integer t_i = 0u, lolita::integer t_j = 0u> (
                auto & set_elements_imp
            )
            mutable
            {
                auto const constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                makeElement<t_element>(module);
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1) {
                    set_elements_imp.template operator()<t_i, t_j + 1u>(set_elements_imp);
                }
                else if constexpr (t_i < DomainTraits<t_domain>::getNumElements() - 1) {
                    set_elements_imp.template operator()<t_i + 1u, 0u>(set_elements_imp);
                }
            };
            set_elements(set_elements);
            auto make_elements = [&] <lolita::integer t_i = 0u, lolita::integer t_j = 0u> (
                auto & self
            )
            mutable
            {
                auto const constexpr t_element = DomainTraits<t_domain>::template getElement<t_i, t_j>();
                buildElementN<t_element>();
                if constexpr (t_j < DomainTraits<t_domain>::template getNumElements<t_i>() - 1) {
                    self.template operator()<t_i, t_j + 1u>(self);
                }
                else if constexpr (t_i < DomainTraits<t_domain>::getNumElements() - 1) {
                    self.template operator()<t_i + 1u, 0u>(self);
                }
            };
            make_elements(make_elements);
        }

        static
        MeshSetBase<t_FiniteElement, t_domain, t_args...>
        makeMesh(
            std::basic_string_view<lolita::character> mesh_file_path
        )
        {
            return MeshParser(mesh_file_path).mesh_data_;
        }

    protected:

        template<Element t_element>
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

        template<Element t_element>
        void
        makeElement(
            std::shared_ptr<t_FiniteElement<t_element, t_domain, t_args...>> & ptr_element,
            std::array<lolita::integer, t_element.num_nodes_> node_tags
        )
        requires(!t_element.isNode())
        {
            using t_ElementDescription = ElementTraits<t_element, t_domain>;
            using t_MeshDescription = DomainTraits<t_domain>;
            auto const constexpr t_element_coordinates = t_MeshDescription::template getElementCoordinates<t_element>();
            auto const constexpr t_node_coordinates = t_ElementDescription::template getComponentCoordinates<Element::node()>();
            auto & elements = mesh_data_.elements_.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>();
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
            ptr_element->tag_ = elements.size();
            ptr_element->domains_.assign(domains.begin(), domains.end());
            elements.insert({getElementHash<t_element>(node_tags), ptr_element});
            // for (auto const & domain : domains)
            // {
            //     auto & set = mesh_data_.sets_[* domain].template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>();
            //     set.insert({getElementHash<t_element>(node_tags), ptr_element});
            // }
            
        }

        template<Element t_element, lolita::integer t_i = 0, lolita::integer t_j = 0>
        void
        makeElement2(
            std::shared_ptr<t_FiniteElement<t_element, t_domain, t_args...>> & ptr_element,
            std::array<lolita::integer, t_element.num_nodes_> node_tags
        )
        requires(!t_element.isNode())
        {
            using t_Element = t_FiniteElement<t_element, t_domain, t_args...>;
            using t_ElementDescription = ElementTraits<t_element, t_domain>;
            using t_ComponentDescription = ElementTraits<t_ElementDescription::template getComponent<t_i, t_j>(), t_domain>;
            using t_MeshDescription = DomainTraits<t_domain>;
            auto const constexpr _is_initialized = t_i == 0 && t_j == 0;
            auto const constexpr _component = t_ElementDescription::template getComponent<t_i, t_j>();
            auto const constexpr _component_coordinates = t_MeshDescription::template getElementCoordinates<_component>();
            auto const constexpr _neighbour_coordinates = t_ComponentDescription::template getNeighbourCoordinates<t_element>();
            using t_Component = t_FiniteElement<_component, t_domain, t_args...>;
            auto & components = mesh_data_.elements_.template getElements<_component_coordinates.dim_, _component_coordinates.tag_>();
            auto & element_component_array = ptr_element->template getComponents<t_i, t_j>();
            for (auto i = 0; i < element_component_array.size(); ++i)
            {
                auto component_hash = std::basic_string<lolita::character>();
                if constexpr(!_component.isNode())
                {
                    auto component_node_tags = std::array<lolita::integer, _component.num_nodes_>();
                    for (lolita::integer j = 0; j < _component.num_nodes_; ++j)
                    {
                        auto const k = t_Element::template getComponentNodeConnection<t_i, t_j>(i, j);
                        component_node_tags[j] = node_tags[k];
                    }
                    component_hash = getElementHash<_component>(component_node_tags);
                    if (!components.contains(component_hash))
                    {
                        auto ptr_component = std::make_shared<t_Component>(t_Component());
                        makeElement2<_component, 0, 0>(ptr_component, component_node_tags);
                    }
                }
                else
                {
                    component_hash = std::to_string(node_tags[t_Element::template getComponentNodeConnection<t_i, t_j>(i, 0)]);
                }
                element_component_array[i] = components[component_hash];
                components[component_hash]->template getNeighbours<_neighbour_coordinates.dim_, _neighbour_coordinates.tag_>().push_back(ptr_element);
            }
            if constexpr (t_j < t_ElementDescription::template getNumComponents<t_i>() - 1)
            {
                makeElement2<t_element, t_i, t_j + 1u>(ptr_element, node_tags);
            }
            else if constexpr (t_i < t_ElementDescription::getNumComponents() - 1)
            {
                makeElement2<t_element, t_i + 1u, 0u>(ptr_element, node_tags);
            }
            if constexpr (_is_initialized)
            {
                makeElement<t_element>(ptr_element, node_tags);
            }
        }

        template<Element t_element>
        void
        makeElement2(
            std::array<lolita::integer, t_element.num_nodes_> node_tags
        )
        requires(!t_element.isNode())
        {
            using t_Element = t_FiniteElement<t_element, t_domain, t_args...>;
            auto ptr_element = std::make_shared<t_Element>(t_Element());
            makeElement2<t_element, 0, 0>(ptr_element, node_tags);
        }

        template<Element t_element>
        void
        makeElement2(
                lolita::natural const & tag,
                Point const & coordinates,
                std::vector<std::shared_ptr<std::basic_string<lolita::character>>> const & domains
        )
        requires(t_element.isNode())
        {
            using t_Element = t_FiniteElement<t_element, t_domain, t_args...>;
            auto ptr_element = std::make_shared<t_Element>(t_Element());
            auto const constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            ptr_element->tag_ = tag;
            ptr_element->domains_ = domains;
            ptr_element->coordinates_ = std::make_shared<Point>(coordinates);
            // auto elem_hash = std::to_string(tag);
            mesh_data_.elements_.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>().insert({std::to_string(tag), ptr_element});
            // for (auto const & domain : domains)
            // {
            //     auto & set = mesh_data_.element_sets_[* domain].template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>();
            //     set.insert({std::to_string(tag), ptr_element});
            // }
        }
    
        void
        makeDomains(
            t_Module const & module
        )
        {
            static_cast<t_Implementation *>(this)->makeDomains(module);
        }
        
        template<Element t_element>
        void
        makeElement(
            t_Module const & module
        )
        {
            static_cast<t_Implementation *>(this)->template makeElement<t_element>(module);
        }

        template<Element t_element, lolita::integer _k = 0>
        void
        buildElementN()
        {
            auto const constexpr t_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
            auto & element_map = mesh_data_.elements_.template getElements<t_coordinates.dim_, t_coordinates.tag_>();
            for (auto & element_map_item : element_map)
            {
                auto & element = element_map_item.second;
                if constexpr (!t_element.isNode())
                {
                    auto const & element_nds = element->template getComponents<t_coordinates.dim_ - 1, 0>();
                    for (lolita::integer i = 0; i < element_nds.size(); ++i)
                    {
                        auto const & nde = element_nds[i];
                        auto const & ngs = nde->template getNeighbours<t_coordinates.dim_ - 1, _k>();
                        for (auto const & neighbour : ngs)
                        {
                            if (((neighbour->tag_ != element->tag_) && _k == t_coordinates.tag_) || (_k != t_coordinates.tag_))
                            {
                                auto & element_ngs = element->template getNeighbours<0, _k>();
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
                    if constexpr (_k < DomainTraits<t_domain>::template getNumElements<t_coordinates.dim_>() - 1)
                    {
                        buildElementN<t_element, _k + 1>();
                    }
                }
            }
        }

    public:
    
        MeshSetBase<t_FiniteElement, t_domain, t_args...> mesh_data_;

    };

    template<template<Element, Domain, auto...> typename t_FiniteElement, MeshFileFormat t_mesh, Domain t_domain, auto... t_args>
    requires(t_mesh.isGmsh())
    struct MeshParserModule<t_FiniteElement, t_mesh, t_domain, t_args...>
    {

        lolita::integer static constexpr node_index_offset = 1;

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

        struct Module
        {

            struct PhysicalEntity
            {

                lolita::boolean
                operator==(
                    PhysicalEntity const &
                    other
                )
                const = default;

                lolita::boolean
                operator!=(
                    PhysicalEntity const &
                    other
                )
                const = default;

                lolita::integer tag_;

                lolita::integer dim_;

                std::basic_string<lolita::character> name_;

            };

            struct GeometricalEntity
            {

                lolita::boolean
                operator==(
                    GeometricalEntity const &
                    other
                )
                const = default;

                lolita::boolean
                operator!=(
                    GeometricalEntity const &
                    other
                )
                const = default;

                lolita::integer tag_;

                lolita::integer dim_;

                std::vector<lolita::integer> physical_entities_tags_;

                std::vector<lolita::integer> bounding_entities_tags_;

            };

            struct PhysicalGroup
            {

                lolita::boolean
                operator==(
                    PhysicalGroup const &
                    other
                )
                const = default;

                lolita::boolean
                operator!=(
                    PhysicalGroup const &
                    other
                )
                const = default;

                std::basic_string<lolita::character> name_;

                lolita::integer dim_;

                std::array<std::vector<lolita::integer>, 4> geometrical_entities_tags_;

            };

            explicit
            Module(
                std::basic_string_view<lolita::character> mesh_file_path
            )
            :
            file_(mesh_file_path)
            {
                setGeometricalEntities();
                setPhysicalEntities();
                setPhysicalGroups();
            }

            std::array<std::vector<GeometricalEntity>, 4>
            makeGeometricalEntities()
            {
                auto geometrical_entities = std::array<std::vector<GeometricalEntity>, 4>();
                auto const & file_lines = file_.lines_;
                lolita::integer line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Entities"));
                lolita::integer offset = 1;
                std::basic_stringstream<lolita::character> line_stream(file_lines[line_start + offset]);
                lolita::integer num_points;
                lolita::integer num_curves;
                lolita::integer num_surfaces;
                lolita::integer num_volumes;
                line_stream >> num_points >> num_curves >> num_surfaces >> num_volumes;
                auto num_domains = std::array<lolita::integer, 4>{num_points, num_curves, num_surfaces, num_volumes};
                offset += 1;
                for (lolita::integer i = 0; i < 4; ++i) {
                    for (lolita::integer j = 0; j < num_domains[i]; ++j) {
                        line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                        lolita::integer tag;
                        line_stream >> tag;
                        if (i == 0) {
                            for (lolita::integer k = 0; k < 3; ++k) {
                                lolita::real a;
                                line_stream >> a;
                            }
                        } else {
                            for (lolita::integer k = 0; k < 6; ++k) {
                                lolita::real a;
                                line_stream >> a;
                            }
                        }
                        lolita::integer num_physical_entities;
                        line_stream >> num_physical_entities;
                        std::vector<lolita::integer> physical_entities_tags;
                        for (lolita::integer k = 0; k < num_physical_entities; ++k) {
                            lolita::integer physical_entity_tag;
                            line_stream >> physical_entity_tag;
                            physical_entities_tags.push_back(physical_entity_tag);
                        }
                        std::vector<lolita::integer> bounding_entities_tags = {};
                        if (i > 0) {
                            lolita::integer num_bounding_entities;
                            line_stream >> num_bounding_entities;
                            for (lolita::integer k = 0; k < num_bounding_entities; ++k) {
                                lolita::integer bounding_entity_tag;
                                line_stream >> bounding_entity_tag;
                                bounding_entities_tags.push_back(std::abs(bounding_entity_tag));
                            }
                        }
                        geometrical_entities[i].push_back(GeometricalEntity{
                                tag,
                                i,
                                physical_entities_tags,
                                bounding_entities_tags
                        });
                        offset += 1;
                    }
                }
                return geometrical_entities;
            }

            std::vector<PhysicalEntity>
            makePhysicalEntities()
            {
                auto physical_entities = std::vector<PhysicalEntity>();
                auto const & file_lines = this->file_.lines_;
                lolita::integer line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$PhysicalNames"));
                lolita::integer offset = 1;
                std::basic_stringstream<lolita::character> line_stream(file_lines[line_start + offset]);
                lolita::integer num_physical_names;
                line_stream >> num_physical_names;
                offset += 1;
                for (lolita::integer i = 0; i < num_physical_names; ++i) {
                    line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                    lolita::integer dim;
                    lolita::integer tag;
                    std::basic_string<lolita::character> name;
                    line_stream >> dim >> tag >> name;
                    lolita::utility::removeCharacter(name, '"');
                    physical_entities.push_back(PhysicalEntity{tag, dim, name});
                    offset += 1;
                }
                return physical_entities;
            }

            std::array<std::set<lolita::integer>, 4>
            getSubGeometricalEntities(
                lolita::integer d,
                lolita::integer t,
                std::array<std::set<lolita::integer>, 4> & a
            )
            {
                a[d].insert(t);
                for (lolita::integer i = 0; i < geometrical_entities_[d][t - 1].bounding_entities_tags_.size(); ++i) {
                    lolita::integer const d2 = d - 1;
                    lolita::integer const t2 = geometrical_entities_[d][t - 1].bounding_entities_tags_[i];
                    a = getSubGeometricalEntities(d2, t2, a);
                }
                return a;
            }

            std::vector<PhysicalGroup>
            makePhysicalGroups()
            {
                auto get_subs = [&] (
                        auto d,
                        auto t,
                        std::array<std::set<lolita::integer>, 4> & a,
                        auto & self
                )
                        mutable
                {
                    a[d].insert(t);
                    for (lolita::integer i = 0; i < geometrical_entities_[d][t - 1].bounding_entities_tags.size(); ++i) {
                        a = self(d - 1, geometrical_entities_[d][t - 1].bounding_entities_tags[i], a, self);
                    }
                    return a;
                };
                /*
                 *
                 */
                auto physical_groups = std::vector<PhysicalGroup>();
                for (lolita::integer i = 0; i < physical_entities_.size(); ++i) {
                    std::array<std::set<lolita::integer>, 4> tags;
                    for (lolita::integer j = 0; j < 4; ++j) {
                        for (lolita::integer k = 0; k < geometrical_entities_[j].size(); ++k) {
                            for (lolita::integer l = 0; l < geometrical_entities_[j][k].physical_entities_tags_.size(); ++l) {
                                lolita::integer const & t0 = geometrical_entities_[j][k].physical_entities_tags_[l];
                                lolita::integer const & t1 = physical_entities_[i].tag_;
                                if (t0 == t1) {
//                                    tags = get_subs(j, k + 1, tags, get_subs);
                                    tags = getSubGeometricalEntities(j, k + 1, tags);
                                }
                            }
                        }
                    }
                    std::array<std::vector<lolita::integer>, 4> group_tags;
                    for (lolita::integer j = 0; j < 4; ++j) {
                        group_tags[j].assign(tags[j].begin(), tags[j].end());
                    }
                    physical_groups.push_back(PhysicalGroup{physical_entities_[i].name_, physical_entities_[i].dim_, group_tags});
                }
                return physical_groups;
            }

            void
            setGeometricalEntities()
            {
                auto const & file_lines = file_.lines_;
                lolita::integer line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Entities"));
                lolita::integer offset = 1;
                std::basic_stringstream<lolita::character> line_stream(file_lines[line_start + offset]);
                lolita::integer num_points;
                lolita::integer num_curves;
                lolita::integer num_surfaces;
                lolita::integer num_volumes;
                line_stream >> num_points >> num_curves >> num_surfaces >> num_volumes;
                auto num_domains = std::array<lolita::integer, 4>{num_points, num_curves, num_surfaces, num_volumes};
                offset += 1;
                for (lolita::integer i = 0; i < 4; ++i) {
                    for (lolita::integer j = 0; j < num_domains[i]; ++j) {
                        line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                        lolita::integer tag;
                        line_stream >> tag;
                        if (i == 0) {
                            for (lolita::integer k = 0; k < 3; ++k) {
                                lolita::real a;
                                line_stream >> a;
                            }
                        } else {
                            for (lolita::integer k = 0; k < 6; ++k) {
                                lolita::real a;
                                line_stream >> a;
                            }
                        }
                        lolita::integer num_physical_entities;
                        line_stream >> num_physical_entities;
                        std::vector<lolita::integer> physical_entities_tags;
                        for (lolita::integer k = 0; k < num_physical_entities; ++k) {
                            lolita::integer physical_entity_tag;
                            line_stream >> physical_entity_tag;
                            physical_entities_tags.push_back(physical_entity_tag);
                        }
                        std::vector<lolita::integer> bounding_entities_tags = {};
                        if (i > 0) {
                            lolita::integer num_bounding_entities;
                            line_stream >> num_bounding_entities;
                            for (lolita::integer k = 0; k < num_bounding_entities; ++k) {
                                lolita::integer bounding_entity_tag;
                                line_stream >> bounding_entity_tag;
                                bounding_entities_tags.push_back(std::abs(bounding_entity_tag));
                            }
                        }
                        geometrical_entities_[i].push_back(GeometricalEntity{
                                tag,
                                i,
                                physical_entities_tags,
                                bounding_entities_tags
                        });
                        offset += 1;
                    }
                }
            }

            void
            setPhysicalEntities()
            {
                auto const & file_lines = file_.lines_;
                lolita::integer line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$PhysicalNames"));
                lolita::integer offset = 1;
                std::basic_stringstream<lolita::character> line_stream(file_lines[line_start + offset]);
                lolita::integer num_physical_names;
                line_stream >> num_physical_names;
                offset += 1;
                for (lolita::integer i = 0; i < num_physical_names; ++i) {
                    line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                    lolita::integer dim;
                    lolita::integer tag;
                    std::basic_string<lolita::character> name;
                    line_stream >> dim >> tag >> name;
                    lolita::utility::removeCharacter(name, '"');
                    physical_entities_.push_back(PhysicalEntity{tag, dim, name});
                    offset += 1;
                }
            }

            void
            setPhysicalGroups()
            {
                for (lolita::integer i = 0; i < physical_entities_.size(); ++i) {
                    std::array<std::set<lolita::integer>, 4> tags;
                    for (lolita::integer j = 0; j < 4; ++j) {
                        for (lolita::integer k = 0; k < geometrical_entities_[j].size(); ++k) {
                            for (lolita::integer l = 0; l < geometrical_entities_[j][k].physical_entities_tags_.size(); ++l) {
                                lolita::integer const & t0 = geometrical_entities_[j][k].physical_entities_tags_[l];
                                lolita::integer const & t1 = physical_entities_[i].tag_;
                                if (t0 == t1) {
                                    tags = getSubGeometricalEntities(j, k + 1, tags);
                                }
                            }
                        }
                    }
                    std::array<std::vector<lolita::integer>, 4> group_tags;
                    for (lolita::integer j = 0; j < 4; ++j) {
                        group_tags[j].assign(tags[j].begin(), tags[j].end());
                    }
                    physical_groups_.push_back(PhysicalGroup{physical_entities_[i].name_, physical_entities_[i].dim_, group_tags});
                }
            }

            lolita::utility::File file_;

            std::array<std::vector<GeometricalEntity>, 4> geometrical_entities_;

            std::vector<PhysicalEntity> physical_entities_;

            std::vector<PhysicalGroup> physical_groups_;

        };

        struct Implementation : public MeshParser<t_FiniteElement, t_mesh, t_domain, t_args...>
        {

        private:

            using t_MeshParser = MeshParser<t_FiniteElement, t_mesh, t_domain, t_args...>;

        public:

            void
            makeDomains(
                    Module const & module
            )
            {
                auto const & physical_groups = module.physical_groups_;
                for (lolita::integer k = 0; k < physical_groups.size(); ++k) {
                    // this->mesh_data_.domains_.push_back(std::make_shared<std::basic_string<lolita::character>>(physical_groups[k].name_));
                    // this->mesh_data_.domains_.[physical_groups[k].name_] = std::make_shared<std::basic_string<lolita::character>>(physical_groups[k].name_);
                    // this->mesh_data.element_sets_[physical_groups[k].name_] = typename Mesh<t_domain>::FiniteElements();
                    // this->mesh_data.element_sets_.insert({physical_groups[k].name_, typename Mesh<t_domain>::FiniteElements()});
                    this->mesh_data_.domains_.insert({physical_groups[k].name_, std::make_shared<std::basic_string<lolita::character>>(physical_groups[k].name_)});
                }
            }

            template<Element t_element>
            void
            makeElement(
                    Module const & module
            )
            requires(t_element.isNode())
            {
                auto const & file_lines = module.file_.lines_;
                auto const & physical_groups = module.physical_groups_;
                auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Nodes"));
                auto offset = 1;
                auto line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                auto num_entity_block = lolita::integer();
                auto num_nodes = lolita::integer();
                auto min_node_tag = lolita::integer();
                auto max_node_tag = lolita::integer();
                line_stream >> num_entity_block >> num_nodes >> min_node_tag >> max_node_tag;
                offset += 1;
                for (lolita::integer i = 0; i < num_entity_block; ++i) {
                    line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                    auto entity_dim = lolita::integer();
                    auto entity_tag = lolita::integer();
                    auto parametric = lolita::integer();
                    auto num_nodes_in_block = lolita::integer();
                    line_stream >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
                    offset += 1;
                    for (lolita::integer j = 0; j < num_nodes_in_block; ++j) {
                        auto tag_ = lolita::natural();
                        auto coordinates_ = Point();
                        auto domains_ = std::vector<std::shared_ptr<std::basic_string<lolita::character>>>();
                        line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                        line_stream >> tag_;
                        tag_ --;
                        line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset + num_nodes_in_block]);
                        coordinates_.setZero();
                        for (lolita::integer k = 0; k < 3; ++k) {
                            line_stream >> coordinates_(k);
                        }
                        for (lolita::integer k = 0; k < physical_groups.size(); ++k) {
                            auto const & node_set_name = physical_groups[k].name_;
                            for (lolita::integer l = 0; l < physical_groups[k].geometrical_entities_tags_[entity_dim].size(); ++l) {
                                lolita::integer tag = physical_groups[k].geometrical_entities_tags_[entity_dim][l];
                                if (entity_tag == tag) {
                                    // auto is_equal = [&] (std::shared_ptr<std::basic_string<lolita::character>> const & ptr_domain) {
                                    //     return * ptr_domain == node_set_name;
                                    // };
                                    // auto const & domains = this->mesh_data_.domains_;
                                    // auto domain_index = std::distance(domains.begin(), std::find_if(domains.begin(), domains.end(), is_equal));
                                    // domains_.push_back(domains[domain_index]);
                                    domains_.push_back(this->mesh_data_.domains_[node_set_name]);
                                }
                            }
                        }
                        t_MeshParser::template makeElement2<t_element>(tag_, coordinates_, domains_);
                        // t_FiniteElement<t_element, t_domain, t_args...>::makeElement(
                        //     this->mesh_data_.elements_,
                        //     std::make_shared<t_FiniteElement<t_element, t_domain, t_args...>>(t_FiniteElement<t_element, t_domain, t_args...>()),
                        //     tag_, coordinates_, domains_
                        // );
                        t_FiniteElement<t_element, t_domain, t_args...>::hello(
                            this->mesh_data_.elements_,
                            std::make_shared<t_FiniteElement<t_element, t_domain, t_args...>>(t_FiniteElement<t_element, t_domain, t_args...>()),
                            tag_, coordinates_, domains_
                        );
                        offset += 1;
                    }
                    offset += num_nodes_in_block;
                }
            }

            template<Element t_element>
            void
            makeElement(
                    Module const &
                    module
            )
            requires(!t_element.isNode())
            {
                auto const constexpr t_element_coordinates = DomainTraits<t_domain>::template getElementCoordinates<t_element>();
                auto const & file_lines = module.file_.lines_;
                auto & elements = this->mesh_data_.elements_.template getElements<t_element_coordinates.dim_, t_element_coordinates.tag_>();
                auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Elements"));
                auto offset = 1;
                auto line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                auto num_entity_blocks = 0;
                auto num_elements = 0;
                auto min_element_tag = 0;
                auto max_element_tag = 0;
                line_stream >> num_entity_blocks >> num_elements >> min_element_tag >> max_element_tag;
                offset += 1;
                for (lolita::integer i = 0; i < num_entity_blocks; ++i) {
                    line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                    auto entity_dim = 0;
                    auto entity_tag = 0;
                    auto element_type_tag = 0;
                    auto num_elements_in_block = 0;
                    line_stream >> entity_dim >> entity_tag >> element_type_tag >> num_elements_in_block;
                    offset += 1;
                    auto const element = getElemType(element_type_tag);
                    if (entity_dim == t_domain.dim_ && element == t_element) {
                        for (lolita::integer j = 0; j < num_elements_in_block; ++j) {
                            auto node_tags = std::array<lolita::integer, t_element.num_nodes_>();
                            line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                            auto tag = lolita::integer();
                            line_stream >> tag;
                            for (lolita::integer k = 0; k < element.num_nodes_; ++k) {
                                auto node_tag = lolita::integer();
                                line_stream >> node_tag;
                                node_tags[k] = node_tag - 1;
                            }
                            t_MeshParser::template makeElement2<t_element>(node_tags);
                            offset += 1;
                        }
                    }
                    else {
                        offset += num_elements_in_block;
                    }
                }
            }

        };

    };

    // template<template<Element, Domain, auto...> typename t_FiniteElement, MeshFileFormat t_mesh, Domain t_domain, auto... t_args>
    // struct MeshParser

    template<MeshFileFormat t_mesh, Domain t_domain, auto... t_args>
    using MyMeshParser = MeshParser<FiniteElementHolder, t_mesh, t_domain, t_args...>;

}

#endif /* CC79BDC5_49DB_4A81_8A93_18ABD6551AF1 */
