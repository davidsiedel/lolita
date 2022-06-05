//
// Created by dsiedel on 04/06/22.
//

#ifndef LOLITA_LOLITA_CORE_4_HXX
#define LOLITA_LOLITA_CORE_4_HXX

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"
#include "lolita/lolita_core_1.hxx"
#include "lolita/lolita_core_2.hxx"
#include "lolita/lolita_core_3.hxx"

namespace lolita::core2::finite_element
{

    template<lolita::core2::geometry::Element _element, lolita::domain::Domain _domain, auto _finite_element>
    struct FEObject;

}

namespace lolita::core2::mesh
{

    /**
     * @brief
     * @tparam _T
     * @tparam _domain
     * @tparam _arg
     */
    template<template<lolita::core2::geometry::Element, lolita::domain::Domain, auto> typename _T, lolita::domain::Domain _domain, auto _args>
    struct ElementCollection
    {

    private:

        using _Elements = lolita::core2::geometry::Elements<_T, _domain, _args>;

    public:

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, _Elements>> const &
        getElements()
        const
        {
            return std::get<_j>(std::get<_i>(elements_));
        }

        /**
         * @brief
         * @tparam _i
         * @tparam _j
         * @return
         */
        template<lolita::index _i, lolita::index _j>
        std::tuple_element_t<_j, std::tuple_element_t<_i, _Elements>> &
        getElements()
        {
            return std::get<_j>(std::get<_i>(elements_));
        }

        /**
         * @brief The list of all available elements in a given domain
         */
        _Elements elements_;

    };

    /**
     * @brief
     */
    struct FiniteElementLinearSystem
    {

        /**
         * @brief
         */
        lolita::natural num_unknowns_ = 0;

        /**
         * @brief
         */
        lolita::natural num_bindings_ = 0;

        /**
         * @brief
         */
        Eigen::SparseMatrix<lolita::real> lhs_;

        /**
         * @brief
         */
        lolita::matrix::Vector<lolita::real> rhs_;

    };

    /**
     * @brief
     * @tparam _domain
     * @tparam _finite_element
     */
    template<lolita::domain::Domain _domain, auto _finite_element>
    struct Mesh
    {

    private:

        /**
         * @brief
         */
        template<lolita::core2::geometry::Element _element, lolita::domain::Domain __domain, auto __finite_element>
        using _ElementPointerMap = std::unordered_map<
                std::basic_string<lolita::character>,
                std::shared_ptr<lolita::core2::finite_element::FEObject<_element, __domain, __finite_element>>
        >;

    public:

        /**
         * @brief
         * @param os
         * @param mesh
         * @return
         */
        friend
        std::ostream &
        operator<<(
                std::ostream & os,
                Mesh const & mesh
        )
        {
            /*
             * print all components of a given element
             */
            auto print_element_components = [&] <lolita::core2::geometry::Element _element, lolita::integer _i = 0, lolita::integer _j = 0> (
                    auto const & elem_arg,
                    auto & self
            )
                    mutable
            {
                if constexpr (!_element.isPoint()) {
                    auto const constexpr _component = lolita::core2::geometry::ElementGeometryTraits<_element, _domain>::template getComponent<_i, _j>();
                    for (auto const & c_ : elem_arg->template getComponents<_i, _j>()) {
                        os << "layer : " << _i << " type : " << _j << " <-- " << _component << " " << c_->hash() << std::endl;
                    }
                    if constexpr (_j < lolita::core2::geometry::ElementGeometryTraits<_element, _domain>::template getNumComponents<_i>() - 1) {
                        self.template operator()<_element, _i, _j + 1>(elem_arg, self);
                    }
                    else if constexpr (_i < lolita::core2::geometry::ElementGeometryTraits<_element, _domain>::getNumComponents() - 1) {
                        self.template operator()<_element, _i + 1, 0>(elem_arg, self);
                    }
                }
            };
            /*
             * print all neighbours of a given element
             */
            auto print_element_neighbours = [&] <lolita::core2::geometry::Element _element, lolita::integer _i = 0, lolita::integer _j = 0> (
                    auto const & elem_arg,
                    auto & self
            )
                    mutable
            {
                auto const constexpr _neighbour = lolita::core2::geometry::ElementGeometryTraits<_element, _domain>::template getNeighbour<_i, _j>();
                for (auto const & c_ : elem_arg->template getNeighbours<_i, _j>()) {
                    if constexpr (!_element.isPoint() && _i == 0) {
                        os << "layer : " << _i << " type : " << _j << " <-> " << _neighbour << " " << c_->hash() << std::endl;
                    }
                    else {
                        os << "layer : " << _i << " type : " << _j << " --> " << _neighbour << " " << c_->hash() << std::endl;
                    }
                }
                if constexpr (_j < lolita::core2::geometry::ElementGeometryTraits<_element, _domain>::template getNumNeighbours<_i>() - 1) {
                    self.template operator()<_element, _i, _j + 1>(elem_arg, self);
                }
                else if constexpr (_i < lolita::core2::geometry::ElementGeometryTraits<_element, _domain>::getNumNeighbours() - 1) {
                    self.template operator()<_element, _i + 1, 0>(elem_arg, self);
                }
            };
            /*
             * print an element, with its neighbours and components
             */
            auto print_elements = [&] <lolita::integer _i = 0, lolita::integer _j = 0> (
                    auto & self
            )
                    mutable
            {
                auto constexpr _element = lolita::core2::geometry::DomainGeometryTraits<_domain>::template getElement<_i, _j>();
                if constexpr (_i == 0 && _j == 0) {
                    os << "*** Elements : " << std::endl;
                }
                for (auto const & element : mesh.elements_.template getElements<_i, _j>()) {
                    os << "* Element : " << _element << " " << element.second->hash() << " is in domains : ";
                    for (auto const & element_domain : element.second->domains_) {
                        os << * element_domain << " ";
                    }
                    os << std::endl;
                    print_element_components.template operator()<_element>(element.second, print_element_components);
                    print_element_neighbours.template operator()<_element>(element.second, print_element_neighbours);
                }
                if constexpr (_j < lolita::core2::geometry::DomainGeometryTraits<_domain>::template getNumElements<_i>() - 1) {
                    self.template operator()<_i, _j + 1>(self);
                }
                else if constexpr (_i < lolita::core2::geometry::DomainGeometryTraits<_domain>::getNumElements() - 1) {
                    self.template operator()<_i + 1, 0>(self);
                }
            };
            print_elements(print_elements);
            return os;
        }

        /**
         * @brief The mesh options, depending on the mesh format
         */
        struct Options
        {

            /**
             * @brief The offset at which starts the indexation of elements. Must be 1 or 0
             */
            lolita::index offset_ = 0;

        };

        /**
         * @brief The list of all elements in the mesh
         */
        lolita::core2::mesh::ElementCollection<_ElementPointerMap, _domain, _finite_element> elements_;

        /**
         * @brief The list of physical domains
         */
        std::vector<std::shared_ptr<std::basic_string<lolita::character>>> domains_;

        /**
         * @brief The list of loads applied on physical domains. Default is the zero natural load
         */
        std::vector<lolita::finite_element::Load> loads_;

        /**
         * @brief The list of behaviours applied on physical domains
         */
        std::vector<lolita::behaviour::MgisBehaviour> behaviours_;

        /**
         * @brief The list of linear problems available
         */
        std::array<lolita::core2::mesh::FiniteElementLinearSystem, _finite_element.count()> systems_;

        /**
         * @brief The mesh options, depending on the mesh format
         */
        Options options_;

    };

    /**
     * @brief
     * @tparam _element
     */
    template<lolita::core2::geometry::Element _element>
    struct ElementInitializationData;

    /**
     * @brief
     * @tparam _element
     */
    template<lolita::core2::geometry::Element _element>
    requires(_element.isPoint())
    struct ElementInitializationData<_element>
    {

        /**
         * @brief
         */
        lolita::natural tag_;

        /**
         * @brief
         */
        lolita::domain::Point coordinates_;

        /**
         * @brief
         */
        std::vector<std::shared_ptr<std::basic_string<lolita::character>>> domains_;

    };

    /**
     * @brief
     * @tparam _element
     */
    template<lolita::core2::geometry::Element _element>
    requires(!_element.isPoint())
    struct ElementInitializationData<_element>
    {

        /**
         * @brief
         */
        lolita::natural tag_;

        /**
         * @brief
         */
        std::array<lolita::index, _element.num_nodes_> node_tags_;

    };

    /**
     * @brief
     * @tparam _mesh
     * @tparam _domain
     * @tparam _finite_element
     */
    template<lolita::mesh::Format _mesh, lolita::domain::Domain _domain, auto _finite_element>
    struct MeshParserModule;

    /**
     * @brief
     * @tparam _mesh
     * @tparam _domain
     * @tparam _finite_element
     */
    template<lolita::mesh::Format _mesh, lolita::domain::Domain _domain, auto _finite_element>
    struct MeshParser
    {

    private:

        using _Module = typename lolita::core2::mesh::MeshParserModule<_mesh, _domain, _finite_element>::Module;

        using _Implementation = typename lolita::core2::mesh::MeshParserModule<_mesh, _domain, _finite_element>::Implementation;

    public:

        /**
         * @brief
         * @param mesh_file_path
         * @param loads
         * @param behaviours
         */
        MeshParser(
                std::basic_string_view<lolita::character> && mesh_file_path,
                std::vector<lolita::finite_element::Load> && loads,
                std::vector<lolita::behaviour::MgisBehaviour> && behaviours
        )
        {
            /*
             * simply load behaviours loads and options
             */
            this->mesh_data_.loads_ = std::forward<std::vector<lolita::finite_element::Load>>(loads);
            this->mesh_data_.options_.offset_ = 1;
            this->mesh_data_.behaviours_ = std::forward<std::vector<lolita::behaviour::MgisBehaviour>>(behaviours);
            /*
             * make the module object
             */
            auto const module = _Module(std::forward<std::basic_string_view<lolita::character>>(mesh_file_path));
            /*
             * read and load physical domains in the mesh
             */
            makeDomains(module);
            /*
             * set elements starting from the node, to the cell
             */
            auto set_elements = [&] <lolita::index _i = 0u, lolita::index _j = 0u> (auto & set_elements_imp)
                    mutable
            {
                auto const constexpr _element = lolita::core2::geometry::DomainGeometryTraits<_domain>::template getElement<_i, _j>();
                makeElement<_element>(module);
                if constexpr (_j < lolita::core2::geometry::DomainGeometryTraits<_domain>::template getNumElements<_i>() - 1) {
                    set_elements_imp.template operator()<_i, _j + 1u>(set_elements_imp);
                }
                else if constexpr (_i < lolita::core2::geometry::DomainGeometryTraits<_domain>::getNumElements() - 1) {
                    set_elements_imp.template operator()<_i + 1u, 0u>(set_elements_imp);
                }
            };
            set_elements(set_elements);
            /*
             * initialize elements, starting from the node, to the cell
             */
            auto make_elements = [&] <lolita::index _i = 0u, lolita::index _j = 0u> (auto & self)
                    mutable
            {
                auto const constexpr _element = lolita::core2::geometry::DomainGeometryTraits<_domain>::template getElement<_i, _j>();
                buildElement<_element>();
                if constexpr (_j < lolita::core2::geometry::DomainGeometryTraits<_domain>::template getNumElements<_i>() - 1) {
                    self.template operator()<_i, _j + 1u>(self);
                }
                else if constexpr (_i < lolita::core2::geometry::DomainGeometryTraits<_domain>::getNumElements() - 1) {
                    self.template operator()<_i + 1u, 0u>(self);
                }
            };
            make_elements(make_elements);
        }

        static
        lolita::core2::mesh::Mesh<_domain, _finite_element>
        makeMesh(
                std::basic_string_view<lolita::character> && mesh_file_path,
                std::vector<lolita::finite_element::Load> && loads,
                std::vector<lolita::behaviour::MgisBehaviour> && behaviours
        )
        {
            return MeshParser(mesh_file_path, loads, behaviours).mesh_data_;
        }

    protected:

        /**
         * @brief
         * @param module
         */
        void
        makeDomains(
                _Module const & module
        )
        {
            static_cast<_Implementation *>(this)->makeDomains(module);
        }

        /**
         * @brief
         * @tparam _element
         * @param module
         */
        template<lolita::core2::geometry::Element _element>
        void
        makeElement(
                _Module const & module
        )
        {
            static_cast<_Implementation *>(this)->template makeElement<_element>(module);
        }

        /**
         * @brief
         * @tparam _element
         */
        template<lolita::core2::geometry::Element _element>
        void
        buildElement()
        {
            /*
             * get element coordinates and element list
             */
//            auto const constexpr _coordinates = lolita::core2::elementPosition<_domain, _element>();
            auto const constexpr _coordinates = lolita::core2::geometry::DomainGeometryTraits<_domain>::template getElementCoordinates<_element>();
            auto & element_map = mesh_data_.elements_.template getElements<_coordinates.dim_, _coordinates.tag_>();
            /*
             * loop over each element
             */
            for (auto & element_map_item : element_map) {
                auto & element = element_map_item.second;
                /*
                 * for each element, read its neighbours its adjacent elements of same type, i.e. the other segments a segment is connected to
                 */
                auto make_element_neighbourhood_ng = [&] <lolita::index _k = 0u> (
                        auto & self
                )
                        mutable
                {
                    if constexpr (!_element.isPoint()) {
                        auto const & element_nds = element->template getComponents<_coordinates.dim_ - 1, 0>();
                        for (int i = 0; i < element_nds.size(); ++i) {
                            auto const & nde = element_nds[i];
                            auto const & ngs = nde->template getNeighbours<_coordinates.dim_ - 1, _k>();
                            for (auto const & neighbour : ngs) {
                                if (((neighbour->tag_ != element->tag_) && _k == _coordinates.tag_) || (_k != _coordinates.tag_)) {
                                    auto & element_ngs = element->template getNeighbours<0, _k>();
                                    auto found = false;
                                    for (auto const & ngb: element_ngs) {
                                        if (ngb->tag_ == neighbour->tag_) {
                                            found = true;
                                            break;
                                        }
                                    }
                                    if (!found) {
                                        element_ngs.push_back(neighbour);
                                    }
                                }
                            }
                        }
//                        if constexpr (_k < lolita::core2::numElements<_domain, _coordinates.dim_>() - 1) {
                        if constexpr (_k < lolita::core2::geometry::DomainGeometryTraits<_domain>::template getNumElements<_coordinates.dim_>() - 1) {
                            self.template operator()<_k + 1>(self);
                        }
                    }
                };
                make_element_neighbourhood_ng(make_element_neighbourhood_ng);
                /*
                 * initialize the element, once all its adjacent elements have been found
                 */
                element->template initialize<_element>(this->mesh_data_);
            }
        }

    public:

        /**
         * @brief
         */
        lolita::core2::mesh::Mesh<_domain, _finite_element> mesh_data_;

    };

    template<lolita::mesh::Format _mesh, lolita::domain::Domain _domain, auto _finite_element>
    requires(_mesh == lolita::mesh::Format::Gmsh)
    struct MeshParserModule<_mesh, _domain, _finite_element>
    {

        std::unordered_map<lolita::index, lolita::core2::geometry::Element> const element_table_ = {
                {15, lolita::core2::geometry::Element::Node()},
                {01, lolita::core2::geometry::Element::LinearSegment()},
                {02, lolita::core2::geometry::Element::LinearTriangle()},
                {03, lolita::core2::geometry::Element::LinearQuadrangle()},
        };

        static constexpr
        lolita::core2::geometry::Element
        getElemType(
                lolita::index tag
        )
        {
            if (tag == 15) return lolita::core2::geometry::Element::Node();
            else if (tag == 01) return lolita::core2::geometry::Element::LinearSegment();
            else if (tag == 02) return lolita::core2::geometry::Element::LinearTriangle();
            else if (tag == 03) return lolita::core2::geometry::Element::LinearQuadrangle();
            else return lolita::core2::geometry::Element::Node();
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

                lolita::index tag_;

                lolita::index dim_;

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

                lolita::index tag_;

                lolita::index dim_;

                std::vector<lolita::index> physical_entities_tags_;

                std::vector<lolita::index> bounding_entities_tags_;

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

                std::array<std::vector<lolita::index>, 4> geometrical_entities_tags_;

            };

            explicit
            Module(
                    std::basic_string_view<lolita::character> const &
                    mesh_file_path
            )
                    :
                    file_(mesh_file_path)
            {
                setGeometricalEntities();
                setPhysicalEntities();
                setPhysicalGroups();
            }

            explicit
            Module(
                    std::basic_string_view<lolita::character> &&
                    mesh_file_path
            )
                    :
                    file_(std::forward<std::basic_string_view<lolita::character>>(mesh_file_path))
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
                lolita::index line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Entities"));
                lolita::index offset = 1;
                std::basic_stringstream<lolita::character> line_stream(file_lines[line_start + offset]);
                lolita::index num_points;
                lolita::index num_curves;
                lolita::index num_surfaces;
                lolita::index num_volumes;
                line_stream >> num_points >> num_curves >> num_surfaces >> num_volumes;
                auto num_domains = std::array<lolita::index, 4>{num_points, num_curves, num_surfaces, num_volumes};
                offset += 1;
                for (lolita::index i = 0; i < 4; ++i) {
                    for (lolita::index j = 0; j < num_domains[i]; ++j) {
                        line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                        lolita::index tag;
                        line_stream >> tag;
                        if (i == 0) {
                            for (lolita::index k = 0; k < 3; ++k) {
                                lolita::real a;
                                line_stream >> a;
                            }
                        } else {
                            for (lolita::index k = 0; k < 6; ++k) {
                                lolita::real a;
                                line_stream >> a;
                            }
                        }
                        lolita::index num_physical_entities;
                        line_stream >> num_physical_entities;
                        std::vector<lolita::index> physical_entities_tags;
                        for (lolita::index k = 0; k < num_physical_entities; ++k) {
                            lolita::integer physical_entity_tag;
                            line_stream >> physical_entity_tag;
                            physical_entities_tags.push_back(physical_entity_tag);
                        }
                        std::vector<lolita::index> bounding_entities_tags = {};
                        if (i > 0) {
                            lolita::index num_bounding_entities;
                            line_stream >> num_bounding_entities;
                            for (lolita::index k = 0; k < num_bounding_entities; ++k) {
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
                lolita::index line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$PhysicalNames"));
                lolita::index offset = 1;
                std::basic_stringstream<lolita::character> line_stream(file_lines[line_start + offset]);
                lolita::index num_physical_names;
                line_stream >> num_physical_names;
                offset += 1;
                for (lolita::index i = 0; i < num_physical_names; ++i) {
                    line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                    lolita::index dim;
                    lolita::index tag;
                    std::basic_string<lolita::character> name;
                    line_stream >> dim >> tag >> name;
                    lolita::utility::removeCharacter(name, '"');
                    physical_entities.push_back(PhysicalEntity{tag, dim, name});
                    offset += 1;
                }
                return physical_entities;
            }

            std::array<std::unordered_set<lolita::index>, 4>
            getSubGeometricalEntities(
                    lolita::index d,
                    lolita::index t,
                    std::array<std::unordered_set<lolita::index>, 4> & a
            )
            {
                a[d].insert(t);
                for (lolita::index i = 0; i < geometrical_entities_[d][t - 1].bounding_entities_tags_.size(); ++i) {
                    lolita::index const d2 = d - 1;
                    lolita::index const t2 = geometrical_entities_[d][t - 1].bounding_entities_tags_[i];
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
                        std::array<std::unordered_set<lolita::index>, 4> & a,
                        auto & self
                )
                        mutable
                {
                    a[d].insert(t);
                    for (lolita::index i = 0; i < geometrical_entities_[d][t - 1].bounding_entities_tags.size(); ++i) {
                        a = self(d - 1, geometrical_entities_[d][t - 1].bounding_entities_tags[i], a, self);
                    }
                    return a;
                };
                /*
                 *
                 */
                auto physical_groups = std::vector<PhysicalGroup>();
                for (lolita::index i = 0; i < physical_entities_.size(); ++i) {
                    std::array<std::unordered_set<lolita::index>, 4> tags;
                    for (lolita::index j = 0; j < 4; ++j) {
                        for (lolita::index k = 0; k < geometrical_entities_[j].size(); ++k) {
                            for (lolita::index l = 0; l < geometrical_entities_[j][k].physical_entities_tags_.size(); ++l) {
                                lolita::index const & t0 = geometrical_entities_[j][k].physical_entities_tags_[l];
                                lolita::index const & t1 = physical_entities_[i].tag_;
                                if (t0 == t1) {
//                                    tags = get_subs(j, k + 1, tags, get_subs);
                                    tags = getSubGeometricalEntities(j, k + 1, tags);
                                }
                            }
                        }
                    }
                    std::array<std::vector<lolita::index>, 4> group_tags;
                    for (lolita::index j = 0; j < 4; ++j) {
                        group_tags[j].assign(tags[j].begin(), tags[j].end());
                    }
                    physical_groups.push_back(PhysicalGroup{physical_entities_[i].name_, group_tags});
                }
                return physical_groups;
            }

            void
            setGeometricalEntities()
            {
                auto const & file_lines = file_.lines_;
                lolita::index line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Entities"));
                lolita::index offset = 1;
                std::basic_stringstream<lolita::character> line_stream(file_lines[line_start + offset]);
                lolita::index num_points;
                lolita::index num_curves;
                lolita::index num_surfaces;
                lolita::index num_volumes;
                line_stream >> num_points >> num_curves >> num_surfaces >> num_volumes;
                auto num_domains = std::array<lolita::index, 4>{num_points, num_curves, num_surfaces, num_volumes};
                offset += 1;
                for (lolita::index i = 0; i < 4; ++i) {
                    for (lolita::index j = 0; j < num_domains[i]; ++j) {
                        line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                        lolita::index tag;
                        line_stream >> tag;
                        if (i == 0) {
                            for (lolita::index k = 0; k < 3; ++k) {
                                lolita::real a;
                                line_stream >> a;
                            }
                        } else {
                            for (lolita::index k = 0; k < 6; ++k) {
                                lolita::real a;
                                line_stream >> a;
                            }
                        }
                        lolita::index num_physical_entities;
                        line_stream >> num_physical_entities;
                        std::vector<lolita::index> physical_entities_tags;
                        for (lolita::index k = 0; k < num_physical_entities; ++k) {
                            lolita::integer physical_entity_tag;
                            line_stream >> physical_entity_tag;
                            physical_entities_tags.push_back(physical_entity_tag);
                        }
                        std::vector<lolita::index> bounding_entities_tags = {};
                        if (i > 0) {
                            lolita::index num_bounding_entities;
                            line_stream >> num_bounding_entities;
                            for (lolita::index k = 0; k < num_bounding_entities; ++k) {
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
                lolita::index line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$PhysicalNames"));
                lolita::index offset = 1;
                std::basic_stringstream<lolita::character> line_stream(file_lines[line_start + offset]);
                lolita::index num_physical_names;
                line_stream >> num_physical_names;
                offset += 1;
                for (lolita::index i = 0; i < num_physical_names; ++i) {
                    line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                    lolita::index dim;
                    lolita::index tag;
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
                auto get_subs = [&] (
                        auto d,
                        auto t,
                        std::array<std::unordered_set<lolita::index>, 4> & a,
                        auto & self
                )
                        mutable
                {
                    a[d].insert(t);
                    for (lolita::index i = 0; i < geometrical_entities_[d][t - 1].bounding_entities_tags.size(); ++i) {
                        a = self(d - 1, geometrical_entities_[d][t - 1].bounding_entities_tags[i], a, self);
                    }
                    return a;
                };
                /*
                 *
                 */
                for (lolita::index i = 0; i < physical_entities_.size(); ++i) {
                    std::array<std::unordered_set<lolita::index>, 4> tags;
                    for (lolita::index j = 0; j < 4; ++j) {
                        for (lolita::index k = 0; k < geometrical_entities_[j].size(); ++k) {
                            for (lolita::index l = 0; l < geometrical_entities_[j][k].physical_entities_tags_.size(); ++l) {
                                lolita::index const & t0 = geometrical_entities_[j][k].physical_entities_tags_[l];
                                lolita::index const & t1 = physical_entities_[i].tag_;
                                if (t0 == t1) {
//                                    tags = get_subs(j, k + 1, tags, get_subs);
                                    tags = getSubGeometricalEntities(j, k + 1, tags);
                                }
                            }
                        }
                    }
                    std::array<std::vector<lolita::index>, 4> group_tags;
                    for (lolita::index j = 0; j < 4; ++j) {
                        group_tags[j].assign(tags[j].begin(), tags[j].end());
                    }
                    physical_groups_.push_back(PhysicalGroup{physical_entities_[i].name_, group_tags});
                }
            }

            lolita::utility::File file_;

            std::array<std::vector<GeometricalEntity>, 4> geometrical_entities_;

            std::vector<PhysicalEntity> physical_entities_;

            std::vector<PhysicalGroup> physical_groups_;

        };

        struct Implementation : public MeshParser<_mesh, _domain, _finite_element>
        {

            void
            makeDomains(
                    Module const & module
            )
            {
                auto const & physical_groups = module.physical_groups_;
                for (int k = 0; k < physical_groups.size(); ++k) {
                    this->mesh_data_.domains_.push_back(std::make_shared<std::basic_string<lolita::character>>(physical_groups[k].name_));
                }
            }

            template<lolita::core2::geometry::Element _element>
            void
            makeElement(
                    Module const & module
            );

            template<lolita::core2::geometry::Element _element>
            void
            makeElement(
                    Module const & module
            )
            requires(_element.isPoint())
            {
//                using __Element = lolita::core2::finite_element::FiniteElementFinal<_element, _domain, _finite_element...>;
                using __Element = lolita::core2::finite_element::FEObject<_element, _domain, _finite_element>;
//                auto const constexpr _element_coordinates2 = lolita::core2::elementPosition<_domain, _element>();
//                auto const constexpr _element_coordinates2 = lolita::core2::MeshDescription<_domain>::template getElementCoordinates<_element>();
                auto const & file_lines = module.file_.lines_;
                auto const & physical_groups = module.physical_groups_;
                auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Nodes"));
                auto offset = 1;
                auto line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                auto num_entity_block = lolita::index();
                auto num_nodes = lolita::index();
                auto min_node_tag = lolita::index();
                auto max_node_tag = lolita::index();
                line_stream >> num_entity_block >> num_nodes >> min_node_tag >> max_node_tag;
                offset += 1;
                for (int i = 0; i < num_entity_block; ++i) {
                    line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                    auto entity_dim = lolita::index();
                    auto entity_tag = lolita::index();
                    auto parametric = lolita::index();
                    auto num_nodes_in_block = lolita::index();
                    line_stream >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
                    offset += 1;
                    for (int j = 0; j < num_nodes_in_block; ++j) {
                        auto element_initialization_data = lolita::core2::mesh::ElementInitializationData<_element>();
                        line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                        line_stream >> element_initialization_data.tag_;
                        line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset + num_nodes_in_block]);
                        element_initialization_data.coordinates_.setZero();
                        for (int k = 0; k < 3; ++k) {
                            line_stream >> element_initialization_data.coordinates_(k);
                        }
                        for (int k = 0; k < physical_groups.size(); ++k) {
                            auto const & node_set_name = physical_groups[k].name_;
                            for (int l = 0; l < physical_groups[k].geometrical_entities_tags_[entity_dim].size(); ++l) {
                                lolita::index tag = physical_groups[k].geometrical_entities_tags_[entity_dim][l];
                                if (entity_tag == tag) {
                                    auto is_equal = [&] (std::shared_ptr<std::basic_string<lolita::character>> const & ptr_domain) {
                                        return * ptr_domain == node_set_name;
                                    };
                                    auto const & domains = this->mesh_data_.domains_;
                                    auto domain_index = std::distance(domains.begin(), std::find_if(domains.begin(), domains.end(), is_equal));
//                                    node_data2.domains_.push_back(node_set_name);
                                    element_initialization_data.domains_.push_back(domains[domain_index]);
                                }
                            }
                        }
                        __Element::template makeElement<_element>(element_initialization_data, this->mesh_data_);
                        offset += 1;
                    }
                    offset += num_nodes_in_block;
                }
            }

            template<lolita::core2::geometry::Element _element>
            void
            makeElement(
                    Module const &
                    module
            )
            requires(!_element.isPoint())
            {
//                using __Element = lolita::core2::finite_element::FiniteElementFinal<_element, _domain, _finite_element...>;
                using __Element = lolita::core2::finite_element::FEObject<_element, _domain, _finite_element>;
//                auto const constexpr _element_coordinates = lolita::core2::elementPosition<_domain, _element>();
                auto const constexpr _element_coordinates = lolita::core2::geometry::DomainGeometryTraits<_domain>::template getElementCoordinates<_element>();
                auto const & file_lines = module.file_.lines_;
                auto & elements = this->mesh_data_.elements_.template getElements<_element_coordinates.dim_, _element_coordinates.tag_>();
                auto line_start = std::distance(file_lines.begin(), std::find(file_lines.begin(), file_lines.end(), "$Elements"));
                auto offset = 1;
                auto line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                auto num_entity_blocks = 0;
                auto num_elements = 0;
                auto min_element_tag = 0;
                auto max_element_tag = 0;
                line_stream >> num_entity_blocks >> num_elements >> min_element_tag >> max_element_tag;
                offset += 1;
                for (int i = 0; i < num_entity_blocks; ++i) {
                    line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                    auto entity_dim = 0;
                    auto entity_tag = 0;
                    auto element_type_tag = 0;
                    auto num_elements_in_block = 0;
                    line_stream >> entity_dim >> entity_tag >> element_type_tag >> num_elements_in_block;
                    offset += 1;
                    auto const element = getElemType(element_type_tag);
                    if (entity_dim == _domain.dim_ && element == _element) {
                        for (int j = 0; j < num_elements_in_block; ++j) {
                            auto element_initialization_data = lolita::core2::mesh::ElementInitializationData<_element>();
                            line_stream = std::basic_stringstream<lolita::character>(file_lines[line_start + offset]);
                            auto tag = lolita::index();
                            line_stream >> tag;
                            for (lolita::index k = 0; k < element.num_nodes_; ++k) {
                                line_stream >> element_initialization_data.node_tags_[k];
                            }
                            element_initialization_data.tag_ = elements.size();
                            __Element::template makeElement<_element>(element_initialization_data, this->mesh_data_);
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

}

#endif //LOLITA_LOLITA_CORE_4_HXX
