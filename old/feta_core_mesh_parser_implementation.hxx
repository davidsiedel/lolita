//
// Created by dsiedel on 14/03/2022.
//

#ifndef FETA_FETA_CORE_MESH_PARSER_IMPLEMENTATION_HXX
#define FETA_FETA_CORE_MESH_PARSER_IMPLEMENTATION_HXX

#include "new/feta_core_mesh_parser_interface.hxx"

namespace feta::core
{

//    enum struct MeshFormatType
//    {
//
//        Gmsh
//
//    };

    // *****************************************************************************************************************
    // ******** GMSH AUXILIARY DATA
    // *****************************************************************************************************************

    template<auto Med>
    struct ParserImplementation<MeshFormatType::Gmsh, Med> : public ParserInterface<MeshFormatType::Gmsh, Med>
    {

        using Base = ParserInterface<MeshFormatType::Gmsh, Med>;

        using Base::Base;

        template<elt::ElementDescription F>
        using EEl = typename Base::template Element<F>;

        static
        elt::ElementDescription
        getElementType(
                Indx
                tag_arg
        )
        {
            if (tag_arg == 15) {
                return elt::pnt_0;
            }
            else if (tag_arg == 1) {
                return elt::seg_2;
            }
            else if (tag_arg == 2) {
                return elt::tri_3;
            }
            else if (tag_arg == 3) {
                return elt::qua_4;
            }
            else {
                assert(false);
            }
        }

        void
        setNodes(
                file::File<file::FileType::Input> const &
                mesh_file
        )
        {
            auto & nodes = this->element_collection.template get<0>().template get<0>();
            Indx line_start = mesh_file.getLineIndex("$Nodes");
            Indx offset = 1;
            StrgStream line_stream(mesh_file.getLine(line_start + offset));
            Indx num_entity_block;
            Indx num_nodes;
            Indx min_node_tag;
            Indx max_node_tag;
            line_stream >> num_entity_block >> num_nodes >> min_node_tag >> max_node_tag;
            offset += 1;
            for (Indx i = 0; i < num_entity_block; ++i) {
                line_stream = StrgStream(mesh_file.getLine(line_start + offset));
                Indx entity_dim;
                Indx entity_tag;
                Indx parametric;
                Indx num_nodes_in_block;
                line_stream >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
                offset += 1;
                for (Indx j = 0; j < num_nodes_in_block; ++j) {
                    line_stream = StrgStream(mesh_file.getLine(line_start + offset));
                    Indx node_tag;
                    line_stream >> node_tag;
                    line_stream = StrgStream(mesh_file.getLine(line_start + offset + num_nodes_in_block));
                    Matrix<Real, Med.dim_euclidean> coordinates;
                    for (int k = 0; k < Med.dim_euclidean; ++k) {
                        line_stream >> coordinates(k);
                    }
                    this->setElement(node_tag, coordinates);
                    offset += 1;
                }
                offset += num_nodes_in_block;
            }
        }

        void
        initializeElementSets(
                file::File<file::FileType::Input> const &
                mesh_file,
                ParserHelper<MeshFormatType::Gmsh> const &
                helper
        )
        {
            auto const & physical_groups = helper.physical_groups;
            for (int i = 0; i < physical_groups.getSize(); ++i) {
                this->element_sets.get(physical_groups(i).name);
            }
        }

        void
        setNodesSet(
                file::File<file::FileType::Input> const &
                mesh_file,
                Strg
                element_set_name,
                ParserHelper<MeshFormatType::Gmsh> const &
                helper
        )
        {
            auto const & physical_groups = helper.physical_groups;
            Indx physical_index;
            for (int i = 0; i < physical_groups.getSize(); ++i) {
                if (physical_groups(i).name == element_set_name) {
                    physical_index = i;
                }
            }
            auto & element_set_hashes = this->element_sets.get(element_set_name).template get<0>().template get<0>();
            auto & nodes = this->element_collection.template get<0>().template get<0>();
            for (Indx i = 0; i < 4; ++i) {
                for (Indx j = 0; j < physical_groups(physical_index).geometrical_entities_tags(i).getSize(); ++j) {
                    Indx tag = physical_groups(physical_index).geometrical_entities_tags(i)(j);
                    Indx line_start = mesh_file.getLineIndex("$Nodes");
                    Indx offset = 1;
                    StrgStream line_stream(mesh_file.getLine(line_start + offset));
                    Indx num_entity_block;
                    Indx num_nodes;
                    Indx min_node_tag;
                    Indx max_node_tag;
                    line_stream >> num_entity_block >> num_nodes >> min_node_tag >> max_node_tag;
                    offset += 1;
                    for (Indx k = 0; k < num_entity_block; ++k) {
                        line_stream = StrgStream(mesh_file.getLine(line_start + offset));
                        Indx entity_dim;
                        Indx entity_tag;
                        Indx parametric;
                        Indx num_nodes_in_block;
                        line_stream >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
                        offset += 1;
                        if (entity_dim == i && entity_tag == tag) {
                            for (Indx l = 0; l < num_nodes_in_block; ++l) {
                                line_stream = StrgStream(mesh_file.getLine(line_start + offset));
//                                hashes.data.insert(line_stream.str());
                                Strg node_hash = line_stream.str();
                                element_set_hashes.data.insert({node_hash, nodes.get(node_hash)});
//                                nodes(node_hash)().domain_names.data.push_back(element_set_name);
                                offset += 1;
                            }
                        } else {
                            offset += num_nodes_in_block;
                        }
                        offset += num_nodes_in_block;
                    }
                }
            }
        }

        template<Indx I = 0>
        void
        setCells(
                file::File<file::FileType::Input> const & mesh_file
        )
        {
            auto const constexpr dim_euclidean = Med.dim_euclidean;
            using Elements = typename core::Elements<dim_euclidean, EEl>;
            using ElementTraits = typename Elements::template Type<dim_euclidean>::template Type<I>;
            auto const constexpr elt = ElementTraits::element_description;
            Indx line_start = mesh_file.getLineIndex("$Elements");
            Indx offset = 1;
            StrgStream line_stream(mesh_file.getLine(line_start + offset));
            Indx num_entity_blocks;
            Indx num_elements;
            Indx min_element_tag;
            Indx max_element_tag;
            line_stream >> num_entity_blocks >> num_elements >> min_element_tag >> max_element_tag;
            offset += 1;
            for (Indx i = 0; i < num_entity_blocks; ++i) {
                line_stream = StrgStream(mesh_file.getLine(line_start + offset));
                Indx entity_dim = 0;
                Indx entity_tag = 0;
                Indx element_type_tag = 0;
                Indx num_elements_in_block = 0;
                line_stream >> entity_dim >> entity_tag >> element_type_tag >> num_elements_in_block;
                offset += 1;
                elt::ElementDescription const constexpr element_type = ElementTraits::element_description;
                elt::ElementDescription const element_type_arg = getElementType(element_type_tag);
                if (entity_dim == dim_euclidean and getElementType(element_type_tag) == element_type) {
                    for (Indx j = 0; j < num_elements_in_block; ++j) {
                        line_stream = StrgStream(mesh_file.getLine(line_start + offset));
                        Indx tag;
                        line_stream >> tag;
                        Array<Indx, ElementTraits::num_nodes> node_tags;
                        for (Indx k = 0; k < ElementTraits::num_nodes; ++k) {
                            line_stream >> node_tags(k);
                        }
                        typename EEl<elt>::Components cell_boundary;
                        SharedPointer<EEl<elt>> ptr_cell = SharedPointer<EEl<elt>>(EEl<elt>());
                        this->template setElement<elt>(node_tags, cell_boundary, ptr_cell);
//                        this->template setElement<CurrentElement::type>(node_tags, cell_boundary);
                        offset += 1;
                    }
                } else {
                    offset += num_elements_in_block;
                }
            }
        }

    };

}

#endif //FETA_FETA_CORE_MESH_PARSER_IMPLEMENTATION_HXX
