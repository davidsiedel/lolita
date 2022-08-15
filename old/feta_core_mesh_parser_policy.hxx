//
// Created by dsiedel on 14/03/2022.
//

#ifndef FETA_FETA_CORE_MESH_PARSER_POLICY_HXX
#define FETA_FETA_CORE_MESH_PARSER_POLICY_HXX

#include <execution>

#include "new/_feta.hxx"
#include "new/_feta_array.hxx"
#include "new/_feta_files.hxx"
#include "new/feta_feta_unordered_set.hxx"
//#include "feta/feta/functions.hxx"

namespace feta::core
{

    enum struct MeshFormatType
    {

        Gmsh

    };

    template<MeshFormatType>
    struct ParserHelper;

    template<>
    struct ParserHelper<MeshFormatType::Gmsh>
    {

        struct PhysicalEntity
        {

            Indx tag;
            Indx dim;
            Strg name;

            Bool
            operator==(
                    PhysicalEntity const &
                    other
            )
            const
            {
                return tag == other.tag && dim == other.dim && name == other.name;
                auto eq_0 = tag == other.tag;
                auto eq_1 = dim == other.dim;
                auto eq_2 = name == other.name;
                return eq_0 && eq_1 && eq_2;
            }

            Bool
            operator!=(
                    PhysicalEntity const &
                    other
            )
            const
            {
                return !(other == * this);
            }

        };

        struct GeometricalEntity
        {

            Indx tag;
            Indx dim;
            Array<Indx> physical_entities_tags;
            Array<Indx> bounding_entities_tags;

            Bool
            operator==(
                    GeometricalEntity const &
                    other
            )
            const
            {
                auto eq_0 = tag == other.tag;
                auto eq_1 = dim == other.dim;
                auto eq_2 = physical_entities_tags == other.physical_entities_tags;
                auto eq_3 = bounding_entities_tags == other.bounding_entities_tags;
                return eq_0 && eq_1 && eq_2 && eq_3;
            }

            Bool
            operator!=(
                    GeometricalEntity const &
                    other
            )
            const
            {
                return !(other == * this);
            }

        };

        struct PhysicalGroup
        {

            Strg name;
            Array<Array<Indx>, 4> geometrical_entities_tags;

            Bool
            operator==(
                    PhysicalGroup const &
                    other
            )
            const
            {
                auto eq_0 = name == other.name;
                auto eq_1 = geometrical_entities_tags == other.geometrical_entities_tags;
                return eq_0 && eq_1;
            }

            Bool
            operator!=(
                    PhysicalGroup const &
                    other
            )
            const
            {
                return !(other == * this);
            }

        };

        using PhysicalGroups = Array<PhysicalGroup>;
        using GeometricalEntities = Array<Array<GeometricalEntity>, 4>;
        using PhysicalEntities = Array<PhysicalEntity>;

        GeometricalEntities geometrical_entities;
        PhysicalEntities physical_entities;
        PhysicalGroups physical_groups;

        explicit
        ParserHelper(
                file::File<file::FileType::Input> const & mesh_file
        )
        :
        geometrical_entities(getGeometricalEntities(mesh_file)),
        physical_entities(getPhysicalEntities(mesh_file)),
        physical_groups(getPhysicalGroups(mesh_file))
        {}

        Bool
        operator==(
                ParserHelper const &
                other
        )
        const
        {
            auto eq_0 = geometrical_entities == other.geometrical_entities;
            auto eq_1 = physical_entities == other.physical_entities;
            auto eq_2 = physical_groups == other.physical_groups;
            return eq_0 && eq_1 && eq_2;
        }

        Bool
        operator!=(
                ParserHelper const &
                other
        )
        const
        {
            return !(other == * this);
        }

    private:

        static
        GeometricalEntities
        getGeometricalEntities(
                file::File<file::FileType::Input> const & mesh_file
        )
        {
            GeometricalEntities geometrical_entities;
            Indx line_start = mesh_file.getLineIndex("$Entities");
            Indx offset = 1;
            StrgStream line_stream(mesh_file.getLine(line_start + offset));
            Indx num_points;
            Indx num_curves;
            Indx num_surfaces;
            Indx num_volumes;
            line_stream >> num_points >> num_curves >> num_surfaces >> num_volumes;
            Array<Indx, 4> num_domains = {
                    num_points,
                    num_curves,
                    num_surfaces,
                    num_volumes
            };
            offset += 1;
            for (Indx i = 0; i < 4; ++i) {
                for (Indx j = 0; j < num_domains(i); ++j) {
                    line_stream = StrgStream(mesh_file.getLine(line_start + offset));
                    Indx tag;
                    line_stream >> tag;
                    if (i == 0) {
                        for (Indx k = 0; k < 3; ++k) {
                            Real a;
                            line_stream >> a;
                        }
                    } else {
                        for (Indx k = 0; k < 6; ++k) {
                            Real a;
                            line_stream >> a;
                        }
                    }
                    Indx num_physical_entities;
                    line_stream >> num_physical_entities;
                    Array<Indx> physical_entities_tags;
                    for (Indx k = 0; k < num_physical_entities; ++k) {
                        Intg physical_entity_tag;
                        line_stream >> physical_entity_tag;
                        physical_entities_tags.data.push_back(physical_entity_tag);
                    }
                    Array<Indx> bounding_entities_tags = {};
                    if (i > 0) {
                        Indx num_bounding_entities;
                        line_stream >> num_bounding_entities;
                        for (Indx k = 0; k < num_bounding_entities; ++k) {
                            Intg bounding_entity_tag;
                            line_stream >> bounding_entity_tag;
                            bounding_entities_tags.data.push_back(std::abs(bounding_entity_tag));
                        }
                    }
                    geometrical_entities(i).data.push_back(GeometricalEntity{
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

        static
        PhysicalEntities
        getPhysicalEntities(
                file::File<file::FileType::Input> const & mesh_file
        )
        {
            PhysicalEntities physical_entities;
            Indx line_start = mesh_file.getLineIndex("$PhysicalNames");
            Indx offset = 1;
            StrgStream line_stream(mesh_file.getLine(line_start + offset));
            Indx num_physical_names;
            line_stream >> num_physical_names;
            offset += 1;
            for (Indx i = 0; i < num_physical_names; ++i) {
                line_stream = StrgStream(mesh_file.getLine(line_start + offset));
                Indx dim;
                Indx tag;
                Strg name;
                line_stream >> dim >> tag >> name;
                file::removeCharacter(name, '"');
                physical_entities.data.push_back(PhysicalEntity{tag, dim, name});
                offset += 1;
            }
            return physical_entities;
        }

        Array<UnorderedSet<Indx>, 4>
        getSubGeometricalEntities(
                Indx const & d,
                Indx const & t,
                Array<UnorderedSet<Indx>, 4> & a
        )
        {
            a(d).data.insert(t);
            for (Indx i = 0; i < geometrical_entities(d)(t - 1).bounding_entities_tags.getSize(); ++i) {
                Indx const d2 = d - 1;
                Indx const t2 = geometrical_entities(d)(t - 1).bounding_entities_tags(i);
                a = getSubGeometricalEntities(d2, t2, a);
            }
            return a;
        }

        PhysicalGroups
        getPhysicalGroups(
                file::File<file::FileType::Input> const & mesh_file
        )
        {
            PhysicalGroups groups;
            for (Indx i = 0; i < physical_entities.getSize(); ++i) {
                Array<UnorderedSet<Indx>, 4> tags;
                for (Indx j = 0; j < 4; ++j) {
                    for (Indx k = 0; k < geometrical_entities(j).getSize(); ++k) {
                        for (Indx l = 0; l < geometrical_entities(j)(k).physical_entities_tags.getSize(); ++l) {
                            Indx const & t0 = geometrical_entities(j)(k).physical_entities_tags(l);
                            Indx const & t1 = physical_entities(i).tag;
                            if (t0 == t1) {
                                tags = getSubGeometricalEntities(j, k + 1, tags);
                            }
                        }
                    }
                }
                Array<Array<Indx>, 4> group_tags;
                for (Indx j = 0; j < 4; ++j) {
                    group_tags(j).data.assign(tags(j).data.begin(), tags(j).data.end());
                }
                groups.data.push_back(PhysicalGroup{physical_entities(i).name, group_tags});
            }
            return groups;
        }

    };

}

#endif //FETA_FETA_CORE_MESH_PARSER_POLICY_HXX
