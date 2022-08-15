//
// Created by dsiedel on 14/03/2022.
//

#ifndef FETA_FETA_CORE_MESH_PARSER_INTERFACE_HXX
#define FETA_FETA_CORE_MESH_PARSER_INTERFACE_HXX

//#include "feta/core/mesh_parser.hxx"
#include "new/feta_core_mesh_parser_base.hxx"
#include "new/feta_core_mesh_parser_policy.hxx"

namespace feta::core
{

//    template<typename, Indx, typename ...>
//    struct ParserInterface;

//    template<MeshFormatType, Indx, typename ...>
//    struct ParserInterface;

    template<MeshFormatType, auto>
    struct ParserImplementation;

    template<MeshFormatType Mft, auto Med>
    struct ParserInterface : public ParserBase<Med>
    {

        ParserInterface(
                file::File<file::FileType::Input> const &
                mesh_file,
                Array<MeshInteriorDomain<Med>> &&
                domains
        )
        {
            ParserHelper<Mft> const auxiliary_data(mesh_file);
            this->setDomainPointers(domains);
            this->setNodes(mesh_file);
            this->initializeElementSets(mesh_file, auxiliary_data);
            for (auto const & [key, val] : this->element_sets.data) {
                this->setNodesSet(mesh_file, key, auxiliary_data);
            }
            this->setCells(mesh_file);
            this->setNeighbours();
        }

    protected:

        void
        initializeElementSets(
                file::File<file::FileType::Input> const & mesh_file,
                ParserHelper<Mft> const & helper
        )
        {
            static_cast<ParserImplementation<Mft, Med> *>(this)->initializeElementSets(mesh_file, helper);
        }

        void
        setNodesSet(
                file::File<file::FileType::Input> const & mesh_file,
                Strg element_set_name,
                ParserHelper<Mft> const & helper
        )
        {
            static_cast<ParserImplementation<Mft, Med> *>(this)->setNodesSet(mesh_file, element_set_name, helper);
        }

        void
        setNodes(
                file::File<file::FileType::Input> const & mesh_file
        )
        {
            static_cast<ParserImplementation<Mft, Med> *>(this)->setNodes(mesh_file);
        }

        template<Indx I = 0>
        void
        setCells(
                file::File<file::FileType::Input> const & mesh_file
        )
        {
            static_cast<ParserImplementation<Mft, Med> *>(this)->template setCells<I>(mesh_file);
            auto const constexpr ddd = ParserBase<Med>::ElementPointers::template Type<Med.dim_euclidean>::getSize();
            if constexpr (I < ddd - 1) setCells<I + 1>(mesh_file);
        }

    };

//    template<typename Implementation, Indx D>
//    struct ParserInterface<Implementation, D>
//            :
//                    public ParserBase<D>
//    {
//
//        using BaseType = ParserBase<D>;
//
//        explicit
//        ParserInterface(
//                file::File<file::FileType::Input> const & mesh_file
//        )
//        {
//            this->setNodes(mesh_file);
//            this->initializeElementSets(mesh_file);
//            for (auto const & [key, val] : this->element_sets) {
//                this->setNodesSet(mesh_file, key);
//            }
//            this->setCells(mesh_file);
//            this->setContiguous();
//            this->removeIsolatedNodes();
//        }
//
//    protected:
//
//        void
//        setNodes(
//                file::File<file::FileType::Input> const & mesh_file
//        )
//        {
//            static_cast<Implementation *>(this)->setNodes(mesh_file);
//        }
//
//        void
//        initializeElementSets(
//                file::File<file::FileType::Input> const & mesh_file
//        )
//        {
//            static_cast<Implementation *>(this)->initializeElementSets(mesh_file);
//        }
//
//        void
//        setNodesSet(
//                file::File<file::FileType::Input> const & mesh_file,
//                Strg element_set_name
//        )
//        {
//            static_cast<Implementation *>(this)->setNodesSet(mesh_file, element_set_name);
//        }
//
//        template<Indx I = 0>
//        void
//        setCells(
//                file::File<file::FileType::Input> const & mesh_file
//        )
//        {
//            static_cast<Implementation *>(this)->template setCells<I>(mesh_file);
//            if constexpr (I < BaseType::Elements::template Type<D>::size - 1) setCells<I + 1>(mesh_file);
//        }
//
//    };

}

#endif //FETA_FETA_CORE_MESH_PARSER_INTERFACE_HXX
