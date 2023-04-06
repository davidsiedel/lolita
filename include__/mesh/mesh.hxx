/**
 * @file mesh.hxx
 * @author David Siedel (davidsiedel@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2023-04-06
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef D678D463_2196_4839_A0B5_7140F2F48756
#define D678D463_2196_4839_A0B5_7140F2F48756

#include "geometry/frame.hxx"
#include "geometry/domain.hxx"
#include "geometry/point.hxx"
#include "geometry/shape.hxx"
#include "mesh/region.hxx"
#include "mesh/element.hxx"
#include "mesh/gmsh.hxx"

namespace lolita::mesh
{

    using MeshFormats = std::tuple<
        Msh
    >;

    struct MeshFormatLibrary
    {
        
        template<Integer i_>
        using MeshFormat = std::tuple_element_t<i_, MeshFormats>;

        static constexpr
        Integer
        getNumComponents()
        {
            return std::tuple_size_v<MeshFormats>;
        }

        static
        void
        apply(
            auto const & fun
        )
        {
            auto apply = [&] <Integer i_ = 0> (
                auto & apply_
            )
            mutable
            {
                fun.template operator()<MeshFormat<i_>>();
                if constexpr (i_ < getNumComponents() - 1)
                {
                    apply_.template operator()<i_ + 1>(apply_);
                }
            };
            apply(apply);
        }

    };
    
    template<geometry::FrameConcept Frame_>
    struct MyMesh
    {        

        using Mesh_ = std::unique_ptr<Mesh<Frame_>>;

        static
        Mesh_
        readMesh(
            std::basic_string<Character> && file_path
        )
        {
            Mesh_ mesh;
            auto set_mesh = [&] <mesh::ReaderConcept MeshReader_> ()
            {
                auto mesh_format_tag = "." + std::basic_string<Character>(MeshReader_::getFileNameExtension());
                if (std::forward<std::basic_string<Character>>(file_path).find(mesh_format_tag) != std::string::npos)
                {
                    mesh = MeshFactory<Frame_, MeshReader_>().makeMesh(std::forward<std::basic_string<Character>>(file_path));
                }
            };
            MeshFormatLibrary::apply(set_mesh);
            if (mesh == nullptr)
            {
                throw std::runtime_error("Mesh format not implemented");
            }
            return mesh;
        }

        explicit
        MyMesh(
            std::basic_string<Character> && file_path
        )
        :
        mesh_(readMesh(std::forward<std::basic_string<Character>>(file_path)))
        {}

        void
        working()
        {
            std::cout << "first node tag : " << mesh_->template getElements<geometry::Node>()[0]->getTag() << std::endl;
        }

    private:

        Mesh_ mesh_;

    };

} // namespace lolita::mesh


#endif /* D678D463_2196_4839_A0B5_7140F2F48756 */
