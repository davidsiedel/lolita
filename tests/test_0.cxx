#include "gtest/gtest.h"
#include "lolita/lolita_core_n_11.hxx"

TEST(t0, t0)
{
    
    auto constexpr domain = lolita2::Domain(2, lolita2::Domain::Cartesian);
    auto constexpr field = lolita2::Field(2);
    auto constexpr basis = lolita2::Basis(lolita2::Basis::Monomial, 2);


    auto constexpr u1 = lolita2::Unknown(lolita2::Field(1), lolita2::Mapping(lolita2::Mapping::Gradient));
    auto constexpr u2 = lolita2::Unknown(lolita2::Field(2), lolita2::Mapping(lolita2::Mapping::Gradient), lolita2::Mapping(lolita2::Mapping::Identity));

    auto constexpr fe1 =  lolita2::HybridDiscontinuousGalerkin(u1, basis, basis);
    auto constexpr fe2 =  lolita2::HybridDiscontinuousGalerkin(u2, basis, basis);

    auto constexpr mesh_fmt = lolita2::MeshFileFormat(lolita2::MeshFileFormat::Gmsh);
    
    auto file_path = "/home/dsiedel/projetcs/lolita/lolita/applications/data/meshes/unit_square_3_cpp.msh";
    
    // auto mshh = lolita2::geometry::MeshFileParser<mesh_fmt, domain>::readMesh(file_path);

    // auto kkk = lolita2::geometry::MeshSetBase2<lolita2::geometry::FiniteElementHolder, domain, fe1, fe2>(mshh);

    // for (auto const & elem : mshh.elements_.getElements<0, 0>())
    // {
    //     std::cout << elem.tag_ << std::endl;
    //     std::cout << * elem.coordinates_ << std::endl;
    //     for (auto const & dom : elem.domains_)
    //     {
    //         std::cout << dom->tag_ << std::endl;
    //     }
        
    // }
    // for (auto const & elem : mshh.elements_.getElements<1, 0>())
    // {
    //     std::cout << "SEGS :" << std::endl;
    //     for (auto val : elem.node_tags_)
    //     {
    //         std::cout << val << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // for (auto const & elem : mshh.elements_.getElements<2, 0>())
    // {
    //     std::cout << "TRIS :" << std::endl;
    //     for (auto val : elem.node_tags_)
    //     {
    //         std::cout << val << " ";
    //     }
    //     std::cout << std::endl;
    // }

    auto gmshfile = lolita2::geometry::NewGmsh(file_path);

    // auto set = std::set<std::string>();

    // for (auto const & [h1, item] : gmshfile.geometrical_entities_)
    // {
    //     std::cout << "---" << std::endl;
    //     for (auto const & physent : item->physical_entities_)
    //     {
    //         std::cout << physent->name_ << std::endl;
    //     }
    // }

    auto lmpo = lolita2::geometry::ICILA<domain, fe1, fe2>::makeIt(file_path);
    std::cout << lmpo << std::endl;
    

    // auto kkk = lolita2::geometry::MeshSetBase2<lolita2::geometry::FiniteElementHolder, domain, fe1, fe2>(mshh);
    // auto kkk = mshh.template make<lolita2::geometry::FiniteElementHolder, domain, fe1, fe2>();

    // std::cout << kkk << std::endl;

    // for (auto const & [hash, elem] : kkk.elements_.getElements<2, 0>())
    // {
    //     auto a = 2;
    //     elem->getElement<0>();
    //     elem->hash();
    //     std::cout << elem->hash() << std::endl;
    //     // std::cout << elem.getElement<0>()hash() << std::endl;
    //     // std::cout << * elem.coordinates_ << std::endl;
    //     // for (auto const & dom : elem.domains_)
    //     // {
    //     //     std::cout << dom->tag_ << std::endl;
    //     // }        
    // }
    
    // auto msh = lolita2::geometry::MeshParser<lolita2::MeshFileFormat::Gmsh, domain, fe1, fe2>::makeMesh(file_path);
    // std::cout << msh << std::endl;
    // auto msh = lolita2::geometry::MyMeshParser<mesh_fmt, domain, fe1, fe2>::makeMesh(file_path);
    // // std::cout << msh << std::endl;
    
}