//
// Created by dsiedel on 13/04/2022.
//

#include <span>
#include "gtest/gtest.h"

#include "lolita/lolita_core_mesh.hxx"

TEST(test_lolita, test_lolita_1)
{

    using namespace lolita;
    using namespace lolita::core;

    auto file_path = "";
//    file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/perforated_strip_huge.msh";
    file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/meshes/unit_square_3_cpp.msh";
    auto msh_file = lolita::file::File(file_path);

    auto const constexpr dom = Domain(EuclideanFrame::Cartesian, 2);
    auto const constexpr qad = Quadrature(QuadratureRule::Gauss, 2);
    auto const constexpr fet = FiniteElement('D', 1, HHO(1, 1), qad, MappingOperator::Gradient, MappingOperator::Identity);
    auto const constexpr met = MixedElement<fet>(Model::Solid);
    auto const constexpr ets = Elements<fet>();
    auto const constexpr cet = CoupledElement<met, met>(QuadratureRule::Gauss);


    auto msh = mesh::MeshBase<MeshFormatType::Gmsh, dom, ets>(msh_file, {});
    auto & elem = msh.elements_.get<2>().get<0>().get("345").get();

    print("ns :");
    print(elem.get<0>().get().ns.get(0,0));

    elem.get<0>().get().operators.get(0);
    print(elem.get<0>().get().operators.get(0));

    auto myvec = Vector<Real>();

    myvec.resize(4);
    myvec.setZero();
    print(myvec);


    print(msh);

}

