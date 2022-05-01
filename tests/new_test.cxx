//
// Created by dsiedel on 01/05/22.
//

#include <span>
#include "gtest/gtest.h"

#include "lolita/lolita.hxx"
#include "lolita/lolita_core_element_geometry.hxx"
#include "lolita/lolita_core.hxx"
#include "lolita/lolita_core_finite_element.hxx"
#include "lolita/lolita_core_mesh.hxx"

TEST(test_lolita2, test_lolita_2)
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

    print(msh);

}