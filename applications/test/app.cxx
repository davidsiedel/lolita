// #include "2/core/traits/_include.hxx"
// #include "2/core/mesh2.hxx"
// #include "geometry/domain.hxx"

#include "mesh/mesh.hxx"
#include "physics/field.hxx"

#include "tensor.hxx"
// #include "quadrature/gauss.hxx"

int
main(int argc, char** argv)
{

    using Frame = lolita::geometry::CartesianFrame<2>;
    using Parser3 = lolita::mesh::MeshFactory<Frame, lolita::mesh::Msh>;
    using MyMesh = lolita::mesh::MyMesh<Frame>;

    auto mesh_file_path = "/home/dsiedel/projetcs/lolita/applications/test/mymesh.msh";

    auto my_mesh = MyMesh(mesh_file_path);
    my_mesh.working();

    // auto elem_label = "0000000017000000001800000000230000000024";
    // auto curve_label = "00000000150000000017";
    // auto const & quad_element = parser3.elements2_.getComponent<lolita::geometry::Quadrangle>()[elem_label];
    // auto const & line_element = parser3.elements2_.getComponent<lolita::geometry::Segment>()[curve_label];
    // //
    // auto telem0 = std::chrono::high_resolution_clock::now();
    // // auto segs0 = quad_element->getNeighbors<lolita::geometry::Segment>();
    // auto telem1 = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elem_dur = telem1 - telem0;
    // std::cout << "neighborhood : " << elem_dur.count() << "s\n";
    // //
    // auto tneigh0 = std::chrono::high_resolution_clock::now();
    // auto segs1 = quad_element->getInnerNeighbors<lolita::geometry::Segment>();
    // auto tneigh1 = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elem_nei = tneigh1 - tneigh0;
    // std::cout << "neighborhood : " << elem_nei.count() << "s\n";
    // std::cout << "sizeof quad : " << sizeof(* quad_element) << std::endl;
    // std::cout << "sizeof line : " << sizeof(* line_element) << std::endl;

    // static_assert(std::same_as<lolita::geometry::Node, typename lolita::geometry::ShapeCollection<Frame>::Shape<0, 0>>);
    // static_assert(lolita::geometry::ShapeCollection<Frame>::getNumComponents<0>() == 1);
    // auto ici = lolita::geometry::ShapeCollection<Frame>();
    // auto cc = ici.getComponent<lolita::geometry::Node>();

    // auto constexpr coord = lolita::geometry::ShapeInnerNeighborhood<lolita::geometry::Triangle>::ShapeTraits<lolita::geometry::Segment>::coordinates_;
    // std::cout << coord[0] << ", " << coord[1] << std::endl;

    // auto constexpr res = lolita::geometry::ShapeInnerNeighborhoodTraits<lolita::geometry::Triangle>::getCoordinate<lolita::geometry::Node>(0);
    // std::cout << res << std::endl;

    // using Node = lolita::geometry::Node;
    // using Quadrangle = lolita::geometry::Quadrangle;
    // using Triangle = lolita::geometry::Triangle;
    // using Tetrahedron = lolita::geometry::Tetrahedron;
    // using Frame = lolita::geometry::CartesianFrame<2>;
    // using NodeNeighbor = lolita::geometry::ShapeOuterNeighborhoodTraits<Node, Frame>::OuterNeighbor<2, 1>;
    // auto constexpr res = lolita::geometry::ShapeOuterNeighborhoodTraits<Node, Frame>::getNumOuterNeighbors<>();
    // std::cout << res << std::endl;
    // static_assert(std::same_as<NodeNeighbor, Quadrangle>);
    // auto constexpr pnt = lolita::geometry::Point<2>();

    // auto constexpr pair = std::pair<int, int>(1, 1);

    // auto constexpr f1 = pair.first;
    // auto constexpr f2 = pair.second;

    // std::cout << sizeof(std::tuple<>) << std::endl;
    // std::cout << std::tuple_size_v<std::tuple<>> << std::endl;

    // ---------------------------------------------------------------------------------------------------------------------------------------------------------

    auto tns = lolita::tensor::StaticTensor<lolita::Real, 3>();
    static_assert(lolita::tensor::TensorConcept<lolita::tensor::StaticTensor<lolita::Real, 3>, lolita::Real, 1>);
    static_assert(lolita::tensor::DynamicTensor<lolita::Real, 3>::NumIndices == 3);

    auto constexpr n_ = 4;

    lolita::tensor::DynamicTensor<lolita::Real, 2> td(n_, n_);
    lolita::tensor::DynamicTensor<lolita::Real, 1> vd(n_);
    for (int i = 0; i < n_; i++)
    {
        vd(i) = 1;
        for (int j = 0; j < n_; j++)
        {
            td(i, j) = 1;
        }
    }
    lolita::tensor::StaticTensor<lolita::Real, n_, n_> ts;
    lolita::tensor::StaticTensor<lolita::Real, n_> vs;
    for (int i = 0; i < n_; i++)
    {
        vs(i) = 1;
        for (int j = 0; j < n_; j++)
        {
            ts(i, j) = 1;
        }
    }

    auto ts0 = std::chrono::high_resolution_clock::now();
    auto ress = ts * vs;
    auto ts1 = std::chrono::high_resolution_clock::now();

    auto td0 = std::chrono::high_resolution_clock::now();
    auto resd = td * vd;
    auto td1 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> fd = td1 - td0;
    std::chrono::duration<double> fs = ts1 - ts0;

    std::cout << "dynamic : " << fd.count() << "s\n";
    std::cout << "static : " << fs.count() << "s\n";
    // lolita::tensor::DynamicTensor<lolita::Real, 1> resss = resd;
    // std::cout << resss << std::endl;

}