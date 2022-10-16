#include "2/core/traits/_include.hxx"
#include "2/core/mesh2.hxx"

int
main(int argc, char** argv)
{

    auto constexpr d = lolita::CartesianMesh(2);
    auto constexpr hdg = lolita::HybridDiscontinuousGalerkinDiscretization(lolita::MonomialBasis(1), lolita::MonomialBasis(1));
    auto constexpr displacement = lolita::UnknownField("Displacement", 2, 1, hdg);
    auto constexpr eps = lolita::SmallStrainOperator(displacement);
    auto constexpr quad = lolita::GaussQuadrature(2);
    auto constexpr pot = lolita::InternalPotential(2, quad, eps);
    auto constexpr lag = lolita::Lagrangian("1", pot);
    //
    auto constexpr young_modulus = lolita::Label("YoungModulus");
    auto constexpr poisson = lolita::Label("PoissonRatio");
    auto constexpr damage_h = lolita::Label("Damage");

    auto file_path = "/home/dsiedel/projetcs/lolita/data/meshes/mesh.msh";
    auto lib_displacement_path = "/home/dsiedel/projetcs/lolita/data/behavior/bhv_micromorphic_displacement/src/libBehaviour.so";
    auto lib_displacement_label = "MicromorphicDisplacement";
    auto hyp = mgis::behaviour::Hypothesis::PLANESTRAIN;

    //
    auto linear_system = lolita::LinearSystem<1>();
    //

    auto elements = lolita::core::MeshFileParser(file_path).template makeFiniteElementSet<d>();
    elements->setDomainLagrangian<displacement.getDimDomain(), lag>("ROD");
    elements->setDomainPotential<displacement.getDimDomain(), lag, pot>("ROD", lib_displacement_path, lib_displacement_label, hyp);
    elements->setLagrangian<displacement.getDimDomain(), lag>("ROD");
    elements->setPotential<displacement.getDimDomain(), lag, pot>("ROD");
    //
    elements->setElementDiscreteField<displacement.getDimDomain(), displacement>("ROD");
    elements->setElementDiscreteField<displacement.getDimDomain() - 1, displacement>("ROD");
    elements->addElementDiscreteFieldToLinearSystem<displacement.getDimDomain() - 1, displacement>("ROD", linear_system);
    //
    elements->setPotentialStrainOperators<displacement.getDimDomain(), lag, pot>("ROD");
    elements->setPotentialMaterialProperty<displacement.getDimDomain(), lag, pot, young_modulus>("ROD", [](lolita::Point const & p) { return 200.0; });
    elements->setPotentialMaterialProperty<displacement.getDimDomain(), lag, pot, poisson>("ROD", [](lolita::Point const & p) { return 200.0; });
    elements->setPotentialExternalVariable<displacement.getDimDomain(), lag, pot, damage_h>("ROD", [](lolita::Point const & p) { return 200.0; });
    elements->setPotentialStrains<displacement.getDimDomain(), lag, pot>("ROD");
    elements->integratePotentialConstitutiveEquation<displacement.getDimDomain(), lag, pot>("ROD");
    elements->setLagrangianJacobianMatrix<displacement.getDimDomain(), lag>("ROD");

    // std::cout << * elements << std::endl;

}