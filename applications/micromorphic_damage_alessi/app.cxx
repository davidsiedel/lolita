#include "core/700_mesh.hxx"

#include <mpi.h>

// class Polygon {
//     public:
//     virtual void show() {
//         std::cout<<"It is a polygon"<<std::endl;
//     }
// };

// class Hexagon : public Polygon {
//     public:
//     void show() {
//         std::cout<<"Hexagon is a 6 sided polygon"<<std::endl;
//     }
// };

// class Pentagon : public Polygon {
//     public:
//     void show() {
//         std::cout<<"Pentagon is a 5 sided polygon"<<std::endl;
//     }
// };

template<int a>
struct Slave;

template<int a>
struct Master
{

    explicit
    Master(
        int val
    )
    :
    val_(val)
    {}
    
    void
    makeSlave()
    {
        slave_ = std::make_unique<Slave<a>>(* this);
    }

    void
    print()
    {
        slave_->print();
    }

    int val_;

    std::unique_ptr<Slave<a>> slave_;
};

template<int a>
struct Slave
{

    explicit
    Slave(
        Master<a> const & master
    )
    :
    master_(master)
    {}

    void
    print()
    {
        std::cout << "val : " << master_.val_ << std::endl;
    }

    Master<a> const & master_;

};

static void
dummy(
    lolita::DenseMatrixConcept<lolita::Real, 2, 2> auto matrix
)
{
    std::cout << matrix << std::endl;
}

static void
dummy2(
    lolita::DenseMatrixConcept<lolita::Real> auto matrix
)
{
    std::cout << matrix << std::endl;
}

template<typename T>
struct Ref
{
    Ref(T const & t) : ref_(t) {}
    T const & get() const { return ref_;};
    T const & ref_;
};

template<lolita::Label a>
struct Ref2
{

    // Ref(T const & t) : ref_(t) {}
    // Ref2() : ref_(a) {}
    lolita::Label const & ref_ = a;
};

template<auto a>
struct RefA
{
    auto static constexpr a_ = a;
};

template<typename... T>
struct RefC
{
    
};

template<int... a>
struct RefB
{
    std::array<int, sizeof...(a)> static constexpr a_ = {a...};
};

template<typename T>
struct Youpi
{
    constexpr explicit
    Youpi(T a) : a_(a) {}
    T a_;
};

template<Youpi y>
struct YoupiT
{

};

template<typename T>
struct RefTest
{

};

template<int...>
struct Data;

struct DataBase
{

    DataBase()
    :
    offset_(-1)
    {}

    explicit
    DataBase(
        int offset
    )
    :
    offset_(offset)
    {}

    virtual
    ~DataBase()
    {}

    virtual
    lolita::Real *
    data()
    =0;

    virtual
    lolita::Real const *
    data()
    const
    =0;

    virtual
    lolita::Integer
    size()
    const
    =0;
    
    // virtual
    // lolita::algebra::View<lolita::DenseVector<lolita::Real> const>
    // get()
    // const
    // =0;

    lolita::algebra::View<lolita::DenseVector<lolita::Real> const>
    view(
        int size
    )
    const
    {
        return lolita::algebra::View<lolita::DenseVector<lolita::Real> const>(this->data(), size);
    }

    lolita::algebra::View<lolita::DenseVector<lolita::Real> const>
    view()
    const
    {
        return lolita::algebra::View<lolita::DenseVector<lolita::Real> const>(this->data(), this->size());
    }

    template<int size>
    lolita::algebra::View<lolita::DenseVector<lolita::Real, size> const>
    view()
    const
    {
        return lolita::algebra::View<lolita::DenseVector<lolita::Real, size> const>(this->data());
    }

    template<int size>
    auto
    get()
    const
    {
        return static_cast<Data<size> const *>(this)->get();
    }

    int offset_;

};

template<int a>
struct Data<a> : DataBase
{

    int static constexpr size_ = a;

    Data(
        auto &&... args
    )
    :
    DataBase(),
    data_(std::move(args)...)
    {}

    lolita::DenseVector<lolita::Real, a> const &
    get()
    const
    {
        return data_;
    }

private:

    lolita::Real *
    data()
    override
    {
        return data_.data();
    }

    lolita::Real const *
    data()
    const override
    {
        return data_.data();
    }

    lolita::Integer
    size()
    const override
    {
        return data_.size();
    }
    
    // lolita::algebra::View<lolita::DenseVector<lolita::Real> const>
    // get()
    // const override
    // {
    //     return lolita::algebra::View<lolita::DenseVector<lolita::Real> const>(data_.data(), data_.size());
    // };

    // template<int size>
    // auto
    // view()
    // const
    // {
    //     return lolita::algebra::View<lolita::DenseVector<lolita::Real, size> const>(data_.data());
    // }

    lolita::DenseVector<lolita::Real, a> data_;

};

template<>
struct Data<> : DataBase
{

    int static constexpr size_ = -1;

    Data(
        auto &&... args
    )
    :
    DataBase(),
    data_(std::move(args)...)
    {}

private:

    lolita::Real *
    data()
    override
    {
        return data_.data();
    }

    lolita::Real const *
    data()
    const override
    {
        return data_.data();
    }

    lolita::Integer
    size()
    const override
    {
        return data_.size();
    }
    
    // lolita::algebra::View<lolita::DenseVector<lolita::Real> const>
    // get()
    // const override
    // {
    //     return lolita::algebra::View<lolita::DenseVector<lolita::Real> const>(data_.data(), data_.size());
    // };

    lolita::DenseVector<lolita::Real> data_;

};

template<lolita::Label lab>
static constexpr
auto
makeIt()
{
    if constexpr (lab == "1")
    {
        return int(1);
    }
    else
    {
        return double(3);
    }
}

template<typename... T>
using AGG = lolita::utility::Aggregate<T...>;

int
main(int argc, char** argv)
{

    auto constexpr cell_field = lolita::UnknownField(2, 1, lolita::Basis("Monomial", 1));
    static_assert(lolita::FieldConcept<decltype(cell_field)>);

    static_assert(lolita::LinearOperatorConcept<decltype(lolita::Gradient(cell_field))>);

    std::cout << lolita::FieldTraits<lolita::Gradient(cell_field)>::template getSize<lolita::Domain("Cartesian", 2)>() << std::endl;

    auto master = Master<2>(5);
    master.makeSlave();
    master.print();

    std::unique_ptr<DataBase> ptr_test = std::make_unique<Data<3>>(lolita::DenseVector<lolita::Real, 3>::Zero());
    // std::cout << ptr_test->get().template segment<2>(0) << std::endl;
    std::cout << ptr_test->view<2>() << std::endl;
    std::cout << ptr_test->get<2>() << std::endl;

    auto constexpr ici_1 = makeIt<lolita::Label("1")>();
    auto constexpr ici_2 = makeIt<lolita::Label("2")>();
    static_assert(std::is_same_v<std::decay_t<decltype(ici_1)>, int>);
    // static_assert(std::is_same_v<decltype(ici_2), double>);

    YoupiT<Youpi(3)>();
    //YoupiT<3>();

    // auto vec_t = std::vector<int>{1, 2, 3, 4};
    // auto spn_t1 = std::span<int>(vec_t.begin(), vec_t.begin() + 2);
    // auto spn_t2 = std::span<int>(vec_t.begin() + 2, vec_t.end());
    // std::cout << "vec_t" << std::endl;
    // for (auto const & val : vec_t)
    // {
    //     std::cout << val << " ";
    // }
    // std::cout << "spn_t1" << std::endl;
    // for (auto const & val : spn_t1)
    // {
    //     std::cout << val << " ";
    // }
    // std::cout << std::endl;
    // std::cout << "spn_t2" << std::endl;
    // for (auto const & val : spn_t2)
    // {
    //     std::cout << val << " ";
    // }
    // std::cout << std::endl;
    // vec_t.insert(vec_t.begin() + 1, 0);
    // std::cout << "vec_t" << std::endl;
    // for (auto const & val : vec_t)
    // {
    //     std::cout << val << " ";
    // }
    // std::cout << "spn_t1" << std::endl;
    // for (auto const & val : spn_t1)
    // {
    //     std::cout << val << " ";
    // }
    // std::cout << std::endl;
    // std::cout << "spn_t2" << std::endl;
    // for (auto const & val : spn_t2)
    // {
    //     std::cout << val << " ";
    // }
    // std::cout << std::endl;


    using LHS_0 = lolita::utility::tuple_expansion_t<lolita::utility::Aggregate, RefA, AGG<int, int, int>(1, 2, 3), AGG<char, int>('A', 1)>;
    // using LHS_ = lolita::utility::array_expansion_t<RefA, std::array<int, 3>{1, 2, 3}>;
    using RHS_0 = lolita::utility::Aggregate<RefA<1>, RefA<2>, RefA<3>, RefA<'A'>, RefA<1>>;
    static_assert(std::is_same_v<LHS_0, RHS_0>);
    using LHS_ = lolita::utility::variadic_type_array_expansion_t<RefB, std::array<int, 3>{1, 2, 3}, std::array<int, 1>{1}>;
    using RHS_ = RefB<1, 2, 3, 1>;
    static_assert(std::is_same_v<LHS_, RHS_>);
    using LHS_1 = lolita::utility::variadic_type_array_expansion_t<RefB, std::array<int, 1>{1}>;
    using RHS_1 = RefB<1>;
    static_assert(std::is_same_v<LHS_1, RHS_1>);
    using LHS_2 = lolita::utility::tuple_cat_t<RefB<1, 2, 3>, RefB<1, 2, 3>, RefB<1, 2, 3>, RefB<1, 2, 3>>;
    using RHS_2 = RefB<1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3>;
    static_assert(std::is_same_v<LHS_2, RHS_2>);
    using LHS_3 = lolita::utility::tuple_cat_t<RefC<int, float, int>, RefC<int, int, int>, RefC<int, int, int>, RefC<float, double, char>>;
    using RHS_3 = RefC<int, float, int, int, int, int, int, int, int, float, double, char>;
    static_assert(std::is_same_v<LHS_3, RHS_3>);


    std::cout << sizeof(std::vector<int>) << std::endl;
    std::cout << sizeof(std::vector<double>) << std::endl;
    std::cout << sizeof(lolita::DenseVector<double>) << std::endl;
    std::cout << sizeof(std::basic_string_view<char>) << std::endl;
    std::cout << sizeof(int) << std::endl;
    std::cout << sizeof(std::vector<lolita::DenseMatrix<lolita::Real>>) << std::endl;
    auto lab = lolita::Label();
    auto const & lab_ref = lab;
    std::cout << sizeof(Ref<lolita::Label>) << std::endl;
    // std::cout << sizeof(Ref2<"Stabilization">) << std::endl;
    std::cout << sizeof(double const &) << std::endl;

    // auto constexpr ref2 = Ref2<"Stabilization">();

    static_assert(lolita::DenseMatrixConcept<lolita::DenseMatrix<lolita::Real, 2, 2>, lolita::Real, 2, 2>);
    // std::cout << std::fixed << std::setprecision(3);
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------
    // declaring behavior
    auto lib_displacement_path = "/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_alessi/data/behavior/bhv_micromorphic_displacement/src/libBehaviour.so";
    auto lib_displacement_label = "MicromorphicDisplacement";
    // declaring behavior
    auto lib_damage_path = "/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_alessi/data/behavior/bhv_micromorphic_damage/src/libBehaviour.so";
    auto lib_damage_label = "MicromorphicDamage";
    //
    auto opts = mgis::behaviour::FiniteStrainBehaviourOptions{
        mgis::behaviour::FiniteStrainBehaviourOptions::PK1,
        mgis::behaviour::FiniteStrainBehaviourOptions::DPK1_DF
    };
    // auto hyp = mgis::behaviour::Hypothesis::PLANESTRAIN;
    auto hyp = mgis::behaviour::Hypothesis::PLANESTRAIN;
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------
    // constants
    // auto constexpr domain = lolita::domain("Cartesian", 2);
    // static_assert(domain.isCartesian());
    // auto constexpr quadrature = lolita::quadrature("Gauss", 2);
    // auto constexpr cells = lolita::ElementType::cells(domain);
    // auto constexpr faces = lolita::ElementType::faces(domain);
    // // discretization
    // auto constexpr hdg = lolita::HybridDiscontinuousGalerkin::hybridDiscontinuousGalerkin(1, 1);
    // // generalized strains
    // auto constexpr displacement_generalized_strain = lolita::GeneralizedStrain(0, lolita::field("DenseVector", "U"), lolita::mapping("SmallStrain", lolita::field("DenseVector", "U")));
    // auto constexpr damage_generalized_strain = lolita::GeneralizedStrain(1, lolita::field("Scalar", "D"), lolita::mapping("Gradient", lolita::field("Scalar", "D")), lolita::mapping("Identity", lolita::field("Scalar", "D")));
    // auto constexpr lag_generalized_strain = lolita::GeneralizedStrain(2, lolita::field("Scalar", "F"), lolita::mapping("Identity", lolita::field("Scalar", "F")));
    // // behaviors
    // auto constexpr displacement_behavior = lolita::Behavior(0, displacement_generalized_strain);
    // auto constexpr damage_behavior = lolita::Behavior(1, damage_generalized_strain);
    // auto constexpr lag_behavior = lolita::Behavior(2, lag_generalized_strain, displacement_generalized_strain);
    // // finite elements
    // auto constexpr displacement_element =  lolita::FiniteElementMethod(displacement_generalized_strain, displacement_behavior, hdg, quadrature);
    // auto constexpr damage_element =  lolita::FiniteElementMethod(damage_generalized_strain, damage_behavior, hdg, quadrature);
    // mesh
    auto file_path = "/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_alessi/mesh.msh";
    auto out_displacement_file = "/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_alessi/.out_u.msh";
    auto out_damage_file = "/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_alessi/.out_d.msh";
    auto force_out_stream = std::basic_ofstream<lolita::Character>();
    auto damage_dissipated_energy_out_stream = std::basic_ofstream<lolita::Character>();
    auto damage_stored_energy_out_stream = std::basic_ofstream<lolita::Character>();
    auto displacement_dissipated_energy_out_stream = std::basic_ofstream<lolita::Character>();
    auto displacement_stored_energy_out_stream = std::basic_ofstream<lolita::Character>();
    force_out_stream.open("/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_alessi/.out_force.txt");
    damage_dissipated_energy_out_stream.open("/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_alessi/.out_damage_dissipated_energy.txt");
    damage_stored_energy_out_stream.open("/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_alessi/.out_damage_stored_energy.txt");
    displacement_dissipated_energy_out_stream.open("/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_alessi/.out_displacement_dissipated_energy.txt");
    displacement_stored_energy_out_stream.open("/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_alessi/.out_displacement_stored_energy.txt");
    // mesh build
    // -> TEST -------------------------------------------------------------------------------------------------------------------------------------------------
    /**
     * Domain
     * 
     */
    auto constexpr domain = lolita::Domain("Cartesian", 2);
    auto constexpr face_dim = lolita::Integer(1);
    auto constexpr cell_dim = lolita::Integer(2);
    /**
     * Labels
     * 
     */
    auto constexpr _yg = lolita::Label("YoungModulus");
    auto constexpr _nu = lolita::Label("PoissonRatio");
    auto constexpr _dm = lolita::Label("Damage");
    auto constexpr _stab = lolita::Label("Stabilization");
    /**
     * Quadrature
     * 
     */
    auto constexpr quadrature = lolita::Quadrature("Gauss", 2);
    /**
     * Discretizations
     * 
     */
    auto constexpr hdg_discretization = lolita::HybridDiscontinuousGalerkin(1, 1);
    /**
     * Fields
     * 
     */
    auto constexpr displacement = lolita::DiscreteField("Displacement", 2, 1, lolita::HybridDiscontinuousGalerkin(1, 1));
    auto constexpr displacement2 = lolita::DiscreteField("Displacement", 2, 1, lolita::HybridDiscontinuousGalerkin(1, 2));
    auto constexpr damage = lolita::DiscreteField("Damage", 2, 0, lolita::HybridDiscontinuousGalerkin(1, 1));
    auto constexpr force = lolita::DiscreteField("Force", 1, 0, lolita::Monomial(1));
    auto ref_f1 = Ref(displacement.getLabel());
    auto ref_f2 = Ref(displacement2.getLabel());
    if (ref_f1.get() != ref_f2.get()) throw std::runtime_error("HIII");
    /**
     * Mappings
     * 
     */
    auto constexpr hdg_displacement_gradient = lolita::Mapping("Gradient", displacement); // G(u)
    auto constexpr hdg_displacement_stabilization = lolita::Mapping("Stabilization", displacement); // Z(u)
    auto constexpr hdg_displacement_x = lolita::Mapping("Identity", 0, 0, displacement); // I_x(u)
    auto constexpr hdg_force_view = lolita::Mapping("Identity", force); // I(f)
    /**
     * Potentials
     * 
     */
    auto constexpr elastic_potential = lolita::Potential("Elasticity", quadrature, hdg_displacement_gradient); // psi(I(u))
    auto constexpr stabilization_potential = lolita::Potential("Stabilization", quadrature, hdg_displacement_stabilization); // psi(Z(u))
    auto constexpr force_potential = lolita::Potential("Lagrangian", quadrature, hdg_displacement_x, hdg_force_view); // psi(I_x(u), I(f))
    auto constexpr lagrangian = lolita::Lagrangian("All", elastic_potential, stabilization_potential); // psi(I_x(u), I(f))
    // auto constexpr pot_size = lolita::PotentialTraits<elastic_potential>::template getSize<lolita::Element::triangle(1), domain>();
    // auto constexpr mat0 = lolita::Mapping("Gradient", displacement);
    // auto constexpr mat1 = lolita::Mapping("Gradient", damage);
    // auto constexpr mat2 = lolita::Mapping("Identity", damage);
    auto constexpr pot0 = lolita::Potential("Elasticity", quadrature, lolita::Mapping("Gradient", displacement));
    auto constexpr pot1 = lolita::Potential("Fracture", quadrature, lolita::Mapping("Gradient", damage), lolita::Mapping("Identity", damage));
    auto constexpr pot2 = lolita::Potential("Stab", quadrature, lolita::Mapping("Identity", displacement));
    auto constexpr pot3 = lolita::Potential("Penalization", quadrature, lolita::Mapping("Identity", force));
    auto constexpr lag0 = lolita::Lagrangian("0", pot0, pot1, pot2, pot3);
    auto constexpr lag1 = lolita::Lagrangian("1", pot0, pot1, pot2);
    auto constexpr lag2 = lolita::Lagrangian("2", pot0, pot1);
    auto constexpr lag3 = lolita::Lagrangian("3", pot0);
    using LHS_PP = std::tuple<std::tuple<std::tuple<int>>, std::tuple<char, double>, std::tuple<float>>;
    using RHS_PP = std::tuple<int, char, float>;
    // static_assert(std::is_same_v<lolita::utility::tuple_cat_t<LHS_PP>, RHS_PP>);
    lolita::utility::TD<lolita::utility::tuple_merge_t<LHS_PP>>();
    // lolita::utility::TD<typename lolita::LagTraits<lag0>::TEST2>();
    // auto constexpr pot_size = lolita::PotentialTraits<pot0, pot1, pot2>::template getSize<lolita::Element::triangle(1), domain>();
    std::cout << "num fields : " << lolita::LagTraits<lag0>::getNumFields() << std::endl;
    
    std::cout << "field index : " << lolita::LagTraits<lag0>::template getIndex<displacement>() << std::endl;
    std::cout << "field index : " << lolita::LagTraits<lag0>::template getIndex<damage>() << std::endl;
    std::cout << "field index : " << lolita::LagTraits<lag0>::template getIndex<force>() << std::endl;
    std::cout << "offset : " << lolita::LagTraits<lag0>::template getOffset<domain, lolita::Mapping("Gradient", displacement)>() << std::endl;
    std::cout << "offset : " << lolita::LagTraits<lag0>::template getOffset<domain, lolita::Mapping("Identity", damage)>() << std::endl;
    std::cout << "offset : " << lolita::LagTraits<lag0>::template getOffset<domain, lolita::Mapping("Identity", force)>() << std::endl;
    //
    std::cout << "pot size : " << lolita::PotentialTraits<pot0, pot1, pot2>::template getSize<lolita::Element::triangle(1), domain>() << std::endl;
    std::cout << "pot size : " << lolita::PotentialTraits<pot0, pot1>::template getSize<lolita::Element::triangle(1), domain>() << std::endl;
    std::cout << "pot size : " << lolita::PotentialTraits<pot0>::template getSize<lolita::Element::triangle(1), domain>() << std::endl;
    //
    std::cout << "pot size : " << lolita::PotentialTraits<pot0, pot1, pot2, pot3>::template getSize<domain>() << std::endl;
    std::cout << "pot size : " << lolita::PotentialTraits<pot0, pot1>::template getSize<domain>() << std::endl;
    std::cout << "pot size : " << lolita::PotentialTraits<pot0>::template getSize<domain>() << std::endl;
    //
    std::cout << "field index : " << lolita::PotentialTraits<pot0, pot1, pot2, pot3>::template getIndex<displacement>() << std::endl;
    std::cout << "field index : " << lolita::PotentialTraits<pot0, pot1, pot2, pot3>::template getIndex<damage>() << std::endl;
    std::cout << "field index : " << lolita::PotentialTraits<pot0, pot1, pot2, pot3>::template getIndex<force>() << std::endl;
    //
    std::cout << "pot size : " << lolita::PotentialTraits<pot1>::template getSize<lolita::Element::triangle(1), domain>() << std::endl;
    //
    std::cout << "pot size : " << lolita::PotentialTraits<pot1>::template getSize<domain>() << std::endl;
    //

    //
    std::cout << "pot size element : " << lolita::PotentialTraits<pot0, pot1, pot2>::template getSize<lolita::Element::triangle(1), domain>() << std::endl;
    std::cout << "pot size point : " << lolita::PotentialTraits<pot0, pot1, pot2>::template getSize<domain>() << std::endl;
    std::cout << "field offset : " << lolita::PotentialTraits<pot0, pot1, pot2>::template getOffset<lolita::Element::triangle(1), domain, displacement>() << std::endl;
    std::cout << "field offset : " << lolita::PotentialTraits<pot0, pot1, pot2>::template getOffset<lolita::Element::triangle(1), domain, damage>() << std::endl;
    /**
     * Potentials
     * 
     */
    using HHHIIIA = lolita::AbstractElementLagrangian<lolita::Element::triangle(1), domain>;
    using HHHIIIC = lolita::ElementLagrangian<lolita::Element::triangle(1), domain, lag3>;
    auto constexpr size_test = lolita::LagTraits<lag3>::template getSize<lolita::Element::triangle(1), domain>();
    auto constexpr index_test = lolita::LagTraits<lag3>::template getField<0>();
    std::cout << "lolita::LagTraits<lag3>::template getSize<lolita::Element::triangle(1), domain>() : " << size_test << std::endl;
    std::cout << "lolita::LagTraits<lag3>::getNumFields() : " << lolita::LagTraits<lag3>::getNumFields() << std::endl;
    std::cout << "lolita::LagTraits<lag0>::getNumFields() : " << lolita::LagTraits<lag0>::getNumFields() << std::endl;
    // auto non_abstract_test = HHHIIIC();
    // std::unique_ptr<HHHIIIA> abstract_test = std::make_unique<HHHIIIC>();
    // abstract_test->template setPotential<lag3, pot0>();
    // abstract_test->template getPotential<lag3, pot0>().hello();

    auto elements = lolita::MeshFileParser(file_path).template makeFiniteElementSet<domain>();
    auto linear_system = lolita::LinearSystem<lolita::Strategy::eigenLU()>::make_unique();
    // std::cout << * elements << std::endl;
    /**
     * Domain stuff
     * 
     */
    /**
     * @brief 
     * 
     */
    /**
     * Dofs
     * 
     */
    elements->addElementDiscreteField<face_dim, displacement>("ROD");
    elements->addElementDiscreteField<cell_dim, displacement>("ROD");
    elements->addElementDiscreteFieldDegreeOfFreedom<face_dim, displacement>("ROD");
    elements->addElementDiscreteFieldDegreeOfFreedom<cell_dim, displacement>("ROD");
    elements->addElementDiscreteFieldDegreeOfFreedomToSystem<face_dim, displacement>("ROD", * linear_system);
    elements->addDomainDiscreteField<cell_dim, displacement>("ROD");
    elements->addDomainDiscreteFieldLoad<cell_dim, displacement>("ROD", 0, 0, [](lolita::Point const & p, lolita::Real const & t) { return t; });
    /**
     * Ops
     * 
     */
    elements->setDomainLagrangian<cell_dim, lag3>("ROD");
    elements->setDomainPotential<cell_dim, lag3, elastic_potential>("ROD", lib_displacement_path, lib_displacement_label, hyp);
    elements->setLagrangian<cell_dim, lag3>("ROD");

    auto jklmpoiu = lolita::HHHIIIMMM3::getSize();

    struct TestRefDom
    {

        explicit
        TestRefDom(
            lolita::DomainPotential<2, domain, lag3, elastic_potential> const & r
        )
        :
        r_(r)
        {}

        lolita::DomainPotential<2, domain, lag3, elastic_potential> const & r_;

    };
    elements->getFiniteDomain<cell_dim>("ROD").getLagrangian<lag3>().getPotential<lag3, elastic_potential>().getMgisBehavior();
    auto ggg = TestRefDom(elements->getFiniteDomain<cell_dim>("ROD").getLagrangian<lag3>().getPotential<lag3, elastic_potential>());
    elements->setPotential<cell_dim, lag3, elastic_potential>("ROD");
    elements->setPotentialStrainOperators<cell_dim, lag3, elastic_potential>("ROD");
    elements->setPotentialStrains<cell_dim, lag3, elastic_potential>("ROD");
    std::cout << linear_system->getSize() << std::endl;
    /**
     * Potentials
     * 
     */
    elements->addFormulation<cell_dim, elastic_potential, quadrature>("ROD", lib_displacement_path, lib_displacement_label, hyp);
    elements->addFormulation<face_dim, elastic_potential>("TOP");
    // elements->addFormulationStrainOperator<cell_dim, elastic_potential, elastic_potential.getStrain<0>()>("ROD");
    //
    elements->setFormulationMaterialProperty<cell_dim, elastic_potential, _yg>("ROD", [](lolita::Point const & p) { return 200.0; });
    elements->setFormulationMaterialProperty<cell_dim, elastic_potential, _nu>("ROD", [](lolita::Point const & p) { return 0.2; });
    elements->setFormulationExternalVariable<cell_dim, elastic_potential, _dm>("ROD", [](lolita::Point const & p) { return 0.0; });
    elements->setStrainOperators<cell_dim, elastic_potential>("ROD");
    elements->setStrains<cell_dim, elastic_potential>("ROD");
    elements->integrateFormulationConstitutiveEquation<cell_dim, elastic_potential>("ROD");
    auto internal_energy = elements->getInternalEnergy<cell_dim, elastic_potential>("ROD");
    std::cout << "internal_energy : " << internal_energy << std::endl;
    auto external_energy = elements->getExternalEnergy<cell_dim, displacement>("ROD", 0.1);
    std::cout << "external_energy : " << external_energy << std::endl;
    elements->CALL<face_dim, displacement>("ROD");
    //
    elements->setFormulationStrain<cell_dim, elastic_potential, elastic_potential.getStrain<0>()>("ROD"); // set G(u) in T
    // elements->setFormulationMaterialProperty<cell_dim, elastic_potential, _yg>("ROD", [](lolita::Point const & p) { return 200.0; });
    // elements->setFormulationExternalVariable<cell_dim, elastic_potential, _dm>("ROD", [](lolita::Point const & p) { return 200.0; });
    elements->integrateFormulationConstitutiveEquation<cell_dim, elastic_potential>("ROD");
    /**
     * Potentials
     * 
     */
    elements->recoverElementDiscreteFieldDegreeOfFreedom<cell_dim, displacement>("ROD");
    elements->reserveElementDiscreteFieldDegreeOfFreedom<cell_dim, displacement>("ROD");
    // <- TEST -------------------------------------------------------------------------------------------------------------------------------------------------

    auto defect = [] (lolita::Point const & point)
    {
        auto offset = 0.0025;
        auto x_c = 0.5;
        auto y_c = 1.0;
        if (point(0) > x_c - offset && point(0) < x_c + offset && point(1) > y_c - offset && point(1) < y_c + offset)
        {
            return true;
        }
        else
        {
            return false;
        }
    };

    auto sane = [] (lolita::Point const & point)
    {
        auto offset = 0.0025;
        auto x_c = 0.5;
        auto y_c = 1.0;
        if (point(0) > x_c - offset && point(0) < x_c + offset && point(1) > y_c - offset && point(1) < y_c + offset)
        {
            return false;
        }
        else
        {
            return true;
        }
    };

    elements->addDomain<cell_dim>("ROD", "DEFECT", defect);
    elements->addDomain<cell_dim>("ROD", "SANE", sane);

//     // dofs
//     auto face_displacement = elements->setDegreeOfFreedom<faces, lolita::field("DenseVector", "A"), hdg.getFaceBasis()>("ROD", "Displacement");
//     auto cell_displacement = elements->setDegreeOfFreedom<cells, lolita::field("DenseVector", "A"), hdg.getCellBasis()>("ROD", "Displacement");
//     auto face_damage = elements->setDegreeOfFreedom<faces, lolita::field("Scalar", "A"), hdg.getFaceBasis()>("ROD", "Damage");
//     auto cell_damage = elements->setDegreeOfFreedom<cells, lolita::field("Scalar", "A"), hdg.getCellBasis()>("ROD", "Damage");
//     //
//     auto top_force = elements->setDegreeOfFreedom<faces, lolita::field("Scalar", "A"), hdg.getFaceBasis()>("TOP", "TopForce");
//     auto left_force = elements->setDegreeOfFreedom<faces, lolita::field("Scalar", "A"), hdg.getFaceBasis()>("LEFT", "LeftForce");
//     auto bottom_force = elements->setDegreeOfFreedom<faces, lolita::field("Scalar", "A"), hdg.getFaceBasis()>("BOTTOM", "BottomForce");
//     // systems
//     auto displacement_system = lolita::System::make();
//     auto damage_system = lolita::System::make();
//     displacement_system->setUnknown("Displacement", face_displacement->size());
//     displacement_system->setBinding("TopForce", top_force->size());
//     displacement_system->setBinding("LeftForce", left_force->size());
//     displacement_system->setBinding("BottomForce", bottom_force->size());
//     displacement_system->initialize();
//     damage_system->setUnknown("Damage", face_damage->size());
//     damage_system->initialize();
//     // load
//     auto load0 = elements->setConstraint<faces>("TOP", "Pull", [](lolita::Point const & p, lolita::Real const & t) { return t; }, 1, 0);
//     auto load1 = elements->setConstraint<faces>("LEFT", "FixedL", [](lolita::Point const & p, lolita::Real const & t) { return 0.0; }, 0, 0);
//     // -> RIGHT
//     auto loadr = elements->setConstraint<faces>("RIGHT", "FixedR", [](lolita::Point const & p, lolita::Real const & t) { return 0.0; }, 0, 0);
//     // <- RIGHT
//     auto load2 = elements->setConstraint<faces>("BOTTOM", "FixedB", [](lolita::Point const & p, lolita::Real const & t) { return 0.0; }, 1, 0);
//     // adding behavior
//     // std::cout << lib_displacement_label << std::endl;
//     // std::cout << lib_displacement_path << std::endl;
//     // std::cout << lib_damage_label << std::endl;
//     // std::cout << lib_damage_path << std::endl;
//     auto micromorphic_displacement = elements->setBehavior<cells, quadrature>("ROD", lib_displacement_path, lib_displacement_label, hyp);
//     auto micromorphic_damage = elements->setBehavior<cells, quadrature>("ROD", lib_damage_path, lib_damage_label, hyp);
//     //making operators
//     elements->setStrainOperators<cells, displacement_element, hdg>("ROD", "MicromorphicDisplacement", "Displacement");
//     elements->setStrainOperators<cells, damage_element, hdg>("ROD", "MicromorphicDamage", "Damage");
//     elements->setElementOperators<cells, displacement_element, hdg>("ROD", "DisplacementStabilization");
//     elements->setElementOperators<cells, damage_element, hdg>("ROD", "DamageStabilization");
//     // setting variable
//     elements->setMaterialProperty<cells>("SANE", "MicromorphicDisplacement", "YoungModulus", [](lolita::Point const & p) { return 200.0; });
//     elements->setMaterialProperty<cells>("DEFECT", "MicromorphicDisplacement", "YoungModulus", [](lolita::Point const & p) { return 199.9; });
//     elements->setMaterialProperty<cells>("ROD", "MicromorphicDisplacement", "PoissonRatio", [](lolita::Point const & p) { return 0.2; });
//     elements->setExternalVariable<cells>("ROD", "MicromorphicDisplacement", "Temperature", [](lolita::Point const & p) { return 293.15; });
//     elements->setExternalVariable<cells>("ROD", "MicromorphicDisplacement", "Damage", [](lolita::Point const & p) { return 0.0; });
//     // setting variable
//     elements->setMaterialProperty<cells>("ROD", "MicromorphicDamage", "FractureEnergy", [](lolita::Point const & p) { return 1.0; });
//     elements->setMaterialProperty<cells>("ROD", "MicromorphicDamage", "CharacteristicLength", [](lolita::Point const & p) { return 0.05; });
//     elements->setMaterialProperty<cells>("ROD", "MicromorphicDamage", "PenalisationFactor", [](lolita::Point const & p) { return 300.; });
//     elements->setExternalVariable<cells>("ROD", "MicromorphicDamage", "Temperature", [](lolita::Point const & p) { return 293.15; });
//     elements->setExternalVariable<cells>("ROD", "MicromorphicDamage", "EnergyReleaseRate", [](lolita::Point const & p) { return 0.0; });
//     // setting parameter
//     elements->setParameter<faces>("TOP", "TopForceLagrange", [](lolita::Point const & p) { return 200.0; });
//     elements->setParameter<faces>("LEFT", "LeftForceLagrange", [](lolita::Point const & p) { return 200.0; });
//     elements->setParameter<faces>("BOTTOM", "BottomForceLagrange", [](lolita::Point const & p) { return 200.0; });
//     // stab
//     elements->setParameter<cells>("ROD", "DisplacementStabilization", [](lolita::Point const & p) { return 200.0 / (1.0 + 0.2); });
//     elements->setParameter<cells>("ROD", "DamageStabilization", [](lolita::Point const & p) { return 1.0 / 0.05; });
//     //
//     lolita::GmshFileParser::setOutput<domain>(out_displacement_file, elements, "MicromorphicDisplacement");
//     lolita::GmshFileParser::setOutput<domain>(out_damage_file, elements, "MicromorphicDamage");
//     // -> TEST
//     // std::cout << "elements->getElements<2, 1>()[0]->quadrature_.at(MicromorphicDamage).ips_[0].ops_.at(Damage)" << std::endl;
//     // std::cout << elements->getElements<2, 1>()[0]->quadrature_.at("MicromorphicDamage").ips_[0].ops_.at("Damage") << std::endl;
//     // std::cout << "elements->getElements<2, 1>()[0]->quadrature_.at(MicromorphicDisplacement).ips_[0].ops_.at(Displacement)" << std::endl;
//     // std::cout << elements->getElements<2, 1>()[0]->quadrature_.at("MicromorphicDisplacement").ips_[0].ops_.at("Displacement") << std::endl;

//     // <- TEST
    
//     auto num_steps = 50;
//     std::cout << "times :" << std::endl;
//     auto times = std::vector<lolita::Real>();
//     for (auto i = 0; i < num_steps + 1; i++)
//     {
//         auto val = i * 3.0e-1 / num_steps;
//         std::cout << i << " : " << val << " / ";
//         times.push_back(val);
//     }
//     std::cout << std::endl;
    
//     auto step = 0;
//     auto time = times[step];
//     // -> DISPLACEMENT
//     auto displacement_newton_step = [&] ()
//     {
//         auto iteration = 0;
//         while(iteration < 10)
//         {
//             // displacement_system->initialize();
//             displacement_system->initializeRhs();
//             displacement_system->initializeLhs();
//             displacement_system->initializeNormalization();
//             elements->setStrainValues<cells, displacement_element, hdg>("ROD", "MicromorphicDisplacement", "Displacement");
//             auto res_eval = elements->integrate<cells>("ROD", "MicromorphicDisplacement");
//             if (res_eval.isFailure())
//             {
//                 std::cout << "!!! displacement integration failure" << std::endl;
//                 return false;
//             }
//             else
//             {
//                 elements->assembleUnknownVector<cells, displacement_element, hdg>("ROD", "MicromorphicDisplacement", "Displacement", displacement_system);
//                 elements->assembleBindingVector<faces, displacement_element, hdg>("TOP", "TopForce", "Displacement", "Pull", displacement_system, time);
//                 elements->assembleBindingVector<faces, displacement_element, hdg>("LEFT", "LeftForce", "Displacement", "FixedL", displacement_system, time);
//                 elements->assembleBindingVector<faces, displacement_element, hdg>("BOTTOM", "BottomForce", "Displacement", "FixedB", displacement_system, time);
//                 auto res = displacement_system->getResidualEvaluation();
//                 if (res < 1.e-6)
//                 {
//                     std::cout << "displacement iteration " << iteration << " | res " << std::setprecision(10) << std::scientific << res << " | convergence" << std::endl;
//                     return true;
//                 }
//                 else
//                 {
//                     std::cout << "displacement iteration " << iteration << " | res " << std::setprecision(10) << std::scientific << res << " | solve" << std::endl;
//                     elements->assembleUnknownBlock<cells, displacement_element, hdg>("ROD", "MicromorphicDisplacement", "Displacement", displacement_system);
//                     //
//                     elements->assembleBindingBlock<faces, displacement_element, hdg>("TOP", "TopForce", "Displacement", "Pull", displacement_system);
//                     elements->assembleBindingBlock<faces, displacement_element, hdg>("LEFT", "LeftForce", "Displacement", "FixedL", displacement_system);
//                     elements->assembleBindingBlock<faces, displacement_element, hdg>("BOTTOM", "BottomForce", "Displacement", "FixedB", displacement_system);
//                     displacement_system->setCorrection();
//                     elements->updateUnknown<cells, displacement_element, hdg>("ROD", "Displacement", displacement_system);
//                     elements->updateUnknown<faces, displacement_element, hdg>("ROD", "Displacement", displacement_system);
//                     elements->updateBinding<faces, lolita::field("Scalar", "A"), hdg.getFaceBasis()>("TOP", "TopForce", displacement_system);
//                     elements->updateBinding<faces, lolita::field("Scalar", "A"), hdg.getFaceBasis()>("LEFT", "LeftForce", displacement_system);
//                     elements->updateBinding<faces, lolita::field("Scalar", "A"), hdg.getFaceBasis()>("BOTTOM", "BottomForce", displacement_system);
//                 }
//             }
//             iteration ++;
//             // total_iteration ++;
//         }
//         std::cout << "!!! displacement max iters" << std::endl;
//         return false;
//     };
//     // <- DISPLACEMENT
//     // -> DAMAGE
//     auto damage_newton_step = [&] ()
//     {
//         auto iteration = 0;
//         while(iteration < 10)
//         {
//             // displacement_system->initialize();
//             damage_system->initializeRhs();
//             damage_system->initializeLhs();
//             damage_system->initializeNormalization(1);
//             elements->setStrainValues<cells, damage_element, hdg>("ROD", "MicromorphicDamage", "Damage");
//             // std::cout << "damage iteration : " << iteration << std::endl;
//             auto res_eval = elements->integrate<cells>("ROD", "MicromorphicDamage");
//             if (res_eval.isFailure())
//             {
//                 std::cout << "!!! damage integration failure" << std::endl;
//                 return false;
//             }
//             else
//             {
//                 elements->assembleUnknownVector<cells, damage_element, hdg>("ROD", "MicromorphicDamage", "Damage", damage_system);
//                 auto res = damage_system->getResidualEvaluation();
//                 if (res < 1.e-6)
//                 {
//                     std::cout << "damage iteration " << iteration << " | res " << std::setprecision(10) << std::scientific << res << " | convergence" << std::endl;
//                     return true;
//                 }
//                 else
//                 {
//                     std::cout << "damage iteration " << iteration << " | res " << std::setprecision(10) << std::scientific << res << " | solve" << std::endl;
//                     elements->assembleUnknownBlock<cells, damage_element, hdg>("ROD", "MicromorphicDamage", "Damage", damage_system);
//                     damage_system->setCorrection();
//                     elements->updateUnknown<cells, damage_element, hdg>("ROD", "Damage", damage_system);
//                     elements->updateUnknown<faces, damage_element, hdg>("ROD", "Damage", damage_system);
//                 }
//             }
//             iteration ++;
//             // total_iteration ++;
//         }
//         std::cout << "!!! damage max iters" << std::endl;
//         return false;
//     };
//     // <- DAMAGE
//     auto displacement_convergence = [&] ()
//     {
//         displacement_system->initializeRhs();
//         displacement_system->initializeNormalization();
//         auto res_eval = elements->integrate<cells>("ROD", "MicromorphicDisplacement");
//         if (res_eval.isFailure())
//         {
//             std::cout << "!!! displacement integration failure" << std::endl;
//             return -1;
//         }
//         else
//         {
//             elements->assembleUnknownVector<cells, displacement_element, hdg>("ROD", "MicromorphicDisplacement", "Displacement", displacement_system);
//             elements->assembleBindingVector<faces, displacement_element, hdg>("TOP", "TopForce", "Displacement", "Pull", displacement_system, time);
//             elements->assembleBindingVector<faces, displacement_element, hdg>("LEFT", "LeftForce", "Displacement", "FixedL", displacement_system, time);
//             elements->assembleBindingVector<faces, displacement_element, hdg>("BOTTOM", "BottomForce", "Displacement", "FixedB", displacement_system, time);
//             // std::cout << "displacement res eval : " << displacement_system->getResidualEvaluation() << std::endl;
//             auto res = displacement_system->getResidualEvaluation();
//             // std::cout << "res : " << lolita::DenseMatrix<lolita::Real, 1, -1>(displacement_system->rhs_values_) << std::endl;
//             if (res < 1.e-6)
//             {
//                 // std::cout << "step convergence" << std::endl;
//                 return 1;
//             }
//             else
//             {
//                 return 0;
//             }
//         }
//     };
//     //
//     auto coupled_newton_step = [&] ()
//     {
//         auto glob_iter = 0;
//         while(glob_iter < 1000)
//         {
//             if (displacement_newton_step())
//             {
//                 auto set_energy = [&] (
//                     auto const & finite_element
//                 )
//                 {
//                     auto cnt = 0;
//                     for (auto & ip : finite_element->quadrature_.at("MicromorphicDamage").ips_)
//                     {
//                         auto const & dip = finite_element->quadrature_.at("MicromorphicDisplacement").ips_[cnt];
//                         auto const * rhs = mgis::behaviour::getInternalStateVariable(dip.behavior_data_->s1, "EnergyReleaseRate");
//                         auto * lhs = mgis::behaviour::getExternalStateVariable(ip.behavior_data_->s1, "EnergyReleaseRate");
//                         * lhs = * rhs;
//                         // mgis::behaviour::setExternalStateVariable(ip.behavior_data_->s1, "EnergyReleaseRate", * rhs);
//                         cnt ++;
//                     }
//                 };
//                 elements->caller<cells>("ROD", set_energy);
//             }
//             else
//             {
//                 return false;
//             }
//             if (damage_newton_step())
//             {
//                 auto set_energy = [&] (
//                     auto const & finite_element
//                 )
//                 {
//                     auto cnt = 0;
//                     for (auto & ip : finite_element->quadrature_.at("MicromorphicDisplacement").ips_)
//                     {
//                         auto const & dip = finite_element->quadrature_.at("MicromorphicDamage").ips_[cnt];
//                         auto const * rhs = mgis::behaviour::getInternalStateVariable(dip.behavior_data_->s1, "Damage");
//                         auto * lhs = mgis::behaviour::getExternalStateVariable(ip.behavior_data_->s1, "Damage");
//                         * lhs = * rhs;
//                         // mgis::behaviour::setExternalStateVariable(ip.behavior_data_->s1, "Damage", * rhs);
//                         cnt ++;
//                     }
//                 };
//                 elements->caller<cells>("ROD", set_energy);
//             }
//             else
//             {
//                 return false;
//             }
//             auto conv_test = displacement_convergence();
//             if (conv_test == 1)
//             {
//                 return true;
//             }
//             else if (conv_test == -1)
//             {
//                 return false;
//             }
//             glob_iter ++;
//         }
//         return false;
//     };
//     //
//     while (step < times.size())
//     {
//         // -> DEBUG
//         // break;
//         // <- DEBUG
//         std::cout << "**** step : " << step << " time : " << time << std::endl;
//         if (coupled_newton_step())
//         {
//             std::cout << "-- time step convergence" << std::endl;
//             auto top_force_value = elements->getBindingIntegral<faces, displacement_element, hdg>("TOP", "TopForce", 0, 0);
//             auto damage_stored_energy_value = elements->getStoredEnergy<cells>("ROD", "MicromorphicDamage");
//             auto damage_dissipated_energy_value = elements->getDissipatedEnergy<cells>("ROD", "MicromorphicDamage");
//             auto displacement_stored_energy_value = elements->getStoredEnergy<cells>("ROD", "MicromorphicDisplacement");
//             auto displacement_dissipated_energy_value = elements->getDissipatedEnergy<cells>("ROD", "MicromorphicDisplacement");
// force_out_stream << step << ", "<< time << ", " << std::setprecision(10) << std::scientific << top_force_value << "\n";
// damage_stored_energy_out_stream << step << ", "<< time << ", " << std::setprecision(10) << std::scientific << damage_stored_energy_value << "\n";
// damage_dissipated_energy_out_stream << step << ", "<< time << ", " << std::setprecision(10) << std::scientific << damage_dissipated_energy_value << "\n";
// displacement_stored_energy_out_stream << step << ", "<< time << ", " << std::setprecision(10) << std::scientific << displacement_stored_energy_value << "\n";
// displacement_dissipated_energy_out_stream << step << ", "<< time << ", " << std::setprecision(10) << std::scientific << displacement_dissipated_energy_value << "\n";
//             elements->reserveBehaviorData<cells>("ROD", "MicromorphicDisplacement");
//             elements->reserveUnknownCoefficients<cells, lolita::field("DenseVector", "A"), hdg.getCellBasis()>("ROD", "Displacement");
//             elements->reserveUnknownCoefficients<faces, lolita::field("DenseVector", "A"), hdg.getFaceBasis()>("ROD", "Displacement");
//             elements->reserveUnknownCoefficients<faces, lolita::field("Scalar", "A"), hdg.getFaceBasis()>("TOP", "TopForce");
//             elements->reserveUnknownCoefficients<faces, lolita::field("Scalar", "A"), hdg.getFaceBasis()>("LEFT", "LeftForce");
//             elements->reserveUnknownCoefficients<faces, lolita::field("Scalar", "A"), hdg.getFaceBasis()>("BOTTOM", "BottomForce");
//             //
//             elements->reserveBehaviorData<cells>("ROD", "MicromorphicDamage");
//             elements->reserveUnknownCoefficients<cells, lolita::field("Scalar", "A"), hdg.getCellBasis()>("ROD", "Damage");
//             elements->reserveUnknownCoefficients<faces, lolita::field("Scalar", "A"), hdg.getFaceBasis()>("ROD", "Damage");
//             #ifndef DEBUG
//                 // std::cout << "writing regular output" << std::endl;
//                 lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_displacement_file, elements, step, time, "MicromorphicDisplacement", 0);
//                 lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_displacement_file, elements, step, time, "MicromorphicDisplacement", 1);
//                 lolita::GmshFileParser::addQuadratureStressOutput<2, domain>(out_displacement_file, elements, step, time, "MicromorphicDisplacement", 0);
//                 lolita::GmshFileParser::addQuadratureStressOutput<2, domain>(out_displacement_file, elements, step, time, "MicromorphicDisplacement", 1);
//                 lolita::GmshFileParser::addQuadratureDofOutput<2, domain, displacement_element, hdg>(out_displacement_file, elements, step, time, "Displacement", "MicromorphicDisplacement", 0, 0);
//                 lolita::GmshFileParser::addQuadratureDofOutput<2, domain, displacement_element, hdg>(out_displacement_file, elements, step, time, "Displacement", "MicromorphicDisplacement", 1, 0);
//                 //
//                 lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_damage_file, elements, step, time, "MicromorphicDamage", 0);
//                 lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_damage_file, elements, step, time, "MicromorphicDamage", 1);
//                 lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_damage_file, elements, step, time, "MicromorphicDamage", 2);
//                 lolita::GmshFileParser::addQuadratureStressOutput<2, domain>(out_damage_file, elements, step, time, "MicromorphicDamage", 0);
//                 lolita::GmshFileParser::addQuadratureStressOutput<2, domain>(out_damage_file, elements, step, time, "MicromorphicDamage", 1);
//                 lolita::GmshFileParser::addQuadratureStressOutput<2, domain>(out_damage_file, elements, step, time, "MicromorphicDamage", 2);
//                 lolita::GmshFileParser::addQuadratureInternalVariableOutput<2, domain>(out_damage_file, elements, step, time, "MicromorphicDamage", 0);
//                 // lolita::GmshFileParser::addQuadratureInternalVariableOutput<2, domain>(out_damage_file, elements, step, time, "MicromorphicDamage", 1);
//                 // lolita::GmshFileParser::addQuadratureInternalVariableOutput<2, domain>(out_damage_file, elements, step, time, "MicromorphicDamage", 2);
//                 // lolita::GmshFileParser::addQuadratureInternalVariableOutput<2, domain>(out_damage_file, elements, step, time, "MicromorphicDamage", 3);
//                 lolita::GmshFileParser::addQuadratureDofOutput<2, domain, displacement_element, hdg>(out_damage_file, elements, step, time, "Damage", "MicromorphicDamage", 0, 0);
//             #endif
//             step ++;
//             time = times[step];
//             /* code */
//         }
//         else
//         {
//             std::cout << "-- time step split" << std::endl;
//             elements->recoverBehaviorData<cells>("ROD", "MicromorphicDisplacement");
//             elements->recoverUnknownCoefficients<cells, lolita::field("DenseVector", "A"), hdg.getCellBasis()>("ROD", "Displacement");
//             elements->recoverUnknownCoefficients<faces, lolita::field("DenseVector", "A"), hdg.getFaceBasis()>("ROD", "Displacement");
//             elements->recoverUnknownCoefficients<faces, lolita::field("Scalar", "A"), hdg.getFaceBasis()>("TOP", "TopForce");
//             elements->recoverUnknownCoefficients<faces, lolita::field("Scalar", "A"), hdg.getFaceBasis()>("LEFT", "LeftForce");
//             elements->recoverUnknownCoefficients<faces, lolita::field("Scalar", "A"), hdg.getFaceBasis()>("BOTTOM", "BottomForce");
//             //
//             elements->recoverBehaviorData<cells>("ROD", "MicromorphicDamage");
//             elements->recoverUnknownCoefficients<cells, lolita::field("Scalar", "A"), hdg.getCellBasis()>("ROD", "Damage");
//             elements->recoverUnknownCoefficients<faces, lolita::field("Scalar", "A"), hdg.getFaceBasis()>("ROD", "Damage");
//             time = times[step - 1] + (1.0 / 2.0) * (time - times[step - 1]);
//             // -> DEBUG
//             // step ++;
//             // time = times[step];
//             // <- DEBUG
//         }
//     }
}