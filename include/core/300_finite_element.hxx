#ifndef B9B48BEA_09F5_41DB_84A1_45A6E708901C
#define B9B48BEA_09F5_41DB_84A1_45A6E708901C

#include "lolita.hxx"
#include "core/000_physics_traits.hxx"
#include "core/100_geometry.hxx"
#include "core/linear_system.hxx"
#include "core/200_quadrature.hxx"
#include "core/201_finite_element_dof.hxx"
#include "core/202_finite_element_frm.hxx"

namespace lolita
{

    template<Integer t_dim, Domain t_domain>
    struct FiniteDomain
    {

        explicit
        FiniteDomain(
            std::basic_string<Character> const & tag
        )
        :
        tag_(tag)
        {}

        explicit
        FiniteDomain(
            std::basic_string<Character> && tag
        )
        :
        tag_(std::move(tag))
        {}

        std::basic_string<Character> const &
        getLabel()
        const
        {
            return tag_;
        }

        template<FieldConcept auto t_field>
        void
        addDiscreteField()
        {
            if (ptr_data_ == nullptr)
            {
                ptr_data_ = std::make_unique<std::vector<MeshDiscreteField<t_domain>>>();
            }
            for (auto const & item : * ptr_data_)
            {
                if (item.getLabel() == t_field.getLabel())
                {
                    return;
                }
            }
            ptr_data_->push_back(MeshDiscreteField<t_domain>(t_field));
        }

        template<FieldConcept auto t_field>
        Boolean
        hasDiscreteField()
        {
            if (ptr_data_ == nullptr)
            {
                return false;
            }
            else
            {
                for (auto const & item : * ptr_data_)
                {
                    if (item.getLabel() == t_field.getLabel())
                    {
                        return true;
                    }
                }
                throw std::runtime_error("No such field data");
            }
            return false;
        }

        template<FieldConcept auto t_field>
        MeshDiscreteField<t_domain> const &
        getDiscreteField()
        const
        {
            if (ptr_data_ == nullptr)
            {
                throw std::runtime_error("Empty");
            }
            else
            {
                for (auto const & item : * ptr_data_)
                {
                    if (item.getLabel() == t_field.getLabel())
                    {
                        return item;
                    }
                }
                throw std::runtime_error("No such field data");
            }
        }

        template<FieldConcept auto t_field>
        MeshDiscreteField<t_domain> &
        getDiscreteField()
        {
            if (ptr_data_ == nullptr)
            {
                throw std::runtime_error("Empty");
            }
            else
            {
                for (auto & item : * ptr_data_)
                {
                    if (item.getLabel() == t_field.getLabel())
                    {
                        return item;
                    }
                }
                throw std::runtime_error("No such field data");
            }
        }

        template<FieldConcept auto t_field, Integer t_size>
        void
        addDiscreteFieldDegreeOfFreedom(
            auto const &... args
        )
        {
            this->template getDiscreteField<t_field>().template addDegreeOfFreedom<t_field, t_size>(args...);
        }

        template<FieldConcept auto t_field>
        void
        addDiscreteFieldLoad(
            Integer row,
            Integer col,
            std::function<Real(Point const &, Real const &)> && function
        )
        {
            this->template getDiscreteField<t_field>().addLoad(row, col, std::forward<std::function<Real(Point const &, Real const &)>>(function));
        }

    private:

        std::unique_ptr<std::vector<MeshDiscreteField<t_domain>>> ptr_data_;

        std::basic_string<Character> tag_;

    };
    
    template<Element t_element, Domain t_domain>
    struct FiniteElementTraits<t_element, t_domain>
    {
        
        template<DiscretizationConcept auto t_discretization>
        static constexpr
        Integer
        getNumElementCoefficients()
        {
            return DiscretizationTraits3<t_discretization>::template getNumElementCoefficients<t_element, t_domain>();
        }

        template<FieldConcept auto t_field>
        static constexpr
        Integer
        getNumElementCoefficients()
        {
            return getNumElementCoefficients<t_field.getDiscretization()>() * FieldTraits<t_field>::template getSize<t_domain>();
            // return FieldDiscretizationTraits<t_field>::template getNumElementCoefficients<t_element, t_domain, t_field>();
        }
        
        template<DiscretizationConcept auto t_discretization>
        static constexpr
        Integer
        getNumUnknownCoefficients1()
        {
            auto num_unknowns = DiscretizationTraits3<t_discretization>::template getNumElementCoefficients<t_element, t_domain>();
            auto set_num_unknowns = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & t_set_num_unknowns
            )
            constexpr mutable
            {
                auto constexpr inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<t_i, t_j>();
                auto constexpr num_inner_neighbors = ElementTraits<t_element>::template getNumInnerNeighbors<t_i, t_j>();
                num_unknowns += DiscretizationTraits3<t_discretization>::template getNumElementCoefficients<inner_neighbor, t_domain>() * num_inner_neighbors;
                if constexpr (t_j < ElementTraits<t_element>::template getNumInnerNeighbors<t_i>() - 1)
                {
                    t_set_num_unknowns.template operator ()<t_i, t_j + 1>(t_set_num_unknowns);
                }
                else if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<>() - 1)
                {
                    t_set_num_unknowns.template operator ()<t_i + 1, 0>(t_set_num_unknowns);
                }
            };
            set_num_unknowns(set_num_unknowns);
            return num_unknowns;
        }

        template<FieldConcept auto t_field>
        static constexpr
        Integer
        getNumUnknownCoefficients()
        {
            return getNumUnknownCoefficients1<t_field.getDiscretization()>() * FieldTraits<t_field>::template getSize<t_domain>();
        }

    };
    
    template<Element t_element, Domain t_domain>
    struct FiniteElement
    {

    private:
    
        using t_ElementTraits = ElementTraits<t_element>;
        
        template<Element t__element, Domain t__domain>
        using t_Neighbor = std::shared_ptr<FiniteElement<t__element, t__domain>>;

        using t_Domains = typename t_ElementTraits::template Domains<FiniteDomain, t_domain>;
    
        using t_InnerNeighbors = typename t_ElementTraits::template InnerConnectivity<t_Neighbor, t_domain>;
        
        using t_OuterNeighbors = typename t_ElementTraits::template OuterConnectivity<t_Neighbor, t_domain>;

        template<Basis t_basis>
        using t_Basis = typename BasisTraits<t_basis>::template Implementation<t_element, t_domain>;

        template<auto t_discretization>
        using t_Disc = typename DiscretizationTraits<t_discretization>::template Implementation<t_element, t_domain>;

        template<auto t_field>
        using t_Disc2 = typename FieldDiscretizationTraits<t_field>::template Implementation<t_element, t_domain, t_field>;
        
        Natural tag_;
        
        Point coordinates_;

        t_Domains domains_;
        
        t_OuterNeighbors outer_neighbors_;
        
        t_InnerNeighbors inner_neighbors_;

        //

        // std::unique_ptr<std::vector<ElementDiscreteField<t_domain>>> ptr_data_;

        std::unique_ptr<std::vector<ElementDiscreteField<t_element, t_domain>>> fields_;

        std::unique_ptr<std::vector<ElementFormulation<t_domain>>> ptr_formulations_;

        //

    public:

        explicit
        FiniteElement(
            Natural const & tag
        )
        :
        tag_(tag),
        coordinates_(),
        domains_(),
        outer_neighbors_(),
        inner_neighbors_()
        {}

        Point const &
        getCoordinates()
        const
        {
            return coordinates_;
        }
        
        template<Basis t_basis>
        DenseVector<Real, BasisTraits<t_basis>::template getSize<t_element>()>
        getBasisEvaluation(
            auto const &... args
        )
        const
        {
            return static_cast<t_Basis<t_basis> const *>(this)->getBasisEvaluation(args...);
        }
        
        template<Basis t_basis>
        DenseVector<Real, BasisTraits<t_basis>::template getSize<t_element>()>
        getBasisDerivative(
            auto const &... args
        )
        const
        {
            return static_cast<t_Basis<t_basis> const *>(this)->getBasisDerivative(args...);
        }

        Natural const &
        getTag()
        const
        {
            return tag_;
        }

        t_Domains const &
        getDomains()
        const
        {
            return domains_;
        }

        void
        setCoordinates(
            Point const & point
        )
        {
            coordinates_ = point;
        }

        void
        addDomain(
            std::shared_ptr<FiniteDomain<t_element.getDim(), t_domain>> const & domain
        )
        {
            for (auto const & finite_element_domain : domains_)
            {
                if (finite_element_domain == domain)
                {
                    return;
                }
            }
            domains_.push_back(domain);
        }

        /**
         * This part is dedicated to the Dof implementation
         * *****************************************************************************************************************************************************
         */
        
        template<FieldConcept auto t_field>
        void
        addDiscreteField()
        {
            if (fields_ == nullptr)
            {
                fields_ = std::make_unique<std::vector<ElementDiscreteField<t_element, t_domain>>>();
            }
            for (auto const & field : * fields_)
            {
                if (field.getLabel() == t_field.getLabel())
                {
                    return;
                }
            }
            fields_->push_back(ElementDiscreteField<t_element, t_domain>(t_field));
        }
        
        template<FieldConcept auto t_field>
        Boolean
        hasDiscreteField()
        const
        {
            if (fields_ == nullptr)
            {
                return false;
            }
            for (auto const & field : * fields_)
            {
                if (field.getLabel() == t_field.getLabel())
                {
                    return true;
                }
            }
            return false;
        }

        template<FieldConcept auto t_field>
        ElementDiscreteField<t_element, t_domain> const &
        getDiscreteField()
        const
        {
            if (fields_ == nullptr)
            {
                throw std::runtime_error("Empty");
            }
            else
            {
                for (auto const & field : * fields_)
                {
                    if (field.getLabel() == t_field.getLabel())
                    {
                        return field;
                    }
                }
                throw std::runtime_error("No such field data");
            }
        }

        template<FieldConcept auto t_field>
        ElementDiscreteField<t_element, t_domain> &
        getDiscreteField()
        {
            if (fields_ == nullptr)
            {
                throw std::runtime_error("Empty");
            }
            else
            {
                for (auto & field : * fields_)
                {
                    if (field.getLabel() == t_field.getLabel())
                    {
                        return field;
                    }
                }
                throw std::runtime_error("No such field data");
            }
        }

        template<FieldConcept auto t_field, Strategy t_s>
        void
        addDiscreteFieldDegreeOfFreedomToSystem(
            LinearSystem<t_s> & linear_system
        )
        {
            auto & dof = this->template getDiscreteField<t_field>().getDegreeOfFreedom();
            dof.link(linear_system.getRhsSize(), linear_system.getLhsSize());
            // linear_system.addRhsSize(FieldDiscretizationTraits<t_field>::template getNumElementCoefficients<t_element, t_domain, t_field>());
            linear_system.addRhsSize(FiniteElementTraits<t_element, t_domain>::template getNumElementCoefficients<t_field>());
            linear_system.addLhsSize(this->template getBandWidth<t_field>());
        }

        template<FieldConcept auto t_field>
        auto
        getDiscreteFieldDegreeOfFreedom()
        const
        {
            return this->template getDiscreteField<t_field>().getDegreeOfFreedom().template getCoefficients<t_field>();
        }

        template<FieldConcept auto t_field>
        auto
        getDiscreteFieldDegreeOfFreedom(
            Integer row,
            Integer col
        )
        const
        {
            return this->template getDiscreteField<t_field>().getDegreeOfFreedom().template getCoefficients<t_field>(row, col);
        }

        template<FieldConcept auto t_field>
        void
        addDiscreteFieldDegreeOfFreedom()
        {
            this->template getDiscreteField<t_field>().template addDegreeOfFreedom<t_field>();
        }

        template<FieldConcept auto t_field>
        void
        upgradeDiscreteFieldDegreeOfFreedom(
            DenseVectorConcept<Real> auto && input
        )
        {
            this->template getDiscreteField<t_field>().getDegreeOfFreedom().template upgradeCoefficients<t_field>(std::forward<decltype(input)>(input));
        }

        template<FieldConcept auto t_field>
        void
        recoverDiscreteFieldDegreeOfFreedom()
        {
            this->template getDiscreteField<t_field>().getDegreeOfFreedom().template recoverCoefficients<t_field>();
        }

        template<FieldConcept auto t_field>
        void
        reserveDiscreteFieldDegreeOfFreedom()
        {
            this->template getDiscreteField<t_field>().getDegreeOfFreedom().template reserveCoefficients<t_field>();
        }

        template<FieldConcept auto t_field>
        Real
        getDiscreteFieldDegreeOfFreedomValue(
            PointConcept auto const & point,
            Integer row,
            Integer col
        )
        const
        {
            auto const coefficients = this->template getUnknownCoefficients<t_field>();
            auto const basis_vector = this->template getFieldDualVector<t_field>(point, row, col);
            return coefficients.dot(basis_vector);
        }

        template<FieldConcept auto t_field>
        DenseVector<Real, FiniteElementTraits<t_element, t_domain>::template getNumUnknownCoefficients<t_field>()>
        getUnknownCoefficients()
        const
        {
            auto constexpr cell_range = FiniteElementTraits<t_element, t_domain>::template getNumElementCoefficients<t_field>();
            auto offset = 0;
            auto unknown = DenseVector<Real, FiniteElementTraits<t_element, t_domain>::template getNumUnknownCoefficients<t_field>()>();
            auto const & cell_dof = this->template getDiscreteField<t_field>().getDegreeOfFreedom();
            auto cell_block = unknown.template segment<cell_range>(offset);
            cell_block = cell_dof.template getCoefficients<t_field>();
            offset += cell_range;
            auto set_faces_unknowns = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & self
            )
            constexpr mutable
            {
                auto constexpr t_inner_neighbor = ElementTraits<t_element>::template getInnerNeighbor<t_i, t_j>();
                auto constexpr range = FiniteElementTraits<t_inner_neighbor, t_domain>::template getNumElementCoefficients<t_field>();
                if constexpr (range > 0)
                {
                    for (auto const & face : this->template getInnerNeighbors<t_i, t_j>())
                    {
                        auto const & face_dof = face->template getDiscreteField<t_field>().getDegreeOfFreedom();
                        auto face_block = unknown.template segment<range>(offset);
                        face_block = face_dof.template getCoefficients<t_field>();
                        offset += range;
                    }
                }
                if constexpr (t_j < ElementTraits<t_element>::template getNumInnerNeighbors<t_i>() - 1)
                {
                    self.template operator ()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < ElementTraits<t_element>::template getNumInnerNeighbors<>() - 1)
                {
                    self.template operator ()<t_i + 1, 0>(self);
                }
            };
            set_faces_unknowns(set_faces_unknowns);
            return unknown;
            // return static_cast<t_Disc2<t_field> const *>(this)->getUnknownCoefficients();
        }
        
        template<FieldConcept auto t_field>
        Integer
        getBandWidth()
        const
        {
            auto band_width = FiniteElementTraits<t_element, t_domain>::template getNumUnknownCoefficients<t_field>();
            auto set_faces_unknowns = [&] <Integer t_i = 0> (
                auto & self
            )
            constexpr mutable
            {
                auto constexpr t_cell = ElementTraits<t_element>::template getOuterNeighbor<t_domain, t_field.getDomainDim() - t_element.getDim(), t_i>();
                if constexpr (t_element.getDim() <  t_field.getDomainDim())
                {
                    for (auto const & cell : this->template getOuterNeighbors<t_field.getDomainDim() - t_element.getDim(), t_i>())
                    {
                        auto set_faces_unknowns2 = [&] <Integer t_j = 0> (
                            auto & self2
                        )
                        constexpr mutable
                        {
                            auto set_faces_unknowns3 = [&] <Integer t_k = 0> (
                                auto & self3
                            )
                            constexpr mutable
                            {
                                auto constexpr t_cell_neighbor = ElementTraits<t_cell>::template getInnerNeighbor<t_j, t_k>();
                                for (auto const & cell_neighbor : cell->template getInnerNeighbors<t_j, t_k>())
                                {
                                    if (!utility::areEqual(* cell_neighbor, * this) && cell_neighbor->template hasDiscreteField<t_field>())
                                    {
                                        // FiniteElementTraits<t_cell_neighbor, t_domain>::template getNumUnknownCoefficients<t_field>()
                                        band_width += FiniteElementTraits<t_cell_neighbor, t_domain>::template getNumUnknownCoefficients<t_field>();
                                    }
                                }
                                if constexpr (t_k < ElementTraits<t_cell>::template getNumInnerNeighbors<t_j>() - 1)
                                {
                                    self3.template operator ()<t_k + 1>(self3);
                                }
                            };
                            set_faces_unknowns3(set_faces_unknowns3);
                            if constexpr (t_j < t_field.getDomainDim() - t_element.getDim() - 1)
                            {
                                self2.template operator ()<t_j + 1>(self2);
                            }
                        };
                        set_faces_unknowns2(set_faces_unknowns2);
                    }
                    if constexpr (t_i < ElementTraits<t_element>::template getNumOuterNeighbors<t_domain, t_field.getDomainDim() - t_element.getDim()>() - 1)
                    {
                        self.template operator ()<t_i + 1>(self);
                    }
                }
            };
            set_faces_unknowns(set_faces_unknowns);
            return band_width;
        }

        /**
         * Formulation
         * *****************************************************************************************************************************************************
         */
        
        template<PotentialConcept auto... t_potential>
        using ElementJacobianMatrix = typename PotentialTraits<t_potential...>::template ElementJacobianMatrix<t_element, t_domain>;
        
        template<PotentialConcept auto... t_potential>
        using ElementResidualVector = typename PotentialTraits<t_potential...>::template ElementResidualVector<t_element, t_domain>;

        template<FieldConcept auto t_field>
        DenseVector<Real, FiniteElementTraits<t_element, t_domain>::template getNumUnknownCoefficients<t_field>()>
        getFieldDualVector(
            PointConcept auto const & point
        )
        const
        {
            return static_cast<t_Disc2<t_field> const *>(this)->getFieldDualVector(point);
        };

        template<FieldConcept auto t_field>
        DenseVector<Real, FiniteElementTraits<t_element, t_domain>::template getNumUnknownCoefficients<t_field>()>
        getFieldDualVector(
            PointConcept auto const & point,
            Integer row,
            Integer col
        )
        const
        {
            return static_cast<t_Disc2<t_field> const *>(this)->getFieldDualVector(point, row, col);
        };

        template<FieldConcept auto t_field>
        DenseVector<Real, FiniteElementTraits<t_element, t_domain>::template getNumUnknownCoefficients<t_field>()>
        getFieldPrimalVector(
            PointConcept auto const & point
        )
        const
        {
            return static_cast<t_Disc2<t_field> const *>(this)->getFieldPrimalVector(point);
        };

        template<FieldConcept auto t_field>
        DenseVector<Real, FiniteElementTraits<t_element, t_domain>::template getNumUnknownCoefficients<t_field>()>
        getFieldPrimalVector(
            PointConcept auto const & point,
            Integer row,
            Integer col
        )
        const
        {
            return static_cast<t_Disc2<t_field> const *>(this)->getFieldPrimalVector(point, row, col);
        };

        template<PotentialConcept auto... t_potential>
        void
        setElementStrainOperators()
        {
            auto set_i = [&] <Integer t_i = 0> (
                auto & t_set_i
            )
            constexpr mutable
            {
                auto constexpr potential = PotentialTraits<t_potential...>::template getPotential<t_i>();
                auto & frm = this->template getFormulation<potential>();
                auto set_j = [&] <Integer t_j = 0> (
                    auto & t_set_j
                )
                constexpr mutable
                {
                    auto constexpr j_mapping = PotentialTraits<potential>::template getStrain<t_j>();
                    for (auto & ip : frm.getIntegrationPoints())
                    {
                        auto strain_operator = this->template getMapping<j_mapping>(ip.getReferenceCoordinates());
                        ip.template addStrainOperator<j_mapping>(strain_operator);
                    }
                    if constexpr (t_j < PotentialTraits<potential>::getNumMappings() - 1)
                    {
                        t_set_j.template operator()<t_j + 1>(t_set_j);
                    }
                };
                set_j(set_j);
                if constexpr (t_i < PotentialTraits<t_potential...>::getNumPotentials() - 1)
                {
                    t_set_i.template operator()<t_i + 1>(t_set_i);
                }
            };
            set_i(set_i);
        }

        template<PotentialConcept auto... t_potential>
        void
        setElementStrains()
        {
            auto set_i = [&] <Integer t_i = 0> (
                auto & t_set_i
            )
            constexpr mutable
            {
                auto constexpr potential = PotentialTraits<t_potential...>::template getPotential<t_i>();
                auto & frm = this->template getFormulation<potential>();
                auto set_j = [&] <Integer t_j = 0> (
                    auto & t_set_j
                )
                constexpr mutable
                {
                    auto constexpr j_mapping = PotentialTraits<potential>::template getStrain<t_j>();
                    auto constexpr j_field = j_mapping.getField();
                    auto const unknown = this->template getUnknownCoefficients<j_field>();
                    for (auto & ip : frm.getIntegrationPoints())
                    {
                        auto const & mat0 = ip.template getStrainOperator<j_mapping>();
                        ip.template setStrainVectorBlock<potential, j_mapping>(mat0 * unknown);
                    }
                    if constexpr (t_j < PotentialTraits<potential>::getNumMappings() - 1)
                    {
                        t_set_j.template operator()<t_j + 1>(t_set_j);
                    }
                };
                set_j(set_j);
                if constexpr (t_i < PotentialTraits<t_potential...>::getNumPotentials() - 1)
                {
                    t_set_i.template operator()<t_i + 1>(t_set_i);
                }
            };
            set_i(set_i);
        }
        
        template<PotentialConcept auto... t_potential>
        Real
        getElementInternalEnergy()
        const
        {
            auto internal_forces = Real(0);
            auto set_i = [&] <Integer t_i = 0> (
                auto & t_set_i
            )
            constexpr mutable
            {
                auto constexpr potential = PotentialTraits<t_potential...>::template getPotential<t_i>();
                auto const & frm = this->template getFormulation<potential>();
                auto set_j = [&] <Integer t_j = 0> (
                    auto & t_set_j
                )
                constexpr mutable
                {
                    auto constexpr j_mapping = PotentialTraits<potential>::template getStrain<t_j>();
                    auto constexpr j_field = j_mapping.getField();
                    for (auto const & ip : frm.getIntegrationPoints())
                    {
                        auto const stress = ip.template getStressVectorBlock<potential, j_mapping>();
                        auto const strain = ip.template getStrainVectorBlock<potential, j_mapping>();
                        internal_forces += ip.getCurrentWeight() * strain.transpose() * stress;
                    }
                    if constexpr (t_j < PotentialTraits<potential>::getNumMappings() - 1)
                    {
                        t_set_j.template operator()<t_j + 1>(t_set_j);
                    }
                };
                set_j(set_j);
                if constexpr (t_i < PotentialTraits<t_potential...>::getNumPotentials() - 1)
                {
                    t_set_i.template operator()<t_i + 1>(t_set_i);
                }
            };
            set_i(set_i);
            return internal_forces;
        }
        
        template<PotentialConcept auto... t_potential>
        ElementResidualVector<t_potential...>
        getElementInternalForces()
        const
        {
            auto internal_forces = ElementResidualVector<t_potential...>();
            internal_forces.setZero();
            auto set_i = [&] <Integer t_i = 0> (
                auto & t_set_i
            )
            constexpr mutable
            {
                auto constexpr potential = PotentialTraits<t_potential...>::template getPotential<t_i>();
                auto const & frm = this->template getFormulation<potential>();
                auto set_j = [&] <Integer t_j = 0> (
                    auto & t_set_j
                )
                constexpr mutable
                {
                    auto constexpr j_mapping = PotentialTraits<potential>::template getStrain<t_j>();
                    auto constexpr j_field = j_mapping.getField();
                    // --> making block
                    auto constexpr size_j = FiniteElementTraits<t_element, t_domain>::template getNumUnknownCoefficients<j_field>();
                    auto constexpr offset_j = PotentialTraits<t_potential...>::template getOffset<t_element, t_domain, j_field>();
                    for (auto const & ip : frm.getIntegrationPoints())
                    {
                        auto const & mat0 = ip.template getStrainOperator<j_mapping>();
                        auto const mat = ip.template getStressVectorBlock<potential, j_mapping>();
                        internal_forces.template segment<size_j>(offset_j) += ip.getCurrentWeight() * mat0.transpose() * mat;
                    }
                    if constexpr (t_j < PotentialTraits<potential>::getNumMappings() - 1)
                    {
                        t_set_j.template operator()<t_j + 1>(t_set_j);
                    }
                };
                set_j(set_j);
                if constexpr (t_i < PotentialTraits<t_potential...>::getNumPotentials() - 1)
                {
                    t_set_i.template operator()<t_i + 1>(t_set_i);
                }
            };
            set_i(set_i);
            return internal_forces;
        }
        
        template<FieldConcept auto... t_fields>
        Real
        getElementExternalEnergy(
            Real const & time
        )
        const
        {
            auto constexpr fields = utility::Aggregate<std::decay_t<decltype(t_fields)>...>(t_fields...);
            auto constexpr quadrature = Quadrature("Gauss", 4);
            auto external_forces = Real(0);
            auto set_i = [&] <Integer t_i = 0> (
                auto & t_set_i
            )
            constexpr mutable
            {
                auto constexpr field = utility::get<t_i>(fields);
                for (auto const & domain : this->getDomains())
                {
                    if (domain->template hasDiscreteField<field>())
                    {
                        if (domain->template getDiscreteField<field>().hasLoads())
                        {
                            for (auto const & load : domain->template getDiscreteField<field>().getLoads())
                            {
                                for (auto q = 0; q < QuadratureTraits<quadrature>::template Rule<t_element>::getSize(); q++)
                                {
                                    auto weight = this->template getCurrentQuadratureWeight<quadrature>(q);
                                    auto reference_point = this->template getReferenceQuadraturePoint<quadrature>(q);
                                    auto current_point = this->template getCurrentQuadraturePoint<quadrature>(q);
                                    auto i = load.getRow();
                                    auto j = load.getCol();
                                    auto field_value = this->template getDiscreteFieldDegreeOfFreedomValue<field>(reference_point, i, j);
                                    external_forces += weight * load.getValue(current_point, time) * field_value;
                                }
                            }
                        }
                    }
                }
                if constexpr (t_i < sizeof...(t_fields) - 1)
                {
                    t_set_i.template operator()<t_i + 1>(t_set_i);
                }
            };
            set_i(set_i);
            return external_forces;
        }
        
        template<PotentialConcept auto... t_potential>
        Real
        getElementExternalEnergy(
            Real const & time
        )
        const
        {
            auto constexpr quadrature = Quadrature("Gauss", 4);
            auto external_forces = Real(0);
            auto set_i = [&] <Integer t_i = 0> (
                auto & t_set_i
            )
            constexpr mutable
            {
                auto constexpr potential = PotentialTraits<t_potential...>::template getPotential<t_i>();
                auto set_j = [&] <Integer t_j = 0> (
                    auto & t_set_j
                )
                constexpr mutable
                {
                    auto constexpr field = PotentialTraits<potential>::template getField<t_j>();
                    for (auto const & domain : this->getDomains())
                    {
                        if (domain->template hasDiscreteField<field>())
                        {
                            if (domain->template getDiscreteField<field>().hasLoads())
                            {
                                for (auto const & load : domain->template getDiscreteField<field>().getLoads())
                                {
                                    for (auto i = 0; i < QuadratureTraits<quadrature>::template Rule<t_element>::getSize(); i++)
                                    {
                                        auto weight = this->template getCurrentQuadratureWeight<quadrature>(i);
                                        auto reference_point = this->template getReferenceQuadraturePoint<quadrature>(i);
                                        auto current_point = this->template getCurrentQuadraturePoint<quadrature>(i);
                                        auto x = load.getRow();
                                        auto y = load.getCol();
                                        // auto vector = this->template getFieldDualVector<field>(reference_point, x, y);
                                        // auto field_value = this->template getDiscreteFieldDegreeOfFreedomValue<field>(reference_point, x, y);
                                        auto field_value = this->template getDiscreteFieldDegreeOfFreedomValue<field>(reference_point, x, y);
                                        external_forces += weight * load.getValue(current_point, time) * field_value;
                                    }
                                }
                            }
                        }
                    }
                    if constexpr (t_j < PotentialTraits<potential>::getNumFields() - 1)
                    {
                        t_set_j.template operator()<t_j + 1>(t_set_j);
                    }
                };
                set_j(set_j);
                if constexpr (t_i < PotentialTraits<t_potential...>::getNumPotentials() - 1)
                {
                    t_set_i.template operator()<t_i + 1>(t_set_i);
                }
            };
            set_i(set_i);
            return external_forces;
        }
        
        template<PotentialConcept auto... t_potential>
        Real
        getElementResidualEnergy(
            Real const & time
        )
        const
        {
            return this->template getElementExternalEnergy<t_potential...>(time) - this->template getElementInternalEnergy<t_potential...>();
        }

        template<PotentialConcept auto... t_potential>
        ElementResidualVector<t_potential...>
        getElementExternalForces(
            Real const & time
        )
        const
        {
            auto constexpr quadrature = Quadrature("Gauss", 4);
            auto external_forces = ElementResidualVector<t_potential...>();
            external_forces.setZero();
            auto set_i = [&] <Integer t_i = 0> (
                auto & t_set_i
            )
            constexpr mutable
            {
                auto constexpr potential = PotentialTraits<t_potential...>::template getPotential<t_i>();
                auto set_j = [&] <Integer t_j = 0> (
                    auto & t_set_j
                )
                constexpr mutable
                {
                    auto constexpr field = PotentialTraits<potential>::template getField<t_j>();
                    auto constexpr size_j = FiniteElementTraits<t_element, t_domain>::template getNumUnknownCoefficients<field>();
                    auto constexpr offset_j = PotentialTraits<t_potential...>::template getOffset<t_element, t_domain, field>();
                    for (auto const & domain : this->getDomains())
                    {
                        if (domain->template hasDiscreteField<field>())
                        {
                            if (domain->template getDiscreteField<field>().hasLoads())
                            {
                                for (auto const & load : domain->template getDiscreteField<field>().getLoads())
                                {
                                    for (auto i = 0; i < QuadratureTraits<quadrature>::template Rule<t_element>::getSize(); i++)
                                    {
                                        auto weight = this->template getCurrentQuadratureWeight<quadrature>(i);
                                        auto reference_point = this->template getReferenceQuadraturePoint<quadrature>(i);
                                        auto current_point = this->template getCurrentQuadraturePoint<quadrature>(i);
                                        auto vector = this->template getFieldDualVector<field>(reference_point, load.getRow(), load.getCol());
                                        external_forces.template segment<size_j>(offset_j) += weight * load.getValue(current_point, time) * vector;
                                    }
                                }
                            }
                        }
                    }
                    if constexpr (t_j < PotentialTraits<potential>::getNumFields() - 1)
                    {
                        t_set_j.template operator()<t_j + 1>(t_set_j);
                    }
                };
                set_j(set_j);
                if constexpr (t_i < PotentialTraits<t_potential...>::getNumPotentials() - 1)
                {
                    t_set_i.template operator()<t_i + 1>(t_set_i);
                }
            };
            set_i(set_i);
            return external_forces;
        }
        
        template<PotentialConcept auto... t_potential>
        ElementResidualVector<t_potential...>
        getElementResidualForces(
            Real const & time
        )
        const
        {
            return this->template getElementExternalForces<t_potential...>(time) - this->template getElementInternalForces<t_potential...>();
        }
        
        template<PotentialConcept auto... t_potential>
        void
        getElementResidualForces(
            Real const & time,
            std::atomic<Real> & i,
            std::atomic<Real> & o
        )
        const
        {
            auto internal_forces = this->template getElementInternalForces<t_potential...>();
            auto external_forces = this->template getElementExternalForces<t_potential...>(time);
            auto residual_forces = external_forces - internal_forces;
            auto residual_forces_max = residual_forces.cwiseAbs().maxCoeff();
            auto external_forces_max = external_forces.cwiseAbs().maxCoeff();
            if (residual_forces_max > i)
            {
                i = residual_forces_max;
            }
            if (external_forces_max > o)
            {
                o = external_forces_max;
            }
            return this->template getElementExternalForces<t_potential...>(time) - this->template getElementInternalForces<t_potential...>();
        }
        
        template<PotentialConcept auto... t_potential>
        ElementJacobianMatrix<t_potential...>
        getElementJacobianMatrix()
        const
        {
            auto jacobian_matrix = ElementJacobianMatrix<t_potential...>();
            jacobian_matrix.setZero();
            auto set_i = [&] <Integer t_i = 0> (
                auto & t_set_i
            )
            constexpr mutable
            {
                auto constexpr potential = PotentialTraits<t_potential...>::template getPotential<t_i>();
                auto const & frm = this->template getFormulation<potential>();
                // --> J
                auto set_j = [&] <Integer t_j = 0> (
                    auto & t_set_j
                )
                constexpr mutable
                {
                    auto constexpr j_mapping = PotentialTraits<potential>::template getStrain<t_j>();
                    auto constexpr j_field = j_mapping.getField();
                    // --> K
                    auto set_k = [&] <Integer t_k = 0> (
                        auto & t_set_k
                    )
                    constexpr mutable
                    {
                        auto constexpr k_mapping = PotentialTraits<potential>::template getStrain<t_k>();
                        auto constexpr k_field = k_mapping.getField();
                        // --> making block
                        // auto constexpr size_j = FieldDiscretizationTraits<j_field>::template getNumUnknownCoefficients<t_element, t_domain, j_field>();
                        // auto constexpr size_k = FieldDiscretizationTraits<k_field>::template getNumUnknownCoefficients<t_element, t_domain, k_field>();
                        auto constexpr size_j = FiniteElementTraits<t_element, t_domain>::template getNumUnknownCoefficients<j_field>();
                        auto constexpr size_k = FiniteElementTraits<t_element, t_domain>::template getNumUnknownCoefficients<k_field>();
                        auto constexpr offset_j = PotentialTraits<t_potential...>::template getOffset<t_element, t_domain, j_field>();
                        auto constexpr offset_k = PotentialTraits<t_potential...>::template getOffset<t_element, t_domain, k_field>();
                        for (auto const & ip : frm.getIntegrationPoints())
                        {
                            auto const & mat0 = ip.template getStrainOperator<j_mapping>();
                            auto const & mat1 = ip.template getStrainOperator<k_mapping>();
                            auto const mat = ip.template getJacobianMatrixBlock<potential, j_mapping, k_mapping>();
                            jacobian_matrix.template block<size_j, size_k>(offset_j, offset_k) += ip.getCurrentWeight() * mat0.transpose() * mat * mat1;
                        }
                        if constexpr (t_k < PotentialTraits<potential>::getNumMappings() - 1)
                        {
                            t_set_k.template operator()<t_k + 1>(t_set_k);
                        }
                    };
                    set_k(set_k);
                    if constexpr (t_j < PotentialTraits<potential>::getNumMappings() - 1)
                    {
                        t_set_j.template operator()<t_j + 1>(t_set_j);
                    }
                };
                set_j(set_j);
                if constexpr (t_i < PotentialTraits<t_potential...>::getNumPotentials() - 1)
                {
                    t_set_i.template operator()<t_i + 1>(t_set_i);
                }
            };
            set_i(set_i);
            return jacobian_matrix;
        }

        template<PotentialConcept auto t_behavior>
        void
        addFormulation()
        {
            if (ptr_formulations_ == nullptr)
            {
                ptr_formulations_ = std::make_unique<std::vector<ElementFormulation<t_domain>>>();
            }
            for (auto const & item : * ptr_formulations_)
            {
                if (item.getLabel() == t_behavior.getLabel())
                {
                    return;
                }
            }
            ptr_formulations_->push_back(ElementFormulation<t_domain>(t_behavior));
        }

        template<PotentialConcept auto t_behavior>
        void
        addFormulation(
            std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
        )
        {
            if (ptr_formulations_ == nullptr)
            {
                ptr_formulations_ = std::make_unique<std::vector<ElementFormulation<t_domain>>>();
            }
            for (auto const & item : * ptr_formulations_)
            {
                if (item.getLabel() == t_behavior.getLabel())
                {
                    return;
                }
            }
            ptr_formulations_->push_back(ElementFormulation<t_domain>(t_behavior, behavior));
        }

        template<PotentialConcept auto t_behavior, Quadrature t_quadrature>
        void
        addFormulation(
            std::shared_ptr<mgis::behaviour::Behaviour> const & behavior
        )
        {
            if (ptr_formulations_ == nullptr)
            {
                ptr_formulations_ = std::make_unique<std::vector<ElementFormulation<t_domain>>>();
            }
            for (auto const & item : * ptr_formulations_)
            {
                if (item.getLabel() == t_behavior.getLabel())
                {
                    return;
                }
            }
            ptr_formulations_->push_back(ElementFormulation<t_domain>(t_behavior, behavior));
            for (auto i = 0; i < QuadratureTraits<t_quadrature>::template Rule<t_element>::getSize(); i++)
            {
                auto point = getCurrentQuadraturePoint<t_quadrature>(i);
                auto r_point = getReferenceQuadraturePoint<t_quadrature>(i);
                auto weight = getCurrentQuadratureWeight<t_quadrature>(i);
                ptr_formulations_->back().addIntegrationPoint(std::move(r_point), std::move(point), weight, behavior);
            }
        }

        template<PotentialConcept auto t_behavior>
        ElementFormulation<t_domain> const &
        getFormulation()
        const
        {
            if (ptr_formulations_ == nullptr)
            {
                throw std::runtime_error("Empty");
            }
            else
            {
                for (auto const & item : * ptr_formulations_)
                {
                    if (item.getLabel() == t_behavior.getLabel())
                    {
                        return item;
                    }
                }
                throw std::runtime_error("No such field data");
            }
        }

        template<PotentialConcept auto t_behavior>
        ElementFormulation<t_domain> &
        getFormulation()
        {
            if (ptr_formulations_ == nullptr)
            {
                throw std::runtime_error("Empty");
            }
            else
            {
                for (auto & item : * ptr_formulations_)
                {
                    if (item.getLabel() == t_behavior.getLabel())
                    {
                        return item;
                    }
                }
                throw std::runtime_error("No such field data");
            }
        }

        template<MappingConcept auto t_mapping>
        auto
        getMapping(
            Point const & point
        )
        const
        {
            return static_cast<t_Disc2<t_mapping.getField()> const *>(this)->template getMapping<t_mapping>(point);
        }

        template<PotentialConcept auto t_behavior, MappingConcept auto t_strain>
        void
        addFormulationStrainOperator()
        {
            auto quadrature_point_count = 0;
            for (auto & ip : this->template getFormulation<t_behavior>().getIntegrationPoints())
            {
                auto strain_operator = this->template getMapping<t_strain>(ip.getReferenceCoordinates());
                ip.template addStrainOperator<t_strain>(strain_operator);
                quadrature_point_count ++;
            }
        }

        template<PotentialConcept auto t_behavior, MappingConcept auto t_strain>
        void
        setFormulationStrain()
        {
            auto const unknown = this->template getUnknownCoefficients<t_strain.getField()>();
            for (auto & ip : this->template getFormulation<t_behavior>().getIntegrationPoints())
            {
                auto rhs = DenseVector<Real>();
                if (ip.template hasStrainOperator<t_strain>())
                {
                    rhs = ip.template getStrainOperator<t_strain>() * unknown;
                }
                else
                {
                    rhs = this->template getMapping<t_strain>(ip.getReferenceCoordinates()) * unknown;
                }
                auto lhs = algebra::View<DenseVector<Real>>(ip.behavior_data_->s1.gradients.data(), ip.behavior_data_->s1.gradients.size());
                lhs = rhs;
            }
        }

        template<PotentialConcept auto... t_behavior>
        void
        integrateConstitutiveEquation(
            std::atomic<Boolean> & output_handler
        )
        {
            (this->template getFormulation<t_behavior>().integrateConstitutiveEquation(output_handler), ...);
        }

        template<PotentialConcept auto... t_behavior>
        void
        reserveBehaviorData()
        {
            (this->template getFormulation<t_behavior>().reserve(), ...);
        }

        template<PotentialConcept auto... t_behavior>
        void
        recoverBehaviorData()
        {
            (this->template getFormulation<t_behavior>().recover(), ...);
        }

        template<PotentialConcept auto t_behavior>
        void
        setMaterialProperty(
            std::basic_string<Character> && label,
            std::function<Real(Point const &)> && function
        )
        {
            this->template getFormulation<t_behavior>().setMaterialProperty(std::forward<std::basic_string<Character>>(label), std::forward<std::function<Real(Point const &)>>(function));
        }

        template<PotentialConcept auto t_behavior>
        void
        setExternalVariable(
            std::basic_string<Character> && label,
            std::function<Real(Point const &)> && function
        )
        {
            this->template getFormulation<t_behavior>().setExternalVariable(std::forward<std::basic_string<Character>>(label), std::forward<std::function<Real(Point const &)>>(function));
        }

        // template<auto t_discretization, BehaviorConcept auto t_behavior, GeneralizedStrainConcept auto t_strain>
        // void
        // addStrainOperators()
        // {
        //     auto strain_operator_num_rows = GeneralizedStrainTraits<t_strain>::template getSize<t_domain>();
        //     auto strain_operator_num_cols = t_Disc<t_discretization>::template getNumElementUnknowns<t_strain.getField()>();
        //     auto quadrature_point_count = 0;
        //     for (auto & ip : getFormulation<t_behavior>().getIntegrationPoints())
        //     {
        //         auto strain_operator = DenseMatrix<Real>(strain_operator_num_rows, strain_operator_num_cols);
        //         strain_operator.setZero();
        //         auto set_mapping_block = [&] <Integer t_i = 0> (
        //             auto & self
        //         )
        //         constexpr mutable
        //         {
        //             auto constexpr mapping = t_strain.template getMapping<t_i>();
        //             auto constexpr mapping_size = GeneralizedStrainTraits<t_strain>::template getMappingSize<t_domain, mapping>();
        //             auto constexpr offset = GeneralizedStrainTraits<t_strain>::template getMappingOffset<t_domain, mapping>();
        //             auto rhs = this->template getMapping<t_strain.getField(), mapping, t_discretization>(ip.getReferenceCoordinates());
        //             auto lhs = strain_operator.block(mapping_size, strain_operator_num_cols, offset, 0);
        //             lhs = rhs;
        //             if constexpr (t_i < t_strain.getNumMappings() - 1)
        //             {
        //                 self.template operator ()<t_i + 1>(self);
        //             }
        //         };
        //         set_mapping_block(set_mapping_block);
        //         ip.template addStrainOperator<t_strain>(strain_operator);
        //         quadrature_point_count ++;
        //     }
        // }

        // template<auto t_discretization, GeneralizedStrainConcept auto t_strain>
        // void
        // addElementOperator()
        // {
        //     getDof<t_strain>().addElementMatrix(this->template getStabilization<t_strain.getField(), t_discretization>());
        // }

        // template<GeneralizedStrainConcept auto t_strain>
        // void
        // setParameter(
        //     Integer tag,
        //     std::function<Real(Point const &)> && function
        // )
        // {
        //     getDof<t_strain>().getRealParameter(tag) = std::forward<std::function<Real(Point const &)>>(function)(* coordinates_);
        // }

        // -----------------------------------------------------------------------------------------------------------------------------------------------------

        template<Basis t_basis>
        static constexpr
        Integer
        getBasisSize()
        {
            return t_Basis<t_basis>::getSize();
        }

        template<FieldConcept auto t_field>
        static constexpr
        Integer
        getFieldSize()
        {
            return FieldTraits<t_field>::template getSize<t_domain>();
        }

        template<FieldConcept auto t_field, Basis t_basis>
        static constexpr
        Integer
        getFieldSize()
        {
            return t_Basis<t_basis>::getSize() * FieldTraits<t_field>::template getSize<t_domain>();
        }
        
        template<Integer t_i, Integer t_j>
        static constexpr
        Integer
        getInnerNeighborNodeConnection(
            Integer component_index,
            Integer node_index
        )
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(ElementTraits<t_element>::node_connectivity_))[component_index][node_index];
        }
        
        Boolean
        operator==(
            FiniteElement const & other
        )
        const = default;
        
        Boolean
        operator!=(
            FiniteElement const & other
        )
        const = default;

        template<Integer t_i>
        Boolean
        isIn(
            std::basic_string<Character> && domain
        )
        const
        requires(t_i == t_element.getDim())
        {
            for (auto const & dom : domains_)
            {
                if (dom->getLabel() == std::forward<std::basic_string<Character>>(domain))
                {
                    return true;
                }
            }
            return false;
        }

        template<Integer t_i>
        Boolean
        isIn(
            std::basic_string<Character> && domain
        )
        const
        requires(t_i > t_element.getDim())
        {
            auto is_in_domain = Boolean(false);
            auto find_domain = [&] <Integer t_j = 0> (
                auto & self
            )
            constexpr mutable
            {
                for (auto const & neighbor : this->template getOuterNeighbors<t_i - t_element.getDim(), t_j>())
                {
                    if (neighbor->template isIn<t_i>(std::forward<std::basic_string<Character>>(domain)))
                    {
                        is_in_domain = true;
                    }
                }
                if constexpr (t_j < t_ElementTraits::template getNumOuterNeighbors<t_domain, t_i - t_element.getDim()>() - 1)
                {
                    self.template operator()<t_j + 1>(self);
                }                
            };
            find_domain(find_domain);
            return is_in_domain;
        }

        template<Integer t_i>
        Boolean
        isIn(
            std::basic_string<Character> && domain
        )
        const
        requires(t_i < t_element.getDim())
        {
            return false;
        }

        Boolean
        isIn(
            std::basic_string<Character> && domain
        )
        const
        {
            auto is_in_domain = Boolean(false);
            auto find_domain = [&] <Integer t_i = t_element.getDim()> (
                auto & self
            )
            constexpr mutable
            {
                if (this->template isIn<t_i>(std::forward<std::basic_string<Character>>(domain)))
                {
                    is_in_domain = true;
                }
                if constexpr (t_i < t_domain.getDim())
                {
                    self.template operator()<t_i + 1>(self);
                }
            };
            find_domain(find_domain);
            return is_in_domain;
        }
        
        std::basic_string<Character>
        getHash()
        const
        {
            auto hash = std::basic_stringstream<Character>();
            for (auto const & node : getInnerNeighbors<t_element.getDim() - 1, 0>())
            {
                hash << node->getHash();
            }
            return hash.str();
        }
        
        std::basic_string<Character>
        getHash()
        const
        requires(t_element.isNode())
        {
            return std::to_string(this->tag_);
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_OuterNeighbors>> &
        getOuterNeighbors()
        {
            return std::get<t_j>(std::get<t_i>(outer_neighbors_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_OuterNeighbors>> const &
        getOuterNeighbors()
        const
        {
            return std::get<t_j>(std::get<t_i>(outer_neighbors_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_InnerNeighbors>> &
        getInnerNeighbors()
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(inner_neighbors_));
        }
        
        template<Integer t_i, Integer t_j>
        std::tuple_element_t<t_j, std::tuple_element_t<t_i, t_InnerNeighbors>> const &
        getInnerNeighbors()
        const
        requires(!t_element.isNode())
        {
            return std::get<t_j>(std::get<t_i>(inner_neighbors_));
        }
        
        // template<Integer t_i, Integer t_j>
        // Integer
        // getInnerNeighborIndex(
        //     Integer i
        // )
        // const
        // requires(!t_element.isNode())
        // {
        //     auto constexpr t_inner_neighbor = t_ElementTraits::template getInnerNeighbor<t_i, t_j>();
        //     auto constexpr t_coordinates = ElementTraits<t_inner_neighbor, t_domain>::template getOuterNeighborCoordinates<t_element>();
        //     auto const & items = getInnerNeighbors<t_i, t_j>()[i]->template getOuterNeighbors<t_coordinates.dim_, t_coordinates.tag_>();
        //     // auto is_equal = [&] (t_ElementPointer<t_element, t_domain> const & ptr_element)
        //     auto is_equal = [&] (std::shared_ptr<FiniteElement<t_element, t_domain>> const & ptr_element)
        //     {
        //         return * ptr_element == * this;
        //     };
        //     return std::distance(items.begin(), std::find_if(items.begin(), items.end(), is_equal));
        // }
        
        // template<Integer t_i, Integer t_j>
        // Integer
        // getInnerNeighborIndex(
        //     std::shared_ptr<FiniteElement<ElementTraits<t_element, t_domain>::template getInnerNeighbor<t_i, t_j>(), t_domain>> const & ptr_neighbor
        // )
        // const
        // requires(!t_element.isNode())
        // {
        //     auto const & inner_neighbors = getInnerNeighbors<t_i, t_j>();
        //     // auto constexpr t_inner_neighbor = ElementTraits<t_element, t_domain>::template getInnerNeighbor<t_i, t_j>();
        //     // auto is_equal = [&] (std::shared_ptr<FiniteElement<t_inner_neighbor, t_domain>> const & neighbor)
        //     // {
        //     //     return * neighbor == * ptr_neighbor;
        //     // };
        //     // auto neighbor_index = std::distance(inner_neighbors.begin(), std::find_if(inner_neighbors.begin(), inner_neighbors.end(), is_equal));
        //     auto neighbor_index = std::distance(inner_neighbors.begin(), std::find(inner_neighbors.begin(), inner_neighbors.end(), ptr_neighbor));
        //     return getInnerNeighborIndex<t_i, t_j>(neighbor_index);
        // }
        
        template<Integer t_i, Integer t_j>
        Integer
        getInnerNeighborSign(
            Integer i
        )
        const
        requires(!t_element.isNode())
        {
            // return getInnerNeighborIndex<t_i, t_j>(i) == 0 ? 1 : -1;
            auto constexpr t_inner_neighbor = t_ElementTraits::template getInnerNeighbor<t_i, t_j>();
            auto constexpr ggg = t_inner_neighbor.getDim() - 1;
            auto constexpr t_num_inner_neighbor_nodes = ElementTraits<t_inner_neighbor>::template getNumInnerNeighbors<t_inner_neighbor.getDim() - 1, 0>();
            // auto constexpr t_inner_neighbor_num_nodes = ElementTraits<t_inner_neighbor, t_domain>::template getNumInnerNeighbors<ggg, 0>();
            auto ori = 1;
            for (auto node_tag = 0; node_tag < t_num_inner_neighbor_nodes; node_tag++)
            {
                // auto mmm = FiniteElement<t_inner_neighbor, t_domain>::template getInnerNeighborNodeConnection<ggg, 0>(node_tag, 0);
                // auto lll = getInnerNeighbors<t_i, t_j>()[i]->getCurrentCoordinates(node_tag);
                // auto kkk = getCurrentCoordinates(getInnerNeighborNodeConnection<t_i, t_j>(i, node_tag));
                // if (!lll.isApprox(kkk))
                // {
                //     ori = -1;
                // }
                // //
                auto lll = getInnerNeighbors<t_i, t_j>()[i]->template getInnerNeighbors<t_inner_neighbor.getDim() - 1, 0>()[node_tag]->getTag();
                auto kkk = getInnerNeighbors<t_element.getDim() - 1, 0>()[getInnerNeighborNodeConnection<t_i, t_j>(i, node_tag)]->getTag();
                if (lll != kkk)
                {
                    ori = -1;
                    break;
                }
            }
            return ori;
            // auto inner_neighbor = getInnerNeighbors<t_i, t_j>()[i];
            // // auto inner_neighbor_reference_centroid = inner_neighbor->getReferenceCentroid();
            // // auto inner_neighbor_current_centroid = inner_neighbor->getCurrentCentroid();
            // // auto current_centroid = this->getCurrentCentroid();
            // auto inner_neighbor_rotation_matrix = inner_neighbor->getRotationMatrix(inner_neighbor->getReferenceCentroid());
            // return (inner_neighbor_rotation_matrix * (inner_neighbor->getCurrentCentroid() - this->getCurrentCentroid()))(t_element.getDim() - 1) > 0 ? 1 : -1;
        }
        
        template<Integer t_i, Integer t_j>
        Integer
        getInnerNeighborOrientation(
            Integer i
        )
        const
        requires(!t_element.isNode())
        {
            // // return getInnerNeighborIndex<t_i, t_j>(i) == 0 ? 1 : -1;
            // auto constexpr t_inner_neighbor = t_ElementTraits::template getInnerNeighbor<t_i, t_j>();
            // auto constexpr ggg = t_inner_neighbor.getDim() - 1;
            // // auto constexpr t_inner_neighbor_num_nodes = ElementTraits<t_inner_neighbor, t_domain>::template getNumInnerNeighbors<ggg, 0>();
            // auto ori = 1;
            // for (auto node_tag = 0; node_tag < ElementTraits<t_inner_neighbor, t_domain>::template getNumInnerNeighbors<ggg, 0>(); node_tag++)
            // {
            //     // auto mmm = FiniteElement<t_inner_neighbor, t_domain>::template getInnerNeighborNodeConnection<ggg, 0>(node_tag, 0);
            //     auto lll = getInnerNeighbors<t_i, t_j>()[i]->getCurrentCoordinates(node_tag);
            //     auto kkk = getCurrentCoordinates(getInnerNeighborNodeConnection<t_i, t_j>(i, node_tag));
            //     if (!lll.isApprox(kkk))
            //     {
            //         ori = -1;
            //     }
            // }
            auto inner_neighbor = getInnerNeighbors<t_i, t_j>()[i];
            // auto inner_neighbor_reference_centroid = inner_neighbor->getReferenceCentroid();
            // auto inner_neighbor_current_centroid = inner_neighbor->getCurrentCentroid();
            // auto current_centroid = this->getCurrentCentroid();
            auto inner_neighbor_rotation_matrix = inner_neighbor->getRotationMatrix(inner_neighbor->getReferenceCentroid());
            return (inner_neighbor_rotation_matrix * (inner_neighbor->getCurrentCentroid() - this->getCurrentCentroid()))(t_element.getDim() - 1) > 0 ? 1 : -1;
        }
        
        // template<Integer t_i, Integer t_j>
        // Integer
        // getInnerNeighborOrientation(
        //     std::shared_ptr<FiniteElement<ElementTraits<t_element, t_domain>::template getInnerNeighbor<t_i, t_j>(), t_domain>> const & ptr_neighbor
        // )
        // const
        // requires(!t_element.isNode())
        // {
        //     return getInnerNeighborIndex<t_i, t_j>(ptr_neighbor) == 0 ? 1 : -1;
        // }
    
        // Point &
        // getCurrentCoordinates()
        // requires(t_element.isNode())
        // {
        //     return this->coordinates_;
        // }
        
        Point const &
        getCurrentCoordinates()
        const
        requires(t_element.isNode())
        {
            return this->coordinates_;
        }
    
        // Point &
        // getCurrentCoordinates(
        //     Integer node_tag
        // )
        // requires(t_element.isNode())
        // {
        //     return this->coordinates_;
        // }
        
        Point const &
        getCurrentCoordinates(
            Integer node_tag
        )
        const
        requires(t_element.isNode())
        {
            return this->coordinates_;
        }
    
        // Point &
        // getCurrentCoordinates(
        //     Integer node_tag
        // )
        // {
        //     return this->template getInnerNeighbors<t_element.getDim() - 1, 0>()[node_tag]->getCurrentCoordinates();
        // }
        
        Point const &
        getCurrentCoordinates(
            Integer node_tag
        )
        const
        {
            return this->template getInnerNeighbors<t_element.getDim() - 1, 0>()[node_tag]->getCurrentCoordinates();
        }
    
        static
        algebra::View<Point const>
        getReferenceCoordinates(
            Integer node_tag
        )
        {
            return algebra::View<Point const>(ElementTraits<t_element>::reference_nodes_[node_tag].data());
        }
        
        DenseMatrix<Real, 3, t_element.getNumNodes()>
        getCurrentCoordinates()
        const
        requires(!t_element.isNode())
        {
            auto current_nodes_coordinates = DenseMatrix<Real, 3, t_element.getNumNodes()>();
            auto count = 0;
            for (auto const & node : this->template getInnerNeighbors<t_element.getDim() - 1, 0>())
            {
                current_nodes_coordinates.col(count) = node->getCurrentCoordinates();
                count ++;
            }
            return current_nodes_coordinates;
        }
        
        static
        lolita::algebra::Span<DenseMatrix<Real, 3, t_element.getNumNodes(), lolita::algebra::colMajor()> const>
        getReferenceCoordinates()
        {
            using t_ReferenceCoordinates = lolita::algebra::Span<DenseMatrix<Real, 3, t_element.getNumNodes(), lolita::algebra::colMajor()> const>;
            return t_ReferenceCoordinates(ElementTraits<t_element>::reference_nodes_.begin()->begin());
        }
        
        static
        Real
        getShapeMappingEvaluation(
            DenseVector<Real, t_element.getNumNodes()> const & nodal_field_values,
            Point const & reference_point
        )
        {
            return t_ElementTraits::getShapeMappingEvaluation(nodal_field_values, reference_point);
        }
        
        static
        Real
        getShapeMappingDerivative(
            DenseVector<Real, t_element.getNumNodes()> const & nodal_field_values,
            Point const & reference_point,
            Integer derivative_direction
        )
        {
            return t_ElementTraits::getShapeMappingDerivative(nodal_field_values, reference_point, derivative_direction);
        }
        
        Real
        getShapeMappingDifferential(
            Point const & point
        )
        const
        requires(t_element.hasDim(3))
        {
            auto const current_coordinates = this->getCurrentCoordinates();
            auto ru = DenseMatrix<Real, 3, 3>();
            auto du = Real(0);
            ru.setZero();
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                for (auto j = 0; j < t_element.getDim(); ++j)
                {
                    ru(i, j) = FiniteElement::getShapeMappingDerivative(current_coordinates.row(i), point, j);
                }
            }
            du = lolita::numerics::abs((ru.col(0).cross(ru.col(1))).dot(ru.col(2)));
            return du;
        }
        
        Real
        getShapeMappingDifferential(
            Point const & point
        )
        const
        requires(t_element.hasDim(2))
        {
            auto const current_coordinates = this->getCurrentCoordinates();
            auto ru = DenseMatrix<Real, 3, 3>();
            auto du = Real(0);
            ru.setZero();
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                for (auto j = 0; j < t_element.getDim(); ++j)
                {
                    ru(i, j) = FiniteElement::getShapeMappingDerivative(current_coordinates.row(i), point, j);
                }
            }
            du = lolita::numerics::abs((ru.col(0).cross(ru.col(1))).norm());
            if constexpr (t_domain.isAxiSymmetric())
            {
                auto r0 = getShapeMappingEvaluation(current_coordinates.row(0), point);
                if (r0 < 1.e-10)
                {
                    r0 = 1.e-10;
                }
                du *= 2.0 * lolita::numerics::pi() * r0;
            }
            return du;
        }
        
        Real
        getShapeMappingDifferential(
            Point const & point
        )
        const
        requires(t_element.hasDim(1))
        {
            auto const current_coordinates = this->getCurrentCoordinates();
            auto ru = DenseMatrix<Real, 3, 3>();
            auto du = Real(0);
            ru.setZero();
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                for (auto j = 0; j < t_element.getDim(); ++j)
                {
                    ru(i, j) = FiniteElement::getShapeMappingDerivative(current_coordinates.row(i), point, j);
                }
            }
            du = lolita::numerics::abs(ru.col(0).norm());
            if constexpr (t_domain.isAxiSymmetric())
            {
                auto r0 = getShapeMappingEvaluation(current_coordinates.row(0), point);
                if (r0 < 1.e-10)
                {
                    r0 = 1.e-10;
                }
                du *= 2.0 * lolita::numerics::pi() * r0;
            }
            return du;
        }
        
        Real
        getShapeMappingDifferential(
            Point const & point
        )
        const
        requires(t_element.isNode())
        {
            return 1.0;
        }
        
        // Real
        // getShapeMappingDifferential(
        //     Point const & point
        // )
        // const
        // {
        //     auto const current_coordinates = this->getCurrentCoordinates();
        //     auto ru = DenseMatrix<Real, 3, 3>();
        //     auto du = Real(0);
        //     ru.setZero();
        //     for (auto i = 0; i < t_domain.dim_; ++i)
        //     {
        //         for (auto j = 0; j < t_element.dim_; ++j)
        //         {
        //             ru(i, j) = FiniteElement::getShapeMappingDerivative(current_coordinates.row(i), point, j);
        //         }
        //     }
        //     if constexpr (t_element.dim_ == 3)
        //     {
        //         du = lolita::numerics::abs((ru.col(0).template cross(ru.col(1))).template dot(ru.col(2)));
        //     }
        //     else if constexpr (t_element.dim_ == 2)
        //     {
        //         du = lolita::numerics::abs((ru.col(0).template cross(ru.col(1))).norm());
        //     }
        //     else
        //     {
        //         du = lolita::numerics::abs(ru.col(0).norm());
        //     }
        //     if constexpr (t_domain.frame_ == Domain::Frame::AxiSymmetric)
        //     {
        //         Real r0 = getShapeMappingEvaluation(current_coordinates.row(0), point);
        //         if (r0 < 1.e-10)
        //         {
        //             r0 = 1.e-10;
        //         }
        //         du *= 2.0 * lolita::numerics::pi() * r0;
        //     }
        //     return du;
        // }
        
        Real
        getRiemannianDistance(
            Point const & first_point,
            Point const & second_point
        )
        const
        requires(t_element.isSub(t_domain, 0))
        {
            auto const & current_coordinates = this->getCurrentCoordinates();
            auto mp0 = Point();
            auto mp1 = Point();
            mp0.setZero();
            mp1.setZero();
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                mp0(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
                mp1(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
            }
            return (mp1 - mp0).norm();
        }
        
        Real
        getRiemannianDistance(
            Point const & first_point,
            Point const & second_point,
            Integer direction
        )
        const
        requires(t_element.isSub(t_domain, 0))
        {
            auto const & current_coordinates = this->getCurrentCoordinates();
            auto mp0 = Point();
            auto mp1 = Point();
            mp0.setZero();
            mp1.setZero();
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                mp0(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
                mp1(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
            }
            return (mp1 - mp0)(direction);
        }
        
        Real
        getRiemannianDistance(
            Point const & first_point,
            Point const & second_point
        )
        const
        requires(!t_element.isSub(t_domain, 0) && t_element.hasDim(2))
        {
            // using SegmentQuadrature = typename QuadratureTraits<Quadrature::gauss(4)>::template Rule<Element::segment(1)>;
            using SegmentQuadrature = typename QuadratureTraits<Quadrature("Gauss", 4)>::template Rule<Element::segment(1)>;
            auto distance = Real(0);
            auto dt = Real();
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            for (auto q = 0; q < SegmentQuadrature::getSize(); ++q)
            {
                auto pq = SegmentQuadrature::reference_points_[q][0];
                auto wq = SegmentQuadrature::reference_weights_[q];
                auto ru = DenseMatrix<Real, 3, 3>();
                auto difference = second_point - first_point;
                auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
                ru.setZero();
                for (auto i = 0; i < t_domain.dim_; ++i)
                {
                    for (auto j = 0; j < t_element.dim_; ++j)
                    {
                        auto du = (1.0 / 2.0) * (second_point(j) - first_point(j));
                        auto dx = FiniteElement::getShapeMappingDerivative(current_nodes_coordinates.row(i), uq, j);
                        ru(i, j) = dx * du;
                    }
                }
                auto Eff = ru.col(0).dot(ru.col(0));
                auto Fff = ru.col(0).dot(ru.col(1));
                auto Gff = ru.col(1).dot(ru.col(1));
                dt = std::sqrt(Eff + 2.0 * Fff + Gff);
                distance += wq * dt;
            }
            return distance;
        }
        
        // DenseMatrix<Real, 2, 2>
        // getMetricTensor(
        //     Point const & first_point,
        //     Point const & second_point,
        //     Integer direction
        // )
        // const
        // requires(!t_element.isSub(t_domain, 0))
        // {
        //     using SegmentQuadrature = QuadratureTraits<Element::segment(1), Quadrature::gauss(4)>;
        //     auto const current_nodes_coordinates = this->getCurrentCoordinates();
        //     for (auto q = 0; q < SegmentQuadrature::getSize(); ++q)
        //     {
        //         auto pq = SegmentQuadrature::reference_points_[q][0];
        //         auto wq = SegmentQuadrature::reference_weights_[q];
        //         auto ru = DenseMatrix<Real, 3, 3>();
        //         auto difference = second_point - first_point;
        //         auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
        //         ru.setZero();
        //         // for (auto i = 0; i < t_domain.dim_; ++i)
        //         // {
        //         for (auto j = 0; j < t_element.dim_; ++j)
        //         {
        //             // if (i == direction)
        //             // {
        //             auto du = (1.0 / 2.0) * (second_point(j) - first_point(j));
        //             auto dx = FiniteElement::getShapeMappingDerivative(current_nodes_coordinates.row(direction), uq, j);
        //             ru(direction, j) = dx * du;
        //             // }
        //         }
        //         // }
        //         auto Eff = ru.col(0).dot(ru.col(0));
        //         auto Fff = ru.col(0).dot(ru.col(1));
        //         auto Gff = ru.col(1).dot(ru.col(1));
        //         dt = std::sqrt(Eff + 2.0 * Fff + Gff);
        //         distance += wq * dt;
        //     }
        //     return sign * distance;
        // }
        
        Real
        getRiemannianDistance(
            Point const & first_point,
            Point const & second_point,
            Integer direction
        )
        const
        requires(!t_element.isSub(t_domain, 0) && t_element.hasDim(2))
        {
            // using SegmentQuadrature = typename QuadratureTraits<Quadrature::gauss(4)>::template Rule<Element::segment(1)>;
            using SegmentQuadrature = typename QuadratureTraits<Quadrature("Gauss", 4)>::template Rule<Element::segment(1)>;
            // using SegmentQuadrature = typename QuadratureTraits<Quadrature(4, Quadrature::Gauss)>::template Rule<Element::segment(1)>;
            auto distance = Real(0);
            auto dt = Real();
            // -> TEST
            auto sign = (second_point - first_point)(direction) > 0 ? 1 : -1;
            // <_ TEST
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            for (auto q = 0; q < SegmentQuadrature::getSize(); ++q)
            {
                auto pq = SegmentQuadrature::reference_points_[q][0];
                auto wq = SegmentQuadrature::reference_weights_[q];
                auto ru = DenseMatrix<Real, 3, 3>();
                auto difference = second_point - first_point;
                auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
                ru.setZero();
                // for (auto i = 0; i < t_domain.dim_; ++i)
                // {
                for (auto j = 0; j < t_element.dim_; ++j)
                {
                    // if (i == direction)
                    // {
                    auto du = (1.0 / 2.0) * (second_point(j) - first_point(j));
                    auto dx = FiniteElement::getShapeMappingDerivative(current_nodes_coordinates.row(direction), uq, j);
                    ru(direction, j) = dx * du;
                    // }
                }
                // }
                auto Eff = ru.col(0).dot(ru.col(0));
                auto Fff = ru.col(0).dot(ru.col(1));
                auto Gff = ru.col(1).dot(ru.col(1));
                dt = std::sqrt(Eff + 2.0 * Fff + Gff);
                distance += wq * dt;
            }
            return sign * distance;
        }
        
        Real
        getRiemannianDistance(
            Point const & first_point,
            Point const & second_point
        )
        const
        requires(!t_element.isSub(t_domain, 0) && t_element.hasDim(1))
        {
            // using SegmentQuadrature = typename QuadratureTraits<Quadrature::gauss(4)>::template Rule<Element::segment(1)>;
            using SegmentQuadrature = typename QuadratureTraits<Quadrature("Gauss", 4)>::template Rule<Element::segment(1)>;
            // using SegmentQuadrature = typename QuadratureTraits<Quadrature(4, Quadrature::Gauss)>::template Rule<Element::segment(1)>;
            auto distance = Real(0);
            auto dt = Real();
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            for (auto q = 0; q < SegmentQuadrature::getSize(); ++q)
            {
                auto pq = SegmentQuadrature::reference_points_[q][0];
                auto wq = SegmentQuadrature::reference_weights_[q];
                auto ru = DenseMatrix<Real, 3, 3>();
                auto difference = second_point - first_point;
                auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
                ru.setZero();
                for (auto i = 0; i < t_domain.dim_; ++i)
                {
                    for (auto j = 0; j < t_element.dim_; ++j)
                    {
                        auto du = (1.0 / 2.0) * (second_point(j) - first_point(j));
                        auto dx = FiniteElement::getShapeMappingDerivative(current_nodes_coordinates.row(i), uq, j);
                        ru(i, j) = dx * du;
                    }
                }
                auto Eff = ru.col(0).dot(ru.col(0));
                dt = std::sqrt(Eff);
                distance += wq * dt;
            }
            return distance;
        }
        
        Real
        getRiemannianDistance(
            Point const & first_point,
            Point const & second_point,
            Integer direction
        )
        const
        requires(!t_element.isSub(t_domain, 0) && t_element.hasDim(1))
        {
            // using SegmentQuadrature = typename QuadratureTraits<Quadrature::gauss(4)>::template Rule<Element::segment(1)>;
            using SegmentQuadrature = typename QuadratureTraits<Quadrature("Gauss", 4)>::template Rule<Element::segment(1)>;
            // using SegmentQuadrature = typename QuadratureTraits<Quadrature(4, Quadrature::Gauss)>::template Rule<Element::segment(1)>;
            auto distance = Real(0);
            auto dt = Real();
            // -> TEST
            auto sign = (second_point - first_point)(direction) > 0 ? 1 : -1;
            // <_ TEST
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            for (auto q = 0; q < SegmentQuadrature::getSize(); ++q)
            {
                auto pq = SegmentQuadrature::reference_points_[q][0];
                auto wq = SegmentQuadrature::reference_weights_[q];
                auto ru = DenseMatrix<Real, 3, 3>();
                auto difference = second_point - first_point;
                auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
                ru.setZero();
                // for (auto i = 0; i < t_domain.dim_; ++i)
                // {
                for (auto j = 0; j < t_element.dim_; ++j)
                {
                    // if (i == direction)
                    // {
                    auto du = (1.0 / 2.0) * (second_point(j) - first_point(j));
                    auto dx = FiniteElement::getShapeMappingDerivative(current_nodes_coordinates.row(direction), uq, j);
                    ru(direction, j) = dx * du;
                    // }
                }
                // }
                auto Eff = ru.col(0).dot(ru.col(0));
                dt = std::sqrt(Eff);
                distance += wq * dt;
            }
            return sign * distance;
        }
        
        Real
        getRiemannianDistance(
            Point const & first_point,
            Point const & second_point
        )
        const
        requires(!t_element.isSub(t_domain, 0) && t_element.isNode())
        {
            return 0.0;
        }
        
        Real
        getRiemannianDistance(
            Point const & first_point,
            Point const & second_point,
            Integer direction
        )
        const
        requires(!t_element.isSub(t_domain, 0) && t_element.isNode())
        {
            return 0.0;
        }
        
        // Real
        // getRiemannianDistance(
        //     Point const & first_point,
        //     Point const & second_point,
        //     Integer direction = -1
        // )
        // const
        // {
        //     if constexpr (t_element.isSub(t_domain, 0))
        //     {
        //         auto const & current_coordinates = this->getCurrentCoordinates();
        //         auto distance = Real();
        //         auto mp0 = Point();
        //         auto mp1 = Point();
        //         mp0.setZero();
        //         mp1.setZero();
        //         for (auto i = 0; i < t_domain.getDim(); ++i)
        //         {
        //             mp0(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
        //             mp1(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
        //         }
        //         direction == -1 ? distance = (mp1 - mp0).norm() : distance = (mp1 - mp0)(direction);
        //         return distance;
        //     }
        //     else
        //     {
        //         using SegmentQuadrature = QuadratureTraits<Element::segment(1), Quadrature::gauss(4)>;
        //         auto distance = Real(0);
        //         auto dt = Real();
        //         auto const current_nodes_coordinates = this->getCurrentCoordinates();
        //         for (auto q = 0; q < SegmentQuadrature::getSize(); ++q)
        //         {
        //             auto pq = SegmentQuadrature::reference_points_[q][0];
        //             auto wq = SegmentQuadrature::reference_weights_[q];
        //             auto ru = DenseMatrix<Real, 3, 3>();
        //             auto difference = second_point - first_point;
        //             auto uq = (1.0 / 2.0) * difference * pq + (1.0 / 2.0) * difference;
        //             ru.setZero();
        //             for (auto i = 0; i < t_domain.dim_; ++i)
        //             {
        //                 for (auto j = 0; j < t_element.dim_; ++j)
        //                 {
        //                     if (direction == -1 || i == direction)
        //                     {
        //                         auto du = (1.0 / 2.0) * (second_point(j) - first_point(j));
        //                         auto dx = FiniteElement::getShapeMappingDerivative(current_nodes_coordinates.row(i), uq, j);
        //                         ru(i, j) = dx * du;
        //                     }
        //                 }
        //             }
        //             if constexpr (t_element.hasDim(2))
        //             {
        //                 auto Eff = ru.col(0).template dot(ru.col(0));
        //                 auto Fff = ru.col(0).template dot(ru.col(1));
        //                 auto Gff = ru.col(1).template dot(ru.col(1));
        //                 dt = std::sqrt(Eff + 2.0 * Fff + Gff);
        //             }
        //             else if constexpr (t_element.hasDim(1))
        //             {
        //                 auto Eff = ru.col(0).template dot(ru.col(0));
        //                 dt = std::sqrt(Eff);
        //             }
        //             else
        //             {
        //                 dt = 0;
        //             }
        //             distance += wq * dt;
        //         }
        //         return distance;
        //     }
        // }
        
        Real
        getLocalFrameDistance(
            Point const & first_point,
            Point const & second_point,
            Integer kkk
        )
        const
        requires(t_element.isSub(t_domain, 0))
        {
            auto first_point_mapping = Point();
            auto second_point_mapping = Point();
            auto const & current_coordinates = this->getCurrentCoordinates();
            first_point_mapping.setZero();
            second_point_mapping.setZero();
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                first_point_mapping(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
                second_point_mapping(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
            }
            return (second_point_mapping - first_point_mapping)(kkk);
            // return (second_point - first_point)(kkk);
        }
        
        Real
        getLocalFrameDistance(
            Point const & first_point,
            Point const & second_point,
            Integer kkk
        )
        const
        requires(t_element.isSub(t_domain, 1))
        {
            auto rotation_matrix = getRotationMatrix(getReferenceCentroid());
            auto first_point_mapping = Point();
            auto second_point_mapping = Point();
            auto const & current_coordinates = this->getCurrentCoordinates();
            first_point_mapping.setZero();
            second_point_mapping.setZero();
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                first_point_mapping(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
                second_point_mapping(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
            }
            return (rotation_matrix * (second_point_mapping - first_point_mapping))(kkk);
            // return (rot * (second_point - first_point))(kkk);
        }
        
        // Real
        // getLocalFrameDistance(
        //     Point const & first_point,
        //     Point const & second_point,
        //     Integer kkk
        // )
        // const
        // {
        //     auto mp_0 = Point();
        //     auto mp_1 = Point();
        //     auto const & current_coordinates = this->getCurrentCoordinates();
        //     mp_0.setZero();
        //     mp_1.setZero();
        //     for (auto i = 0; i < t_domain.getDim(); ++i)
        //     {
        //         mp_0(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), first_point);
        //         mp_1(i) = FiniteElement::getShapeMappingEvaluation(current_coordinates.row(i), second_point);
        //     }
        //     if constexpr (t_element.isSub(t_domain, 0))
        //     {
        //         return (mp_1 - mp_0)(kkk);
        //     }
        //     else if constexpr (t_element.isSub(t_domain, 1))
        //     {
        //         auto rot = getRotationMatrix(getReferenceCentroid());
        //         return (rot * (mp_1 - mp_0))(kkk);
        //     }
        // }
        
        Point
        getLocalFrameDiameters()
        const
        requires(t_element.isSub(t_domain, 0))
        {
            return getCurrentDiameters();
        }
        
        Point
        getLocalFrameDiameters()
        const
        requires(t_element.isSub(t_domain, 1))
        {
            auto rotation_matrix = getRotationMatrix(getReferenceCentroid());
            auto coordinates = getCurrentCoordinates();
            auto projected_coordinates = rotation_matrix * coordinates;
            auto current_diameters = Point();
            current_diameters.setZero();
            for (auto i = 0; i < t_element.getNumNodes(); ++i)
            {
                for (auto j = i + 1; j < t_element.getNumNodes(); ++j)
                {
                    for (auto k = 0; k < t_element.getDim(); ++k)
                    {
                        auto new_value = projected_coordinates(k, j) - projected_coordinates(k, i);
                        auto & current_value = current_diameters(k);
                        if (new_value > current_value)
                        {
                            current_value = new_value;
                        }
                    }
                }
            }
            return current_diameters;
            // auto reference_coordinates = FiniteElement::getReferenceCoordinates();
            // auto current_diameters = Point();
            // current_diameters.setZero();
            // for (auto i = 0; i < t_element.getNumNodes(); ++i)
            // {
            //     for (auto j = i + 1; j < t_element.getNumNodes(); ++j)
            //     {
            //         auto const & pt0 = reference_coordinates.col(i);
            //         auto const & pt1 = reference_coordinates.col(j);
            //         for (auto k = 0; k < t_element.getDim(); ++k)
            //         {
            //             auto new_value = lolita::numerics::abs(getRiemannianDistance(pt0, pt1, k));
            //             auto & current_value = current_diameters(k);
            //             if (new_value > current_value)
            //             {
            //                 current_value = new_value;
            //             }
            //         }
            //     }
            // }
            // return current_diameters;
        }
        
        // Point
        // getLocalFrameDiameters()
        // const
        // {
        //     if constexpr (t_element.isSub(t_domain, 0))
        //     {
        //         return getCurrentDiameters();
        //     }
        //     else if constexpr (t_element.isSub(t_domain, 1))
        //     {
        //         auto proj_v = getRotationMatrix(getReferenceCentroid()) * getCurrentCoordinates();
        //         auto current_diameters = Point();
        //         current_diameters.setZero();
        //         for (auto i = 0; i < t_element.getNumNodes(); ++i)
        //         {
        //             for (auto j = i + 1; j < t_element.getNumNodes(); ++j)
        //             {
        //                 auto pt0 = proj_v.col(i);
        //                 auto pt1 = proj_v.col(j);
        //                 for (auto k = 0; k < t_element.getDim(); ++k)
        //                 {
        //                     auto new_value = (pt1 - pt0)(k);
        //                     auto & current_value = current_diameters(k);
        //                     if (new_value > current_value)
        //                     {
        //                         current_value = new_value;
        //                     }
        //                 }
        //             }
        //         }
        //         return current_diameters;
        //     }
        // }
        
        Point
        getCurrentDiameters()
        const
        {
            auto reference_coordinates = FiniteElement::getReferenceCoordinates();
            auto current_diameters = Point();
            current_diameters.setZero();
            for (auto i = 0; i < t_element.getNumNodes(); ++i)
            {
                for (auto j = i + 1; j < t_element.getNumNodes(); ++j)
                {
                    auto const & pt0 = reference_coordinates.col(i);
                    auto const & pt1 = reference_coordinates.col(j);
                    for (auto k = 0; k < 3; ++k)
                    {
                        auto new_value = lolita::numerics::abs(getRiemannianDistance(pt0, pt1, k));
                        auto & current_value = current_diameters(k);
                        if (new_value > current_value)
                        {
                            current_value = new_value;
                        }
                    }
                }
            }
            return current_diameters;
        }
        
        static
        Point
        getReferenceCentroid()
        {
            auto reference_nodes_coordinates = FiniteElement::getReferenceCoordinates();
            auto barycenter = Point();
            barycenter.setZero();
            for (auto i = 0; i < t_element.getNumNodes(); ++i)
            {
                barycenter += reference_nodes_coordinates.col(i);
            }
            barycenter /= Real(t_element.getNumNodes());
            return barycenter;
        }
        
        Point
        getCurrentCentroid()
        const
        {
            // auto const current_nodes_coordinates = this->getCurrentCoordinates();
            // auto barycenter = Point();
            // barycenter.setZero();
            // for (auto i = 0; i < t_element.getNumNodes(); ++i) {
            //     barycenter += current_nodes_coordinates.col(i);
            // }
            // barycenter /= Real(t_element.getNumNodes());
            // return barycenter;
            auto const reference_centroid = FiniteElement::getReferenceCentroid();
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto barycenter = Point();
            barycenter.setZero();
            for (auto i = 0; i < t_domain.getDim(); i++)
            {
                barycenter(i) = FiniteElement::getShapeMappingEvaluation(current_nodes_coordinates.row(i), reference_centroid);
            }
            return barycenter;
        }
        
        static
        Point
        getReferenceDiameters()
        {
            auto dts = Point();
            auto nds = FiniteElement::getReferenceCoordinates();
            dts.setZero();
            for (auto i = 0; i < t_element.getNumNodes(); ++i)
            {
                for (auto j = i + 1; j < t_element.getNumNodes(); ++j)
                {
                    for (auto k = 0; k < 3; ++k)
                    {
                        auto new_value = lolita::numerics::abs(nds(k, i) - nds(k, j));
                        auto & current_value = dts(k);
                        if (new_value > current_value)
                        {
                            current_value = new_value;
                        }
                    }
                }
            }
            return dts;
        }
        
        // Point
        // getNormalVector(
        //     Point const & point
        // )
        // const
        // requires(t_element.isSub(t_domain, 1))
        // {
        //     auto const current_nodes_coordinates = this->getCurrentCoordinates();
        //     auto ru = DenseMatrix<Real, 3, t_element.getDim()>();
        //     ru.setZero();
        //     for (auto i = 0; i < t_domain.getDim(); ++i)
        //     {
        //         for (auto j = 0; j < t_element.getDim(); ++j)
        //         {
        //             ru(i, j) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, j);
        //         }
        //     }
        //     if constexpr (t_element.isNode()) {
        //         return Point{0, 0, 0};
        //     }
        //     else if constexpr (t_element.hasDim(1))
        //     {
        //         return Point{ru(1)/ru.norm(), -ru(0)/ru.norm(), 0};
        //     }
        //     else
        //     {
        //         return (ru.col(0) / ru.col(0).norm()).cross((ru.col(1) / ru.col(1).norm()));
        //     }
        // }
        
        Point
        getNormalVector(
            Point const & point
        )
        const
        requires(t_element.isSub(t_domain, 1) && t_domain.hasDim(2))
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto ru = DenseMatrix<Real, 3, t_element.getDim()>();
            auto normal_vector = Point();
            ru.setZero();
            normal_vector.setZero();
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                for (auto j = 0; j < t_element.getDim(); ++j)
                {
                    ru(i, j) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, j);
                }
            }
            normal_vector(0) = + ru(1)/ru.norm();
            normal_vector(1) = - ru(0)/ru.norm();
            return normal_vector;
        }
        
        Point
        getNormalVector(
            Point const & point
        )
        const
        requires(t_element.isSub(t_domain, 1) && t_domain.hasDim(3))
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto ru = DenseMatrix<Real, 3, t_element.getDim()>();
            auto normal_vector = Point();
            ru.setZero();
            normal_vector.setZero();
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                for (auto j = 0; j < t_element.getDim(); ++j)
                {
                    ru(i, j) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, j);
                }
            }
            // normal_vector = (ru.col(0) / ru.col(0).norm()).cross((ru.col(1) / ru.col(1).norm()));
            normal_vector = ru.col(0).cross(ru.col(1));
            normal_vector.normalize();
            return normal_vector;
        }

        DenseMatrix<Real, 3, 3>
        getRotationMatrix(
            Point const & point
        )
        const
        requires(t_element.isSub(t_domain, 1) && t_domain.hasDim(3))
        {
            // getting the element centroid and the first node coordinates
            auto centroid = getCurrentCentroid();
            auto node_coordinates = getCurrentCoordinates(0);
            // setting basis vectors
            auto basis_vector_0 = getNormalVector(point);
            auto basis_vector_1 = (node_coordinates - centroid) / (node_coordinates - centroid).norm();
            auto basis_vector_2 = basis_vector_0.cross(basis_vector_1);
            // setting rotation matrix
            auto rotation_matrix = DenseMatrix<Real, 3, 3>();
            rotation_matrix.setZero();
            rotation_matrix(0, 0) = + basis_vector_0(0);
            rotation_matrix(0, 1) = + basis_vector_0(1);
            rotation_matrix(0, 2) = + basis_vector_0(2);
            rotation_matrix(1, 0) = + basis_vector_1(0);
            rotation_matrix(1, 1) = + basis_vector_1(1);
            rotation_matrix(1, 2) = + basis_vector_1(2);
            rotation_matrix(2, 0) = + basis_vector_2(0);
            rotation_matrix(2, 1) = + basis_vector_2(1);
            rotation_matrix(2, 2) = + basis_vector_2(2);
            return rotation_matrix;
        }

        DenseMatrix<Real, 3, 3>
        getRotationMatrix(
            Point const & point
        )
        const
        requires(t_element.isSub(t_domain, 1) && t_domain.hasDim(2))
        {
            // setting the only basis vector needed
            auto basis_vector_0 = getNormalVector(point);
            // setting rotation matrix
            auto rotation_matrix = DenseMatrix<Real, 3, 3>();
            rotation_matrix.setZero();
            rotation_matrix(0, 0) = - basis_vector_0(1);
            rotation_matrix(0, 1) = + basis_vector_0(0);
            rotation_matrix(1, 0) = + basis_vector_0(0);
            rotation_matrix(1, 1) = + basis_vector_0(1);
            return rotation_matrix;
        }

        DenseMatrix<Real, 3, 3>
        getRotationMatrix(
            Point const & point
        )
        const
        requires(t_element.isSub(t_domain, 1) && t_domain.hasDim(1))
        {
            // setting rotation matrix
            auto rotation_matrix = DenseMatrix<Real, 3, 3>();
            rotation_matrix.setZero();
            rotation_matrix(0, 0) = 1.0;
            return rotation_matrix;
        }
        
        Point
        getTangentVector(
            Point const & point,
            Integer direction
        )
        const
        requires(t_element.hasDim(1))
        {
            auto const current_nodes_coordinates = this->getCurrentCoordinates();
            auto tangent_vector = Point();
            tangent_vector.setZero();
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                tangent_vector(i) = this->getShapeMappingDerivative(current_nodes_coordinates.row(i), point, direction);
            }
            return tangent_vector;
        }
        
        template<Quadrature t_quadrature>
        static
        Real
        getReferenceQuadratureWeight(
            Integer index
        )
        {
            return QuadratureTraits<t_quadrature>::template Rule<t_element>::reference_weights_[index];
        }
        
        template<Quadrature t_quadrature>
        static
        lolita::algebra::Span<Point const>
        getReferenceQuadraturePoint(
            Integer index
        )
        {
            return lolita::algebra::Span<Point const>(QuadratureTraits<t_quadrature>::template Rule<t_element>::reference_points_[index].begin());
        }
        
        template<Quadrature t_quadrature>
        Real
        getCurrentQuadratureWeight(
            Integer index
        )
        const
        {
            auto w = getReferenceQuadratureWeight<t_quadrature>(index);
            return w * getShapeMappingDifferential(getReferenceQuadraturePoint<t_quadrature>(index));
        }
        
        template<Quadrature t_quadrature>
        Point
        getCurrentQuadraturePoint(
            Integer index
        )
        const
        {
            auto p = Point();
            auto const nds = this->getCurrentCoordinates();
            p.setZero();
            for (auto j = 0; j < t_domain.getDim(); ++j)
            {
                p(j) = FiniteElement::getShapeMappingEvaluation(nds.row(j), getReferenceQuadraturePoint<t_quadrature>(index));
            }
            return p;
        }
        
        template<Quadrature t_quadrature, Integer _i, Integer _j>
        static
        Real
        getInnerNeighborReferenceQuadratureWeight(
            Integer index
        )
        requires(!t_element.isNode())
        {
            auto const constexpr t_component = t_ElementTraits::template getInnerNeighbor<_i, _j>();
            return FiniteElement<t_component, t_domain>::template getReferenceQuadratureWeight<t_quadrature>(index);
        }
        
        template<Quadrature t_quadrature, Integer _i, Integer _j>
        static
        Point
        getInnerNeighborReferenceQuadraturePoint(
            Integer component_index,
            Integer index
        )
        requires(!t_element.isNode())
        {
            auto constexpr t_component = t_ElementTraits ::template getInnerNeighbor<_i, _j>();
            auto p = Point();
            p.setZero();
            auto const & elt_reference_nodes = ElementTraits<t_element>::reference_nodes_;
            for (auto i = 0; i < t_domain.getDim(); ++i)
            {
                auto cpt_coordinates = DenseVector<Real, t_component.getNumNodes()>();
                for (auto j = 0; j < t_component.getNumNodes(); ++j)
                {
                    auto const node_tag = getInnerNeighborNodeConnection<_i, _j>(component_index, j);//.get(component_index).get(j);
                    cpt_coordinates(j) = elt_reference_nodes[node_tag][i];
                }
                auto cpt_reference_point = FiniteElement<t_component, t_domain>::template getReferenceQuadraturePoint<t_quadrature>(index);
                p(i) = FiniteElement<t_component, t_domain>::getShapeMappingEvaluation(cpt_coordinates, cpt_reference_point);
            }
            return p;
        }
        
        template<Quadrature t_quadrature, Integer t_i, Integer t_j>
        Real
        getInnerNeighborCurrentQuadratureWeight(
            Integer component_index,
            Integer index
        )
        const
        requires(!t_element.isNode())
        {
            auto const & cmp =  this->template getInnerNeighbors<t_i, t_j>()[component_index];//.template get<I>().template get<J>().get(component_index).get();
            return cmp->template getCurrentQuadratureWeight<t_quadrature>(index);
        }
        
        template<Quadrature t_quadrature, Integer t_i, Integer t_j>
        Point
        getInnerNeighborCurrentQuadraturePoint(
            Integer component_index,
            Integer index
        )
        const
        requires(!t_element.isNode())
        {
            auto p = Point();
            auto const cpt_ref_pnt = getInnerNeighborReferenceQuadraturePoint<t_quadrature, t_i, t_j>(component_index, index);
            auto const nds = this->getCurrentCoordinates();
            p.setZero();
            for (auto j = 0; j < t_domain.getDim(); ++j)
            {
                p(j) = FiniteElement::getShapeMappingEvaluation(nds.row(j), cpt_ref_pnt);
            }
            return p;
        }

    };   
    
} // namespace lolita

#endif /* B9B48BEA_09F5_41DB_84A1_45A6E708901C */
