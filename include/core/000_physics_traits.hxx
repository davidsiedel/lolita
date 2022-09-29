#ifndef DA9C5D5A_13CE_4129_A31C_3550B18DAB24
#define DA9C5D5A_13CE_4129_A31C_3550B18DAB24

#include "lolita.hxx"

namespace lolita
{
        
    template<Basis t_basis>
    struct BasisTraits;

    template<auto t_discretization>
    struct DiscretizationTraits;

    template<DiscreteFieldConcept auto t_field>
    struct DiscretizationTraits2;

    template<FieldConcept auto t_field>
    struct FieldTraits
    {

        template<Domain t_domain>
        static constexpr
        Integer
        getCols()
        requires(t_field.isTensor(0))
        {
            return lolita::numerics::pow(t_domain.getDim(), 0);
        }

        template<Domain t_domain>
        static constexpr
        Integer
        getCols()
        requires(t_field.isTensor(1))
        {
            return lolita::numerics::pow(t_domain.getDim(), 0);
        }

        template<Domain t_domain>
        static constexpr
        Integer
        getCols()
        requires(t_field.isTensor(2))
        {
            return lolita::numerics::pow(t_domain.getDim(), 1);
        }

        template<Domain t_domain>
        static constexpr
        Integer
        getCols()
        requires(t_field.isTensor(3))
        {
            return lolita::numerics::pow(t_domain.getDim(), 1);
        }

        template<Domain t_domain>
        static constexpr
        Integer
        getCols()
        requires(t_field.isTensor(4))
        {
            return lolita::numerics::pow(t_domain.getDim(), 2);
        }

        // ---

        template<Domain t_domain>
        static constexpr
        Integer
        getRows()
        requires(t_field.isTensor(0))
        {
            return lolita::numerics::pow(t_domain.getDim(), 0);
        }

        template<Domain t_domain>
        static constexpr
        Integer
        getRows()
        requires(t_field.isTensor(1))
        {
            return lolita::numerics::pow(t_domain.getDim(), 1);
        }

        template<Domain t_domain>
        static constexpr
        Integer
        getRows()
        requires(t_field.isTensor(2))
        {
            return lolita::numerics::pow(t_domain.getDim(), 1);
        }

        template<Domain t_domain>
        static constexpr
        Integer
        getRows()
        requires(t_field.isTensor(3))
        {
            return lolita::numerics::pow(t_domain.getDim(), 2);
        }

        template<Domain t_domain>
        static constexpr
        Integer
        getRows()
        requires(t_field.isTensor(4))
        {
            return lolita::numerics::pow(t_domain.getDim(), 2);
        }

        // ---

        template<Domain t_domain>
        static constexpr
        Integer
        getSize()
        {
            return getRows<t_domain>() * getCols<t_domain>();
        }

    };

    struct MappingValues
    {

        constexpr
        MappingValues(
            Integer row,
            Integer col,
            Integer position,
            Real coefficient
        )
        :
        row_(row),
        col_(col),
        position_(position),
        coefficient_(coefficient)
        {}

        constexpr
        Integer
        row()
        const
        {
            return row_;
        }

        constexpr
        Integer
        col()
        const
        {
            return col_;
        }

        constexpr
        Integer
        rank()
        const
        {
            return position_;
        }

        constexpr
        Real
        value()
        const
        {
            return coefficient_;
        }

        Integer row_;

        Integer col_;

        Integer position_;

        Real coefficient_;

    };

    template<MappingConcept auto t_mapping>
    struct MappingTraits;

    template<MappingConcept auto t_mapping>
    requires(t_mapping.isIdentity())
    struct MappingTraits<t_mapping>
    {

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getRows()
        {
            return FieldTraits<t_field>::template getRows<t_domain>();
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getCols()
        {
            return FieldTraits<t_field>::template getCols<t_domain>();
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getSize()
        {
            return getRows<t_domain, t_field>() * getCols<t_domain, t_field>();
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        std::array<MappingValues, getSize<t_domain, t_field>()>
        getValues()
        requires(t_field.isTensor(0))
        {
            return {
                MappingValues{0, 0, 0, 1},
            };
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        std::array<MappingValues, getSize<t_domain, t_field>()>
        getValues()
        requires(t_field.isTensor(1) && t_domain.hasDim(2))
        {
            return {
                MappingValues{0, 0, 0, 1},
                MappingValues{0, 1, 1, 1},
            };
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        std::array<MappingValues, getSize<t_domain, t_field>()>
        getValues()
        requires(t_field.isTensor(1) && t_domain.hasDim(3))
        {
            return {
                MappingValues{0, 0, 0, 1},
                MappingValues{0, 1, 1, 1},
                MappingValues{0, 2, 2, 1},
            };
        }
        
        static
        void
        non_linear(
                auto & gradient
        )
        {}
        
        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static
        std::array<Real, getSize<t_domain, t_field>()>
        getStressValues(
            auto const & thermodynamic_forces,
            Integer offset
        )
        {
            auto stress = std::array<Real, getSize<t_domain, t_field>()>();
            for (auto i = 0; i < getSize<t_domain, t_field>(); i++)
            {
                stress[i] = thermodynamic_forces[i + offset];
            }
            return stress;
        }
        
        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static
        std::array<Real, getSize<t_domain, t_field>()>
        getStrainValues(
            auto const & gradients,
            Integer offset
        )
        {
            auto strain = std::array<Real, getSize<t_domain, t_field>()>();
            for (auto i = 0; i < getSize<t_domain, t_field>(); i++)
            {
                strain[i] = gradients[i + offset];
            }
            return strain;
        }
        
    };

    template<MappingConcept auto t_mapping>
    requires(t_mapping.isGradient())
    struct MappingTraits<t_mapping>
    {

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getRows()
        {
            return FieldTraits<Field(t_field.getDim() + 1)>::template getRows<t_domain>();
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getCols()
        {
            return FieldTraits<Field(t_field.getDim() + 1)>::template getCols<t_domain>();
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getSize()
        {
            return getRows<t_domain, t_field>() * getCols<t_domain, t_field>();
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        std::array<MappingValues, getSize<t_domain, t_field>()>
        getValues()
        requires(t_field.isTensor(0) && t_domain.hasDim(1))
        {
            return {
                MappingValues{0, 0, 0, 1},
            };
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        std::array<MappingValues, getSize<t_domain, t_field>()>
        getValues()
        requires(t_field.isTensor(0) && t_domain.hasDim(2))
        {
            return {
                MappingValues{0, 0, 0, 1},
                MappingValues{0, 1, 1, 1},
            };
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        std::array<MappingValues, getSize<t_domain, t_field>()>
        getValues()
        requires(t_field.isTensor(0) && t_domain.hasDim(3))
        {
            return {
                MappingValues{0, 0, 0, 1},
                MappingValues{0, 1, 1, 1},
                MappingValues{0, 2, 2, 1},
            };
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        std::array<MappingValues, getSize<t_domain, t_field>()>
        getValues()
        requires(t_field.isTensor(1) && t_domain.hasDim(2))
        {
            return {
                MappingValues{0, 0, 0, 1},
                MappingValues{0, 1, 1, 1},
                MappingValues{1, 0, 2, 1},
                MappingValues{1, 1, 3, 1},
            };
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        std::array<MappingValues, getSize<t_domain, t_field>()>
        getValues()
        requires(t_field.isTensor(1) && t_domain.hasDim(3))
        {
            return {
                MappingValues{0, 0, 0, 1},
                MappingValues{0, 1, 1, 1},
                MappingValues{0, 2, 2, 1},
                MappingValues{1, 0, 3, 1},
                MappingValues{1, 1, 4, 1},
                MappingValues{1, 2, 5, 1},
                MappingValues{2, 0, 6, 1},
                MappingValues{2, 1, 7, 1},
                MappingValues{2, 2, 8, 1},
            };
        }
        
        static
        void
        non_linear(
                auto & gradient
        )
        {}
        
        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static
        std::array<Real, getSize<t_domain, t_field>()>
        getStressValues(
            auto const & thermodynamic_forces,
            Integer offset
        )
        {
            auto stress = std::array<Real, getSize<t_domain, t_field>()>();
            for (auto i = 0; i < getSize<t_domain, t_field>(); i++)
            {
                stress[i] = thermodynamic_forces[i + offset];
            }
            return stress;
        }
        
        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static
        std::array<Real, getSize<t_domain, t_field>()>
        getStrainValues(
            auto const & gradients,
            Integer offset
        )
        {
            auto strain = std::array<Real, getSize<t_domain, t_field>()>();
            for (auto i = 0; i < getSize<t_domain, t_field>(); i++)
            {
                strain[i] = gradients[i + offset];
            }
            return strain;
        }
        
    };

    template<MappingConcept auto t_mapping>
    requires(t_mapping.isSmallStrain())
    struct MappingTraits<t_mapping>
    {

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getRows()
        requires(t_field.isTensor(1) && t_domain.hasDim(2))
        {
            return 4;
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getRows()
        requires(t_field.isTensor(1) && t_domain.hasDim(3))
        {
            return 6;
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getCols()
        {
            return 1;
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getSize()
        {
            return getRows<t_domain, t_field>() * getCols<t_domain, t_field>();
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        std::array<MappingValues, getSize<t_domain, t_field>() - 1>
        getValues()
        requires(t_field.isTensor(1) && t_domain.hasDim(2))
        {
            return {
                MappingValues{0, 0, 0, 1},
                MappingValues{1, 1, 1, 1},
                // MappingValues{2, 2, 2, 0},
                MappingValues{0, 1, 3, lolita::numerics::sqrt(2)},
            };
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        std::array<MappingValues, getSize<t_domain, t_field>()>
        getValues()
        requires(t_field.isTensor(1) && t_domain.hasDim(3))
        {
            return {
                MappingValues{0, 0, 0, 1},
                MappingValues{1, 1, 1, 1},
                MappingValues{2, 2, 2, 1},
                MappingValues{0, 1, 3, lolita::numerics::sqrt(2)},
                MappingValues{0, 2, 4, lolita::numerics::sqrt(2)},
                MappingValues{1, 2, 5, lolita::numerics::sqrt(2)},
            };
        }
        
        static
        void
        non_linear(
            auto & gradient
        )
        {}
        
        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static
        std::array<Real, 9>
        getStressValues(
            auto const & thermodynamic_forces,
            Integer offset
        )
        requires(t_field.isTensor(1) && t_domain.hasDim(2))
        {
            auto stress = std::array<Real, 9>();
            stress[0] = thermodynamic_forces[0 + offset];
            stress[1] = thermodynamic_forces[1 + offset];
            stress[2] = thermodynamic_forces[2 + offset];
            stress[3] = (1.0 / std::sqrt(2.0)) * thermodynamic_forces[3 + offset];
            stress[4] = (1.0 / std::sqrt(2.0)) * thermodynamic_forces[3 + offset];
            stress[5] = 0;
            stress[6] = 0;
            stress[7] = 0;
            stress[8] = 0;
            return stress;
        }
        
        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static
        std::array<Real, 9>
        getStrainValues(
            auto const & gradients,
            Integer offset
        )
        requires(t_field.isTensor(1) && t_domain.hasDim(2))
        {
            auto strain = std::array<Real, 9>();
            strain[0] = gradients[0 + offset];
            strain[1] = gradients[1 + offset];
            strain[2] = gradients[2 + offset];
            strain[3] = (1.0 / std::sqrt(2.0)) * gradients[3 + offset];
            strain[4] = (1.0 / std::sqrt(2.0)) * gradients[3 + offset];
            strain[5] = 0;
            strain[6] = 0;
            strain[7] = 0;
            strain[8] = 0;
            return strain;
        }
        
    };

    template<MappingConcept auto t_mapping>
    requires(t_mapping.isLargeStrain())
    struct MappingTraits<t_mapping>
    {

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getCols()
        {
            return 1;
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getRows()
        requires(t_field.isTensor(1) && t_domain.hasDim(2))
        {
            return 5;
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getRows()
        requires(t_field.isTensor(1) && t_domain.hasDim(3))
        {
            return 9;
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getSize()
        {
            return getRows<t_domain, t_field>() * getCols<t_domain, t_field>();
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        std::array<MappingValues, getSize<t_domain, t_field>() - 1>
        getValues()
        requires(t_field.isTensor(1) && t_domain.hasDim(2))
        {
            return {
                MappingValues{0, 0, 0, 1},
                MappingValues{1, 1, 1, 1},
                // MappingValues{2, 2, 2, 1},
                MappingValues{0, 1, 3, 1},
                MappingValues{1, 0, 4, 1},
            };
        }

        template<Domain t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        std::array<MappingValues, getSize<t_domain, t_field>()>
        getValues()
        requires(t_field.isTensor(1) && t_domain.hasDim(3))
        {
            return {
                MappingValues{0, 0, 0, 1},
                MappingValues{1, 1, 1, 1},
                MappingValues{2, 2, 2, 1},
                MappingValues{0, 1, 3, 1},
                MappingValues{1, 0, 4, 1},
                MappingValues{0, 2, 5, 1},
                MappingValues{2, 0, 6, 1},
                MappingValues{1, 2, 7, 1},
                MappingValues{2, 1, 8, 1},
            };
        }

        static
        void
        non_linear(
                auto & gradient
        )
        {
            gradient(0) += 1.0;
            gradient(1) += 1.0;
            gradient(2) += 1.0;
        }
        
    };

    template<PotentialConcept auto... t_potentials>
    struct PotentialTraits
    {

    private:

        template<PotentialConcept auto t_potential>
        struct PotentialView
        {
            
            static constexpr
            auto
            getPotential()
            {
                return t_potential;
            }
        
        };

        template<MappingConcept auto t_mapping>
        struct MappingView
        {
            
            static constexpr
            auto
            getMapping()
            {
                return t_mapping;
            }
        
        };

        template<FieldConcept auto t_field>
        struct FieldView
        {
            
            static constexpr
            auto
            getField()
            {
                return t_field;
            }
        
        };

        using Potentials = utility::tuple_unique_t<std::tuple<PotentialView<t_potentials>...>>;

        using Strains = utility::tuple_unique_t<utility::tuple_expansion_t<std::tuple, MappingView, t_potentials.getStrains()...>>;

        using Fields = utility::tuple_unique_t<utility::tuple_expansion_t<std::tuple, FieldView, t_potentials.getFields()...>>;

    public:
        
        template<Integer t_i>
        static constexpr
        auto
        getPotential()
        {
            return std::tuple_element_t<t_i, Potentials>::getPotential();
        }
        
        template<Integer t_i>
        static constexpr
        auto
        getStrain()
        {
            return std::tuple_element_t<t_i, Strains>::getMapping();
        }
        
        template<Integer t_i>
        static constexpr
        auto
        getField()
        {
            return std::tuple_element_t<t_i, Fields>::getField();
        }
        
        static constexpr
        Integer
        getNumPotentials()
        {
            return std::tuple_size_v<Potentials>;
        }
        
        static constexpr
        Integer
        getNumMappings()
        {
            return std::tuple_size_v<Strains>;
        }
        
        static constexpr
        Integer
        getNumFields()
        {
            return std::tuple_size_v<Fields>;
        }

        template<Domain t_domain>
        static constexpr
        Integer
        getSize()
        {
            auto size = Integer(0);
            auto set_size = [&] <Integer t_i = 0> (
                auto & t_set_size
            )
            constexpr mutable
            {
                size += MappingTraits<getStrain<t_i>()>::template getSize<t_domain>();
                if constexpr (t_i < getNumMappings() - 1)
                {
                    t_set_size.template operator ()<t_i + 1>(t_set_size);
                }
            };
            set_size(set_size);
            return size;
        }

        template<auto t_element, Domain t_domain>
        static constexpr
        Integer
        getSize()
        {
            auto size = Integer(0);
            auto set_size = [&] <Integer t_i = 0> (
                auto & t_set_size
            )
            constexpr mutable
            {
                size += DiscretizationTraits2<getField<t_i>()>::template getStaticSize<t_element, t_domain>();
                if constexpr (t_i < getNumFields() - 1)
                {
                    t_set_size.template operator ()<t_i + 1>(t_set_size);
                }
            };
            set_size(set_size);
            return size;
        }

        template<PotentialConcept auto t_potential>
        static constexpr
        Integer
        getIndex()
        {
            auto index = Integer(-1);
            auto set_index = [&] <Integer t_i = 0> (
                auto & t_set_index
            )
            constexpr mutable
            {
                if (utility::areEqual(t_potential, getPotential<t_i>()))
                {
                    index = t_i;
                }
                if constexpr (t_i < getNumPotentials() - 1)
                {
                    t_set_index.template operator ()<t_i + 1>(t_set_index);
                }
            };
            set_index(set_index);
            return index;
        }

        template<MappingConcept auto t_mapping>
        static constexpr
        Integer
        getIndex()
        {
            auto index = Integer(-1);
            auto set_index = [&] <Integer t_i = 0> (
                auto & t_set_index
            )
            constexpr mutable
            {
                if (utility::areEqual(t_mapping, getStrain<t_i>()))
                {
                    index = t_i;
                }
                if constexpr (t_i < getNumMappings() - 1)
                {
                    t_set_index.template operator ()<t_i + 1>(t_set_index);
                }
            };
            set_index(set_index);
            return index;
        }

        template<FieldConcept auto t_field>
        static constexpr
        Integer
        getIndex()
        {
            auto index = Integer(-1);
            auto set_index = [&] <Integer t_i = 0> (
                auto & t_set_index
            )
            constexpr mutable
            {
                if (utility::areEqual(t_field, getField<t_i>()))
                {
                    index = t_i;
                }
                if constexpr (t_i < getNumFields() - 1)
                {
                    t_set_index.template operator ()<t_i + 1>(t_set_index);
                }
            };
            set_index(set_index);
            return index;
        }

        template<Domain t_domain, MappingConcept auto t_mapping>
        static constexpr
        Integer
        getOffset()
        {
            auto offset = Integer(0);
            auto set_offset = [&] <Integer t_i = 0> (
                auto & t_set_offset
            )
            constexpr mutable
            {
                if constexpr (t_i > 0)
                {
                    offset += MappingTraits<getStrain<t_i - 1>()>::template getSize<t_domain>();
                }
                if constexpr (utility::areEqual(getStrain<t_i>(), t_mapping))
                {
                    return;
                }
                else if constexpr (t_i < getNumMappings() - 1)
                {
                    t_set_offset.template operator ()<t_i + 1>(t_set_offset);
                }
            };
            set_offset(set_offset);
            return offset;
        }

        template<auto t_element, Domain t_domain, FieldConcept auto t_field>
        static constexpr
        Integer
        getOffset()
        {
            auto offset = Integer(0);
            auto set_offset = [&] <Integer t_i = 0> (
                auto & t_set_offset
            )
            constexpr mutable
            {
                if constexpr (t_i > 0)
                {
                    offset += DiscretizationTraits2<getField<t_i - 1>()>::template getStaticSize<t_element, t_domain>();
                }
                if constexpr (utility::areEqual(getField<t_i>(), t_field))
                {
                    return;
                }
                else if constexpr (t_i < getNumFields() - 1)
                {
                    t_set_offset.template operator ()<t_i + 1>(t_set_offset);
                }
            };
            set_offset(set_offset);
            return offset;
        }

    };

} // namespace lolita

#endif /* DA9C5D5A_13CE_4129_A31C_3550B18DAB24 */
