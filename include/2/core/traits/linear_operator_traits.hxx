#ifndef D45B3361_1CF7_4F26_BC33_AAD54510770A
#define D45B3361_1CF7_4F26_BC33_AAD54510770A

#include "2/core/traits/_include.hxx"
#include "2/core/traits/field_traits.hxx"
#include "2/core/traits/basis_traits.hxx"

namespace lolita::core
{
    
    template<LinearOperatorConcept auto linear_operator_>
    struct LinearOperatorView
    {
    
        static constexpr
        LinearOperatorConcept auto const &
        getLinearOperator()
        {
            return linear_operator_;
        }

        using Type = decltype(linear_operator_);

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

    template<LinearOperatorConcept auto t_mapping>
    struct LinearOperatorTraits;

    template<LinearOperatorConcept auto t_mapping>
    requires(TraceOperatorConcept<decltype(t_mapping)>)
    struct LinearOperatorTraits<t_mapping>
    {

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getNumRows()
        {
            return FieldTraits<t_field>::template getNumRows<t_domain>();
        }

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getNumCols()
        {
            return FieldTraits<t_field>::template getNumCols<t_domain>();
        }

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getSize()
        {
            return getNumRows<t_domain, t_field>() * getNumCols<t_domain, t_field>();
        }

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        std::array<MappingValues, getSize<t_domain, t_field>()>
        getValues()
        requires(t_field.isTensor(0))
        {
            return {
                MappingValues{0, 0, 0, 1},
            };
        }

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
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

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
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
        
        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
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
        
        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
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

    template<LinearOperatorConcept auto t_mapping>
    requires(GradientOperatorConcept<decltype(t_mapping)>)
    struct LinearOperatorTraits<t_mapping>
    {

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getNumRows()
        {
            // return FieldTraits<Field(t_field.getDimDomain(), t_field.getDimTensor() + 1)>::template getNumRows<t_domain>();
            return FieldTraits<t_mapping>::template getNumRows<t_domain>();
        }

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getNumCols()
        {
            // return FieldTraits<Field(t_field.getDimDomain(), t_field.getDimTensor() + 1)>::template getNumCols<t_domain>();
            return FieldTraits<t_mapping>::template getNumCols<t_domain>();
        }

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getSize()
        {
            return getNumRows<t_domain, t_field>() * getNumCols<t_domain, t_field>();
        }

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        std::array<MappingValues, getSize<t_domain, t_field>()>
        getValues()
        requires(t_field.isTensor(0) && t_domain.hasDim(1))
        {
            return {
                MappingValues{0, 0, 0, 1},
            };
        }

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
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

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
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

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
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

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
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
        
        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
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
        
        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
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

    template<LinearOperatorConcept auto t_mapping>
    requires(SmallStrainOperatorConcept<decltype(t_mapping)>)
    struct LinearOperatorTraits<t_mapping>
    {

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getNumRows()
        requires(t_field.isTensor(1) && t_domain.hasDim(2))
        {
            return 4;
        }

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getNumRows()
        requires(t_field.isTensor(1) && t_domain.hasDim(3))
        {
            return 6;
        }

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getNumCols()
        {
            return 1;
        }

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getSize()
        {
            return getNumRows<t_domain, t_field>() * getNumCols<t_domain, t_field>();
        }

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
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

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
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
        
        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
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
        
        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
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

    template<LinearOperatorConcept auto t_mapping>
    requires(LargeStrainOperatorConcept<decltype(t_mapping)>)
    struct LinearOperatorTraits<t_mapping>
    {

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getNumCols()
        {
            return 1;
        }

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getNumRows()
        requires(t_field.isTensor(1) && t_domain.hasDim(2))
        {
            return 5;
        }

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getNumRows()
        requires(t_field.isTensor(1) && t_domain.hasDim(3))
        {
            return 9;
        }

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
        static constexpr
        Integer
        getSize()
        {
            return getNumRows<t_domain, t_field>() * getNumCols<t_domain, t_field>();
        }

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
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

        template<MeshConcept auto t_domain, FieldConcept auto t_field = t_mapping.getField()>
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

} // namespace lolita::core


#endif /* D45B3361_1CF7_4F26_BC33_AAD54510770A */
