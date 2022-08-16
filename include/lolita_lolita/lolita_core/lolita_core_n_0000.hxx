#ifndef DA9C5D5A_13CE_4129_A31C_3550B18DAB24
#define DA9C5D5A_13CE_4129_A31C_3550B18DAB24

#include "lolita_lolita/lolita_core/lolita.hxx"

namespace lolita
{
    template<Field t_field>
    struct FieldTraits
    {

        template<Domain t_domain>
        static constexpr
        Integer
        getRows()
        requires(t_field.isTensor(0))
        {
            return lolita::numerics::pow(t_domain.dim_, 0);
        }

        template<Domain t_domain>
        static constexpr
        Integer
        getRows()
        requires(t_field.isTensor(1))
        {
            return lolita::numerics::pow(t_domain.dim_, 0);
        }

        template<Domain t_domain>
        static constexpr
        Integer
        getRows()
        requires(t_field.isTensor(2))
        {
            return lolita::numerics::pow(t_domain.dim_, 1);
        }

        template<Domain t_domain>
        static constexpr
        Integer
        getRows()
        requires(t_field.isTensor(3))
        {
            return lolita::numerics::pow(t_domain.dim_, 1);
        }

        template<Domain t_domain>
        static constexpr
        Integer
        getRows()
        requires(t_field.isTensor(4))
        {
            return lolita::numerics::pow(t_domain.dim_, 2);
        }

        // ---

        template<Domain t_domain>
        static constexpr
        Integer
        getCols()
        requires(t_field.isTensor(0))
        {
            return lolita::numerics::pow(t_domain.dim_, 0);
        }

        template<Domain t_domain>
        static constexpr
        Integer
        getCols()
        requires(t_field.isTensor(1))
        {
            return lolita::numerics::pow(t_domain.dim_, 1);
        }

        template<Domain t_domain>
        static constexpr
        Integer
        getCols()
        requires(t_field.isTensor(2))
        {
            return lolita::numerics::pow(t_domain.dim_, 1);
        }

        template<Domain t_domain>
        static constexpr
        Integer
        getCols()
        requires(t_field.isTensor(3))
        {
            return lolita::numerics::pow(t_domain.dim_, 2);
        }

        template<Domain t_domain>
        static constexpr
        Integer
        getCols()
        requires(t_field.isTensor(4))
        {
            return lolita::numerics::pow(t_domain.dim_, 2);
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

    template<Mapping t_mapping>
    struct MappingTraits;

    template<Mapping t_mapping>
    requires(t_mapping.isIdentity())
    struct MappingTraits<t_mapping>
    {

        template<Domain t_domain, Field t_field>
        static constexpr
        Integer
        getRows()
        {
            return FieldTraits<t_field>::template getRows<t_domain>();
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        Integer
        getCols()
        {
            return FieldTraits<t_field>::template getCols<t_domain>();
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        Integer
        getSize()
        {
            return getRows<t_domain, t_field>() * getCols<t_domain, t_field>();
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        std::array<MappingValues, getSize<t_domain, t_field>()>
        getValues()
        requires(t_field.isTensor(0))
        {
            return {
                MappingValues{0, 0, 0, 1},
            };
        }

        template<Domain t_domain, Field t_field>
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

        template<Domain t_domain, Field t_field>
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
        
    };

    template<Mapping t_mapping>
    requires(t_mapping.isGradient())
    struct MappingTraits<t_mapping>
    {

        template<Domain t_domain, Field t_field>
        static constexpr
        Integer
        getRows()
        {
            return FieldTraits<Field(t_field.getDim() + 1)>::template getRows<t_domain>();
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        Integer
        getCols()
        {
            return FieldTraits<Field(t_field.getDim() + 1)>::template getCols<t_domain>();
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        Integer
        getSize()
        {
            return getRows<t_domain, t_field>() * getCols<t_domain, t_field>();
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        std::array<MappingValues, getSize<t_domain, t_field>()>
        getValues()
        requires(t_field.isTensor(0) && t_domain.hasDim(1))
        {
            return {
                MappingValues{0, 0, 0, 1},
            };
        }

        template<Domain t_domain, Field t_field>
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

        template<Domain t_domain, Field t_field>
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

        template<Domain t_domain, Field t_field>
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

        template<Domain t_domain, Field t_field>
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
        
    };

    template<Mapping t_mapping>
    requires(t_mapping.isSmallStrain())
    struct MappingTraits<t_mapping>
    {

        template<Domain t_domain, Field t_field>
        static constexpr
        Integer
        getRows()
        requires(t_field.isTensor(1) && t_domain.hasDim(2))
        {
            return 4;
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        Integer
        getRows()
        requires(t_field.isTensor(1) && t_domain.hasDim(3))
        {
            return 6;
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        Integer
        getCols()
        {
            return 1;
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        Integer
        getSize()
        {
            return getRows<t_domain, t_field>() * getCols<t_domain, t_field>();
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        std::array<MappingValues, getSize<t_domain, t_field>()>
        getValues()
        requires(t_field.isTensor(1) && t_domain.hasDim(2))
        {
            return {
                MappingValues{0, 0, 0, 1},
                MappingValues{1, 1, 1, 1},
                MappingValues{2, 2, 2, 1},
                MappingValues{0, 1, 3, lolita::numerics::sqrt(2)},
            };
        }

        template<Domain t_domain, Field t_field>
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
        
    };

    template<Mapping t_mapping>
    requires(t_mapping.isLargeStrain())
    struct MappingTraits<t_mapping>
    {

        template<Domain t_domain, Field t_field>
        static constexpr
        Integer
        getCols()
        {
            return 1;
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        Integer
        getRows()
        requires(t_field.isTensor(1) && t_domain.hasDim(2))
        {
            return 5;
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        Integer
        getRows()
        requires(t_field.isTensor(1) && t_domain.hasDim(3))
        {
            return 9;
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        Integer
        getSize()
        {
            return getRows<t_domain, t_field>() * getCols<t_domain, t_field>();
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        std::array<MappingValues, getSize<t_domain, t_field>()>
        getValues()
        requires(t_field.isTensor(1) && t_domain.hasDim(2))
        {
            return {
                MappingValues{0, 0, 0, 1},
                MappingValues{1, 1, 1, 1},
                MappingValues{2, 2, 2, 1},
                MappingValues{0, 1, 3, 1},
                MappingValues{1, 0, 4, 1},
            };
        }

        template<Domain t_domain, Field t_field>
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

    template<GeneralizedStrainConcept auto t_generalized_strain>
    struct GeneralizedStrainTraits
    {

        template<Domain t_domain>
        static constexpr
        Integer
        getSize()
        {
            auto constexpr t_field = t_generalized_strain.getField();
            auto size = Integer(0);
            auto set_size = [&] <Integer t_i = 0> (
                auto & self
            )
            constexpr mutable
            {
                auto constexpr t_mapping = t_generalized_strain.template getMapping<t_i>();
                size += MappingTraits<t_mapping>::template getSize<t_domain, t_field>();
                if constexpr (t_i < t_generalized_strain.getNumMappings() - 1)
                {
                    self.template operator ()<t_i + 1>(self);
                }
            };
            set_size(set_size);
            return size;
        }

    };

    template<BehaviorConcept auto t_behavior>
    struct BehaviorTraits
    {

        template<Domain t_domain>
        static constexpr
        Integer
        getGeneralizedStrainSize()
        {
            auto p = Integer(0);
            auto cnt = [&] <Integer t_i = 0> (
                auto & self
            )
            constexpr mutable
            {
                p += GeneralizedStrainTraits<t_behavior.template getGeneralizedStrain<t_i>()>::template getSize<t_domain>();
                if constexpr (t_i < t_behavior.getNumGeneralizedStrains() - 1)
                {
                    self.template operator ()<t_i + 1>(self);
                }
            };
            cnt(cnt);
            return p;
        }

    };

    template<FiniteElementMethodConcept auto t_finite_element_method>
    struct FiniteElementMethodTraits
    {

        template<Domain t_domain>
        static constexpr
        Integer
        getGeneralizedStrainSize()
        {
            return GeneralizedStrainTraits<t_finite_element_method.getGeneralizedStrain()>::template getSize<t_domain>();
        }

        template<Domain t_domain>
        static constexpr
        Integer
        getGeneralizedStrainOffset()
        {
            auto constexpr t_finite_element_generalized_strain = t_finite_element_method.getGeneralizedStrain();
            auto offset = Integer(0);
            auto is_set = false;
            auto set_offset = [&] <Integer t_i = 0> (
                auto & self
            )
            constexpr mutable
            {
                auto constexpr t_generalized_strain = t_finite_element_method.getBehavior().template getGeneralizedStrain<t_i>();
                if constexpr (std::is_same_v<std::decay_t<decltype(t_generalized_strain)>, std::decay_t<decltype(t_finite_element_generalized_strain)>>)
                {
                    if constexpr (t_generalized_strain == t_finite_element_generalized_strain)
                    {
                        is_set = true;
                    }
                }
                if (!is_set)
                {
                    offset += GeneralizedStrainTraits<t_generalized_strain>::template getSize<t_domain>();
                }
                if constexpr (t_i < t_finite_element_method.getBehavior().getNumGeneralizedStrains() - 1)
                {
                    self.template operator ()<t_i + 1>(self);
                }
            };
            set_offset(set_offset);
            return offset;
        }

        template<Domain t_domain, Mapping t_mapping>
        static constexpr
        Integer
        getMappingSize()
        {
            return MappingTraits<t_mapping>::template getSize<t_domain, t_finite_element_method.getField()>();
        }

    };

} // namespace lolita

#endif /* DA9C5D5A_13CE_4129_A31C_3550B18DAB24 */
