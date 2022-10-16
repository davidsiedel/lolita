#ifndef D6DB7C02_99A9_4252_850A_30F636C3824D
#define D6DB7C02_99A9_4252_850A_30F636C3824D

#include "2/core/traits/_include.hxx"
#include "2/core/traits/potential_traits.hxx"
#include "2/core/traits/field_traits.hxx"
#include "2/core/traits/linear_operator_traits.hxx"

namespace lolita::core
{
    
    template<LagrangianConcept auto t_lag>
    struct LagTraits
    {

    private:

        template<PotentialConcept auto t_potential>
        using StrainViews_ = typename PotentialTraits2<t_potential>::Strains;

        template<PotentialConcept auto t_potential>
        using FieldViews_ = typename PotentialTraits2<t_potential>::Fields;

    public:

        using Potentials = utility::tuple_unique_t<utility::tuple_expansion_t<std::tuple, PotentialView, t_lag.getPotentials()>>;

        using Strains = utility::tuple_unique_t<utility::tuple_merge_t<utility::tuple_expansion_t<std::tuple, StrainViews_, t_lag.getPotentials()>>>;

        using Fields = utility::tuple_unique_t<utility::tuple_merge_t<utility::tuple_expansion_t<std::tuple, FieldViews_, t_lag.getPotentials()>>>;
    
        static constexpr
        auto const &
        getLagrangian()
        {
            return t_lag;
        }
        
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
        getLinearOperator()
        {
            return std::tuple_element_t<t_i, Strains>::getLinearOperator();
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
        getNumLinearOperators()
        {
            return std::tuple_size_v<Strains>;
        }
        
        static constexpr
        Integer
        getNumFields()
        {
            return std::tuple_size_v<Fields>;
        }

        template<MeshConcept auto t_domain>
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
                size += FieldTraits<getLinearOperator<t_i>()>::template getSize<t_domain>();
                if constexpr (t_i < getNumLinearOperators() - 1)
                {
                    t_set_size.template operator ()<t_i + 1>(t_set_size);
                }
            };
            set_size(set_size);
            return size;
        }

        template<auto t_element, MeshConcept auto t_domain>
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
                size += DiscretizationTraits<getField<t_i>()>::template getTensorSpaceSize<t_element, t_domain>();
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
        getPotentialIndex()
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

        template<LinearOperatorConcept auto t_mapping>
        static constexpr
        Integer
        getLinearOperatorIndex()
        {
            auto index = Integer(-1);
            auto set_index = [&] <Integer t_i = 0> (
                auto & t_set_index
            )
            constexpr mutable
            {
                if (utility::areEqual(t_mapping, getLinearOperator<t_i>()))
                {
                    index = t_i;
                }
                if constexpr (t_i < getNumLinearOperators() - 1)
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
        getFieldIndex()
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

        template<MeshConcept auto t_domain, LinearOperatorConcept auto t_mapping>
        static constexpr
        Integer
        getLinearOperatorOffset()
        {
            auto offset = Integer(0);
            auto set_offset = [&] <Integer t_i = 0> (
                auto & t_set_offset
            )
            constexpr mutable
            {
                if constexpr (t_i > 0)
                {
                    offset += FieldTraits<getLinearOperator<t_i - 1>()>::template getSize<t_domain>();
                }
                if constexpr (utility::areEqual(getLinearOperator<t_i>(), t_mapping))
                {
                    return;
                }
                else if constexpr (t_i < getNumLinearOperators() - 1)
                {
                    t_set_offset.template operator ()<t_i + 1>(t_set_offset);
                }
            };
            set_offset(set_offset);
            return offset;
        }

        template<auto t_element, MeshConcept auto t_domain, FieldConcept auto t_field>
        static constexpr
        Integer
        getFieldOffset()
        {
            auto offset = Integer(0);
            auto set_offset = [&] <Integer t_i = 0> (
                auto & t_set_offset
            )
            constexpr mutable
            {
                if constexpr (t_i > 0)
                {
                    offset += DiscretizationTraits<getField<t_i - 1>()>::template getTensorSpaceSize<t_element, t_domain>();
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

} // namespace lolita::core

#endif /* D6DB7C02_99A9_4252_850A_30F636C3824D */
