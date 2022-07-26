#ifndef DA9C5D5A_13CE_4129_A31C_3550B18DAB24
#define DA9C5D5A_13CE_4129_A31C_3550B18DAB24

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_defs.hxx"

namespace lolita2::geometry
{

    template<Field t_field>
    struct FieldTraits
    {

        template<Domain t_domain>
        static constexpr
        lolita::matrix::Shape
        shape()
        requires(t_field.isTensor(0))
        {
            return {
                lolita::numerics::pow(t_domain.dim_, 0),
                lolita::numerics::pow(t_domain.dim_, 0)
            };
        }

        template<Domain t_domain>
        static constexpr
        lolita::matrix::Shape
        shape()
        requires(t_field.isTensor(1))
        {
            return {
                lolita::numerics::pow(t_domain.dim_, 0),
                lolita::numerics::pow(t_domain.dim_, 1)
            };
        }

        template<Domain t_domain>
        static constexpr
        lolita::matrix::Shape
        shape()
        requires(t_field.isTensor(2))
        {
            return {
                lolita::numerics::pow(t_domain.dim_, 1),
                lolita::numerics::pow(t_domain.dim_, 1)
            };
        }

        template<Domain t_domain>
        static constexpr
        lolita::matrix::Shape
        shape()
        requires(t_field.isTensor(3))
        {
            return {
                lolita::numerics::pow(t_domain.dim_, 1),
                lolita::numerics::pow(t_domain.dim_, 2)
            };
        }

        template<Domain t_domain>
        static constexpr
        lolita::matrix::Shape
        shape()
        requires(t_field.isTensor(4))
        {
            return {
                lolita::numerics::pow(t_domain.dim_, 2),
                lolita::numerics::pow(t_domain.dim_, 2)
            };
        }

        template<Domain t_domain>
        static constexpr
        lolita::integer
        size()
        {
            return shape<t_domain>().size_;
        }

    };

    struct MappingValues
    {

        constexpr
        MappingValues(
            lolita::index row,
            lolita::index col,
            lolita::index position,
            lolita::real coefficient
        )
        :
        row_(row),
        col_(col),
        position_(position),
        coefficient_(coefficient)
        {}

        constexpr
        lolita::integer
        row()
        const
        {
            return row_;
        }

        constexpr
        lolita::integer
        col()
        const
        {
            return col_;
        }

        constexpr
        lolita::integer
        rank()
        const
        {
            return position_;
        }

        constexpr
        lolita::real
        value()
        const
        {
            return coefficient_;
        }

        lolita::index row_;

        lolita::index col_;

        lolita::index position_;

        lolita::real coefficient_;

    };

    template<Mapping t_mapping>
    struct MappingTraits;

    template<Mapping t_mapping>
    requires(t_mapping.isIdentity())
    struct MappingTraits<t_mapping>
    {

        template<Domain t_domain, Field t_field>
        static constexpr
        lolita::matrix::Shape
        shape()
        {
            return FieldTraits<t_field>::template shape<t_domain>();
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        lolita::integer
        size()
        {
            return shape<t_domain, t_field>().size_;
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        std::array<MappingValues, size<t_domain, t_field>()>
        getValues()
        requires(t_field.isTensor(0))
        {
            return {
                MappingValues{0, 0, 0, 1},
            };
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        std::array<MappingValues, size<t_domain, t_field>()>
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
        std::array<MappingValues, size<t_domain, t_field>()>
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
        lolita::matrix::Shape
        shape()
        {
            return FieldTraits<Field(t_field.dim_ + 1)>::template shape<t_domain>();
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        lolita::integer
        size()
        {
            return shape<t_domain, t_field>().size_;
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        std::array<MappingValues, size<t_domain, t_field>()>
        getValues()
        requires(t_field.isTensor(0) && t_domain.hasDim(1))
        {
            return {
                MappingValues{0, 0, 0, 1},
            };
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        std::array<MappingValues, size<t_domain, t_field>()>
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
        std::array<MappingValues, size<t_domain, t_field>()>
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
        std::array<MappingValues, size<t_domain, t_field>()>
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
        std::array<MappingValues, size<t_domain, t_field>()>
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
        lolita::matrix::Shape
        shape()
        requires(t_field.isTensor(1) && t_domain.hasDim(2))
        {
            return {1, 4};
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        lolita::matrix::Shape
        shape()
        requires(t_field.isTensor(1) && t_domain.hasDim(3))
        {
            return {1, 6};
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        lolita::integer
        size()
        {
            return shape<t_domain, t_field>().size_;
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        std::array<MappingValues, size<t_domain, t_field>()>
        getValues()
        requires(t_field.isTensor(1) && t_domain.hasDim(2))
        {
            return {
                MappingValues{0, 0, 0, 1},
                MappingValues{1, 1, 1, 1},
                MappingValues{2, 2, 2, 1},
                MappingValues{0, 1, 3, lolita::numerics::sqrt_2},
            };
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        std::array<MappingValues, size<t_domain, t_field>()>
        getValues()
        requires(t_field.isTensor(1) && t_domain.hasDim(3))
        {
            return {
                MappingValues{0, 0, 0, 1},
                MappingValues{1, 1, 1, 1},
                MappingValues{2, 2, 2, 1},
                MappingValues{0, 1, 3, lolita::numerics::sqrt_2},
                MappingValues{0, 2, 4, lolita::numerics::sqrt_2},
                MappingValues{1, 2, 5, lolita::numerics::sqrt_2},
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
        lolita::matrix::Shape
        shape()
        requires(t_field.isTensor(1) && t_domain.hasDim(2))
        {
            return {1, 5};
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        lolita::matrix::Shape
        shape()
        requires(t_field.isTensor(1) && t_domain.hasDim(3))
        {
            return {1, 9};
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        lolita::integer
        size()
        {
            return shape<t_domain, t_field>().size_;
        }

        template<Domain t_domain, Field t_field>
        static constexpr
        std::array<MappingValues, size<t_domain, t_field>()>
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
        std::array<MappingValues, size<t_domain, t_field>()>
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

    template<BehaviorConcept auto t_behavior>
    struct BehaviorTraits
    {

        template<Domain t_domain>
        static constexpr
        lolita::integer
        strainSize()
        {
            auto p = lolita::integer(0);
            auto cnt = [&] <lolita::integer t_i = 0, lolita::integer t_j = 0> (
                auto & self
            )
            constexpr mutable
            {
                auto constexpr t_field = t_behavior.template getUnknown<t_i>().getField();
                auto constexpr t_mapping = t_behavior.template getUnknown<t_i>().template getMapping<t_j>();
                p += MappingTraits<t_mapping>::template size<t_domain, t_field>();
                if constexpr (t_j < t_behavior.template getUnknown<t_i>().getNumMappings() - 1)
                {
                    self.template operator ()<t_i, t_j + 1>(self);
                }
                else if constexpr (t_i < t_behavior.getNumUnknowns() - 1)
                {
                    self.template operator ()<t_i + 1, 0>(self);
                }
            };
            cnt(cnt);
            return p;
        }

    };

}

#endif /* DA9C5D5A_13CE_4129_A31C_3550B18DAB24 */
