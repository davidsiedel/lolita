#ifndef A2996F51_A92E_4F9C_9645_40D250D9C8BB
#define A2996F51_A92E_4F9C_9645_40D250D9C8BB

#include "2/core/_include.hxx"

#include "2/core/region_potential.hxx"
#include "2/core/region_lagrangian_implementation.hxx"
#include "2/core/region_lagrangian_interface.hxx"

#include "2/core/dof.hxx"

namespace lolita::core
{

    template<DomainConcept auto t_dim, MeshConcept auto t_domain>
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
        setDiscreteField()
        {
            if (ptr_data_ == nullptr)
            {
                ptr_data_ = std::make_unique<std::vector<std::unique_ptr<MeshDiscreteField<t_domain>>>>();
            }
            for (auto const & item : * ptr_data_)
            {
                if (item->getLabel() == t_field.getLabel())
                {
                    return;
                }
            }
            ptr_data_->push_back(std::make_unique<MeshDiscreteField<t_domain>>(t_field));
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
                    if (item->getLabel() == t_field.getLabel())
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
                    if (item->getLabel() == t_field.getLabel())
                    {
                        return * item;
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

        using t_Lag = AbstractDomainLagrangian<t_dim, t_domain>;

        template<LagrangianConcept auto t_lag>
        using t_LagImpl = DomainLagrangian<t_dim, t_domain, t_lag>;

        template<LagrangianConcept auto t_lag>
        void
        setLagrangian()
        {
            if (lags_ == nullptr)
            {
                lags_ = std::make_unique<std::vector<std::unique_ptr<t_Lag>>>();
            }
            for (auto & kk : * lags_)
            {
                if (kk->getLabel() == t_lag.getLabel())
                {
                    kk = std::make_unique<t_LagImpl<t_lag>>(* this);
                    return;
                }
            }
            lags_->push_back(std::make_unique<t_LagImpl<t_lag>>(* this));
        }

        template<LagrangianConcept auto t_lag>
        Boolean
        hasLagrangian()
        const
        {
            if (lags_ == nullptr)
            {
                return false;
            }
            for (auto const & kk : * lags_)
            {
                if (kk->getLabel() == t_lag.getLabel())
                {
                    return true;
                }
            }
            return false;
        }

        template<LagrangianConcept auto t_lag>
        t_Lag const &
        getLagrangian()
        const
        {
            if (lags_ == nullptr)
            {
                throw std::runtime_error("NO");
            }
            for (auto const & kk : * lags_)
            {
                if (kk->getLabel() == t_lag.getLabel())
                {
                    return * kk;
                }
            }
            throw std::runtime_error("NO");
        }

        template<LagrangianConcept auto t_lag>
        t_Lag &
        getLagrangian()
        {
            if (lags_ == nullptr)
            {
                throw std::runtime_error("NO");
            }
            for (auto & kk : * lags_)
            {
                if (kk->getLabel() == t_lag.getLabel())
                {
                    return * kk;
                }
            }
            throw std::runtime_error("NO");
        }

        // template<LagrangianConcept auto t_lag, PotentialConcept auto t_potential>
        // void
        // setPotential(
        //     auto const &... args
        // )
        // {
        //     if (!this->template hasLagrangian<t_lag>())
        //     {
        //         this->template setLagrangian<t_lag>();
        //     }
        //     this->template getLagrangian<t_lag>().template setPotential<t_lag, t_potential>(args...);
        // }

    private:

        std::unique_ptr<std::vector<std::unique_ptr<t_Lag>>> lags_;

        std::unique_ptr<std::vector<MeshDiscreteField<t_domain>>> ptr_data_;

        std::basic_string<Character> tag_;

    };
    
} // namespace lolita


#endif /* A2996F51_A92E_4F9C_9645_40D250D9C8BB */
