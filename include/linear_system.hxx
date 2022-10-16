#ifndef E835FAB2_0F70_40F4_B4AE_8420BF21D16F
#define E835FAB2_0F70_40F4_B4AE_8420BF21D16F

#include "config.hxx"
#include "algebra.hxx"

namespace lolita
{

    template<auto t_s>
    struct LinearSystem
    {

        static
        std::unique_ptr<LinearSystem>
        make_unique()
        {
            return std::make_unique<LinearSystem>();
        }

        LinearSystem()
        :
        rhs_size_(0),
        lhs_size_(0)
        {}

        inline
        Boolean
        operator==(
            LinearSystem const & other
        )
        const = default;

        inline
        Boolean
        operator!=(
            LinearSystem const & other
        )
        const = default;

        std::atomic<Natural> const &
        getSize()
        const
        {
            return rhs_size_;
        }

        std::atomic<Natural> &
        getSize()
        {
            return rhs_size_;
        }

        std::atomic<Natural> const &
        getLhsSize()
        const
        {
            return lhs_size_;
        }

        std::atomic<Natural> const &
        getRhsSize()
        const
        {
            return rhs_size_;
        }

        void
        addRhsSize(
            Natural && size
        )
        {
            rhs_size_ += std::forward<Natural>(size);
        }

        void
        addLhsSize(
            Natural && size
        )
        {
            lhs_size_ += std::forward<Natural>(size);
        }

    private:

        std::atomic<Natural> rhs_size_;

        std::atomic<Natural> lhs_size_;

    };

} // namespace lolita

#endif /* E835FAB2_0F70_40F4_B4AE_8420BF21D16F */
