#ifndef DC705F9F_5889_448A_AE5F_A4973293370D
#define DC705F9F_5889_448A_AE5F_A4973293370D

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"
#include "lolita/lolita_core_1.hxx"
#include "lolita/lolita_core_2.hxx"
#include "lolita/lolita_core_3.hxx"

namespace lolita::core::system
{
    
    /**
     * @brief 
     * 
     */
    struct FiniteElementLinearSystem
    {

        /**
         * @brief
         */
        lolita::natural num_unknowns_ = 0;

        /**
         * @brief
         */
        lolita::natural num_bindings_ = 0;

        /**
         * @brief
         */
        Eigen::SparseMatrix<lolita::real> lhs_;

        /**
         * @brief
         */
        lolita::matrix::Vector<lolita::real> rhs_;

    };

}

#endif /* DC705F9F_5889_448A_AE5F_A4973293370D */
