/**
 * @file point.hxx
 * @author David Siedel (davidsiedel@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2023-03-29
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef DA6BDF79_2135_42C7_A76F_92B4C1BD50EF
#define DA6BDF79_2135_42C7_A76F_92B4C1BD50EF

#include "config.hxx"
// #include "numerics.hxx"
// #include "geometry/frame.hxx"
#include "tensor.hxx"

namespace lolita::geometry
{
    
    /**
     * @brief 
     * 
     * @tparam dim_euclidean_ 
     */
    template<Integer dim_euclidean_>
    struct Point
    {
        
        /**
         * @brief Default constructor with zero values in one dimension.
         * 
         */
        constexpr
        Point()
        requires(dim_euclidean_ == 1)
        :
        coordinates_({0})
        {}
        
        /**
         * @brief Default constructor with zero values in two dimensions.
         * 
         */
        constexpr
        Point()
        requires(dim_euclidean_ == 2)
        :
        coordinates_({0, 0})
        {}
        
        /**
         * @brief Default constructor with zero values in three dimensions.
         * 
         */
        constexpr
        Point()
        requires(dim_euclidean_ == 3)
        :
        coordinates_({0, 0, 0})
        {}

        /**
         * @brief Copy constructor using known coordinates.
         * 
         * @tparam Reals_ The coordinates types.
         */
        template<typename... Reals_>
        constexpr
        Point(
            Reals_ const &... x
        )
        requires(sizeof...(Reals_) == dim_euclidean_)
        :
        coordinates_({x...})
        {}

        /**
         * @brief Move constructor using known coordinates.
         * 
         * @tparam Reals_ The coordinates types.
         */
        template<typename... Reals_>
        constexpr
        Point(
            Reals_ &&... x
        )
        requires(sizeof...(Reals_) == dim_euclidean_)
        :
        coordinates_({std::move(x)...})
        {}

        /**
         * @brief Get the ith coordinate.
         * 
         * @param i The ith coordinate to get.
         * @return constexpr Real const& 
         */
        constexpr
        Real const &
        getCoordinate(
            Integer i
        )
        const
        {
            return coordinates_[i];
        }

        /**
         * @brief Set the ith coordinate.
         * 
         * @param i 
         * @param value 
         */
        constexpr
        void
        setCoordinate(
            Integer i,
            Real const & value
        )
        {
            coordinates_[i] = value;
        }

        /**
         * @brief 
         * 
         * @return tensor::StaticTensor<lolita::Real, 3> 
         */
        tensor::StaticTensor<lolita::Real, 3>
        asTensor()
        const
        {
            return tensor::StaticTensor<lolita::Real, 3>();
        }

    private:

        /**
         * @brief 
         * 
         */
        std::array<Real, dim_euclidean_> coordinates_;

    };

} // namespace lolita::geometry


#endif /* DA6BDF79_2135_42C7_A76F_92B4C1BD50EF */
