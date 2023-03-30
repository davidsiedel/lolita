/**
 * @file frame.hxx
 * @author David Siedel (davidsiedel@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2023-03-29
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef FBA75E4B_A209_4011_87A7_FF975994C267
#define FBA75E4B_A209_4011_87A7_FF975994C267

namespace lolita::geometry
{

    /**
     * @brief 
     * 
     * @tparam T_ 
     */
    template<typename T_>
    concept FrameConcept = requires
    {

        { T_::getDimEuclidean() } -> std::same_as<lolita::Integer>;
        
    };
    
    /**
     * @brief 
     * 
     * @tparam dim_euclidean_ 
     */
    template<Integer dim_euclidean_>
    struct CartesianFrame
    {

        /**
         * @brief Get the Dim Euclidean object
         * 
         * @return constexpr Integer 
         */
        static constexpr
        Integer
        getDimEuclidean()
        {
            return dim_euclidean_;
        }

    };

    /**
     * @brief 
     * 
     */
    struct AxiSymmetricFrame
    {

        /**
         * @brief Get the Dim Euclidean object
         * 
         * @return constexpr Integer 
         */
        static constexpr
        Integer
        getDimEuclidean()
        {
            return 2;
        }

    };
    
} // namespace lolita::geometry


#endif /* FBA75E4B_A209_4011_87A7_FF975994C267 */
