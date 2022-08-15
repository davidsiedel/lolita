//
// Created by dsiedel on 16/02/2022.
//

#ifndef FETA__FETA_SHARED_POINTER_HXX
#define FETA__FETA_SHARED_POINTER_HXX

#include "new/_feta.hxx"

namespace feta
{

    /**
     * @brief
     */
    template<
            typename T
    >
    struct SharedPointer
    {

        /**
         * @brief
         */
        using Type = T;

        /**
         * @brief
         */
        using Data = std::shared_ptr<T>;

        struct Hash
        {

            Indx
            operator() (
                    SharedPointer const & other
            )
            const
            {
                return std::hash<Data>()(other.data);
            }

        };

        /**
         * @brief
         */
        SharedPointer()
        :
        data(nullptr)
        {}

        /**
         * @brief
         * @param t
         */
        explicit
        SharedPointer(
                T const &
                t
        )
        :
        data(setData(t))
        {}

        /**
         * @brief
         * @param t
         */
        explicit
        SharedPointer(
                T &&
                t
        )
        :
        data(setData(t))
        {}

        Bool
        operator==(
                SharedPointer const &
                other
        )
        const
        {
            auto eq_0 = data == other.data;
            return eq_0;
        }

        Bool
        operator!=(
                SharedPointer const &
                other
        )
        const
        {
            return !(other == * this);
        }

//        Bool
//        operator ==(
//                SharedPointer const &
//                other
//        )
//        const
//        {
//            return data == other.data;
//        }

//        /**
//         * @brief
//         * @return
//         */
//        T &
//        operator ()()
//        {
//            return * data;
//        }
//
//        /**
//         * @brief
//         * @return
//         */
//        T const &
//        operator ()()
//        const
//        {
//            return * data;
//        }

        /**
         * @brief
         * @return
         */
        T &
        get()
        {
            return * data;
        }

        /**
         * @brief
         * @return
         */
        T const &
        get()
        const
        {
            return * data;
        }

    private:

        /**
         * @brief
         * @param t
         * @return
         */
        static
        Data
        setData(
                T const &
                t
        )
        {
            return std::make_shared<T>(t);
        }

        /**
         * @brief
         * @param t
         * @return
         */
        static
        Data
        setData(
                T &&
                t
        )
        {
            return std::make_shared<T>(t);
        }

    public:

        /**
         * @brief
         */
        Data data;

    };

}

#endif //FETA__FETA_SHARED_POINTER_HXX
