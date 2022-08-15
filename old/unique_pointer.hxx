//
// Created by dsiedel on 16/02/2022.
//

#ifndef FETA_UNIQUE_POINTER_HXX
#define FETA_UNIQUE_POINTER_HXX

namespace feta
{

    /**
     * @brief
     */
    template<typename T>
    struct UniquePointer
    {

        /**
         * @brief
         */
        using Type = T;

        /**
         * @brief
         */
        using Data = std::unique_ptr<
                T
        >;

        struct Hash
        {

            Indx
            operator() (
                    UniquePointer const & other
            )
            const
            {
                return std::hash<Data>()(other.data);
            }

        };

        /**
         * @brief
         */
        UniquePointer()
        :
        data(nullptr)
        {}

        /**
         * @brief
         * @param t
         */
        explicit
        UniquePointer(
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
        UniquePointer(
                T &&
                t
        )
        :
        data(setData(t))
        {}

        /**
         * @brief
         * @param other
         */
        UniquePointer(
                UniquePointer const &
                other
        )
        :
        data(copyData(other))
        {}

        /**
         * @brief
         * @param other
         */
        UniquePointer(
                UniquePointer &&
                other
        )
        :
        data(std::move(other.data))
        {}

        Bool
        operator==(
                UniquePointer const &
                other
        )
        const
        {

            auto eq_0 = data == other.data;
            return eq_0;
        }

        Bool
        operator!=(
                UniquePointer const &
                other
        )
        const
        {
            return !(other == * this);
        }

        /**
         * @brief
         * @param other
         * @return
         */
        UniquePointer &
        operator =(
                UniquePointer const &
                other
        )
        {
            data = copyData(other);
            return * this;
        }

        /**
         * @brief
         * @param other
         * @return
         */
        UniquePointer &
        operator =(
                UniquePointer &&
                other
        )
        noexcept
        {
            data = std::move(other.data);
            return * this;
        }

        /**
         * @brief
         * @return
         */
        T &
        operator () ()
        {
            return * data;
        }

        /**
         * @brief
         * @return
         */
        T const &
        operator () ()
        const
        {
            return * data;
        }

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
            return std::make_unique<T>(t);
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
            return std::make_unique<T>(t);
        }

        /**
         * @brief
         * @param other
         * @return
         */
        static
        Data
        copyData(
                UniquePointer const &
                other
        )
        {
            Data data_arg = nullptr;
            if (other.data != nullptr) {
                data_arg = std::make_unique<T>(other.get());
            }
            return data_arg;
        }

        /**
         * @brief
         * @param other
         * @return
         */
        static
        Data
        moveData(
                UniquePointer const &
                other
        )
        {
            Data data_arg = nullptr;
            if (other.data != nullptr) {
                data_arg = std::move(other.data);
            }
            return data_arg;
        }

    public:

        /**
         * @brief
         */
        Data data;

    };

}

#endif //FETA_UNIQUE_POINTER_HXX
