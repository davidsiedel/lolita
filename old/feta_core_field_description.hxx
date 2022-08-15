//
// Created by dsiedel on 23/03/2022.
//

#ifndef FETA_FETA_CORE_FIELD_DESCRIPTION_HXX
#define FETA_FETA_CORE_FIELD_DESCRIPTION_HXX

#include "new/_feta.hxx"
#include "new/_feta_matrix.hxx"
#include "new/_feta_tensor.hxx"
#include "new/feta_core_mapping_description.hxx"
//#include "feta/core/00_tenor.hxx"

namespace feta::core
{

    enum struct Field
    {

        Scalar,
        Vector,
        Tensor,

    };

    struct FieldComponent
    {

        constexpr
        FieldComponent()
        :
        row(0),
        col(0)
        {}

        constexpr
        FieldComponent(
                auto
                row_arg,
                auto
                col_arg
        )
        :
        row(row_arg),
        col(col_arg)
        {}

        Indx row;

        Indx col;

    };

//    auto const static constexpr fld_scalar = Field::Scalar;
//    auto const static constexpr fld_vector = Field::Vector;

    struct TensorDescription
    {

        constexpr
        TensorDescription()
        :
        rows(-1),
        cols(-1),
        size(-1),
        storage_option(matrix::StorageOption::RowMajor)
        {}

        constexpr explicit
        TensorDescription(
                Intg
                rows_arg,
                matrix::StorageOption
                storage_option_arg = matrix::StorageOption::RowMajor
        )
        :
        rows(rows_arg),
        cols(1),
        size(rows_arg),
        storage_option(storage_option_arg)
        {
            assert(rows_arg >= -1);
        }

        constexpr
        TensorDescription(
                Intg
                rows_arg,
                Intg
                cols_arg,
                matrix::StorageOption
                storage_option_arg = matrix::StorageOption::RowMajor
        )
        :
        rows(rows_arg),
        cols(cols_arg),
        size(rows_arg * cols_arg),
        storage_option(storage_option_arg)
        {
            assert(rows_arg >= -1);
            assert(cols_arg >= -1);
        }

        constexpr
        Bool
        operator==(
                TensorDescription const &
                other
        )
        const
        {
            auto eq_0 = rows == other.rows;
            auto eq_1 = cols == other.cols;
            auto eq_2 = size == other.size;
            auto eq_3 = storage_option == other.storage_option;
            return eq_0 && eq_1 && eq_2 && eq_3;
        }

        constexpr
        Bool
        operator!=(
                TensorDescription const &
                other
        )
        const
        {
            return !(other == * this);
        }

        Intg rows;

        Intg cols;

        Intg size;

        matrix::StorageOption storage_option;

    };

    struct FieldDescription
    {

        constexpr
        FieldDescription()
        :
        field(Field::Scalar),
        dim_field(0),
        tensor_description()
        {}

        constexpr inline
        FieldDescription(
                Field
                field_arg,
                Indx
                dim_field_arg
        )
        :
        field(field_arg),
        dim_field(dim_field_arg),
        tensor_description(setTensorDescription(field_arg, dim_field_arg))
        {}

        constexpr
        Bool
        operator==(
                FieldDescription const &
                other
        )
        const
        {
            auto eq_0 = field == other.field;
            auto eq_1 = dim_field == other.dim_field;
            auto eq_2 = tensor_description == other.tensor_description;
            return eq_0 && eq_1 && eq_2;
        }

        constexpr
        Bool
        operator!=(
                FieldDescription const &
                other
        )
        const
        {
            return !(other == * this);
        }

    private:

        static constexpr inline
        TensorDescription
        setTensorDescription(
                Field
                field_arg,
                Indx
                dim_field_arg
        )
        {
            if (field_arg == Field::Scalar) {
                return TensorDescription(1);
            }
            else if (field_arg == Field::Vector) {
                return TensorDescription(dim_field_arg);
            }
            else if (field_arg == Field::Tensor) {
                return TensorDescription(dim_field_arg, dim_field_arg);
            }
            else {
                return TensorDescription();
            }
        }

    public:

        static constexpr inline
        FieldDescription
        fromMapping(
                FieldDescription const &
                field_description_arg,
                Mapping
                operator_type_arg
        )
        {
            if (operator_type_arg == Mapping::Gradient || operator_type_arg == Mapping::SymmetricGradient) {
                if (field_description_arg.field == Field::Scalar) {
                    return FieldDescription(Field::Vector, field_description_arg.dim_field);
                }
                else if (field_description_arg.field == Field::Vector) {
                    return FieldDescription(Field::Tensor, field_description_arg.dim_field);
                }
                else {
                    assert(false);
                }
            }
            else if (operator_type_arg == Mapping::LargeStrain || operator_type_arg == Mapping::SmallStrain) {
                if (field_description_arg.field == Field::Scalar) {
                    return FieldDescription(Field::Vector, 3);
                }
                else if (field_description_arg.field == Field::Vector) {
                    return FieldDescription(Field::Tensor, 3);
                }
                else {
                    assert(false);
                }
            }
            else {
                return field_description_arg;
            }
        }

        Field field;

        Indx dim_field;

        TensorDescription tensor_description;

    };

}

#endif //FETA_FETA_CORE_FIELD_DESCRIPTION_HXX
