/** \file
 * Parameters implementation
 */

#pragma once

#include "utils/num.h"
#include "utils/str.h"
#include "utils/vec.h"


namespace ql {

// enum class parameter_type_t {
//     PBOOL, PINT, PREAL, PCOMPLEX
// };

/**
 * cparam_
 */
class cparam {
public:
    cparam(utils::Str typeStr);
    cparam(utils::Str typeStrStr, utils::Str name);
    cparam(utils::Str typeStrStr, utils::Str name, utils::UInt value);
    cparam(utils::Str typeStr, utils::Str name, utils::Real value);
    cparam(utils::Str typeStr, utils::Str name, utils::Complex value);
    cparam(utils::Str typeStr, utils::UInt value);
    cparam(utils::Str typeStr, utils::Real value);
    cparam(utils::Str typeStr, utils::Complex value);
    utils::UInt id;
    // parameter_type_t type() const;
    void print() const;
    utils::Str parameter_name;

private:
    utils::Str typeStr_;
};
} // namespace ql