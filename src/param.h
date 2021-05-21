/** \file
 * Parameters implementation
 */

#pragma once

#include "utils/num.h"
#include "utils/str.h"
#include "utils/vec.h"


namespace ql {

enum class parameter_type_t {
    PERROR, PINT, PREAL, PCOMPLEX, PBOOL
};


/**
 * cparam_
 */
class cparam {
public:
    cparam();
    cparam(const utils::Str typeStr);
    cparam(const utils::Str typeStr, const utils::Str name);
    // cparam(utils::Str typeStr, utils::Str name, utils::UInt value);
    cparam(const utils::Str typeStr, const utils::Str name, const utils::Real value);
    cparam(const utils::Str typeStr, const utils::Str name, const utils::Complex value);
    // cparam(utils::Str typeStr, utils::UInt value);
    cparam(const utils::Str typeStr, const utils::Real value);
    cparam(const utils::Str typeStr, const utils::Complex value);
    utils::Str typeStr;
    utils::Str name;
    utils::Complex value; 
    utils::Bool assigned = false;

    utils::Int bool_value;
    utils::UInt int_value;
    utils::Real real_value; 
    utils::Complex complex_value; 


    parameter_type_t type();
    utils::Str qasm() const;
    void print() const;
private:
    const utils::Real rand(const utils::Real val_max, utils::UInt * p_rng);
    const char* paramnamerand(); 
    parameter_type_t type_;
    parameter_type_t set_type(const utils::Str typeStr);
};
} // namespace ql