/** \file
 * Parameters implementation
 */

#pragma once

#include "utils/num.h"
#include "utils/str.h"
#include "utils/vec.h"


namespace ql {

enum class parameter_type_t {
    PBOOL, PINT, PREAL, PCOMPLEX
};

/**
 * cparam_
 */
class cparam {
public:
    cparam(utils::Str typeStr);
    cparam(utils::Str typeStr, utils::Str name);
    // cparam(utils::Str typeStr, utils::Str name, utils::UInt value);
    cparam(utils::Str typeStr, utils::Str name, utils::Real value);
    cparam(utils::Str typeStr, utils::Str name, utils::Complex value);
    // cparam(utils::Str typeStr, utils::UInt value);
    cparam(utils::Str typeStr, utils::Real value);
    cparam(utils::Str typeStr, utils::Complex value);

    utils::Str name;
    utils::Int bool_value = 0;
    utils::UInt int_value = 0;
    utils::Real real_value = 0.0;  

    utils::Str typeStr;
    parameter_type_t type();
    utils::Str qasm() const;
    void print() const;
private:
    const utils::Real rand(const utils::Real val_max, utils::UInt * p_rng);
    const char* paramnamerand(); 
};

template<typename Type>
class ParamTemp{
    protected:
        Type value;
        utils::Str name;
        utils::Str paramnamerand();
    public:
    ParamTemp(){ };
    ParamTemp(const Type& initValue) :  value(initValue){ };
    // Type Param(const utils::Str name, const Type& initValue);
    // virtual utils::Str toString() const;
};

class PInt : public ParamTemp<utils::UInt>{
//     private:
//         // utils::Str paramnamerand();
//         // utils::Str name;
//         // utils::UInt value;
    public:
        PInt();
        PInt(const utils::UInt &value);

        utils::Str toString();
};
} // namespace ql