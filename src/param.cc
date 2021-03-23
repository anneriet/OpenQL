/** \file
 * Parameters implementation
 */

#include "param.h"

namespace ql {

/**
 * @brief   (Classical) parameters constructor
 * @param   typeStr    typeStr of the parameter
 */
cparam::cparam(utils::Str typeStr){
    typeStr_ = typeStr;
    
}
cparam::cparam(utils::Str typeStr, utils::Str name){
    typeStr_ = typeStr;
    parameter_name = name;
}
cparam::cparam(utils::Str typeStr, utils::Str name, utils::UInt value){
    typeStr_ = typeStr;
    parameter_name = name;
}
cparam::cparam(utils::Str typeStr, utils::Str name, utils::Real value){
    typeStr_ = typeStr;
    parameter_name = name;
}
cparam::cparam(utils::Str typeStr, utils::Str name, utils::Complex value){
    typeStr_ = typeStr;
    parameter_name = name;
}
cparam::cparam(utils::Str typeStr, utils::UInt value){
    typeStr_ = typeStr;
}
cparam::cparam(utils::Str typeStr, utils::Real value){
    typeStr_ = typeStr;
}
cparam::cparam(utils::Str typeStr, utils::Complex value){
    typeStr_ = typeStr;
}

utils::UInt id;
// parameter_typeStr_t typeStr();
void print();
utils::Str parameter_name;


} // namespace ql