/** \file
 * Parameters implementation
 */

#include "param.h"

#include "utils/logger.h"
#include "utils/str.h"

namespace ql {

/**
 * @brief   (Classical) parameters constructor
 * @param   typeStr    typeStr of the parameter
 */

cparam::cparam() : typeStr("")
{}

cparam::cparam(const utils::Str typeStr) : typeStr(typeStr){
    name = paramnamerand();
    set_type(typeStr);
    QL_DOUT("Param of typeStr: "<< typeStr<< " name: " << name);
}
cparam::cparam(const utils::Str typeStr, const utils::Str name)  : typeStr(typeStr) , name(name){
    set_type(typeStr);
    QL_DOUT("Param of typeStr: "<< typeStr<< " name: " << name);
}

cparam::cparam(const utils::Str typeStr, const utils::Str name, const utils::Real value) : typeStr(typeStr), name(name), value(value), assigned(true){
    set_type(typeStr);
    if(type_ == parameter_type_t::PINT)
    {
        int_value = (int) value;
    }
    else if (type_ == parameter_type_t::PREAL)
    {
        real_value = value;
    }
    else
    {
        complex_value = value;
    }
    QL_DOUT("Param of typeStr: "<< typeStr<< " name: " << name << " value: " << value);
}

cparam::cparam(const utils::Str typeStr, const utils::Str name, const utils::Complex value) : typeStr(typeStr) , name(name), value(value), assigned(true){
    complex_value = value;
    QL_DOUT("Param of typeStr: "<< typeStr<< " name: " << name << " value: " << value);
}

cparam::cparam(const utils::Str typeStr, const utils::Real value) : typeStr(typeStr), value(value), assigned(true){
    set_type(typeStr);
    if(type_ == parameter_type_t::PINT)
    {
        int_value = (int) value;
    }
    else if (type_ == parameter_type_t::PREAL)
    {
        real_value = value;
    }
    else
    {
        complex_value = value;
        }
    name = paramnamerand();    
    QL_DOUT("Param of type: "<< typeStr << ", name: " << name << " and with value: " << value);
}
cparam::cparam(const utils::Str typeStr, const utils::Complex value) : typeStr(typeStr), value(value), assigned(true){
    complex_value = value;
    name = paramnamerand();
    QL_DOUT("Param of typeStr: "<< typeStr<< " (complex) name: " << name << " value: " << value);
}


utils::Str cparam::qasm() const
{
    print();
    return name;
};

parameter_type_t cparam::set_type(utils::Str typeStr)
{
    if(typeStr == "INT")
    {
        type_ = parameter_type_t::PINT;
    }
    else if (typeStr == "REAL" || typeStr == "ANGLE")
    {
        type_ = parameter_type_t::PREAL;
    }
    else if (typeStr == "COMPLEX" || typeStr == "STATE")
    {
        type_ = parameter_type_t::PCOMPLEX;
    }
    else
    {
        QL_FATAL("Type " << typeStr << " is not known, please use 'INT', 'REAL' or 'COMPLEX'");
    }
    return type_;
}
parameter_type_t cparam::type()
{
    return type_;
}

void cparam::print() const
{
    QL_DOUT("Parameter:\n\t |- type\t\t: " << typeStr << "\n\t |- name\t\t: " << name << "\n\t |- assigned\t\t: "<< assigned << "\n\t |- int value\t\t: " << int_value);
}

// Return a random parametername.
const char* cparam::paramnamerand() {
      static char randomid[9];
      for (unsigned int k = 0; k<8; ++k) {
        const int v = (int)std::rand()%3;
        randomid[k] = (char)(v==0?('0' + ((int)std::rand()%10)):
                             (v==1?('a' + ((int)std::rand()%26)):
                              ('A' + ((int)std::rand()%26))));
      }
      return randomid;
    }

} // namespace ql