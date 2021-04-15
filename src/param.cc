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

cparam::cparam(const utils::Str typeStr) : typeStr(typeStr){
    name = paramnamerand();
    QL_DOUT("Param of typeStr: "<< typeStr<< " name: " << name);
}
cparam::cparam(const utils::Str typeStr, const utils::Str name)  : typeStr(typeStr) , name(name){
    QL_DOUT("Param of typeStr: "<< typeStr<< " name: " << name);
}
// cparam::cparam(utils::Str typeStr, utils::Str name, utils::UInt value){
//     typeStr_ = typeStr;
//     name = name;
// }
cparam::cparam(const utils::Str typeStr, const utils::Str name, const utils::Real value) : typeStr(typeStr) , name(name), assigned(true){
    QL_DOUT("Param of typeStr: "<< typeStr<< " name: " << name << " value: " << value);

    if(typeStr == "INT")
    {
        int_value = (int) value;
    }
    else
    {
        real_value = value;
    }
}
cparam::cparam(const utils::Str typeStr, const utils::Str name, const utils::Complex value) : typeStr(typeStr) , name(name), assigned(true){
       QL_DOUT("Param of typeStr: "<< typeStr<< " name: " << name << " value: " << value);
    complex_value = value;
}
// cparam::cparam(utils::Str typeStr, utils::UInt value) : typeStr(typeStr){
// }
cparam::cparam(const utils::Str typeStr, const utils::Real value) : typeStr(typeStr), assigned(true){
    name = paramnamerand();
    if(typeStr == "INT")
    {
    QL_DOUT("Param of typeStr: "<< typeStr<< " (Int) name: " << name << " value: " << value);
        int_value = (int) value;
    }
    else
    {
        real_value = value;
    }
}
cparam::cparam(const utils::Str typeStr, const utils::Complex value) : typeStr(typeStr), assigned(true){
    name = paramnamerand();
    QL_DOUT("Param of typeStr: "<< typeStr<< " name: " << name << " value: " << value);
    complex_value = value;
}

utils::UInt id;
parameter_type_t type_;
utils::Str cparam::qasm() const
{
    print();
    return name;
};


parameter_type_t type()
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
// template<typename Type>
// class Param{
//     protected:
//         Type value;
//         utils::Str name;
//     public:
//     // Type Param(){}
//     Type Param(const Type& initValue) : value(initValue){ }
//     // Type Param(const utils::Str name, const Type& initValue) : value(initValue), name(name){ }
//     virtual utils::Str toString()const;
// };



// PInt::PInt() : ParamTemp(){
//     name = PInt::paramnamerand();
//     QL_DOUT("Parameter: |- type     : " << "Int" << "\n\t |- name     : " << name);
//  }

// PInt::PInt(const utils::UInt &value) : ParamTemp(){ };

// utils::Str PInt::toString()  {
//     std::stringstream strm;
//     strm << value;
//     return strm.str();
// };

// // Return a random parametername.
// utils::Str ParamTemp::paramnamerand() {
//       static char randomid[9];
//       for (unsigned int k = 0; k<8; ++k) {
//         const int v = (int)std::rand()%3;
//         randomid[k] = (char)(v==0?('0' + ((int)std::rand()%10)):
//                              (v==1?('a' + ((int)std::rand()%26)):
//                               ('A' + ((int)std::rand()%26))));
//       }
//       return randomid;
//     }

} // namespace ql