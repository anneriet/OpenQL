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

cparam::cparam(utils::Str typeStr) : typeStr(typeStr){
    name = paramnamerand();
    QL_DOUT("Param of typeStr: "<< typeStr<< " name: " << name);
}
cparam::cparam(utils::Str typeStr, utils::Str name)  : typeStr(typeStr) , name(name){
    QL_DOUT("Param of typeStr: "<< typeStr<< " name: " << name);
}
// cparam::cparam(utils::Str typeStr, utils::Str name, utils::UInt value){
//     typeStr_ = typeStr;
//     name = name;
// }
cparam::cparam(utils::Str typeStr, utils::Str name, utils::Real value) : typeStr(typeStr) , name(name){
    QL_DOUT("Param of typeStr: "<< typeStr<< " name: " << name);
}
cparam::cparam(utils::Str typeStr, utils::Str name, utils::Complex value) : typeStr(typeStr) , name(name){
    QL_DOUT("Param of typeStr: "<< typeStr<< " name: " << name);
}
// cparam::cparam(utils::Str typeStr, utils::UInt value) : typeStr(typeStr){
// }
cparam::cparam(utils::Str typeStr, utils::Real value) : typeStr(typeStr){
    name = paramnamerand();
    QL_DOUT("Param of typeStr: "<< typeStr<< " name: " << name);
}
cparam::cparam(utils::Str typeStr, utils::Complex value) : typeStr(typeStr){
    name = paramnamerand();
    QL_DOUT("Param of typeStr: "<< typeStr<< " name: " << name);
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
    QL_DOUT("Parameter: |- type     : " << typeStr << "\n\t |- name     : " << name);
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



PInt::PInt() : ParamTemp(){
    name = PInt::paramnamerand();
    QL_DOUT("Parameter: |- type     : " << "Int" << "\n\t |- name     : " << name);
 }

PInt::PInt(const utils::UInt &value) : ParamTemp(){ };

utils::Str PInt::toString()  {
    std::stringstream strm;
    strm << value;
    return strm.str();
};

// Return a random parametername.
utils::Str ParamTemp::paramnamerand() {
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