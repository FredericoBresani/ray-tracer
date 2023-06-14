#ifndef __SHADER__
#define __SHADER__

#include "RGBColor.h"

class Shader {
    public:
        Shader() {}
        ~Shader() {}
        virtual RGBColor shade() = 0;
};


#endif