#ifndef __MATERIAL__
#define __MATERIAL__

#include "RGBColor.h"

typedef struct ObjectMaterial {
    RGBColor color;
    double difuseK, specularK, ambientalK, reflectiveK, transmissionK, roughK;
} Material;


#endif