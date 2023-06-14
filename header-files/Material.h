#ifndef __MATERIAL__
#define __MATERIAL__

#include "RGBColor.h"
#include "HitInfo.h"

class Material {
    public:
        RGBColor color;
        double difuseK, specularK, ambientalK, reflectiveK, transmissionK, roughK;
        Material(RGBColor c, double d, double s, double a, double r, double t, double p): 
            color(c), difuseK(d), specularK(s), ambientalK(a), reflectiveK(r), transmissionK(t), roughK(p) {}
        ~Material() {}
        RGBColor shade(HitInfo &hit);
        RGBColor area_light_shade(HitInfo &hit);
        RGBColor path_shade(HitInfo &hit);
};


#endif