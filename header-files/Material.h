#ifndef __MATERIAL__
#define __MATERIAL__

#include "RGBColor.h"
#include "HitInfo.h"

class Material {
    public:
        RGBColor color;
        double difuseK, specularK, ambientalK, reflectiveK, transmissionK, roughK;
        bool getShadows;
        Material(RGBColor c, double d, double s, double a, double r, double t, double p, bool ms): 
            color(c), difuseK(d), specularK(s), ambientalK(a), reflectiveK(r), transmissionK(t), roughK(p), getShadows(ms) {}
        ~Material() {}
        RGBColor shade(HitInfo &hit);
        RGBColor area_light_shade(HitInfo &hit);
        RGBColor path_shade(HitInfo &hit);
};


#endif