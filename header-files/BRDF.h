#ifndef __BRDF__
#define __BRDF__


#include "RGBColor.h"
#include "Vectors.h"
#include "HitInfo.h"
#include "Sampler.h"

class BRDF {
    public:
        virtual RGBColor f(const HitInfo &hit, const Vec3D &l, const Vec3D &r) const = 0;
        virtual RGBColor sample_f(const HitInfo &hit, Vec3D &l, const Vec3D &r) const = 0;
        virtual RGBColor rho(const HitInfo &hit, const Vec3D &r) const = 0;
    protected:
        Sampler *sample_p;
};


#endif