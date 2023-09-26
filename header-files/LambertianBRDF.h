#ifndef __LAMBERTIAN__
#define __LAMBERTIAN__

#include "BRDF.h"
#include "RGBColor.h"
#include "Vectors.h"
#include "HitInfo.h"
#include "Definitions.h"

class LambertianBRDF: public BRDF {
    public:
        LambertianBRDF() {}
        ~LambertianBRDF() {}
        virtual RGBColor f(const HitInfo &hit, const Vec3D &l, const Vec3D &r) const = 0;
        virtual RGBColor sample_f(const HitInfo &hit, Vec3D &l, const Vec3D &r) const = 0;
        virtual RGBColor rho(const HitInfo &hit, const Vec3D &r) const = 0;
    private:
        double kd;
        RGBColor dc;

};

RGBColor LambertianBRDF::f(const HitInfo &hit, const Vec3D &l, const Vec3D &r) const 
{
    return (dc * kd)/M_PI;
}

RGBColor LambertianBRDF::rho(const HitInfo &hit, const Vec3D &r) const 
{
    return (dc * kd);
}

#endif