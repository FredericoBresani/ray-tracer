#ifndef __GEOMETRICOBJECT__
#define __GEOMETRICOBJECT__

#include "Vectors.h"
#include "Ray.h"
#include "HitInfo.h"
#include "RGBColor.h"
#include "Material.h"


class Object 
{
    public:
        RGBColor color;
        Material *material;
        bool castShadows;
        Object() {}
        virtual ~Object() {}
        virtual bool rayObjectIntersect(const Ray &ray, double *tmin, HitInfo& info) = 0;
        virtual RGBColor getColor() const = 0;
        virtual double getKd() const = 0;
        virtual double getKs() const = 0;
        virtual double getKa() const = 0;
        virtual double getKr() const = 0;
        virtual double getKt() const = 0;
        virtual double getIor() const = 0;
        virtual double getPhongExp() const = 0;
        virtual Vec3D getNormal(const Point3D &hit, const Ray &ray) const = 0;
        virtual bool getShadows() const = 0;
        virtual bool getCastShadows() const = 0;
};

#endif