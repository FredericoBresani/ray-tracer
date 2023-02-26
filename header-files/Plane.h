#ifndef __PLANE__
#define __PLANE__

#include "GeometricObject.h"
#include "Points.h"
#include "Vectors.h"
#include "Definitions.h"


class Plane: public Object 
{
    public:
        Vec3D normal;
        Point3D pp;
        Vec3D color;
        double difuseK, specularK, ambientK, reflectionK, transmissionK, phongExp;
        Plane(const Vec3D &n, const Point3D &p, const Vec3D &RGB, double difuse, double specular, double ambient, double reflection, double transmission, double phong): normal(n), pp(p), color(RGB), difuseK(difuse), specularK(specular), ambientK(ambient), reflectionK(reflection), transmissionK(transmission), phongExp(phong) {}
        ~Plane() {}
        bool rayObjectIntersect(const Ray &ray, double *tmin, const HitInfo& info) const 
        {
            double t = ((pp - ray.origin) * this->normal) / (ray.direction * this->normal);
            Point3D location = ray.origin + ray.direction*t;
            if (t > kEpsilon && t < (*tmin))
            {
                (*tmin) = t;
                return true;

            } else {
                return false;
            }
        }
        Vec3D getColor() const
        {
            return this->color;
        }
        double getKd() const
        {
            return this->difuseK;
        }
        double getKs() const
        {
            return this->specularK;
        }
        double getKa() const
        {
            return this->ambientK;
        }
        double getKr() const
        {
            return this->reflectionK;
        }
        double getKt() const
        {
            return this->transmissionK;
        }
        double getPhongExp() const
        {
            return this->phongExp;
        }
        Vec3D getNormal(const Point3D &hit, const Ray &ray) const
        {
            return normal.normalize(normal);
        }
};

#endif