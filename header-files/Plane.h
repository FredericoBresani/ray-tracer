#ifndef __PLANE__
#define __PLANE__

#include "GeometricObject.h"
#include "Points.h"
#include "Vectors.h"
#include "Definitions.h"
#include "RGBColor.h"
#include "Material.h"


class Plane: public Object 
{
    public:
        Vec3D normal;
        Point3D pp;
        Material *material;
        Plane(const Vec3D &n, const Point3D &p, Material *m): normal(n), pp(p), material(m) {}
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
        RGBColor getColor() const
        {
            return this->material->color;
        }
        double getKd() const
        {
            return this->material->difuseK;
        }
        double getKs() const
        {
            return this->material->specularK;
        }
        double getKa() const
        {
            return this->material->ambientalK;
        }
        double getKr() const
        {
            return this->material->reflectiveK;
        }
        double getKt() const
        {
            return this->material->transmissionK;
        }
        double getPhongExp() const
        {
            return this->material->roughK;
        }
        Vec3D getNormal(const Point3D &hit, const Ray &ray) const
        {
            return normal.normalize(normal);
        }
};

#endif