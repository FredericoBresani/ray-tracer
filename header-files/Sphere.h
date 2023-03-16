#ifndef __SPHERE__
#define __SPHERE__

#include "Vectors.h"
#include "Points.h"
#include "Definitions.h"
#include "GeometricObject.h"
#include "HitInfo.h"
#include "RGBColor.h"
#include "Material.h"
#include <math.h>


class Sphere: public Object 
{
    public:
        Point3D center;
        double radius;
        Material *material;
        Sphere(const Point3D &c, double r, Material *m): center(c), radius(r), material(m) {}
        ~Sphere() {}
        bool rayObjectIntersect(const Ray &ray, double *tmin, const HitInfo& info)
        {
            double a = pow(ray.direction.norma(ray.direction), 2.0);
            double b = ((ray.origin - this->center) * ray.direction) * 2.0;
            double c = ((this->center ^ this->center) + (ray.origin ^ ray.origin)) + (-2.0)*(ray.origin ^ this->center) - (this->radius*this->radius);
            double delta = (b*b) - 4.0*a*c; 
            if (delta == 0.0)
            {
                double t = (b*b)/(2.0*a);
                if (t > kEpsilon && t < (*tmin))
                {
                    (*tmin) = t;
                    return true;
                } else {
                    return false;
                }
            }
            else if (delta > 0.0)
            {
                double sqrtDelta = pow(delta, 0.5);
                double t1 = (((-1)*b) + sqrtDelta)/(2.0*a);
                double t2 = (((-1)*b) - sqrtDelta)/(2.0*a);
                if ((t1 < t2) && (t1 > kEpsilon) && (t1 < (*tmin)))
                {
                    (*tmin) = t1;
                    return true;
                }
                else if (t2 > kEpsilon && t2 < (*tmin))
                {
                    (*tmin) = t2;
                    return true;
                }
                return false;
            }
            return false;
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
            Vec3D normal = hit - this->center; 
            return normal.normalize(normal);
        }
};

#endif