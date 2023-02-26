#ifndef __SPHERE__
#define __SPHERE__

#include "Vectors.h"
#include "Points.h"
#include "Definitions.h"
#include "GeometricObject.h"
#include "HitInfo.h"
#include <math.h>


class Sphere: public Object 
{
    public:
        Point3D center;
        Vec3D color;
        double radius, difuseK, specularK, ambientK, reflectionK, transmissionK, phongExp;
        Sphere(const Point3D &c, const Vec3D &RGB, double r, double difuse, double specular, double ambient, double reflection, double transmission, double phong): center(c), color(RGB), radius(r), difuseK(difuse), specularK(specular), ambientK(ambient), reflectionK(reflection), transmissionK(transmission), phongExp(phong) {}
        ~Sphere() {}
        bool rayObjectIntersect(const Ray &ray, double *tmin, const HitInfo& info) const
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
            Vec3D normal = hit - this->center; 
            return normal.normalize(normal);
        }
};

#endif