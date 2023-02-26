#ifndef __LINES__
#define __LINES__

#include "GeometricObject.h"
#include "Vectors.h"
#include "Points.h"
#include "Definitions.h"
#include "Ray.h"

class Line: public Object 
{
    public:
        Point3D origin;
        Vec3D direction, normal;
        Vec3D color;
        double difuseK, specularK, ambientK, reflectionK, transmissionK, phongExp;
        Line(const Point3D &o, const Vec3D &d, const Vec3D &RGB): origin(o), direction(d), color(RGB) {}
        ~Line() {}
        bool rayObjectIntersect(const Ray &ray, double *tmin, const HitInfo& hit) const
        { 
            double t = ((ray.origin - origin)*normal)/((normal*ray.direction)*(-1.0));
            Point3D location = ray.origin + (ray.direction*t);
            if (location.x > 7 || location.x < -7) return false;
            if (location.y > 7 || location.y < -7) return false;
            if (location.z > 7 || location.z < -7) return false;
            Vec3D toLoc = location - origin;
            toLoc = toLoc.normalize(toLoc);
            double cos = toLoc*(direction.normalize(direction));
            if (t > kEpsilon && t < (*tmin) && (cos >= -1.0 - 0.000002 && cos <= -1.0 + 0.000002 || cos >= 1.0 - 0.000002 && cos <= 1.0 + 0.000002))
            {
                Vec3D color = this->getColor();
                (*tmin) = t;
                // hit.hit_location = Point3D(color.x, color.y, color.z);
                // hit.normal = direction.normalize(direction);
                return true;
            }
            return false;
        }
        void setNormal(const Point3D& C)
        {
            double closerT = ((origin - C)*direction*(-1.0))/(direction*direction);
            Point3D closerP = origin + (direction*closerT);
            this->normal = C - closerP;
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
            return Vec3D();
        }
};

#endif