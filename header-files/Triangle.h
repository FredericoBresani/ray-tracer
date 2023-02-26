#ifndef __TRIANGLE__
#define __TRIANGLE__

#include "GeometricObject.h"
#include "Points.h"
#include "Vectors.h"
#include "Definitions.h"
#include "Plane.h"

class Triangle: public Object 
{
    public:
        Point3D A;
        Point3D B;
        Point3D C;
        Vec3D color;
        double difuseK, specularK, ambientK, reflectionK, transmissionK, phongExp;
        Triangle(const Point3D &a, const Point3D &b, const Point3D &c, const Vec3D &RGB, double difuse, double specular, double ambient, double reflection, double transmission, double phong): A(a), B(b), C(c), color(RGB), difuseK(difuse), specularK(specular), ambientK(ambient), reflectionK(reflection), transmissionK(transmission), phongExp(phong) {}
        ~Triangle() {}
        bool rayObjectIntersect(const Ray& ray, double *tmin, const HitInfo& info) const
        {
            Vec3D tPlaneNormal = (this->A - this->B) ^ (this->A - this->C);
            if (((ray.direction * (-1)) * tPlaneNormal) <= 0)
            {
                tPlaneNormal = (this->A - this->C) ^ (this->A - this->B);
            }
            Plane *tPlane = new Plane(tPlaneNormal, this->A, color, 0, 0, 0, 0, 0, 0);
            Point3D pHit;
            if (tPlane->rayObjectIntersect(ray, tmin, info)) 
            {   
                pHit = ray.origin + ray.direction*(*tmin);      
                Vec3D temp;                     // _      _  _     _
                Vec3D v0 = Vec3D(A.x, B.x, C.x);//|a1 b1 c1||a|   |X|
                Vec3D v1 = Vec3D(A.y, B.y, C.y);//|a2 b2 c2||b| = |Y|
                Vec3D v2 = Vec3D(A.z, B.z, C.z);//|a3 b3 c3||g|   |Z|
                double point[] = { pHit.x, pHit.y, pHit.z };
                double tempP;
                if (v1.x != 0.0 && v0.x == 0.0)
                {
                    temp = v0;
                    v0 = v1;
                    v1 = temp;
                    tempP = point[0];
                    point[0] = point[1];
                    point[1] = tempP;
                }
                else if(v2.x != 0.0 && v0.x == 0.0) 
                {
                    temp = v0;
                    v0 = v2;
                    v2 = temp;
                    tempP = point[0];
                    point[0] = point[2];
                    point[2] = tempP;
                }
                if (v0.x != 0.0)
                {
                    point[1] = point[1] + (point[0]*(-1*(v1.x/v0.x)));
                    point[2] = point[2] + (point[0]*(-1*(v2.x/v0.x)));
                    v1 = v1 + v0*(-1*(v1.x/v0.x));
                    v2 = v2 + v0*(-1*(v2.x/v0.x));
                }
                if (v2.y != 0.0 && v1.y == 0.0)
                {
                    temp = v1;
                    v1 = v2;
                    v2 = temp;
                    tempP = point[1];
                    point[1] = point[2];
                    point[2] = tempP;
                }
                else if (v0.y != 0.0 && v1.y == 0.0 && v0.x == 0.0)
                {
                    temp = v1;
                    v1 = v0;
                    v0 = temp;
                    tempP = point[1];
                    point[1] = point[0];
                    point[0] = tempP;
                }
                if (v1.y != 0.0)
                {
                    point[0] = point[0] + (point[1]*(-1*(v0.y/v1.y)));
                    point[2] = point[2] + (point[1]*(-1*(v2.y/v1.y)));
                    v0 = v0 + v1*(-1*(v0.y/v1.y));
                    v2 = v2 + v1*(-1*(v2.y/v1.y));
                }
                if (v0.x == 0.0 && v0.y == 0.0 && v0.z != 0.0)
                {
                    temp = v2;
                    v2 = v0;
                    v0 = temp;
                    tempP = point[2];
                    point[2] = point[0];
                    point[0] = tempP;
                }
                else if (v1.x == 0.0 && v1.y == 0.0 && v1.z != 0.0)
                {
                    temp = v2;
                    v2 = v1;
                    v1 = temp;
                    tempP = point[2];
                    point[2] = point[1];
                    point[1] = tempP;
                }
                if (v2.z != 0.0)
                {
                    point[0] = point[0] + (point[2]*(-1*(v0.z/v2.z)));
                    point[1] = point[1] + (point[2]*(-1*(v1.z/v2.z)));
                    v0 = v0 + v2*(-1*(v0.z/v2.z));
                    v1 = v1 + v2*(-1*(v1.z/v2.z));
                }
                double alpha = 0.0, beta = 0.0, gama = 0.0;
                if (v2.z == 0.0 && (point[2] < -kEpsilon || point[2] > kEpsilon)) return false;
                if (v2.z != 0.0)
                {
                    gama = point[2]/v2.z; 
                }
                if (v1.y == 0.0 && (point[1] < -kEpsilon || point[1] > kEpsilon)) return false;
                if (v1.y != 0.0)
                {
                    beta = point[1]/v1.y;
                }
                if (v0.x == 0.0 && (point[0] < -kEpsilon || point[0] > kEpsilon)) return false;
                if (v0.x != 0.0)
                {
                    alpha = point[0]/v0.x;
                }
                if (alpha == 0.0) alpha = 1.0 - (beta + gama);
                if (beta == 0.0) beta = 1.0 - (alpha + gama);
                if (gama == 0.0) gama = 1.0 - (alpha + beta);
                double ABGsum = alpha + beta + gama;
                if (ABGsum <= 1.0 + kEpsilon && ABGsum >= 1.0 - kEpsilon)
                {   
                    if (alpha > 1.0 || beta > 1.0 || gama > 1.0) return (false);
                    if (alpha < 0.0 || beta < 0.0 || gama < 0.0) return (false);
                    return true;
                } else {
                    return false;
                }
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
            Vec3D tPlaneNormal = (this->A - this->B) ^ (this->A - this->C);
            if (((ray.direction * (-1)) * tPlaneNormal) <= 0)
            {
                tPlaneNormal = (this->A - this->C) ^ (this->A - this->B);
            }
            return tPlaneNormal.normalize(tPlaneNormal);
        }
};

#endif