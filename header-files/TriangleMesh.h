#ifndef __TMESH__
#define __TMESH__

#include <vector>
#include "GeometricObject.h"
#include "RGBColor.h"
#include "Material.h"
#include "Plane.h"
#include "Points.h"
#include "Ray.h"
#include "HitInfo.h"
#include "Vectors.h"

class TriangleMesh: public Object {
    public:
        Material *material;
        bool castShadows;
        int nTriangles, nVertices, triangleIndice = 0;
        std::vector<Point3D> vertices;
        std::vector<Point3I> triangles;
        std::vector<Vec3D> triangleNormals;
        std::vector<Vec3D> verticesNormals;
        TriangleMesh(int n, int v, Material *m, bool s): nTriangles(n), nVertices(n), material(m), castShadows(s) {}
        ~TriangleMesh() {}
        bool rayObjectIntersect(const Ray &ray, double *tmin, HitInfo &info) 
        {
            double min = infinity;
            bool hit = false;
            for (int i = 0; i < triangles.size(); i++)
            {
                Point3D A, B, C;
                A = vertices[triangles[i].x];
                B = vertices[triangles[i].y];
                C = vertices[triangles[i].z]; 
                Vec3D tPlaneNormal = (A - B) ^ (A - C);
                if (((ray.direction * (-1)) * tPlaneNormal) <= 0)
                {
                    tPlaneNormal = (A - C) ^ (A - B);
                }
                Material *tempMaterial = new Material{
                    color, 0, 0, 0, 0, 0, 0, false
                };
                Plane *tPlane = new Plane(tPlaneNormal, A, tempMaterial, true);
                Point3D pHit;
                if (tPlane->rayObjectIntersect(ray, tmin, info)) // to-do: the intersecion fails when 3 respective coordinates on
                {   // diferent points, equals to 0
                    free(tPlane);
                    pHit = ray.origin + ray.direction*(*tmin);      
                    Vec3D temp;                     // _      _  _     _
                    Vec3D v0 = Vec3D(A.x, B.x, C.x);//|a1 b1 c1||a|   |X|
                    Vec3D v1 = Vec3D(A.y, B.y, C.y);//|a2 b2 c2||b| = |Y|
                    Vec3D v2 = Vec3D(A.z, B.z, C.z);//|a3 b3 c3||g|   |Z|
                    double point[] = { pHit.x, pHit.y, pHit.z };
                    double tempP;
                    if (v1.x != 0.0 && v0.x == 0.0) // try to escalonate, reordering the vectors and values to avoid 0's
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

                    if (v2.y != 0.0 && v1.y == 0.0)
                    {
                        temp = v1;
                        v1 = v2;
                        v2 = temp;
                        tempP = point[1];
                        point[1] = point[2];
                        point[2] = tempP;
                    }

                    if (v0.x != 0.0)
                    {
                        point[1] = point[1] + (point[0]*(-1*(v1.x/v0.x)));
                        point[2] = point[2] + (point[0]*(-1*(v2.x/v0.x)));
                        v1 = v1 + v0*(-1*(v1.x/v0.x));
                        v2 = v2 + v0*(-1*(v2.x/v0.x));
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
                    if (v2.z == 0.0 && (point[2] < -kEpsilon || point[2] > kEpsilon)) 
                    {
                        (*tmin) = infinity;
                    }
                    if (v2.z != 0.0)
                    {
                        gama = point[2]/v2.z; 
                    }
                    if (v1.y == 0.0 && (point[1] < -kEpsilon || point[1] > kEpsilon))
                    { 
                        (*tmin) = infinity;
                    }
                    if (v1.y != 0.0)
                    {
                        beta = point[1]/v1.y;
                    }
                    if (v0.x == 0.0 && (point[0] < -kEpsilon || point[0] > kEpsilon))
                    {
                        (*tmin) = infinity;
                    }
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
                        if (alpha > 1.0 || beta > 1.0 || gama > 1.0) 
                        {
                            (*tmin) = infinity;
                        }
                        else if (alpha < 0.0 || beta < 0.0 || gama < 0.0)
                        {
                            (*tmin) = infinity;
                        } else {
                            hit = true;
                            if ((*tmin) < min)
                            {
                                info.hit_object = true;
                                this->triangleIndice = i;
                                min = (*tmin);
                            }
                        }
                    } else {
                        (*tmin) = infinity;
                    }
                } else {
                    (*tmin) = infinity;
                }
            }
            if (hit)
            {
                (*tmin) = min;
                return true;
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
            Point3D A, B, C;
            A = vertices[triangles[triangleIndice].x];
            B = vertices[triangles[triangleIndice].y];
            C = vertices[triangles[triangleIndice].z];
            Vec3D tPlaneNormal = (A - B) ^ (A - C);
            if (((ray.direction * (-1)) * tPlaneNormal) <= 0)
            {
                tPlaneNormal = (A - C) ^ (A - B);
            }
            return Vec3D::normalize(tPlaneNormal);
        }
        bool getShadows() const
        {
            return this->material->getShadows;
        }
        bool getCastShadows() const
        {
            return this->castShadows;
        }
};



#endif