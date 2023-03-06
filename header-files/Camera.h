#ifndef __CAMERA__
#define __CAMERA__

#include "Vectors.h"
#include "Points.h"
#include "Matrices.h"
#include "RegularSampler.h"
#include "JitteredSampler.h"

class Camera 
{
    public:
        int hr, vr, samples;
        double distance, pixelSize, pixelQtnH, pixelQtnV;
        Vec3D up, u, v, w, right, iup;
        Point3D cameraPos, lookAt;
        Sampler *sampler_ptr;
        Camera(int _hr, int _vr, double d, const Vec3D& _up, const Point3D& pos, const Point3D &_lookAt, float p, int s): hr(_hr), vr(_vr), distance(d), up(_up), cameraPos(pos), lookAt(_lookAt), pixelSize(p), samples(s) {} 
        void makeCamera()
        {
            // Tha camera base should follow this order {w = z, v = y, u = x}
            pixelQtnH = (double)hr/pixelSize;
            pixelQtnV = (double)vr/pixelSize;
            Vec3D toScreen = Vec3D();
            toScreen = toScreen.normalize(lookAt - cameraPos);  
            w = toScreen;
            Vec3D WUP = w ^ up;
            if (WUP.x == 0.0 && WUP.y == 0.0 && WUP.x == 0.0)
            {
                up = Vec3D(1.0, 0.0, 0.0);
            }
            v = up - (w*((up*w)/(w*w)));
            v = v.normalize(v);
            u = v ^ w;
            if (pixelQtnH <= pixelQtnV) {
                right =  u*(2.0/pixelQtnH);
                iup = v*(2.0/pixelQtnH); 
            } else {
                right =  u*(2.0/pixelQtnV);
                iup = v*(2.0/pixelQtnV); 
            }
            this->set_sampler();
            
        }
        ~Camera() {}

        void transformCamera(const Matrix4D &matrix) 
        {
            cameraPos = cameraPos*matrix;
        }

        Point3D worldToCameraCoordinates(const Point3D &point, const Point3D &cameraPosition)
        {
            Vec3D worldPoint = point - cameraPos; // litle cheating here, using a vector as a point
            double x = (w.x*(worldPoint.x)) + (w.y*(worldPoint.y)) + (w.z*(worldPoint.z));
            double y = (v.x*(worldPoint.x)) + (v.y*(worldPoint.y)) + (v.z*(worldPoint.z));
            double z = (u.x*(worldPoint.x)) + (u.y*(worldPoint.y)) + (u.z*(worldPoint.z));
            return Point3D(z, y, x); // Remenber that the base order is like that
        }

        Point2D worldToScreenCoordinates(const Point3D &point, const Point3D &cameraLocation) //this is a point on the screen already, althoug in world coordinates
        {                                               //that is why no projection is needed
            Point3D cameraCoordinates = Camera::worldToCameraCoordinates(point, cameraLocation);
            return Point2D(cameraCoordinates.x, cameraCoordinates.y);
        }

        void set_sampler()
        {
            sampler_ptr = new RegularSampler(samples);
        }
};


#endif