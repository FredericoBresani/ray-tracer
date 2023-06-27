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
        void makeCamera();
        void transformCamera(const Matrix4D &matrix);
        Point3D worldToCameraCoordinates(const Point3D &point, const Point3D &cameraPosition);
        Point2D worldToScreenCoordinates(const Point3D &point, const Point3D &cameraLocation);
        Point3D getPos();
        virtual void set_sampler() = 0;
        virtual void render(std::vector<Object*> objetos, std::vector<Light*>& lights, Ambient& ambient) = 0;
    protected:
        Sampler *sampler_ptr;
        double distance, pixelSize, pixelQtnH, pixelQtnV, focalPlaneDistance, fishAngle;
        int hr, vr, samples;
        Vec3D up, u, v, w, right, iup;
        Point3D cameraPos, lookAt;
};

Point3D Camera::getPos()
{
    return cameraPos;
}


void Camera::makeCamera()
{
    // Tha camera base should follow this order {w = z, v = y, u = x}
    pixelQtnH = (double)hr/pixelSize;
    pixelQtnV = (double)vr/pixelSize;
    Vec3D toScreen = Vec3D::normalize(lookAt - cameraPos);
    w = toScreen;
    Vec3D WUP = w ^ up;
    if (WUP.x == 0.0 && WUP.y == 0.0 && WUP.x == 0.0)
    {
        up = Vec3D(1.0, 0.0, 0.0);
    }
    v = Vec3D::normalize(up - (w*((up*w)/(w*w))));
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


void Camera::transformCamera(const Matrix4D &matrix) 
{
    cameraPos = cameraPos*matrix;
}

Point3D Camera::worldToCameraCoordinates(const Point3D &point, const Point3D &cameraPosition)
{
    Vec3D worldPoint = point - cameraPos; // litle cheating here, using a vector as a point
    double x = (w.x*(worldPoint.x)) + (w.y*(worldPoint.y)) + (w.z*(worldPoint.z));
    double y = (v.x*(worldPoint.x)) + (v.y*(worldPoint.y)) + (v.z*(worldPoint.z));
    double z = (u.x*(worldPoint.x)) + (u.y*(worldPoint.y)) + (u.z*(worldPoint.z));
    return Point3D(z, y, x); // Remenber that the base order is like that
}

Point2D Camera::worldToScreenCoordinates(const Point3D &point, const Point3D &cameraLocation) //this is a point on the screen already, althoug in world coordinates
{                                               //that is why no projection is needed
    Point3D cameraCoordinates = Camera::worldToCameraCoordinates(point, cameraLocation);
    return Point2D(cameraCoordinates.x, cameraCoordinates.y);
}

void Camera::set_sampler()
{
    sampler_ptr = new JitteredSampler(samples);
}


#endif