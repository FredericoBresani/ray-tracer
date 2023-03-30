#ifndef __POINTS__
#define __POINTS__

#include "Vectors.h"
#include "Matrices.h"

template <typename T> class Point2
{
    public:
        T x, y;
        Point2(): x(T(0)), y(T(0)) {}
        Point2(T _x): x(_x), y(_x) {}
        Point2(T _x, T _y): x(_x), y(_y) {} 
        Point2<T> operator + (const Vec2<T> &v) const { return Point2<T>(x + v.x, y + v.y); }
        Point2<T> operator / (const T& t) const { return Point2<T>(x/t, y/t); }
        Point2<T> operator * (const T& t) const { return Point2<T>(x*t, y*t); }
        Point2<T> operator - (const T& t) const { return Point2<T>(x - t, y - t); }
};

template <typename T> class Point3
{
    public:
        T x, y, z;
        Point3(): x(T(0)), y(T(0)), z(T(0)) {}
        Point3(T _x): x(T(_x)), y(T(_x)), z(T(_x)) {}
        Point3(T _x, T _y, T _z): x(T(_x)), y(T(_y)), z(T(_z)) {}
        Point3<T> operator + (const Vec3<T>& v) const { return Point3<T>(x + v.x, y + v.y, z + v.z); }
        Point3<T> operator - (const Vec3<T>& v) const { return Point3<T>(x - v.x, y - v.y, z - v.z); }
        Point3<T> operator * (const T& t) const { return Point3<T>(x*t, y*t, z*t); }
        Point3<T> operator / (const T& t) const { return Point3<T>(x/t, y/t, z/t); }
        T operator ^ (const Point3<T>& p) const { return (x*p.x) + (y*p.y) + (z*p.z); }
        Point3<T> operator * (const Matrix4<T>& m) const
        {
            T _x, _y, _z;
            _x = x*m.matrix[0][0] + y*m.matrix[0][1] + z*m.matrix[0][2] + m.matrix[0][3]; 
            _y = x*m.matrix[1][0] + y*m.matrix[1][1] + z*m.matrix[1][2] + m.matrix[1][3];
            _z = x*m.matrix[2][0] + y*m.matrix[2][1] + z*m.matrix[2][2] + m.matrix[2][3];
            return Point3<T>(_x, _y, _z); 
        }
        Vec3<T> operator - (const Point3<T>& p) const { return Vec3<T>(x - p.x, y - p.y, z - p.z); }
        T distance(const Point3<T>& p) const
        {
            Vec3<T> v = Vec3<T>(x - p.x, y - p.y, z - p.z);
            return v.norma(v);
        }
};

typedef Point3<double> Point3D;
typedef Point2<double> Point2D;
typedef Point3<int> Point3I;

#endif