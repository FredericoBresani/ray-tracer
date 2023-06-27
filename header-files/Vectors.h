#ifndef __VECTORS__
#define __VECTORS__

#include "Matrices.h"

template <typename T> class Vec2
{
    public:
        T x, y;
        Vec2(): x(T(0)), y(T(0)) {}
        Vec2(T _x): x(T(_x)), y(T(_x)) {}
        Vec2(T _x, T _y): x(T(_x)), y(T(_y)) {}
        Vec2<T> operator + (const Vec2<T>& v) const { return Vec2<T>(x + v.x, y + v.y); }
        Vec2<T> operator - (const Vec2<T>& v) const { return Vec2<T>(x - v.x, y - v.y); }
        Vec2<T> operator * (const T& t) const { return Vec2<T>(x*t, y*t); }
        Vec2<T> operator / (const T& t) const { return Vec2<T>(x/t, y/t); }
        T operator * (const Vec2<T>& v) const { return (x*v.x) + (y*v.y); }
        static T norma(const Vec2<T>& v)
        {
            return pow((v.x*v.x) + (v.y*v.y), 0.5);
        }
        static Vec2<T> normalize(const Vec2<T>& v)
        {
            const T norma = v->norma();
            return Vec2<T>(v / norma);
        }
        // multiplicacao por matriz 3x3 a esquerda
};

template <typename T> class Vec3
{
    public:
        T x, y, z;
        Vec3(): x(T(0)), y(T(0)), z(T(0)) {} //build all the vec3 operation inside this public section
        Vec3(T _x): x(T(_x)), y(T(_x)), z(T(_x)) {}
        Vec3(T _x, T _y, T _z): x(T(_x)), y(T(_y)), z(T(_z)) {}
        Vec3<T> operator + (const Vec3<T>& v) const { return Vec3<T>(x + v.x, y + v.y, z + v.z); }
        Vec3<T> operator - (const Vec3<T>& v) const { return Vec3<T>(x - v.x, y - v.y, z - v.z); }
        Vec3<T> operator * (const T& t) const { return Vec3<T>(x*t, y*t, z*t); }
        Vec3<T> operator / (const T& t) const { return Vec3<T>(x/t, y/t, z/t); }
        T operator * (const Vec3<T>& v) const { return (x*v.x) + (y*v.y) + (z*v.z); }
        Vec3<T> operator ^ (const Vec3<T>& v) const { return Vec3<T>(y*v.z - (z*v.y), z*v.x - (x*v.z), x*v.y - (y*v.x)); }
        static T norma(const Vec3<T>& v)
        {
            return pow((v.x*v.x) + (v.y*v.y) + (v.z*v.z), 0.5);
        }
        static Vec3<T> normalize(const Vec3<T>& v)
        {
            const T norma = v.norma(v);
            return Vec3<T>(v / norma);
        }
        Vec3<T> operator * (const Matrix4<T>& m) const 
        { 
            T _x, _y, _z;
            for(size_t column = 0; column < m->matrix.size(); column++)
            {
                _x += x*m->matrix[0][column]; 
                _y += y*m->matrix[1][column];
                _z += z*m->matrix[2][column];
            }
            return Vec3<T>(_x, _y, _z); 
        }
};

typedef Vec3<double> Vec3D;
typedef Vec2<double> Vec2D;


#endif