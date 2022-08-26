#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <cstdio>
#include <math.h>

// lets use doubles for object-ray intersection and floats for shading calculations

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
        T norma(const Vec2<T>& v) const
        {
            return pow((v.x*v.x) + (v.y*v.y), 0.5);
        }
        Vec2<T> normalize(const Vec2<T>& v) const
        {
            return Vec2<T>(v.x/v->norma(), v.y/v->norma());
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
        // normalização
        // norma
        // multiplicacao por matriz 4x4 a esquerda
        // ponto - ponto
};

template <typename T> class Point2
{
    
};

template <typename T> class Point3
{
    public:
        T x, y, z;
        Point3(): x(T(0)), y(T(0)), z(T(0)) {}
        Point3(T _x): x(T(_x)), y(T(_x)), z(T(_x)) {}
        Point3(T _x, T _y, T _z): x(T(_x)), y(T(_y)), z(T(_z)) {}
        // ponto + vetor, ponto - vetor
        // ponto multiplicado por escalar pela esquerda ou direita
        // soma varicentrica
        // interpolação
        // ponto multiplicado por matriz 4x4 pela esquerda
        // distancia entre dois pontos
        // ponto - ponto
};

template <typename T> class Matrix4
{  
    public:
        std::vector<std::vector<T>> matrix;
        Matrix4(): matrix(std::vector<std::vector<T>>(4, std::vector<T>(4))) {}
};

typedef Vec2<double> Vec2D;
typedef Vec3<double> Vec3D;
typedef Point3<double> Point3D;
typedef Matrix4<double> Matrix4D;
typedef Vec3<double> RGBColor;

class Ray {
    public:
        Point3D origin;
        Vec3D direction;
        Ray(void);                              // default constructor
        Ray(const Point3D& o, const Vec3D& dir);// constructor
        Ray(const Ray& ray);                    // copy constructor
        Ray& operator= (const Ray& rhs);        // assignment operator
        ~Ray() {}                               // destructor
};

class HitInfo {
    public:
        bool hit_object;
        Point3D hit_location;
        Vec3D normal;
        RGBColor surface_color;
        HitInfo(void);
        HitInfo(const HitInfo& hitInfo);
        ~HitInfo(void);
        HitInfo& operator= (const HitInfo& rhs);
};

HitInfo::HitInfo():hit_object(false),hit_location(),normal(),surface_color() {} // constructor

class Object 
{
    public:
        Object() {}
        virtual ~Object() {}
        virtual bool rayObjectIntersect(const Ray &ray, double& tmin) const = 0; //consider t from e = [10ˆ-6, +infiniyt)
};

class Sphere: public Object 
{
    public:
        Sphere() {}
        ~Sphere() {}
        bool rayObjectIntersect(const Ray &ray, double& tmin) const
        {
            return 1;
        }
};

class Triangle: public Object 
{
    public:
        Triangle() {}
        ~Triangle() {}
        bool rayObjectIntersect(const Ray &ray, double& tmin) const
        {
            return 1;
        }
};

class Plane: public Object 
{
    public:
        Plane() {}
        ~Plane() {}
        bool rayObjectIntersect(const Ray &ray, double& tmin) const 
        {
            return 1;
        }
};

class Line: public Object 
{
    public:
        Line() {}
        ~Line() {}
        bool rayObjectIntersect(const Ray &ray, double& tmin) const
        {
            return 1;
        }
};

class Light
{
    public:
        Point3D lightPos;
        Vec3D lightColor;
        Light(Point3D pos, Vec3D color): lightPos(pos), lightColor(color) {}
        ~Light() {}
};

class Camera 
{
    public:
        int hr, vr;
        float distance;
        Vec3D up, toScreen;
        Point3D cameraPos;
        Camera(int _hr, int _vr, float d, Vec3D _up, Vec3D _toScreen, Point3D pos): hr(_hr), vr(_vr), distance(d), up(_up), toScreen(_toScreen), cameraPos(pos) {} 
        ~Camera() {}
};

int main()
{
    Matrix4D matriz = Matrix4D();
    std::vector<Object*> objetos;
    std::vector<Light*> lights;
    Camera *camera;
    char objectType;
    int _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12;
    while (scanf("%c %i %i %i %i %i %i %i %i %i %i %i %i\n", &objectType, &_1, &_2, &_3, &_4, &_5, &_6, &_7, &_8, &_9, &_10, &_11, &_12) != EOF) 
    {
        switch(objectType)
        {
            case 's': 
            {
                Sphere *e = new Sphere();
                objetos.push_back(e);
                break;    
            }
            case 'p':
            {
                Plane *p = new Plane();
                objetos.push_back(p);
                break;
            }
            case 't':
            {
                Triangle *t = new Triangle();
                objetos.push_back(t);
                break;
            }
            case 'l':
            {
                Light *l = new Light(Point3D(), Vec3D());
                lights.push_back(l);
                break;
            }
            case 'c':
            {
                camera = new Camera(1920, 1080, 10.0, Vec3D(), Vec3D(), Point3D());
                break;
            }
            case 'x' || 'y' || 'z':
            {
                Line *line = new Line();
                objetos.push_back(line);
                break;
            }   
        }
    }
    return 0;
}