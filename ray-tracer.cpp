#include <iostream>
#include <math.h>
#include <vector>

#define infinity 1e8
#define kEpsilon 10e-6

// lets use doubles for object-ray intersection and floats for shading calculations

template <typename T> class Matrix4
{  
    public:
        std::vector<std::vector<T>> matrix;
        Matrix4(): matrix(std::vector<std::vector<T>>(4, std::vector<T>(4))) {}
};
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
        T norma(const Vec3<T>& v) const
        {
            return pow((v.x*v.x) + (v.y*v.y) + (v.z*v.z), 0.5);
        }
        Vec3<T> normalize(const Vec3<T>& v) const
        {
            const T norma = v->norma();
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
        Point3<T> operator + (const Vec3<T>& v) const { return Point3<T>(x + v.x, y + v.y, z + v.z); }
        Point3<T> operator - (const Vec3<T>& v) const { return Point3<T>(x - v.x, y - v.y, z - z.v); }
        Point3<T> operator * (const T& t) const { return Point3<T>(x*t, y*t, z*t); }
        T operator ^ (const Point3<T>& p) const { return (x*p.x) + (y*p.y) + (z*p.z); }
        Point3<T> operator * (const Matrix4<T>& m) const 
        {
            T _x, _y, _z;
            for(size_t column = 0; column < m->matrix.size(); column++)
            {
                _x += x*m->matrix[0][column]; 
                _y += y*m->matrix[1][column];
                _z += z*m->matrix[2][column];
            }
            return Point3<T>(_x, _y, _z); 
        }
        Vec3<T> operator - (const Point3<T>& p) const { return Vec3<T>(x - p.x, y - p.y, z - p.z); }
        T distance(const Point3<T>& p) const
        {
            Vec3<T> v = Vec3<T>(x - p.x, y - p.y, z - p.z);
            return v.norma(v);
        }
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
        virtual bool rayObjectIntersect(const Ray &ray, double& tmin) const = 0;
};

class Sphere: public Object 
{
    public:
        Point3D center;
        double radius;
        Sphere(Point3D c, double r): center(c), radius(r) {}
        ~Sphere() {}
        bool rayObjectIntersect(const Ray &ray, double& tmin) const
        {
            double a = pow(ray.direction.norma(ray.direction), 2);
            double b = ((ray.origin - this->center) * ray.direction) * 2;
            double c = ((this->center ^ this->center) + (ray.origin ^ ray.origin)) + (-2)*(ray.origin ^ this->center) - (this->radius*this->radius);
            double delta = (b*b) - 4*a*c; 
            if (delta == 0)
            {
                double t = (b*b)/2*a;
                if (t > kEpsilon)
                {
                    tmin = t;
                    return (true);
                } else {
                    return (false);
                }
            }
            else if (delta > 0)
            {
                double sqrtDelta = pow(delta, 0.5);
                double t1 = ((b*b) + sqrtDelta)/2*a;
                double t2 = ((b*b) - sqrtDelta)/2*a;
                if (t1 > t2 && t1 > kEpsilon)
                {
                    tmin = t1;
                    return (true);
                }
                else if (t2 > kEpsilon)
                {
                    tmin = t2;
                    return (true);
                }
                return false;
            }
            return (false);
        }
};

class Plane: public Object 
{
    public:
        Vec3D normal;
        Point3D pp;
        Plane(Vec3D n, Point3D p): normal(n), pp(p) {}
        ~Plane() {}
        bool rayObjectIntersect(const Ray &ray, double& tmin) const 
        {
            double t = ((pp - ray.origin) * this->normal) / (ray.direction * this->normal);
            if (t > kEpsilon)
            {
                tmin = tmin;
                return (true);

            } else {
                return (false);
            }
        }
};

class Triangle: public Object 
{
    public:
        Point3D A;
        Point3D B;
        Point3D C;
        Triangle(Point3D a, Point3D b, Point3D c): A(a), B(b), C(c) {}
        ~Triangle() {}
        bool rayObjectIntersect(const Ray& ray, double& tmin) const
        {
            Vec3D tPlaneNormal = (this->A - this->B) ^ (this->A - this->C);
            if (((ray.direction * (-1)) * tPlaneNormal) <= 0)
            {
                tPlaneNormal = (this->A - this->C) ^ (this->A - this->B);
            }
            Plane *tPlane = new Plane(tPlaneNormal, this->A);
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
                Sphere *e = new Sphere(Point3D(), 10.0);
                objetos.push_back(e);
                break;    
            }
            case 'p':
            {
                Plane *p = new Plane(Vec3D(), Point3D());
                objetos.push_back(p);
                break;
            }
            case 't':
            {
                Triangle *t = new Triangle(Point3D(), Point3D(), Point3D());
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
            case 'r':
            {
                Line *line = new Line();
                objetos.push_back(line);
                break;
            }  
        }
    }
    return 0;
}