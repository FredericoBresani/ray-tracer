#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <algorithm>

#define infinity 1e8
#define kEpsilon 10e-6

// lets use doubles for object-ray intersection and floats for shading calculations

template <typename T> class Matrix4
{  
    public:
        std::vector<std::vector<T>> matrix;
        Matrix4(): matrix(std::vector<std::vector<T>>(4, std::vector<T>(4))) {}
};

template <typename T> class Matrix3
{
    public:
        std::vector<std::vector<T>> matrix;
        Matrix3(): matrix(std::vector<std::vector<T>>(3, std::vector<T>(3))) {}
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

template <typename T> class Point2
{
    public:
        T x, y;
        Point2(): x(T(0)), y(T(0)) {}
        Point2(T _x): x(_x), y(_x) {}
        Point2(T _x, T _y): x(_x), y(_y) {} 
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
typedef Point2<double> Point2D;
typedef Matrix4<double> Matrix4D;
typedef Matrix3<double> Matrix3D;
typedef Vec3<double> RGBColor;


class Ray {
    public:
        Point3D origin;
        Vec3D direction;
        Ray(void);                              // default constructor
        Ray(const Point3D& o, const Vec3D& dir): origin(o), direction(dir) {}// constructor
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
        Vec3D color;
        Object() {}
        virtual ~Object() {}
        virtual bool rayObjectIntersect(const Ray &ray, double& tmin, HitInfo& info) const = 0;
        virtual Vec3D getColor() const = 0;
};

class Sphere: public Object 
{
    public:
        Point3D center;
        Vec3D color;
        double radius;
        Sphere(Point3D c, Vec3D RGB, double r): center(c), color(RGB), radius(r) {}
        ~Sphere() {}
        bool rayObjectIntersect(const Ray &ray, double& tmin, HitInfo& info) const
        {
            double a = pow(ray.direction.norma(ray.direction), 2.0);
            double b = ((ray.origin - this->center) * ray.direction) * 2.0;
            double c = ((this->center ^ this->center) + (ray.origin ^ ray.origin)) + (-2.0)*(ray.origin ^ this->center) - (this->radius*this->radius);
            double delta = (b*b) - 4.0*a*c; 
            if (delta == 0.0)
            {
                double t = (b*b)/(2.0*a);
                if (t > kEpsilon)
                {
                    tmin = t;
                    info.hit_location = ray.origin + ray.direction*tmin;
                    info.normal = info.hit_location - center;
                    info.normal = info.normal.normalize(info.normal);
                    return (true);
                } else {
                    return (false);
                }
            }
            else if (delta > 0.0)
            {
                double sqrtDelta = pow(delta, 0.5);
                double t1 = (((-1)*b) + sqrtDelta)/(2.0*a);
                double t2 = (((-1)*b) - sqrtDelta)/(2.0*a);
                if (t1 < t2 && t1 > kEpsilon)
                {
                    tmin = t1;
                    info.hit_location = ray.origin + ray.direction*tmin;
                    info.normal = info.hit_location - center;
                    info.normal = info.normal.normalize(info.normal);
                    return (true);
                }
                else if (t2 > kEpsilon)
                {
                    tmin = t2;
                    info.hit_location = ray.origin + ray.direction*tmin;
                    info.normal = info.hit_location - center;
                    info.normal = info.normal.normalize(info.normal);
                    return (true);
                }
                return false;
            }
            return (false);
        }
        Vec3D getColor() const
        {
            return this->color;
        }
};

class Plane: public Object 
{
    public:
        Vec3D normal;
        Point3D pp;
        Vec3D color;
        Plane(Vec3D n, Point3D p, Vec3D RGB): normal(n), pp(p), color(RGB) {}
        ~Plane() {}
        bool rayObjectIntersect(const Ray &ray, double& tmin, HitInfo& info) const 
        {
            double t = ((pp - ray.origin) * this->normal) / (ray.direction * this->normal);
            Point3D location = ray.origin + ray.direction*t;
            if (location.x > 5 || location.x < -5) return (false);
            if (location.y > 5 || location.y < -5) return (false);
            if (location.z > 5 || location.z < -5) return (false);
            if (t > kEpsilon)
            {
                tmin = t;
                info.hit_location = ray.origin + ray.direction*tmin;
                info.normal = normal.normalize(normal);
                return (true);

            } else {
                return (false);
            }
        }
        Vec3D getColor() const
        {
            return this->color;
        }
};

class Triangle: public Object 
{
    public:
        Point3D A;
        Point3D B;
        Point3D C;
        Vec3D color;
        Triangle(Point3D a, Point3D b, Point3D c, Vec3D RGB): A(a), B(b), C(c), color(RGB) {}
        ~Triangle() {}
        bool rayObjectIntersect(const Ray& ray, double& tmin, HitInfo& info) const
        {
            Vec3D tPlaneNormal = (this->A - this->B) ^ (this->A - this->C);
            if (((ray.direction * (-1)) * tPlaneNormal) <= 0)
            {
                tPlaneNormal = (this->A - this->C) ^ (this->A - this->B);
            }
            Plane *tPlane = new Plane(tPlaneNormal, this->A, color);
            Point3D pHit;
            if (tPlane->rayObjectIntersect(ray, tmin, info)) 
            {   
                pHit = info.hit_location;       
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
                if (v2.z == 0.0 && (point[2] < -kEpsilon || point[2] > kEpsilon)) return (false);
                if (v2.z != 0.0)
                {
                    gama = point[2]/v2.z; 
                }
                if (v1.y == 0.0 && (point[1] < -kEpsilon || point[1] > kEpsilon)) return (false);
                if (v1.y != 0.0)
                {
                    beta = point[1]/v1.y;
                }
                if (v0.x == 0.0 && (point[0] < -kEpsilon || point[0] > kEpsilon)) return (false);
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
                    return (true);
                } else {
                    return (false);
                }
            } else {
                return 0;
            }
        }
        Vec3D getColor() const
        {
            return this->color;
        }
};

class Line: public Object 
{
    public:
        Point3D origin;
        Vec3D direction, normal;
        Vec3D color;
        Line(Point3D o, Vec3D d, Vec3D RGB): origin(o), direction(d), color(RGB) {}
        ~Line() {}
        bool rayObjectIntersect(const Ray &ray, double& tmin, HitInfo& hit) const
        { 
            double t = ((ray.origin - origin)*normal)/((normal*ray.direction)*(-1.0));
            Point3D location = ray.origin + (ray.direction*t);
            if (location.x > 7 || location.x < -7) return (false);
            if (location.y > 7 || location.y < -7) return (false);
            if (location.z > 7 || location.z < -7) return (false);
            Vec3D toLoc = location - origin;
            toLoc = toLoc.normalize(toLoc);
            double cos = toLoc*(direction.normalize(direction));
            if (t > kEpsilon && (cos >= -1.0 - 0.000002 && cos <= -1.0 + 0.000002 || cos >= 1.0 - 0.000002 && cos <= 1.0 + 0.000002))
            {
                Vec3D color = this->getColor();
                tmin = t;
                hit.hit_location = Point3D(color.x, color.y, color.z);
                hit.normal = direction.normalize(direction);
                return (true);
            }
            return (false);
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
        double distance, pixelsize, pixelQtnH, pixelQtnV;
        Vec3D up, u, v, w, right, iup;
        Point3D cameraPos, lookAt;
        Camera(int _hr, int _vr, double d, Vec3D _up, Point3D pos, Point3D _lookAt): hr(_hr), vr(_vr), distance(d), up(_up), cameraPos(pos), lookAt(_lookAt) {} 
        void makeCamera(double pixel)
        {
            // Tha camera base should follow this order {w = z, v = y, u = x}
            pixelsize = pixel;
            pixelQtnH = (double)hr/pixelsize;
            pixelQtnV = (double)vr/pixelsize;
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
            
        }
        ~Camera() {}

        Point3D worldToCameraCoordinates(Point3D point)
        {
            double x = (w.x*(point.x)) + (w.y*(point.y)) + (w.z*(point.z));
            double y = (v.x*(point.x)) + (v.y*(point.y)) + (v.z*(point.z));
            double z = (u.x*(point.x)) + (u.y*(point.y)) + (u.z*(point.z));
            return Point3D(z, y, x); // Remenber that the base order is like that
        }

        Point2D worldToScreenCoordinates(Point3D point) //this is a point on the screen already, althoug in world coordinates
        {                                               //that is why no projection is needed
            Point3D cameraCoordinates = Camera::worldToCameraCoordinates(point);
            return Point2D(cameraCoordinates.x, cameraCoordinates.y);
        }
};

Vec3D setPixelColorNormal(Vec3D &normal)
{
    return Vec3D(std::min(double(1), normal.x), std::min(double(1), normal.y), std::min(double(1), -normal.z))*255.0;
}

Vec3D setPixelColorCoordinates(Point3D &location)
{   
    //x = red
    //y = green
    //z = blue
    double red = 0.0, green = 0.0, blue = 0.0;
    Vec3D aux = Vec3D(location.x, location.y, location.z)/255.0;
    red += aux.x < 0.0 ? 0.0 : aux.x;
    green += aux.y < 0.0 ? 0.0 : aux.y;
    blue += aux.z < 0.0 ? 0.0 : aux.z;
    aux = Vec3D(std::min(double(1), red*40.0), std::min(double(1), green*40.0), std::min(double(1), blue*40.0))*255.0;
    return aux;
}

Vec3D trace(const Point3D& origin, const Point3D& pixel, std::vector<Object*>& objetos, HitInfo& hit, Camera *camera)
{
    double t = infinity;
    double tmin = infinity;
    Ray *ray = new Ray(origin, pixel - origin);
    Vec3D color;
    for (int i = 0; i < objetos.size(); i++)
    {
        if (objetos[i]->rayObjectIntersect(*ray, t, hit))
        {
            if (t < tmin)
            {
                tmin = t;
                // color = objetos[i]->getColor();
                color = setPixelColorNormal(hit.normal);
                // color = setPixelColorCoordinates(hit.hit_location);
            }
        }
    }
    if (tmin == infinity)
    {
        Point2D screenCoordinates = camera->worldToScreenCoordinates(pixel);
        double maxScreen = 0;
        if (camera->hr >= camera->vr) {
            maxScreen = double(camera->hr)/double(camera->vr);
        } else {
            maxScreen = double(camera->vr)/double(camera->hr);
        }
        double yNormalization = std::abs(screenCoordinates.y)/maxScreen;
        double xNormalization = std::abs(screenCoordinates.x)/maxScreen;
        double atenuation = 2.0*(1.0 - yNormalization);
        return Vec3D(121.0, 100.0, 138.0)/(atenuation >= 1.0 ? atenuation : 1.0);
    }
    return color;
}

void render(std::vector<Object*>& objetos, std::vector<Light*>& lights, Camera& camera)
{
    Vec3D toPixel = camera.w*camera.distance + camera.right*(-camera.pixelQtnH/2.0) + camera.iup*(camera.pixelQtnV/2.0) - (camera.iup/2.0) + (camera.right/2.0);
    Point3D screenP = camera.cameraPos + toPixel;
    Vec3D down;
    HitInfo *hInfo = new HitInfo();
    std::vector<Vec3D> pixels;
    for (int i = 0; i < camera.pixelQtnH*camera.pixelQtnV; i++)
    {
        if (i == 400*195 - 180)
        {
            int fr = 10;
        }
        if ((i) % (int)camera.pixelQtnH == 0)
        {
            down = down - camera.iup;
            screenP = camera.cameraPos + toPixel;
            screenP = screenP + down;
        } else {
            screenP = screenP + camera.right;
        }
        pixels.push_back(trace(camera.cameraPos, screenP, objetos, *hInfo, &camera));
    }
    std::ofstream pixelOutput("./image.ppm", std::ios::out | std::ios::binary);
    pixelOutput << "P6\n" << camera.pixelQtnH << " " << camera.pixelQtnV << "\n255\n";
    for (int i = 0; i < camera.pixelQtnH*camera.pixelQtnV; i++)
    {
        pixelOutput <<(unsigned char)(std::max(double(1), pixels[i].x)) <<
            (unsigned char)(std::max(double(1), pixels[i].y)) <<
            (unsigned char)(std::max(double(1), pixels[i].z));
    }
    pixelOutput.close();   
}

int main()
{
    // Matrix4D matriz = Matrix4D();
    std::vector<Object*> objetos;
    std::vector<Light*> lights;
    Camera *camera;
    char objectType;
    bool read = true;
    float _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12;
    while (scanf("%c %f %f %f %f %f %f %f %f %f %f %f %f\n", &objectType, &_1, &_2, &_3, &_4, &_5, &_6, &_7, &_8, &_9, &_10, &_11, &_12) != EOF && read) 
    {
        switch(objectType)
        {
            case 's': 
            {
                Sphere *e = new Sphere(Point3D(_1, _2, _3), Vec3D(_5, _6, _7), _4);
                objetos.push_back(e);
                break;    
            }
            case 'p':
            {
                Plane *p = new Plane(Vec3D(_4, _5, _6), Point3D(_1, _2, _3), Vec3D(_7, _8, _9));
                objetos.push_back(p);
                break;
            }
            case 't':
            {
                Triangle *t = new Triangle(Point3D(_1, _2, _3), Point3D(_4, _5, _6), Point3D(_7, _8, _9), Vec3D(_10, _11, _12));
                objetos.push_back(t);
                break;
            }
            case 'l':
            {
                Light *l = new Light(Point3D(_1, _2, _3), Vec3D(_4, _5, _6));
                lights.push_back(l);
                break;
            }
            case 'c':
            {
                camera = new Camera(_1, _2, _3, Vec3D(_4, _5, _6), Point3D(_7, _8, _9), Point3D(_10, _11, _12));
                camera->makeCamera(1.0);
                break;
            }
            case 'r':
            {
                Line *line = new Line(Point3D(_1, _2, _3), Vec3D(_4, _5, _6), Vec3D(_7, _8, _9));
                line->setNormal(camera->cameraPos);
                objetos.push_back(line);
                break;
            }  
            default:
            {
                read = false;
                break;
            }
        }
    }
    render(objetos, lights, *camera);
    return 0;
}