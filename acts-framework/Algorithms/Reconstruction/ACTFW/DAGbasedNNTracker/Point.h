#ifndef _POINT_H_
#define _POINT_H_
// Neural Network based tracker
// M.Kunze, Heidelberg University, 2018

#include <cmath>

//Basics for 3d coordinate representation

struct point {
    double x, y, z;
    point() {}
    point(double x, double y, double z) : x(x),y(y),z(z) {}
    inline point operator-(const point&p) {
        return point(x-p.x, y-p.y, z-p.z);
    }
    inline point operator+(const point&p) {
        return point(x+p.x, y+p.y, z+p.z);
    }
    inline double operator*(const point&p) {
        return x*p.x+y*p.y+z*p.z;
    }
    inline point operator*(double f) {
        return point(x*f, y*f, z*f);
    }
    static double norm(const point &a);
    static double norm2(const point &a);
    static double dist(const point&p);
    static double dist3(point &a,point &b,point &c);
    static double dot(const point &a,const point &b);
    static point cross(const point &a,const point &b);
    static bool sortByRadius(const point &a,const point &b);
    static bool sortByRz(const point &a,const point &b);
};

// Calculate the scalar product between two vectors a and b
inline
double point::dot(const point &a,const point &b)
{
    double r1 = sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
    double r2 = sqrt(b.x*b.x + b.y*b.y + b.z*b.z);
    double product = a.x*b.x + a.y*b.y + a.z*b.z;
    double result = product / (r1 * r2);
    if (result>1.0) result = 1.0;
    if (result<-1.0) result = -1.0;
    return result;
}

// Calculate the vector product between two vectors a and b
inline point point::cross(const point &a,const point &b) {
    return point(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}

inline
double point::norm(const point &a)
{
    return sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
}

inline
double point::norm2(const point &a)
{
    return a.x*a.x+a.y*a.y+a.z*a.z;
}

inline
bool point::sortByRadius(const point &a,const point &b)
{
    return (a.x*a.x+a.y*a.y+a.z*a.z < b.x*b.x+b.y*b.y+b.z*b.z);
}

inline
bool point::sortByRz(const point &a,const point &b)
{
    return (a.x*a.x+a.y*a.y < b.x*b.x+b.y*b.y);
}

inline double point::dist(const point &p) { return sqrt(p.x*p.x+p.y*p.y+p.z*p.z); }

inline double point::dist3(point &a,point &b,point &c) {
    const point x = a-b;
    const point y = a-c;
    const point z = c-b;
    double d = dist( cross(x,y)) / dist(z);
    return d;
    
}
#endif
