#include <iostream>
#include <vector>
#include <math.h>

struct PointD
{
  PointD(float _x, float _y, float _z, float _d) : x(_x), y(_y), z(_z), d(_d) {}
  float x, y, z, d;
};

struct Point
{
  Point(float *p) { xyz[0] = p[0]; xyz[1] = p[1]; xyz[2] = p[2]; }
  Point(float x, float y, float z) { xyz[0] = x; xyz[1] = y; xyz[2] = z; }
  Point(const Point& p) { for (int i = 0; i < 3; i++) xyz[i] = p.xyz[i]; }
  Point() { xyz[0] = xyz[1] = xyz[2] = 0.0; }
  Point(PointD& p) { xyz[0] = p.x; xyz[1] = p.y; xyz[2] = p.z; }
  float& operator[](int i) { return xyz[i]; }
  float x() { return xyz[0]; }
  float y() { return xyz[1]; }
  float z() { return xyz[2]; }
  float xyz[3];
  void print(std::ostream& o) { o << x() << " " << y() << " " << z(); }
  static float dist(Point a, Point b);
  static Point cross(struct Point a, struct Point b)
  {
    return Point(a.y()*b.z() - a.z()*b.y(), a.z()*b.x() - a.x()*b.z(), a.x()*b.y() - a.y()*b.x());
  }
  static float dot(struct Point a, struct Point b)
  {
    return a.x()*b.x() + a.y()*b.y() + a.z()*b.z();
  }
};

typedef std::vector<Point> PointList;

Point operator-(Point a, Point b);
Point operator+(Point a, Point b);
Point operator*(Point a, float r);
Point operator/(Point a, float r);


