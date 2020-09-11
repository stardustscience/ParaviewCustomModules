#include "point.h"

Point operator-(Point a, Point b)
{
  return Point(a.x() - b.x(), a.y() - b.y(), a.z() - b.z());
}

Point operator+(Point a, Point b)
{
  return Point(a.x() + b.x(), a.y() + b.y(), a.z() + b.z());
}

Point operator*(Point a, float r)
{
  return Point(a.x()*r, a.y()*r, a.z()*r);
}

Point operator/(Point a, float r)
{
  return Point(a.x()/r, a.y()/r, a.z()/r);
}


float 
Point::dist(Point a, Point b)
{
  Point d = a - b;
  return sqrt(d.x()*d.x() + d.y()*d.y() + d.z()*d.z());
}


