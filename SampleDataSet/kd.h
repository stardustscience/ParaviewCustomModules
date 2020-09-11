#include <iostream>
#include <vector>
#include <math.h>

#include "point.h"

struct KDNode
{
  KDNode()
  {
    axis = -2;
    left = NULL; 
    right = NULL;
  }

  KDNode(Point p)
  {
    point = p;
    axis = -1;
  }
    
  KDNode(PointList pl)
  { 
    if (pl.size() == 1) {point = pl[0]; left = NULL; right = NULL; axis = -1; }
    else split(pl);
  }

  ~KDNode()
  {
    if (left) delete left;
    if (right) delete right;
  }

  void split(PointList pl)
  {
    Point pmin = pl[0];
    Point pmax = pl[0];

    for (auto pp : pl)
      for (int j = 0; j < 3; j++)
      {
        if (pmin[j] > pp[j]) pmin[j] = pp[j];
        if (pmax[j] < pp[j]) pmax[j] = pp[j];
      }

    Point diff = pmax - pmin;

    axis = ((diff.x() > diff.y()) && (diff.x() > diff.z())) ? 0 : (diff.y() > diff.z()) ? 1 : 2;
    min = pmin[axis];
    max = pmax[axis];
      
    sep = (max + min) / 2.0;

    int l = 0, r = 0;
    for (auto pp : pl)
      if (pp[axis] < sep) l++;
      else r++;

    PointList left_points, right_points;
    left_points.resize(l); right_points.resize(r);
    int ll = 0, rr = 0;
    for (auto pp : pl)
      if (pp[axis] < sep) left_points[ll++] = pp;
      else right_points[rr++] = pp;

    left = new KDNode(left_points);
    right = new KDNode(right_points);
  }

  void insert(Point p)
  {
    if (axis == -2)
    {
      axis = -1;
      point = p;
    }
    else if (axis == -1)
    {
      PointList pl = {point, p};
      split(pl);
    }
    else
    {
      if (p[axis] < sep)
        left->insert(p);
      else
        right->insert(p);
    }
  }

  int count() 
  {
    if (axis == -1) return 1;
    else return left->count() + right->count();
  }

  bool search(Point p, float radius)
  {
    if (axis == -2)
      return false;
    if (axis == -1)
      return Point::dist(p, point) < radius;
    else
    {
      bool l =  p[axis] + radius > sep ? right->search(p, radius) : false;
      bool r =  p[axis] - radius < sep ? left->search(p, radius) : false;
      return l || r;
    }
  }

  void print(std::ostream& o, int indent)
  {
    for (int i = 0; i < indent; i++) o << " ";
    if (axis == -1) { o << "leaf: "; point.print(o); o << "\n";}
    else
    {
      o << "axis " << axis << " min " << min << " max " << max << " sep " << sep << "\n";
      for (int i = 0; i < indent; i++) o << " ";
      o << "left:\n";
      left->print(o, indent + 1);
      for (int i = 0; i < indent; i++) o << " ";
      o << "right:\n";
      right->print(o, indent + 1);
    }
  }

  Point  point;
  float  min;
  float  max;
  int    axis;
  float  sep;
  KDNode *left, *right;
};
