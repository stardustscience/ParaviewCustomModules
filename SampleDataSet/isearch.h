#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>

namespace isearch
{
  struct interval
  {
    interval() : value(0), indx(0) {}
    interval(float v, int i) : value(v), indx(i) {}

    float value;
    int   indx;
  };

  struct IntervalNode
  {
    IntervalNode *left_subtree, *right_subtree;
    float left_value, value;
    int   indx;
    
    IntervalNode(float v, int i) :
        value(v),
        indx(i), 
        left_value(0),
        left_subtree(NULL),
        right_subtree(NULL)
    {
    }

    ~IntervalNode()
    {
      if (left_subtree) delete left_subtree;
      if (right_subtree) delete right_subtree;
    }

    IntervalNode(interval *ilist, int s, int e)
    {
      int n = (e - s) + 1;
      int m = s + (n / 2);
      value = ilist[m].value;
      indx = ilist[m].indx;
      if (m > s)
      {
        left_subtree = new IntervalNode(ilist, s, m-1);
        left_value = ilist[m-1].value;
      }
      else
        left_subtree = NULL;
      if (m < e)
        right_subtree = new IntervalNode(ilist, m+1, e);
      else
        right_subtree = NULL;
    }

    int
    Search(float v)
    {
      if (value > v)
      {
        if (left_subtree && left_value > v)
          return left_subtree->Search(v);
        else
          return indx;
      }
      else
      {
        if (! right_subtree)
        {
          std::cerr << "error!\n";
          std::cerr << "search value: " << v << "\n";
          std::cerr << "node value: " << value << "\n";
          exit(1);
        }
        else
          return right_subtree->Search(v);
      }
    }
  };

  struct IntervalSearchTree
  {
    IntervalSearchTree(float *values, int *indices, int nValues)
    {
      interval *intervals = new interval[nValues];
      intervals[0].value = values[0];
      intervals[0].indx = indices[0];
      for (int i = 1; i < nValues; i++)
      {
        intervals[i].value = intervals[i-1].value + values[i];
        intervals[i].indx = indices[i];
      }
      root = new IntervalNode(intervals, 0, nValues-1);
      max = intervals[nValues-1].value;
      delete[] intervals;
    }

    int Search(float v) { return root->Search(v); }

    IntervalNode *root;
    float max;
  };
}

