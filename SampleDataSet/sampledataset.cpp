#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>

#include "kd.h"
#include "isearch.h"

#define RAND (((float)rand())/RAND_MAX)

PointList samples;

extern "C" void
SampleCartesian(int ns, float *origin, float *spacing, int *counts, float R0, float R1, int pdep, float *data, int nRetained, float *retained_points, float *retained_data)
{
  float *point_centered_data = NULL;
  float *cell_centered_data = NULL;

  std::cerr << "Sample! ns = " << ns << "\n";
  srand(12345);

  // Vertex counts

  int nx = counts[0];
  int ny = counts[1];
  int nz = counts[2];

  int np = nx*ny*nz;
  int nc = (nx-1)*(ny-1)*(nz-1);

  int xp = nz;
  int yxp = ny*nz;

  int xc  = (nz-1);
  int yxc = (ny-1)*(nz-1);

#define PINDEX(i,j,k) ((k)*yxp + (j)*xp + (i))

  if (pdep) 
  {
    point_centered_data = data;

    cell_centered_data = new float[(nx-1)*(ny-1)*(nz-1)];

    // For each cell, sum the 8 vertices and divide

    float *ccd = cell_centered_data;
    for (int x = 0; x < (nx-1); x++)
      for (int y = 0; y < (ny-1); y++)
        for (int z = 0; z < (nz-1); z++)
          *ccd++ = (point_centered_data[PINDEX(x+0,y+0,z+0)] +
                    point_centered_data[PINDEX(x+0,y+0,z+1)] +
                    point_centered_data[PINDEX(x+0,y+1,z+0)] +
                    point_centered_data[PINDEX(x+0,y+1,z+1)] +
                    point_centered_data[PINDEX(x+1,y+0,z+0)] +
                    point_centered_data[PINDEX(x+1,y+0,z+1)] +
                    point_centered_data[PINDEX(x+1,y+1,z+0)] +
                    point_centered_data[PINDEX(x+1,y+1,z+1)]) / 8.0;
          
  }
  else
  {
    cell_centered_data = data;

    // For each cell, add to the adjacent vertices and bump counts
    // then divide through

    point_centered_data = new float[nx*ny*nz];
    int* counts = new int[nx*ny*nz];
    for (int i = 0; i < nx*ny*nz; i++)
    {
      point_centered_data[i] = 0;
      counts[i] = 0;
    }

    float* ccd = cell_centered_data;
    for (int x = 0; x < (nx-1); x++)
      for (int y = 0; y < (ny-1); y++)
        for (int z = 0; z < (nz-1); z++)
        {
          int i;
          i = PINDEX(x+0,y+0,z+0); counts[i]++; point_centered_data[i] += *ccd;
          i = PINDEX(x+0,y+0,z+1); counts[i]++; point_centered_data[i] += *ccd;
          i = PINDEX(x+0,y+1,z+0); counts[i]++; point_centered_data[i] += *ccd;
          i = PINDEX(x+0,y+1,z+1); counts[i]++; point_centered_data[i] += *ccd;
          i = PINDEX(x+1,y+0,z+0); counts[i]++; point_centered_data[i] += *ccd;
          i = PINDEX(x+1,y+0,z+1); counts[i]++; point_centered_data[i] += *ccd;
          i = PINDEX(x+1,y+1,z+0); counts[i]++; point_centered_data[i] += *ccd;
          i = PINDEX(x+1,y+1,z+1); counts[i]++; point_centered_data[i] += *ccd;
          ccd++;
        }

    for (int i = 0; i < nx*ny*nz; i++)
      point_centered_data[i] /= counts[i];

    delete[] counts;
  }

  KDNode kdtree;

  int kept = 0;
  for (int i = 0; i < nRetained; i++)
  {
    float d = retained_data[i];
    if (d > 0)
    {
      Point p(retained_points+3*i);
      float r = R0 + d*(R1 - R0);
      if (! kdtree.search(p, r))
      {
        samples.push_back(p);
        kdtree.insert(p);
        kept ++;
      }
    }
  }

  int nNonZeroCells = 0;
  for (int i = 0; i < nc; i++)
  {
    float d = cell_centered_data[i];
    if (d > 0)
      nNonZeroCells ++;
  }

  std::cerr << nNonZeroCells << " of " << nc << " are non-zero\n";

  if (nNonZeroCells > 0)
  {
    float *weights = new float[nNonZeroCells];
    int   *indices = new int[nNonZeroCells];

    for (int k = 0, i = 0; i < nc; i++)
    {
      float d = cell_centered_data[i];
      if (d > 0)
      {
        weights[k] = d;
        indices[k++] = i;
      }
    }

    isearch::IntervalSearchTree itree(weights, indices, nNonZeroCells);

    delete[] weights;
    delete[] indices;

    int limit = 1000;

    int clast = 0;
    int i = 0, created = 0;
    while (created + nRetained <= ns)
    {
      if (created > clast && (created % 1000) == 0)
      {
        std::cerr << std::dec << created << "(" << i << ")\n";
        clast = created;
      }

      if (i > limit*created)
      {
        std::cerr << "looping stopped when the proportion of accepts fell below 1/" << limit << " of the attempts\n";
        break;
      }

      int indx = itree.Search((rand() / float(RAND_MAX)) * itree.max);

      int ci = indx % xc;
      int cj = (indx % yxc) / xc;
      int ck = indx / yxc;

      float dx = rand() / float(RAND_MAX);
      float dy = rand() / float(RAND_MAX);
      float dz = rand() / float(RAND_MAX);

      Point p(origin[0] + (ci+dx)*spacing[0], origin[1] + (cj+dy)*spacing[1], origin[2] + (ck+dz)*spacing[2]);

      float w111 = (  dx)*(  dy)*(  dz);
      float w110 = (  dx)*(  dy)*(1-dz);
      float w101 = (  dx)*(1-dy)*(  dz);
      float w100 = (  dx)*(1-dy)*(1-dz);
      float w011 = (1-dx)*(  dy)*(  dz);
      float w010 = (1-dx)*(  dy)*(1-dz);
      float w001 = (1-dx)*(1-dy)*(  dz);
      float w000 = (1-dx)*(1-dy)*(1-dz);

      float d = w000*point_centered_data[PINDEX(ci+0,cj+0,ck+0)] +
                w001*point_centered_data[PINDEX(ci+0,cj+0,ck+1)] +
                w010*point_centered_data[PINDEX(ci+0,cj+1,ck+0)] +
                w011*point_centered_data[PINDEX(ci+0,cj+1,ck+1)] +
                w100*point_centered_data[PINDEX(ci+1,cj+0,ck+0)] +
                w101*point_centered_data[PINDEX(ci+1,cj+0,ck+1)] +
                w110*point_centered_data[PINDEX(ci+1,cj+1,ck+0)] +
                w111*point_centered_data[PINDEX(ci+1,cj+1,ck+1)];

      float r = R0 + d*(R1 - R0);

      if (! kdtree.search(p, r))
      {
        samples.push_back(p);
        kdtree.insert(p);
        created++;
      }

      i++;
    }
  }

  if (point_centered_data)
    delete[] cell_centered_data;
  else
    delete[] point_centered_data;

  return;
}

PointD
SampleTet(int* cells, float* points, float *pcd, int i)
{
  int *tet = cells + (i*5);

  struct Point v0(points + 3*tet[1]);
  struct Point v1(points + 3*tet[2]);
  struct Point v2(points + 3*tet[3]);
  struct Point v3(points + 3*tet[4]);

  double s = RAND;
  double t = RAND;
  double u = RAND;

  if (s+t > 1.0)
  {
    s = 1.0 - s;
    t = 1.0 - t;
  }

  if (t+u > 1.0)
  {
    double tmp = u;
    u = 1.0 - s - t;
    t = 1.0 - tmp;
  }
  else if (s+t+u > 1.0)
  {
    double tmp = u;
    u = s + t + u - 1.0;
    s = 1 - t - tmp;
  }

  double a = 1 - (s+t+u);

  float x = v0.x()*a + v1.x()*s + v2.x()*t + v3.x()*u;
  float y = v0.y()*a + v1.y()*s + v2.y()*t + v3.y()*u;
  float z = v0.z()*a + v1.z()*s + v2.z()*t + v3.z()*u;
  float d = pcd[tet[1]]*a + pcd[tet[2]]*s + pcd[tet[3]]*t + pcd[tet[4]]*u;

  return PointD(x, y, z, d);
}

static void
ComputeTetCellData(int nCells, float *point_centered_data, int *cells, float *points, float*& cell_centered_data, float*& volumes)
{
  volumes = new float[nCells];

  if (point_centered_data)
    cell_centered_data = new float[nCells];

  float *v    = volumes;
  float *ccd  = cell_centered_data;
  int   *cell = cells;
  for (int nv, i = 0; i < nCells; i++, cell += 5)
  {
    int p = cell[1];
    int q = cell[2];
    int r = cell[3];
    int s = cell[4];

    struct Point pp(points + 3*p);
    struct Point pq(points + 3*q);
    struct Point pr(points + 3*r);
    struct Point ps(points + 3*s);
    
    if (point_centered_data)
      *ccd++ = (point_centered_data[p] + point_centered_data[q] + point_centered_data[r] + point_centered_data[s]) / 4;
    
    struct Point a = pq - pp;
    struct Point b = pr - pp;
    struct Point c = ps - pp;
    
    *v++ = Point::dot(a, Point::cross(b, c));
  }
}

extern "C" void
SampleTetrahedra(int ns, int nCells, float R0, float R1, float *points, int *cells, int pdep, float *data, int nRetained, float *retained_points, float *retained_data)
{
  float  *point_centered_data = NULL;
  float  *cell_centered_data = NULL;
  float  *volumes   = NULL;

  std::cerr << "SampleTets! ns = " << ns << "\n";
  srand(12345);

  if (pdep) 
    point_centered_data = data;
  else
    cell_centered_data = data;

  // Compute cell volume and, if necessary, cell-centered data

  ComputeTetCellData(nCells, point_centered_data, cells, points, cell_centered_data, volumes);

  KDNode kdtree;

  int kept = 0;
  for (int i = 0; i < nRetained; i++)
  {
    float d = retained_data[i];
    if (d > 0)
    {
    
      Point p(retained_points+3*i);
      float r = R0 + d*(R1 - R0);
      if (! kdtree.search(p, r))
      {
        samples.push_back(p);
        kdtree.insert(p);
        kept ++;
      }
    }
  }

  int nNonZeroCells = 0;
  for (int i = 0; i < nCells; i++)
  {
    float d = cell_centered_data[i];
    if (d > 0)
      nNonZeroCells ++;
  }

  std::cerr << nNonZeroCells << " of " << nCells << " are non-zero\n";

  if (nNonZeroCells > 0)
  {
    float *weights = new float[nNonZeroCells];
    int   *indices = new int[nNonZeroCells];

    for (int k = 0, i = 0; i < nCells; i++)
    {
      float d = cell_centered_data[i];
      if (d > 0)
      {
        weights[k] = d * volumes[i];
        indices[k++] = i;
      }
    }

    isearch::IntervalSearchTree itree(weights, indices, nNonZeroCells);

    delete[] weights;
    delete[] indices;

    int limit = 1000;

    int clast = 0;
    int i = 0, created = 0;
    while (created + nRetained <= ns)
    {
      if (created > clast && (created % 1000) == 0)
      {
        std::cerr << std::dec << created << "(" << i << ")\n";
        clast = created;
      }

      if (i > limit*created)
      {
        std::cerr << "looping stopped when the proportion of accepts fell below 1/" << limit << " of the attempts\n";
        break;
      }

      int indx = itree.Search((rand() / float(RAND_MAX)) * itree.max);
      float d  = cell_centered_data[indx];

      float r = R0 + d*(R1 - R0);
      PointD pd = SampleTet(cells, points, point_centered_data, indx);

      if (! kdtree.search(pd, r))
      {
        samples.push_back(pd);
        kdtree.insert(pd);
        created++;
      }

      i++;
    }
  }

  delete[] volumes;
  if (point_centered_data)
    delete[] cell_centered_data;

  return;
}

extern "C" int
GetNumberOfSamples()
{
  return samples.size();
}

extern "C" void
GetSamples(float *d)
{
  memcpy(d, samples.data(), samples.size()*sizeof(Point));
  samples.clear();
}
