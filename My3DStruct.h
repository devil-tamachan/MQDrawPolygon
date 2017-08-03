#ifndef _MY3DSTRUCT_
#define _MY3DSTRUCT_

#include <stdio.h>

class MyCoordinate
{
public:
  MyCoordinate()
  {
    u = 0.0;
    v = 0.0;
  }
  MyCoordinate(float _u, float _v)
  {
    u = _u;
    v = _v;
  }
  
  float u;
  float v;
};

class MyPoint
{
public:
  MyPoint()
  {
  }
  MyPoint(float _x, float _y, float _z)
  {
    x = _x;
    y = _y;
    z = _z;
  }
  MyPoint operator+(const MyPoint& r)
  {
    return MyPoint(x+r.x, y+r.y, z+r.z);
  }
  MyPoint operator-(const MyPoint& r)
  {
    return MyPoint(x-r.x, y-r.y, z-r.z);
  }
  MyPoint operator*(const float& r)
  {
    return MyPoint(x*r, y*r, z*r);
  }
  MyPoint operator/(const float& r)
  {
    return MyPoint(x/r, y/r, z/r);
  }
  
  float x;
  float y;
  float z;
  
  float length()
  {
    return sqrt(x*x + y*y + z*z);
  }
};

MyPoint CrossProduct(MyPoint a, MyPoint b)
{
  return MyPoint(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

MyPoint Normalize(MyPoint a)
{
  float len = a.length();
  if(len==0.0f)return MyPoint(0.0f, 0.0f, 0.0f);
  return a/len;
}

class _MyObject
{
public:
  
  int AddVertex(MyPoint p)
  {
    verts.push_back(p);
    int newidx = verts.size()-1;
    return newidx;
  }
  
  int AddFace(int num, int *idx)
  {
    std::vector<int> vidx;
    for(int i=0;i<num;i++)vidx.push_back(idx[i]);
    if(vidx.size()>0)
    {
      faces.push_back(vidx);
      coords.push_back(std::vector<MyCoordinate>());
      materials.push_back(-1);
      int newfidx = faces.size()-1;
      return newfidx;
    } else {
      return -1;
    }
  }
  MyPoint GetVertex(int vi)
  {
    return verts[vi];
  }
  void GetFacePointArray(int fi, int *vidx)
  {
    if(fi>=faces.size())return;
    std::vector<int> &fidx = faces[fi];
    int numV = fidx.size();
    for(int i=0;i<numV;i++)
    {
      vidx[i] = fidx[i];
    }
  }
  void GetFacePointArray(int fi, std::vector<int> &vidx)
  {
    if(fi>=faces.size())return;
    std::vector<int> &fidx = faces[fi];
    vidx = fidx;
  }
  void GetFaceMyPointArray(int fi, std::vector<MyPoint> &vp)
  {
    if(fi>=faces.size())return;
    std::vector<int> &fidx = faces[fi];
    int numV = fidx.size();
    for(int vi=0;vi<numV;vi++)
    {
      vp.push_back(GetVertex(fidx[vi]));
    }
  }
  void GetFaceCoordinateArray(int fi, MyCoordinate *_coords)
  {
    if(fi>=coords.size())return;
    std::vector<MyCoordinate> &fc = coords[fi];
    int numV = fc.size();
    for(int i=0;i<numV;i++)
    {
      _coords[i] = fc[i];
    }
  }
  void GetFaceCoordinateArray(int fi, std::vector<MyCoordinate> &_coords)
  {
    if(fi>=coords.size())return;
    std::vector<MyCoordinate> &fc = coords[fi];
    _coords = fc;
  }
  
  void SetFaceCoordinateArray(int fi, MyCoordinate *c)
  {
    if(fi>=coords.size())return;
    int num = faces[fi].size();
    std::vector<MyCoordinate> uv;
    for(int i=0;i<num;i++)
    {
      uv.push_back(c[i]);
    }
    coords[fi] = uv;
  }
  
  void SetFaceMaterial(int fi, int matId)
  {
    if(fi>=materials.size())return;
    materials[fi] = matId;
  }
  
  
  int GetFaceCount()
  {
    return faces.size();
  }
  
  int GetVertexCount()
  {
    return verts.size();
  }
  
  int GetFacePointCount(int fi)
  {
    if(fi>=faces.size() || fi<0)return -1;
    return faces[fi].size();
  }
  int GetFaceMaterial(int fi)
  {
    if(fi>=materials.size())return -1;
    return materials[fi];
  }
  
  std::vector<MyPoint> verts;
  std::vector<std::vector<MyCoordinate> > coords;
  std::vector<std::vector<int> > faces;
  std::vector<int> materials;
  
  bool write(const char *path)
  {
    FILE *fp = fopen(path, "wb");
    if(fp==NULL)return false;
    bool bRet = write(fp);
    fclose(fp);
    return bRet;
  }
  bool write(const wchar_t *path)
  {
    FILE *fp = _wfopen(path, L"wb");
    if(fp==NULL)return false;
    bool bRet = write(fp);
    fclose(fp);
    return bRet;
  }
  bool write(FILE *fp)
  {
    if(coords.size() != faces.size() || materials.size() != faces.size())return false;
    unsigned int numV = verts.size();
    unsigned int numF = faces.size();
    
    if(fputc('P', fp)==EOF)return false;
    if(fputc('O', fp)==EOF)return false;
    if(fputc('L', fp)==EOF)return false;
    if(fputc('Y', fp)==EOF)return false;
    if(fputc(0x0, fp)==EOF)return false;
    
    if(fwrite(&numV, 4, 1, fp)!=1)return false;
    if(fwrite(&numF, 4, 1, fp)!=1)return false;
    
    for(int i=0;i<numV;i++)
    {
      MyPoint &p = verts[i];
      float x = p.x, y = p.y, z = p.z;
      if(fwrite(&x, sizeof(float), 1, fp)!=1)return false;
      if(fwrite(&y, sizeof(float), 1, fp)!=1)return false;
      if(fwrite(&z, sizeof(float), 1, fp)!=1)return false;
    }
    for(int fi=0;fi<numF;fi++)
    {
      std::vector<int> &f = faces[fi];
      unsigned int n = f.size();
      if(fwrite(&n, 4, 1, fp)!=1)return false;
      
      for(int k=0;k<n;k++)
      {
        int idx = f[k];
        if(fwrite(&idx, sizeof(int), 1, fp)!=1)return false;
      }
      
      std::vector<MyCoordinate> &c = coords[fi];
      
      if(f.size()!=c.size())
      {
        n = 0;
        if(fwrite(&n, 4, 1, fp)!=1)return false;
        continue;
      }
      n = c.size();
      if(fwrite(&n, 4, 1, fp)!=1)return false;
      
      for(int k=0;k<n;k++)
      {
        MyCoordinate &uv = c[k];
        if(fwrite(&(uv.u), sizeof(float), 1, fp)!=1)return false;
        if(fwrite(&(uv.v), sizeof(float), 1, fp)!=1)return false;
      }
      int matid = materials[fi];
      if(fwrite(&matid, 4, 1, fp)!=1)return false;
    }
    return true;
  }
  
  bool read(const char *path)
  {
    FILE *fp = fopen(path, "rb");
    if(fp==NULL)return false;
    bool bRet = read(fp);
    fclose(fp);
    return bRet;
  }
  bool read(const wchar_t *path)
  {
    FILE *fp = _wfopen(path, L"rb");
    if(fp==NULL)return false;
    bool bRet = read(fp);
    fclose(fp);
    return bRet;
  }
  
  bool read(FILE *fp)
  {
    unsigned int numV = 0;
    unsigned int numF = 0;
    
    if(fgetc(fp)!='P')return false;
    if(fgetc(fp)!='O')return false;
    if(fgetc(fp)!='L')return false;
    if(fgetc(fp)!='Y')return false;
    if(fgetc(fp)!=0x0)return false;
    
    if(fread(&numV, 4, 1, fp)!=1)return false;
    if(fread(&numF, 4, 1, fp)!=1)return false;
    
    for(int i=0;i<numV;i++)
    {
      float x, y, z;
      if(fread(&x, sizeof(float), 1, fp)!=1)return false;
      if(fread(&y, sizeof(float), 1, fp)!=1)return false;
      if(fread(&z, sizeof(float), 1, fp)!=1)return false;
      verts.push_back(MyPoint(x,y,z));
    }
    for(int i=0;i<numF;i++)
    {
      std::vector<int> f;
      unsigned int numVI = 0;
      if(fread(&numVI, 4, 1, fp)!=1)return false;
      
      for(int k=0;k<numVI;k++)
      {
        int idx;
        if(fread(&idx, sizeof(int), 1, fp)!=1)return false;
        f.push_back(idx);
      }
      faces.push_back(f);
      
      std::vector<MyCoordinate> c;
      
      unsigned int numCoords = 0;
      if(fread(&numCoords, 4, 1, fp)!=1)return false;
      if(numVI!=numCoords)return false;
      
      for(int k=0;k<numCoords;k++)
      {
        MyCoordinate uv;
        if(fread(&(uv.u), sizeof(float), 1, fp)!=1)return false;
        if(fread(&(uv.v), sizeof(float), 1, fp)!=1)return false;
        c.push_back(uv);
      }
      coords.push_back(c);
      
      int matid = -1;
      if(fread(&matid, 4, 1, fp)!=1)return false;
      materials.push_back(matid);
    }
    return true;
  }
};

typedef _MyObject* MyObject;


#endif //_MY3DSTRUCT_