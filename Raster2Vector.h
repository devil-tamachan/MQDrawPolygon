
#ifndef _TAMARASTER2VECTOR_
#define _TAMARASTER2VECTOR_

#include "MyPolygon.h"
//#include "MQ3DLib.h"
#include <Eigen/Dense>
#include <math.h>
#include "My3DStruct.h"

Vector3 MyPointToVector3(MyPoint mp)
{
  return Vector3(mp.x, mp.y, mp.z);
}
void MyPointToVector3(MyPoint *mpArr, Vector3 *vArr, int num)
{
  for(int i=0;i<num;i++)
  {
    MyPoint &mp = mpArr[i];
    vArr[i] = Vector3(mp.x, mp.y, mp.z);
  }
}
MyPoint Vector3ToMyPoint(Vector3 v)
{
  double x = CGAL::to_double(v.x());
  double y = CGAL::to_double(v.y());
  double z = CGAL::to_double(v.z());
  return MyPoint(x, y, z);
}
MyPoint Point3ToMyPoint(Point3 p)
{
  double x = CGAL::to_double(p.x());
  double y = CGAL::to_double(p.y());
  double z = CGAL::to_double(p.z());
  return MyPoint(x, y, z);
}
/*
Point2 MyPointToPoint2XY(MyPoint mp)
{
  return Point2(mp.x, mp.y);
}
Point3 MyPointToPoint3(MyPoint mp)
{
  return Point3(mp.x, mp.y, mp.z);
}*/

double GetRatioStart2Pt(Point2 &start, Point2 &end, Point2 &pt)
{
  double distance = sqrt(CGAL::squared_distance(start, end));
  double distance2= sqrt(CGAL::squared_distance(start, pt));
  return distance2 / distance;
}

double lerp(double start, double end, double ratio)
{
  double diff = end - start;
  return start+diff*ratio;
}
Point2 Point3_Point2XY(Point3 &p)
{
  return Point2(p.x(), p.y());
}
MyPoint Point2_MyPoint(Point3 &pS, Point3 &pE, Point2 &pt)
{
  Point2 ptStart = Point3_Point2XY(pS);
  Point2 ptEnd   = Point3_Point2XY(pE);
  double ratioPt = GetRatioStart2Pt(ptStart, ptEnd, pt);
  double z = lerp(pS.z(), pE.z(), ratioPt);
  return MyPoint(pt.x(), pt.y(), z);
}
MyPoint Point2_MyPoint(Point2 &pt, Eigen::Matrix<double, 1, 3> &matXYToZ)
{
  Eigen::Matrix<double, 3, 1> XY;
  XY << pt.x(),
        pt.y(),
        1.0;
  Eigen::Matrix<double, 1, 1> Z = matXYToZ * XY;
  return MyPoint(pt.x(), pt.y(), Z(0,0));
}

//#define DEBUGCGALINTERSECT

class CRaster2Vector
{
private:
  int opt_edgeproc;
  int opt_vectorconv;
  double opt_threshold;
  int opt_np;
  int opt_thresholdMennuki;
#ifdef DEBUGCGALINTERSECT
  FILE *fp;
  bool bFirst;
#endif
  
  
  std::vector<Point2> otr2_v;
  std::vector<Segment2> otr2_line;
  cv::Mat chkFace;
  
  CRaster2Vector() { }
  
public:
  CRaster2Vector(DrawPolygonOptions &opt)
  {
    opt_edgeproc = opt.edgeproc;
    opt_vectorconv = opt.vectorconv;
    opt_threshold = opt.threshold;
    opt_np = opt.np;
    opt_thresholdMennuki = opt.thresholdMennuki * 255.0;
#ifdef DEBUGCGALINTERSECT
    fp = fopen("c:\\tmp\\cgal.bin", "wb");
    bFirst = true;
#endif
  }
  
  void reset()
  {
    std::vector<Point2>().swap(otr2_v); //clear
    std::vector<Segment2>().swap(otr2_line); //clear
    chkFace.release();
  }

  bool Raster2Vector(cv::Mat &mat)
  {
    reset();
    
    //cv::Mat matf;
    //mat.convertTo(matf, CV_32F, 1.0/255);

    cv::Mat dstEdge;
    Edge(mat, dstEdge);

    //  std::vector<Point2> otr2_v;
    //  std::vector<Segment2> otr2_line;

    Point_property_map point_pmap;
    Mass_property_map  mass_pmap;
    PointMassList pointsMass;
    std::vector<Point2> points;
    chkFace = cv::Mat(dstEdge.rows, dstEdge.cols, CV_8U, cv::Scalar(255));
    cv::Mat srcGray;
    cvtColor(mat, srcGray, CV_RGB2GRAY);

    switch(opt_vectorconv)
    {
    case 0://モノクロ値を重みとして使用
      {
        if(opt_edgeproc==3)chkFace = srcGray;
        else CheckFaceThreshold(254, chkFace, srcGray); //カラーだと色によって判定が異なってしまう。Texがモノクロだと逆にノイズの原因になっている
        ToPolygonMass(dstEdge, pointsMass);
        if(pointsMass.size()==0)return false;
        Otr_2 otr2(pointsMass, point_pmap, mass_pmap, 15);
        otr2.set_random_sample_size(15);
#ifdef MYOUTPUTDEBUG
char s[1025];
int dbg = MIN(opt_np, pointsMass.size());
sprintf(s, "%d\n", dbg);
  MyOutputDebugStringA(s);
#endif
        otr2.run_until(MIN(opt_np, pointsMass.size()));
        otr2.list_output(std::back_inserter(otr2_v), std::back_inserter(otr2_line));
      }
      break;
    case 1://２値化
      {
        CheckFaceThreshold(opt_threshold*255.0, chkFace, srcGray);
        ToPolygonThreshold(dstEdge, points);
        if(points.size()==0)return false;
        Otr otr(points);
        otr.set_random_sample_size(15);
        otr.run_until(MIN(opt_np, points.size()));
        otr.list_output(std::back_inserter(otr2_v), std::back_inserter(otr2_line));
      }
      break;
    default:
      return false;
    }

    return true;
  }
  
  void AddPoint(Vector3 &p)
  {
    otr2_v.push_back(Point2(p.x(), p.y()));
  }
  void AddPoint(Point2 &p)
  {
    otr2_v.push_back(p);
  }
  void AddPoints(Vector3 *tri2d, int num)
  {
    for(int i=0;i<num;i++)
    {
      AddPoint(tri2d[i]);
    }
  }
  void AddPoints(Point2 *tri2d, int num)
  {
    for(int i=0;i<num;i++)
    {
      AddPoint(tri2d[i]);
    }
  }
  
  bool isContainedAll(Point2 *pUvRect, std::vector<Point2> &arrPt2)
  {
    int iFound = 0;
    for(int i=0;i<4;i++)
    {
      Point2 &p = pUvRect[i];
      for(int j=0;j<4;j++)
      {
        Point2 &p2 = arrPt2[j];
        if(p==p2)iFound++;
      }
      if(iFound!=i+1)return false;
    }
    return iFound==4;
  }
  
  void MakeTexRepeatList(std::vector<MyCoordinate> &coordTri, CacheTexInfo &info, std::vector<Point2> &retClip, std::vector<Point2> &retThru)
  {
    float umin, umax, vmin, vmax;
    umin = vmin = FLT_MAX;
    umax = vmax = FLT_MIN;
    float uminTex, vminTex, umaxTex, vmaxTex;
    uminTex = info.minx / info.w;
    vminTex = 1.0-(info.maxy / info.h);
    umaxTex = info.maxx / info.w;
    vmaxTex = 1.0-(info.miny / info.h);
    for(int i=0;i<3;i++)
    {
      MyCoordinate &c = coordTri[i];
      float u = c.u;
      float v = c.v;
      umin = MIN(u, umin);
      vmin = MIN(v, vmin);
      umax = MAX(u, umax);
      vmax = MAX(v, vmax);
    }
    
    int ustart, uend, vstart, vend;
    ustart = floor(floor(umin) + uminTex);
    vstart = floor(floor(vmin) + vminTex);
    uend = ceil(floor(umax) + umaxTex);
    vend = ceil(floor(vmax) + vmaxTex);
    Triangle2 ClipTri = _MakeClipTri(coordTri);
    
    for(int u=ustart;u<=uend;u++)
    {
      for(int v=vstart;v<=vend;v++)
      {
        float u2, u3, u4, v2, v3, v4;
        u2 = u3 = u4 = u;
        v2 = v3 = v4 = v;
        u2 += uminTex;
        u3 += umaxTex;
        v2 += vminTex;
        v3 += vmaxTex;
        u4 *= info.w;
        v4 *= info.h;
        Point2 pUvRect[4];
        pUvRect[0] = Point2(u2, v2);
        pUvRect[1] = Point2(u3, v3);
        pUvRect[2] = Point2(u2, v3);
        pUvRect[3] = Point2(u3, v2);
        Iso_rectangle2 uvRect(pUvRect[0], pUvRect[1]);
        //if(CGAL::do_intersection(uvRect, ClipTri)) ret.push_back(Point2(u4,v4));
        auto result = CGAL::intersection(uvRect, ClipTri);
        if(result)
        {
          if(std::vector<Point2> *arrPt2 = boost::get<std::vector<Point2> >(&*result))
          {
            if(arrPt2->size()==4)
            {
              if(isContainedAll(pUvRect, *arrPt2))
              {
                retThru.push_back(Point2(u4,v4));
                continue;
              }
            }
          }
          retClip.push_back(Point2(u4,v4));
        }
      }
    }
  }
  
  MyCoordinate CalcNewUV(MyPoint &newMqp, Eigen::Matrix<double, 2, 4> &matPtToUV)
  {
    Eigen::Matrix<double, 4, 1> matNewMqp;
    matNewMqp << newMqp.x,
                 newMqp.y,
                 newMqp.z,
                 1.0;
    Eigen::Matrix<double, 2, 1> uv = matPtToUV * matNewMqp;
    return MyCoordinate(uv(0,0), uv(1,0));
  }
  /*
  Eigen::Matrix<double, 3, 4> __GetAffineTri2UV(std::vector<MyPoint> &pts, std::vector<MyCoordinate> &coord)
  {
    MyPoint p0 = pts[0], p1 = pts[1], p2 = pts[2];
    MyCoordinate c0 = coord[0], c1 = coord[1], c2 = coord[2];
    Eigen::Matrix<double, 9, 12> A;
    A << 
         p0.x, p0.y, p0.z, 1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  p0.x, p0.y, p0.z, 1.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  p0.x, p0.y, p0.z, 1.0,
         p1.x, p1.y, p1.z, 1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  p1.x, p1.y, p1.z, 1.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  p1.x, p1.y, p1.z, 1.0,
         p2.x, p2.y, p2.z, 1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  p2.x, p2.y, p2.z, 1.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  p2.x, p2.y, p2.z, 1.0;
    Eigen::Matrix<double, 9, 1> B;
    B << 
         c0.u,
         c0.v,
         0.0,
         c1.u,
         c1.v,
         0.0,
         c2.u,
         c2.v,
         0.0;
    Eigen::Matrix<double, 12, 1> X = A.ldlt().solve(B);
    Eigen::Matrix<double, 3, 4> M = Eigen::Map<Eigen::Matrix<double, 3, 4, Eigen::RowMajor>>(X.data());
    return M;
  }*/
  
  Eigen::Matrix<double, 2, 4> __GetAffineTri2UV(std::vector<MyPoint> &pts, std::vector<MyCoordinate> &coord)
  {
    MyPoint p0 = pts[0], p1 = pts[1], p2 = pts[2];
    MyCoordinate c0 = coord[0], c1 = coord[1], c2 = coord[2];
    Eigen::Matrix<double, 6, 8> A;
    A << 
         p0.x, p0.y, p0.z, 1.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  p0.x, p0.y, p0.z, 1.0,
         p1.x, p1.y, p1.z, 1.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  p1.x, p1.y, p1.z, 1.0,
         p2.x, p2.y, p2.z, 1.0,  0.0,  0.0,  0.0,  0.0,
         0.0,  0.0,  0.0,  0.0,  p2.x, p2.y, p2.z, 1.0;
    Eigen::Matrix<double, 6, 1> B;
    B << 
         c0.u,
         c0.v,
         c1.u,
         c1.v,
         c2.u,
         c2.v;
    Eigen::Matrix<double, 8, 1> X = A.fullPivLu().solve(B);
    Eigen::Matrix<double, 2, 4> M = Eigen::Map<Eigen::Matrix<double, 2, 4, Eigen::RowMajor> >(X.data());
    return M;
  }
  
  //入力ベクトルZが潰れるので注意。すべての入力Zが0だろうがなんだろうが0扱いになる。
  Eigen::Matrix<double, 3, 4> __GetAffineUV2Tri(std::vector<MyCoordinate> &coord, std::vector<MyPoint> &pts)
  {
    MyPoint p0 = pts[0], p1 = pts[1], p2 = pts[2];
    MyCoordinate c0 = coord[0], c1 = coord[1], c2 = coord[2];
    Eigen::Matrix<double, 9, 12> A;
    A << 
         c0.u, c0.v, 0.0, 1.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,
         0.0,  0.0,  0.0, 0.0,  c0.u, c0.v, 0.0, 1.0,  0.0,  0.0,  0.0, 0.0,
         0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,  c0.u, c0.v, 0.0, 1.0,
         c1.u, c1.v, 0.0, 1.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,
         0.0,  0.0,  0.0, 0.0,  c1.u, c1.v, 0.0, 1.0,  0.0,  0.0,  0.0, 0.0,
         0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,  c1.u, c1.v, 0.0, 1.0,
         c2.u, c2.v, 0.0, 1.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,
         0.0,  0.0,  0.0, 0.0,  c2.u, c2.v, 0.0, 1.0,  0.0,  0.0,  0.0, 0.0,
         0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,  c2.u, c2.v, 0.0, 1.0;
    Eigen::Matrix<double, 9, 1> B;
    B << 
         p0.x,
         p0.y,
         p0.z,
         p1.x,
         p1.y,
         p1.z,
         p2.x,
         p2.y,
         p2.z;
    Eigen::Matrix<double, 12, 1> X = A.fullPivLu().solve(B);
    Eigen::Matrix<double, 3, 4> M = Eigen::Map<Eigen::Matrix<double, 3, 4, Eigen::RowMajor> >(X.data());
    return M;
  }
  Eigen::Matrix<double, 3, 4> __GetAffine34UV2Tri_NoScaleZ(std::vector<MyCoordinate> &coord, std::vector<MyPoint> &pts)
  {
    MyPoint p0 = pts[0], p1 = pts[1], p2 = pts[2];
    MyCoordinate c0 = coord[0], c1 = coord[1], c2 = coord[2];
    Eigen::Matrix<double, 12, 12> A;
    A << 
         c0.u, c0.v, 0.0, 1.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,
         0.0,  0.0,  0.0, 0.0,  c0.u, c0.v, 0.0, 1.0,  0.0,  0.0,  0.0, 0.0,
         0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,  c0.u, c0.v, 0.0, 1.0,
         c1.u, c1.v, 0.0, 1.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,
         0.0,  0.0,  0.0, 0.0,  c1.u, c1.v, 0.0, 1.0,  0.0,  0.0,  0.0, 0.0,
         0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,  c1.u, c1.v, 0.0, 1.0,
         c2.u, c2.v, 0.0, 1.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,
         0.0,  0.0,  0.0, 0.0,  c2.u, c2.v, 0.0, 1.0,  0.0,  0.0,  0.0, 0.0,
         0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,  c2.u, c2.v, 0.0, 1.0,
         
         c1.u, c1.v, 10000.0, 1.0,  0.0,  0.0,  0.0,     0.0,  0.0,  0.0,  0.0,     0.0,
         0.0,  0.0,  0.0,     0.0,  c1.u, c1.v, 10000.0, 1.0,  0.0,  0.0,  0.0,     0.0,
         0.0,  0.0,  0.0,     0.0,  0.0,  0.0,  0.0,     0.0,  c1.u, c1.v, 10000.0, 1.0;
    
    MyPoint p3 = Normalize(CrossProduct(p0-p1, p2-p1))*10000.0 + p1;
    Eigen::Matrix<double, 12, 1> B;
    B << 
         p0.x,
         p0.y,
         p0.z,
         p1.x,
         p1.y,
         p1.z,
         p2.x,
         p2.y,
         p2.z,
         p3.x,
         p3.y,
         p3.z;
    Eigen::Matrix<double, 12, 1> X = A.fullPivLu().solve(B);
    Eigen::Matrix<double, 3, 4> M = Eigen::Map<Eigen::Matrix<double, 3, 4, Eigen::RowMajor> >(X.data());
    return M;
  }
  Eigen::Matrix<double, 4, 4> __GetAffine44UV2Tri_NoScaleZ(std::vector<MyCoordinate> &coord, std::vector<MyPoint> &pts)
  {
    Eigen::MatrixXd M = __GetAffine34UV2Tri_NoScaleZ(coord, pts);
    M.conservativeResize(4,4);
    M(3,3) = 1.0;
    return M;
  }
  Eigen::Matrix<double, 3, 4> __GetAffine34UV2Tri_ScaleZ(std::vector<MyCoordinate> &coord, std::vector<MyPoint> &pts, double scaleZ)
  {
    MyPoint p0 = pts[0], p1 = pts[1], p2 = pts[2];
    MyCoordinate c0 = coord[0], c1 = coord[1], c2 = coord[2];
    Eigen::Matrix<double, 12, 12> A;
    A << 
         c0.u, c0.v, 0.0, 1.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,
         0.0,  0.0,  0.0, 0.0,  c0.u, c0.v, 0.0, 1.0,  0.0,  0.0,  0.0, 0.0,
         0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,  c0.u, c0.v, 0.0, 1.0,
         c1.u, c1.v, 0.0, 1.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,
         0.0,  0.0,  0.0, 0.0,  c1.u, c1.v, 0.0, 1.0,  0.0,  0.0,  0.0, 0.0,
         0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,  c1.u, c1.v, 0.0, 1.0,
         c2.u, c2.v, 0.0, 1.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,
         0.0,  0.0,  0.0, 0.0,  c2.u, c2.v, 0.0, 1.0,  0.0,  0.0,  0.0, 0.0,
         0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,  c2.u, c2.v, 0.0, 1.0,
         
         c1.u, c1.v, 10000.0, 1.0,  0.0,  0.0,  0.0,     0.0,  0.0,  0.0,  0.0,     0.0,
         0.0,  0.0,  0.0,     0.0,  c1.u, c1.v, 10000.0, 1.0,  0.0,  0.0,  0.0,     0.0,
         0.0,  0.0,  0.0,     0.0,  0.0,  0.0,  0.0,     0.0,  c1.u, c1.v, 10000.0, 1.0;
    
    MyPoint p3 = Normalize(CrossProduct(p0-p1, p2-p1))*(10000.0*scaleZ) + p1;
    Eigen::Matrix<double, 12, 1> B;
    B << 
         p0.x,
         p0.y,
         p0.z,
         p1.x,
         p1.y,
         p1.z,
         p2.x,
         p2.y,
         p2.z,
         p3.x,
         p3.y,
         p3.z;
    Eigen::Matrix<double, 12, 1> X = A.fullPivLu().solve(B);
    Eigen::Matrix<double, 3, 4> M = Eigen::Map<Eigen::Matrix<double, 3, 4, Eigen::RowMajor> >(X.data());
    return M;
  }
  Eigen::Matrix<double, 4, 4> __GetAffine44UV2Tri_ScaleZ(std::vector<MyCoordinate> &coord, std::vector<MyPoint> &pts, double scaleZ)
  {
    Eigen::MatrixXd M = __GetAffine34UV2Tri_ScaleZ(coord, pts, scaleZ);
    M.conservativeResize(4,4);
    M(3,3) = 1.0;
    return M;
  }

  
  // XYが一緒だとZが潰れるので注意
  Eigen::Matrix<double, 1, 3> __GetAffineXY2Z(std::vector<MyPoint> &pts)
  {
    MyPoint p0 = pts[0], p1 = pts[1], p2 = pts[2];
    Eigen::Matrix<double, 3, 3> A;
    A << 
         p0.x, p0.y, 1.0,
         p1.x, p1.y, 1.0,
         p2.x, p2.y, 1.0;
    Eigen::Matrix<double, 3, 1> B;
    B << 
         p0.z,
         p1.z,
         p2.z;
    Eigen::Matrix<double, 3, 1> X = A.fullPivLu().solve(B);
    Eigen::Matrix<double, 1, 3> M = Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor> >(X.data());
    return M;
  }
  
  void ShiftTri(std::vector<MyPoint> &tri, Point2 &shift, std::vector<MyPoint> &shiftTri)
  {
    for(int i=0;i<3;i++)
    {
      MyPoint &p = tri[i];
      shiftTri.push_back(MyPoint(p.x+shift.x(), p.y+shift.y(), p.z));
    }
  }
  /*
  void OutputDebugTri(std::vector<MyPoint> &tri)
  {
    char buf[1025];
    sprintf(buf, "Tri = ");
    OutputDebugStringA(buf);
    for(int i=0;i<3;i++)
    {
      if(i!=0)sprintf(buf, ", ");
      MyPoint p = tri[i];
      sprintf(buf, "(%f,%f)", p.x, p.y);
      OutputDebugStringA(buf);
    }
    OutputDebugStringA("\n");
  }*/
  
  bool __SegmentTriangleIntersect25(Triangle2 &ClipTri, std::vector<MyPoint> &tri, std::vector<MyCoordinate> &coord, std::vector<MyPoint> &ret, std::vector<MyCoordinate> &retCoord)
  {
        // p0(p1) ------ p2
    Point2 p01 = Point2(tri[0].x, tri[0].y);
    Point2 p2  = Point2(tri[2].x, tri[2].y);
    CGAL::cpp11::result_of<K::Intersect_2(Segment2, Triangle2)>::type result = CGAL::intersection(Segment2(p01, p2), ClipTri);
    if(result)
    {
      if(Segment2 *s2 = boost::get<Segment2>(&*result))
      {
        Eigen::Matrix<double, 2, 4> matPtToUV = __GetAffineTri2UV(tri, coord);
        
        Point3 q0 = Point3(tri[0].x, tri[0].y, tri[0].z);
        Point3 q1 = Point3(tri[1].x, tri[1].y, tri[1].z);
        Point3 q2 = Point3(tri[2].x, tri[2].y, tri[2].z);
        
        Point2 s2s = s2->source();
        Point2 s2t = s2->target();
        bool bIn01 = CGAL::do_intersect(ClipTri, p01);
        bool bIn2  = CGAL::do_intersect(ClipTri, p2);
        
        MyPoint newMqp = Point2_MyPoint(q0, q2, s2s);
        ret.push_back(newMqp);
        MyCoordinate newCoord = CalcNewUV(newMqp, matPtToUV);
        retCoord.push_back(newCoord);
        
        if(CGAL::do_intersect(ClipTri, p2))
        {
          ret.push_back(tri[2]);
          retCoord.push_back(coord[2]);
        } else {
          MyPoint newMqp2 = Point2_MyPoint(q0, q2, s2t);
          ret.push_back(newMqp2);
          MyCoordinate newCoord2 = CalcNewUV(newMqp2, matPtToUV);
          retCoord.push_back(newCoord2);
          
          MyPoint newMqp3 = Point2_MyPoint(q1, q2, s2t);
          ret.push_back(newMqp3);
          MyCoordinate newCoord3 = CalcNewUV(newMqp3, matPtToUV);
          retCoord.push_back(newCoord3);
        }
        
        MyPoint newMqp4 = Point2_MyPoint(q1, q2, s2s);
        ret.push_back(newMqp4);
        MyCoordinate newCoord4 = CalcNewUV(newMqp4, matPtToUV);
        retCoord.push_back(newCoord4);
        return true;
      }
    }
    return false;
  }
  
  bool _TriangleTriangleZeroAreaIntersect25(Triangle2 &ClipTri, std::vector<MyPoint> &tri, std::vector<MyCoordinate> &coord, std::vector<MyPoint> &ret, std::vector<MyCoordinate> &retCoord)
  {
    bool bEq01 = tri[0].x==tri[1].x && tri[0].y==tri[1].y;
    bool bEq12 = tri[1].x==tri[2].x && tri[1].y==tri[2].y;
    if(bEq01 && bEq12)
    {
      // p0.xy==p1.xy==p2.xy
      if(CGAL::do_intersect(ClipTri, Point2(tri[0].x, tri[0].y)))
      {
        ret = tri;
        retCoord = coord;
        return true;
      }
      return false;
    }
    bool bEq02 = tri[0].x==tri[2].x && tri[0].y==tri[2].y;
    if(bEq01 && !bEq12 && !bEq02)
    {
      // p0(p1) ------ p2
      std::vector<MyPoint> tri2;
      std::vector<MyCoordinate> coord2;
      int idxShuffle[] = {1, 0, 2};
      for(int i=0;i<3;i++)
      {
        int k = idxShuffle[i];
        tri2.push_back(tri[k]);
        coord2.push_back(coord[k]);
      }
      return __SegmentTriangleIntersect25(ClipTri, tri2, coord2, ret, retCoord);
    } else if(!bEq01 && bEq12 && !bEq02) {
      // p1(p2) ------ p0
      std::vector<MyPoint> tri2;
      std::vector<MyCoordinate> coord2;
      int idxShuffle[] = {2, 1, 0};
      for(int i=0;i<3;i++)
      {
        int k = idxShuffle[i];
        tri2.push_back(tri[k]);
        coord2.push_back(coord[k]);
      }
      return __SegmentTriangleIntersect25(ClipTri, tri2, coord2, ret, retCoord);
      
    } else if(!bEq01 && !bEq12 && bEq02) {
      // p0(p2) ------ p1
      std::vector<MyPoint> tri2;
      std::vector<MyCoordinate> coord2;
      int idxShuffle[] = {0, 2, 1};
      for(int i=0;i<3;i++)
      {
        int k = idxShuffle[i];
        tri2.push_back(tri[k]);
        coord2.push_back(coord[k]);
      }
      return __SegmentTriangleIntersect25(ClipTri, tri2, coord2, ret, retCoord);
      
    } else {
      return false;
    }
  }
  
  bool TriangleTriangleIntersect25(Triangle2 &ClipTri, std::vector<MyPoint> &tri, std::vector<MyCoordinate> &coord, std::vector<MyPoint> &ret, std::vector<MyCoordinate> &retCoord)
  {
    bool bReverse = false;
    Triangle2 Tri2 = _MakeClipTri(tri, &bReverse);
    //if(ClipTri.area()==0.0)return false; OutputToMetaseqFast_LineTriに移動
    if(Tri2.area()==0.0)
    {
      return _TriangleTriangleZeroAreaIntersect25(ClipTri, tri, coord, ret, retCoord);
    }

    CGAL::cpp11::result_of<K::Intersect_2(Triangle2, Triangle2)>::type result = CGAL::intersection(Tri2, ClipTri);
    if(result)
    {
      if(Triangle2 *tri3 = boost::get<Triangle2>(&*result))
      {
        Eigen::Matrix<double, 2, 4> matPtToUV = __GetAffineTri2UV(tri, coord);
        Eigen::Matrix<double, 1, 3> matXYToZ = __GetAffineXY2Z(tri);
        
        for(int i=0;i<3;i++)
        {
          Point2 p = (*tri3)[i];
          MyPoint newMqp = Point2_MyPoint(p, matXYToZ);
          ret.push_back(newMqp);
          MyCoordinate newCoord = CalcNewUV(newMqp, matPtToUV);
          retCoord.push_back(newCoord);
        }
        if(!bReverse)
        {
          std::reverse(ret.begin(), ret.end());
          std::reverse(retCoord.begin(), retCoord.end());
        }
        return true;
      } else if(std::vector<Point2> *p2 = boost::get<std::vector<Point2> >(&*result))
      {
        Eigen::Matrix<double, 2, 4> matPtToUV = __GetAffineTri2UV(tri, coord);
        Eigen::Matrix<double, 1, 3> matXYToZ = __GetAffineXY2Z(tri);
        
        int numP2 = p2->size();
        for(int i=0;i<numP2;i++)
        {
          MyPoint newMqp = Point2_MyPoint((*p2)[i], matXYToZ);
          ret.push_back(newMqp);
          MyCoordinate newCoord = CalcNewUV(newMqp, matPtToUV);
          retCoord.push_back(newCoord);
        }
        if(!bReverse)
        {
          std::reverse(ret.begin(), ret.end());
          std::reverse(retCoord.begin(), retCoord.end());
        }
        return true;
      }
    }
    return false;
  }
  
  /*バグあり: クリップ範囲内(ClipTri内)にtriの角が１～２入っているとそこが抜けておかしくなる
  bool TriangleTriangleIntersect2_5(Triangle2 &ClipTri, std::vector<MyPoint> &tri, std::vector<MyCoordinate> &coord, std::vector<MyPoint> &ret, std::vector<MyCoordinate> &retCoord)
  {
    Point2 pt2[3];
    Segment2 seg2[3];
    CGAL::Bounded_side PtIsIn[3];
    //OutputDebugTri(tri);
    for(int i=0;i<3;i++)
    {
      pt2[i] = MyPointToPoint2XY(tri[i]);
      PtIsIn[i] = ClipTri.bounded_side(pt2[i]);
      Point2 &pt2Next = i+1>=3 ? pt2[0] : pt2[i+1];
      seg2[i] = Segment2(pt2[i], pt2Next);
    }
    int iIn = 0;
    for(int i=0;i<3;i++)
    {
      if(PtIsIn[i]!=CGAL::ON_UNBOUNDED_SIDE)iIn++;
    }
    if(iIn==3)
    {
      for(int i=0;i<3;i++)
      {
        MyPoint &mqp = tri[i];
        ret.push_back(mqp);
        retCoord.push_back(coord[i]);
      }
      return true;
    }
    for(int i=0;i<3;i++)
    {
      int k = i+1>=3 ? 0:i+1;
      if(PtIsIn[i]!=CGAL::ON_UNBOUNDED_SIDE && PtIsIn[k]!=CGAL::ON_UNBOUNDED_SIDE)
      {
        ret.push_back(tri[k]);
        retCoord.push_back(coord[k]);
        continue;
      }
      Segment2 seg(pt2[i], pt2[k]);
      CGAL::cpp11::result_of<K::Intersect_2(Segment2, Triangle2)>::type result = CGAL::intersection(seg, ClipTri);
      if(result)
      {
        if(Segment2 *s2 = boost::get<Segment2>(&*result))
        {
          Point3 ptI = MyPointToPoint3(tri[i]);
          Point3 ptK = MyPointToPoint3(tri[k]);
          Eigen::Matrix<double, 2, 4> matPtToUV = __GetAffineTri2UV(tri, coord);
          if(PtIsIn[i]==CGAL::ON_UNBOUNDED_SIDE)
          {
            Point2 s2s = s2->source();
            MyPoint newMqp = Point2_MyPoint(ptI, ptK, s2s);
            ret.push_back(newMqp);
            MyCoordinate newCoord = CalcNewUV(newMqp, matPtToUV);
            retCoord.push_back(newCoord);
          }
          Point2 s2t = s2->target();
          MyPoint newMqp = Point2_MyPoint(ptI, ptK, s2t);
          ret.push_back(newMqp);
          MyCoordinate newCoord = CalcNewUV(newMqp, matPtToUV);
          retCoord.push_back(newCoord);
          continue;
        } else if(Point2 *p2 = boost::get<Point2>(&*result))
        {
          if(PtIsIn[k]!=CGAL::ON_UNBOUNDED_SIDE)
          {
            ret.push_back(tri[k]);
            retCoord.push_back(coord[k]);
          }
          continue;
        }
      }
    }
    if(ret.size()<3)
    {
      OutputDebugStringA("Intersect: p<3");
      ret.clear();
      retCoord.clear();
    } else OutputDebugStringA("Intersect: p>=3");
    return true;
  }
  */
  
  void CoordResize(std::vector<MyCoordinate> &coord, double uScale, double vScale, std::vector<MyCoordinate> &ret)
  {
    int num = coord.size();
    for(int i=0;i<num;i++)
    {
      MyCoordinate &c = coord[i];
      ret.push_back(MyCoordinate(c.u * uScale, c.v * vScale));
    }
  }
  
  /*
  void OutputDebugClipTri(Triangle2 &ClipTri)
  {
    char buf[1025];
    sprintf(buf, "ClipTri = ");
    OutputDebugStringA(buf);
    for(int i=0;i<3;i++)
    {
      if(i!=0)sprintf(buf, ", ");
      Point2 p = ClipTri[i];
      sprintf(buf, "(%f,%f)", (double)(p.x()), (double)(p.y()));
      OutputDebugStringA(buf);
    }
    OutputDebugStringA("\n");
  }
  */
  void FlipShiftY(std::vector<MyCoordinate> &coordTri, float shiftY, std::vector<MyCoordinate> &ret)
  {
    int num = coordTri.size();
    for(int i=0;i<num;i++)
    {
      MyCoordinate &c = coordTri[i];
      ret.push_back(MyCoordinate(c.u, -c.v+shiftY));
    }
  }
  
  Eigen::Matrix<double, 4, 4> flipY_Eigen44()
  {
    Eigen::Matrix<double, 4, 4> M;
    M <<
        1.0,  0.0, 0.0, 0.0,
        0.0, -1.0, 0.0, 0.0,
        0.0,  0.0, 1.0, 0.0,
        0.0,  0.0, 0.0, 1.0;
    return M;
  }
  Eigen::Matrix<double, 4, 4> translateY_Eigen44(double shiftY)
  {
    Eigen::Matrix<double, 4, 4> M;
    M <<
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, shiftY,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0;
    return M;
  }

  void OutputToMetaseqFast_LineTri(MyObject mqoOut, CacheTexInfo &info, std::vector<MyPoint> &ptsTri, std::vector<MyCoordinate> &coordTri, int modeZScale, double zscale)
  {
    MyObject mqoTex = &(info.o);
    //TriangulateObj(doc, mqoTex);MQTexManagerに移動
    
    std::vector<Point2> uvRepeatClip;
    std::vector<Point2> uvRepeatThru;
    MakeTexRepeatList(coordTri, info, uvRepeatClip, uvRepeatThru);
    int numRepeatClip = uvRepeatClip.size();
    int numRepeatThru = uvRepeatThru.size();
    std::vector<MyCoordinate> coordTriScale;
    CoordResize(coordTri, info.w, info.h, coordTriScale);
    std::vector<MyCoordinate> coordTriScaleFlipShiftY;
    FlipShiftY(coordTriScale, info.h, coordTriScaleFlipShiftY);
    Triangle2 ClipTri = _MakeClipTri(coordTriScaleFlipShiftY);
    //OutputDebugClipTri(ClipTri);
    
    if(ClipTri.area()==0.0)return;
    
    int numTexF = mqoTex->GetFaceCount();
    for(int fiTex=0;fiTex<numTexF;fiTex++)
    {
      int numTexFV = mqoTex->GetFacePointCount(fiTex);
      if(numTexFV!=3)continue;
      std::vector<MyPoint> ptsTex(numTexFV);
      std::vector<MyCoordinate> coordTex(numTexFV);
      GetPointAndCoord(mqoTex, fiTex, numTexFV, ptsTex, coordTex);
      
      int mqTexMatIdx = mqoTex->GetFaceMaterial(fiTex);
      
      Eigen::Matrix<double, 4, 4> invMat44;
      switch(modeZScale)
      {
      case 0:
        invMat44 = __GetAffine44UV2Tri_NoScaleZ(coordTriScale, ptsTri) * translateY_Eigen44(info.h) * flipY_Eigen44();
        break;
      case 1:
        {
          double a1 = Triangle2(Point2(coordTri[0].u, coordTri[0].v), Point2(coordTri[1].u, coordTri[1].v), Point2(coordTri[2].u, coordTri[2].v)).area();
          double a2 = 0.5;
          double b1 = sqrt(Triangle3(Point3(ptsTri[0].x, ptsTri[0].y, ptsTri[0].z), Point3(ptsTri[1].x, ptsTri[1].y, ptsTri[1].z), Point3(ptsTri[2].x, ptsTri[2].y, ptsTri[2].z)).squared_area());
          double b2 = info.w*info.h*0.5;
          
          double r1 = sqrt(a2/a1);
          double r2 = sqrt(b1/b2);
          /*{
            char buf[1025];
            sprintf(buf, "a1 = %f, a2 = %f, b1 = %f, b2 = %f, r1 = %f, r2 = %f\n", a1, a2, b1, b2, r1, r2);
            OutputDebugStringA(buf);
          }*/
          invMat44 = __GetAffine44UV2Tri_ScaleZ(coordTriScale, ptsTri, r1*r2)  * translateY_Eigen44(info.h) * flipY_Eigen44();
        }
        break;
      case 2:
        invMat44 = __GetAffine44UV2Tri_ScaleZ(coordTriScale, ptsTri, zscale)  * translateY_Eigen44(info.h) * flipY_Eigen44();
        break;
      }
      Eigen::Matrix<double, 3, 4> invMat = invMat44.block(0,0,3,4);
      
      for(int rep=0;rep<numRepeatThru;rep++)
      {
        Point2 &_shiftUV = uvRepeatThru[rep];
        Point2 shiftUV = Point2((_shiftUV.x()), -(_shiftUV.y()));
        std::vector<MyPoint> ptsTexShifted;
        ShiftTri(ptsTex, shiftUV, ptsTexShifted);
        //OutputDebugStringA("thru!\n");

        _OutputPoly2Metaseq(mqoOut, ptsTexShifted, coordTex, invMat, mqTexMatIdx);
      }
      
      
      for(int rep=0;rep<numRepeatClip;rep++)
      {
        Point2 &_shiftUV = uvRepeatClip[rep];
        Point2 shiftUV = Point2((_shiftUV.x()), -(_shiftUV.y()));
        std::vector<MyPoint> triShifted;
        ShiftTri(ptsTex, shiftUV, triShifted);
        std::vector<MyPoint> retPoly;
        std::vector<MyCoordinate> retPolyCoord;
        
        TriangleTriangleIntersect25(ClipTri, triShifted, coordTex, retPoly, retPolyCoord);
        if(retPoly.size()>=3)
        {
          _OutputPoly2Metaseq(mqoOut, retPoly, retPolyCoord, invMat, mqTexMatIdx);
        }
      }
      
    }
  }

  void OutputToMetaseq_LineTri(MyObject o, cv::Mat &matUVCalc, bool bThresholdMennuki = true, transform3 *invMat = NULL, bool bInvFace = false, int mqMatIdx = -1, Vector3 *cliptri = NULL, Point2 *shift2d = NULL, double padding = 0.0)
  {
    bool bOutputAllFace = opt_edgeproc==3 || !bThresholdMennuki;

    bool bPerFace = false;
    if(cliptri!=NULL && invMat!=NULL && shift2d!=NULL)
      bPerFace = true;

    Polygon2_K2 ClipTri;
    Triangle2 ClipTri2;
    if(bPerFace)
    {
      if(cliptri!=NULL)_MakeClipTri(cliptri, ClipTri, ClipTri2, padding, *shift2d);
      if(ClipTri.area()==0.0)return;
    }

#ifdef MYOUTPUTDEBUG
    if(cliptri!=NULL)
    {
      if(ClipTri.is_counterclockwise_oriented())OutputDebugStringA("ClipTri CounterClockWise!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      else OutputDebugStringA("ClipTri cw\n");
    }
#endif
    CDT cdt;
    cdt.insert(otr2_v.begin(), otr2_v.end());
    int lnum = otr2_line.size();
    
    if(bOutputAllFace)
    {
      double wf = chkFace.cols, hf = chkFace.rows;
      cdt.insert(Point2(0.0, 0.0));
      cdt.insert(Point2(wf,  0.0));
      cdt.insert(Point2(0.0, hf));
      cdt.insert(Point2(wf,  hf));
    }
    
#ifdef MYOUTPUTDEBUG
    for(int i=0;i<3;i++)
    {
      char buf[1025];
      sprintf(buf, "ClipTri2: x=%f, y=%f\n", ClipTri2[i].x(), ClipTri2[i].y());
      MyOutputDebugStringA(buf);
    }
    for(int i=0;i<otr2_v.size();i++)
    {
      char buf[1025];
      sprintf(buf, "ov: x=%f, y=%f\n", otr2_v[i].x(), otr2_v[i].y());
      MyOutputDebugStringA(buf);
    }
#endif

    for(int i=0;i<lnum;i++)
    {
      Segment2 *l = &(otr2_line[i]);
      /*
#ifdef MYOUTPUTDEBUG
    {
      char buf[1025];
      sprintf(buf, "ov: x=%f, y=%f, x=%f, y=%f\n", l->source().x(), l->source().y(), l->target().x(), l->target().y());
      MyOutputDebugStringA(buf);
    }
#endif
*/
      
      if(bPerFace)
      {
        //CGAL::cpp11::result_of<K::Intersect_2(Segment_2, Triangle2)>::type result = CGAL::intersection(*l, ClipTri2);
        auto result = CGAL::intersection(*l, ClipTri2);
        if(result)
        {
          if(Segment2 *s2 = boost::get<Segment2>(&*result))
          {
            cdt.insert_constraint(s2->source(), s2->target());
          } else if(Point2 *p2 = boost::get<Point2>(&*result))
          {
            cdt.insert(*p2);
          }
        }
      } else {
        cdt.insert_constraint(l->source(), l->target());
      }
    }
    for (CDT::Finite_faces_iterator fit=cdt.finite_faces_begin(); fit!=cdt.finite_faces_end(); fit++)
    {
      int idx[6];
      MyCoordinate newcoord[6];
      CDT::Face_handle face = fit;
      Point2 &p = face->vertex(0)->point();
      Point2 &p2 = face->vertex(2)->point();
      Point2 &p3 = face->vertex(1)->point();
      bool bFace = true;
      if(bPerFace)
      {
        if(!bOutputAllFace)bFace = isDrawFace(p, p2, p3, chkFace, opt_thresholdMennuki);
        if(!bFace)continue;

        Polygon2_K2 P = _TriangleIntersectionAfter2(p, p2, p3, ClipTri, *shift2d);

        _Output2Metaseq(o, P, *invMat, mqMatIdx, matUVCalc);
      } else {
        if(!bOutputAllFace)bFace = isDrawFace(p, p2, p3, chkFace, opt_thresholdMennuki);
        if(bFace)
        {
          if(invMat==NULL)
          {
            idx[0] = o->AddVertex(MyPoint(CGAL::to_double(p.x()), CGAL::to_double(p.y()), 0.0));
            idx[bInvFace ? 2:1] = o->AddVertex(MyPoint(CGAL::to_double(p2.x()), CGAL::to_double(p2.y()), 0.0));
            idx[bInvFace ? 1:2] = o->AddVertex(MyPoint(CGAL::to_double(p3.x()), CGAL::to_double(p3.y()), 0.0));
          } else {
            Point3 pz = invMat->transform(Point3(CGAL::to_double(p.x()), CGAL::to_double(p.y()), 0.0));
            idx[0] = o->AddVertex(Point3ToMyPoint(pz));
            pz = invMat->transform(Point3(CGAL::to_double(p2.x()), CGAL::to_double(p2.y()), 0.0));
            idx[bInvFace ? 2:1] = o->AddVertex(Point3ToMyPoint(pz));
            pz = invMat->transform(Point3(CGAL::to_double(p3.x()), CGAL::to_double(p3.y()), 0.0));
            idx[bInvFace ? 1:2] = o->AddVertex(Point3ToMyPoint(pz));
          }
          newcoord[0] = CalcUVByAffineMat(matUVCalc, p);
          newcoord[bInvFace ? 2:1] = CalcUVByAffineMat(matUVCalc, p2);
          newcoord[bInvFace ? 1:2] = CalcUVByAffineMat(matUVCalc, p3);

          int fi = o->AddFace(3, idx);
          if(fi!=-1)
          {
            o->SetFaceCoordinateArray(fi, newcoord);
            if(mqMatIdx>=0)o->SetFaceMaterial(fi, mqMatIdx);
          }
        }
      }
    }
  }


  void OutputToMetaseq_LineOnly(MyObject o, transform3 *invMat = NULL, Vector3 *cliptri = NULL, Point2 *shift2d = NULL, double padding = 0.0)
  {
    int idx[2];

    bool bPerFace = false;
    if(cliptri!=NULL && invMat!=NULL && shift2d!=NULL)
      bPerFace = true;

    Polygon2_K2 ClipTri;
    Triangle2 ClipTri2;
    if(bPerFace)
    {
      if(cliptri!=NULL)_MakeClipTri(cliptri, ClipTri, ClipTri2, padding, *shift2d);
      if(ClipTri.area()==0.0)return;
    }
    
    int lnum = otr2_line.size();
    for(int i=0;i<lnum;i++)
    {
      Segment2 *l = &(otr2_line[i]);
      if(bPerFace)
      {
        auto result = CGAL::intersection(*l, ClipTri2);
        if(result)
        {
          if(Segment2 *s2 = boost::get<Segment2>(&*result))
          {
            l = s2;
          }
        }
        const Point2 &s = l->source();
        const Point2 &t = l->target();
        Point3 pz = invMat->transform(Point3(CGAL::to_double(s.x())+shift2d->x(), CGAL::to_double(s.y())+shift2d->y(), 0.0));
        idx[0] = o->AddVertex(Point3ToMyPoint(pz));
        pz = invMat->transform(Point3(CGAL::to_double(t.x())+shift2d->x(), CGAL::to_double(t.y())+shift2d->y(), 0.0));
        idx[1] = o->AddVertex(Point3ToMyPoint(pz));
        o->AddFace(2, idx);
      } else {
        const Point2 &s = l->source();
        const Point2 &t = l->target();
        if(invMat==NULL)
        {
          float x2 = CGAL::to_double(s.x());
          float y2 = CGAL::to_double(s.y());
          idx[0] = o->AddVertex(MyPoint(x2, y2, 0.0));
          x2 = CGAL::to_double(t.x());
          y2 = CGAL::to_double(t.y());
          idx[1] = o->AddVertex(MyPoint(x2, y2, 0.0));
        } else {
          Point3 pz = invMat->transform(Point3(CGAL::to_double(s.x()), CGAL::to_double(s.y()), 0.0));
          idx[0] = o->AddVertex(Point3ToMyPoint(pz));
          pz = invMat->transform(Point3(CGAL::to_double(t.x()), CGAL::to_double(t.y()), 0.0));
          idx[1] = o->AddVertex(Point3ToMyPoint(pz));
        }
        o->AddFace(2, idx);
      }
    }
  }

private:

  void _Output2Metaseq(MyObject o, const Polygon2_K2 &P, transform3 &invMat, int mqMatIdx, cv::Mat &matUVCalc)
  {
    int idx[6];
    MyCoordinate newcoord[6];
    int numIdx=0;
    int numOutVertex = P.size();

    if(numOutVertex>6 || numOutVertex<3)return;

    for(int m=numOutVertex-1;m>=0;m--)
    {
      const Point2_K2 &px = P[m];
      newcoord[numIdx] = CalcUVByAffineMat_K2(matUVCalc, px);
      Point3 vz = invMat.transform(Point3(CGAL::to_double(px.x()), CGAL::to_double(px.y()), 0.0));
      idx[numIdx] = o->AddVertex(Point3ToMyPoint(vz));
      numIdx++;
    }
    if(numIdx>0)
    {
      int fi = o->AddFace(numIdx, idx);
      if(fi!=-1)o->SetFaceCoordinateArray(fi, newcoord);
      if(mqMatIdx!=-1)o->SetFaceMaterial(fi, mqMatIdx);
    }
  }
  Eigen::Matrix<double, 4, 1> MyPoint2EigenMat41(MyPoint p)
  {
    Eigen::Matrix<double, 4, 1> ret;
    ret << p.x,
                 p.y,
                 p.z,
                 1.0;
    return ret;
  }
  
  void OutputDebugStringA_MQP(MyPoint p)
  {
    char buf[1025];
    sprintf(buf, "  (%f, %f, %f)\n", p.x, p.y, p.z);
    MyOutputDebugStringA(buf);
  }
  
  MyPoint transformMQP(MyPoint p, Eigen::Matrix<double, 3, 4> &invMat)
  {
    //OutputDebugStringA("transformMQP:\n");
    //OutputDebugStringA_MQP(p);
    Eigen::Matrix<double, 4, 1> p2 = MyPoint2EigenMat41(p);
    Eigen::Matrix<double, 3, 1> X = invMat * p2;
    MyPoint ret(X(0,0), X(1,0), X(2,0));
    //OutputDebugStringA_MQP(ret);
    return ret;
  }
  void _OutputPoly2Metaseq(MyObject o, std::vector<MyPoint> &poly, std::vector<MyCoordinate> &coord, Eigen::Matrix<double, 3, 4> &invMat, int mqMatIdx)
  {
    int idx[6];
    int numIdx=0;
    int numOutVertex = poly.size();

    if(numOutVertex>6 || numOutVertex<3)return;

    //for(int m=numOutVertex-1;m>=0;m--)
    for(int m=0;m<numOutVertex;m++)
    {
      MyPoint &px = poly[m];
      MyPoint newp = transformMQP(px, invMat);
      idx[numIdx] = o->AddVertex(newp);
      numIdx++;
    }
    if(numIdx>0)
    {
      int fi = o->AddFace(numIdx, idx);
      if(fi!=-1)o->SetFaceCoordinateArray(fi, &(*coord.begin()));
      if(mqMatIdx!=-1)o->SetFaceMaterial(fi, mqMatIdx);
    }
  }
  void _Output2Vector(std::vector<MyPolygon> &polygons, const Polygon2_K2 &P, transform3 &invMat, int mqMatIdx, cv::Mat &matUVCalc)
  {
    int idx[6];
    MyCoordinate newcoord[6];
    int numIdx=0;
    int numOutVertex = P.size();
    MyPolygon polygon;
    polygon.matid = mqMatIdx;

    if(numOutVertex>6 || numOutVertex<3)return;

    for(int m=numOutVertex-1;m>=0;m--)
    {
      const Point2_K2 &px = P[m];
      polygon.coord.push_back(CalcUVByAffineMat_K2(matUVCalc, px));
      Point3 pz = invMat.transform(Point3(CGAL::to_double(px.x()), CGAL::to_double(px.y()), 0.0));
      polygon.vert.push_back(pz);
      numIdx++;
    }
    if(numIdx>0)
    {
      polygons.push_back(polygon);
    }
  }
  
  void Edge(cv::Mat &src, cv::Mat &dst)
  {
    cv::Mat tmp;
    switch(opt_edgeproc)
    {
    case 2:
    {
      //Sobel
      cv::Mat srcf;
      src.convertTo(srcf, CV_32F, 1.0/255);
      cv::Sobel(srcf, tmp, srcf.depth(),1,0,3);
      cvtColor(tmp, dst, CV_RGB2GRAY);
    }
      break;
    case 0:
      cv::Canny(src, tmp, 100, 200);
      tmp.convertTo(dst, CV_32F, 1.0/255);
      break;
    case 1:
    {
      cv::Mat srcf;
      src.convertTo(srcf, CV_32F, 1.0/255);
      cv::Laplacian(srcf, tmp, srcf.depth());
      cvtColor(tmp, dst, CV_RGB2GRAY);
    }
      break;
    default:
    {
      cv::Mat srcf;
      src.convertTo(srcf, CV_32F, 1.0/255);
      cvtColor(srcf, tmp, CV_RGB2GRAY);
      tmp.convertTo(tmp, CV_8U, 255.0);
      tmp = ~tmp;
      tmp.convertTo(dst, CV_32F, 1.0/255.0);
    }
      break;
    }
  }


  void ToPolygonMass(cv::Mat &src, PointMassList& points)
  {
    int w = src.cols, h = src.rows;
    if(src.type() !=5 || src.depth() !=5 || src.channels()!=1)return;

    float *s = (float*)(src.data);
    float *p = s;

    Point2 pt;
    FT mass;
    for(int y = 0;y<h;y++)
    {
      for(int x=0;x<w;x++)
      {
        float v = *p;
        if(v>0.0f)
        {
          pt = Point2(x, y);
          mass = FT(v);
          points.push_back(std::make_pair(pt, mass));
        }
        p++;
      }
    }
  }

  void CheckFaceThreshold(int threshold, cv::Mat &chkFace, cv::Mat &srcGray)
  {
    int w = srcGray.cols, h = srcGray.rows;

    unsigned char *pChk = chkFace.data;
    unsigned char *pGray = srcGray.data;

    for(int y = 0;y<h;y++)
    {
      for(int x=0;x<w;x++)
      {
        *pChk = (*pGray>=threshold) ? 255:0;
        pChk++;
        pGray++;
      }
    }
  }

  void ToPolygonThreshold(cv::Mat &src, std::vector<Point2>& points)
  {
    int w = src.cols, h = src.rows;
    if(src.type() !=5 || src.depth() !=5 || src.channels()!=1)return;

    float *s = (float*)(src.data);
    float *p = s;

    float wf = w, hf = h;

    Point2 pt;
    for(int y = 0;y<h;y++)
    {
      for(int x=0;x<w;x++)
      {
        float v = *p;
        if(v>=opt_threshold)
        {
          pt = Point2(x, y);
          points.push_back(pt);
        }
        p++;
      }
    }
  }




  void _Triangulation(std::vector<Segment2> &otr2_line, std::vector<Point2> &otr2_v, cv::Mat &chkFace, bool bOutputAllFace, int thresholdMennuki, std::vector<std::vector<Point2> > &triArr)
  {
    CDT cdt;
    cdt.insert(otr2_v.begin(), otr2_v.end());
    int lnum = otr2_line.size();

    for(int i=0;i<lnum;i++)
    {
      Segment2 *l = &(otr2_line[i]);
      cdt.insert_constraint(l->source(), l->target());
    }
    triArr.clear();
    triArr.resize(cdt.number_of_faces());
    int i=0;
    for (CDT::Finite_faces_iterator fit=cdt.finite_faces_begin(); fit!=cdt.finite_faces_end(); fit++)
    {
      CDT::Face_handle face = fit;
      Point2 &p = face->vertex(0)->point();
      Point2 &p2 = face->vertex(2)->point();
      Point2 &p3 = face->vertex(1)->point();
      bool bFace = true;

      if(!bOutputAllFace && !isDrawFace(p, p2, p3, chkFace, thresholdMennuki))continue;
      std::vector<Point2> &arr = triArr[i];
      arr.resize(3);
      arr[0] = p;
      arr[1] = p2;
      arr[2] = p3;

      i++;
    }
    triArr.resize(i);
  }

  Polygon2_K2 _TriangleIntersectionAfter2(Point2 &p, Point2 &p2, Point2 &p3, Polygon2_K2 &ClipTri, Point2 &shift2d)
  {
    Point2_K2 ps(p.x()+shift2d.x(), p.y()+shift2d.y());
    Point2_K2 p2s(p2.x()+shift2d.x(), p2.y()+shift2d.y());
    Point2_K2 p3s(p3.x()+shift2d.x(), p3.y()+shift2d.y());

    Polygon2_K2 triPoly;
    triPoly.push_back(ps);
    triPoly.push_back(p2s);
    triPoly.push_back(p3s);
#ifdef DEBUGCGALINTERSECT
      double x = CGAL::to_double(p.x()+shift2d.x()), y = CGAL::to_double(p.y()+shift2d.y());
      fwrite(&x, sizeof(double), 1, fp);
      fwrite(&y, sizeof(double), 1, fp);
      x = CGAL::to_double(p2.x()+shift2d.x()), y = CGAL::to_double(p2.y()+shift2d.y());
      fwrite(&x, sizeof(double), 1, fp);
      fwrite(&y, sizeof(double), 1, fp);
      x = CGAL::to_double(p3.x()+shift2d.x()), y = CGAL::to_double(p3.y()+shift2d.y());
      fwrite(&x, sizeof(double), 1, fp);
      fwrite(&y, sizeof(double), 1, fp);
      fflush(fp);
#endif

    if(triPoly.area()!=0.0)
    {
      if(triPoly.is_clockwise_oriented())
      {
        triPoly.reverse_orientation();
      }
      Pwh_list_2 intR;
      CGAL::intersection(ClipTri, triPoly, std::back_inserter(intR));
  
      for (Pwh_list_2::const_iterator it = intR.begin(); it != intR.end(); it++)
      {
        const Polygon_with_holes_2 &pwh = *it;
        if(!pwh.is_unbounded())
        {
          const Polygon2_K2 &P = pwh.outer_boundary();
          return P;
        }
      }
    }

    Polygon2_K2 empty;
    return empty;
  }

  Polygon2_K2 _TriangleIntersectionAfter2(std::vector<Point2> &tri, Polygon2_K2 &ClipTri, Point2 &shift2d)
  {
    Point2 &p = tri[0];
    Point2 &p2 = tri[2];
    Point2 &p3 = tri[1];
    return _TriangleIntersectionAfter2(p, p2, p3, ClipTri, shift2d);
  }
  
  Triangle2 _MakeClipTri(std::vector<MyCoordinate> &coordTri)
  {
    Triangle2 ret;
    if(coordTri.size()!=3)return ret;
    Polygon2 p;//回転方向判定用
    
    for(int i=0;i<3;i++)
    {
      MyCoordinate &c = coordTri[i];
      p.push_back(Point2(c.u, c.v));
    }
    MyCoordinate &c0 = coordTri[0];
    MyCoordinate &c1 = coordTri[1];
    MyCoordinate &c2 = coordTri[2];
    Point2 p0(c0.u, c0.v);
    Point2 p1(c1.u, c1.v);
    Point2 p2(c2.u, c2.v);
    
    if(p.is_clockwise_oriented())
    {
      ret = Triangle2(p2, p1, p0);
    } else {
      ret = Triangle2(p0, p1, p2);
    }
    return ret;
  }
  Triangle2 _MakeClipTri(std::vector<MyPoint> &tri, bool *bReverse = NULL)
  {
    Triangle2 ret;
    if(tri.size()!=3)return ret;
    Polygon2 poly;//回転方向判定用
    MyPoint &pt0 = tri[0];
    MyPoint &pt1 = tri[1];
    MyPoint &pt2 = tri[2];
    
    for(int i=0;i<3;i++)
    {
      MyPoint &pt = tri[i];
      poly.push_back(Point2(pt.x, pt.y));
    }
    if(poly.area()==0.0)
    {
      return ret;
    }
    Point2 p0(pt0.x, pt0.y);
    Point2 p1(pt1.x, pt1.y);
    Point2 p2(pt2.x, pt2.y);
    
    if(poly.is_clockwise_oriented())
    {
      ret = Triangle2(p2, p1, p0);
      if(bReverse!=NULL)*bReverse = true;
    } else {
      ret = Triangle2(p0, p1, p2);
      if(bReverse!=NULL)*bReverse = false;
    }
    return ret;
  }

  void _MakeClipTri(Vector3 *cliptri, Polygon2_K2 &ClipTri, Triangle2 &ClipTri2, double padding, Point2 &shift2d)
  {
    if(cliptri==NULL)return;
    
    for(int i=0;i<3;i++)
    {
      Vector3 &v = cliptri[i];
      ClipTri.push_back(Point2_K2(v.x(), v.y()));
#ifdef DEBUGCGALINTERSECT
      if(bFirst)
      {
        double x = CGAL::to_double(v.x()), y = CGAL::to_double(v.y());
        fwrite(&x, sizeof(double), 1, fp);
        fwrite(&y, sizeof(double), 1, fp);
        bFirst = false;
#endif
      }
    }
    
    Vector3 &v = cliptri[0];
    Vector3 &v2 = cliptri[1];
    Vector3 &v3 = cliptri[2];
    Point2 p1(v.x()-shift2d.x(), v.y()-shift2d.y());
    Point2 p2(v2.x()-shift2d.x(), v2.y()-shift2d.y());
    Point2 p3(v3.x()-shift2d.x(), v3.y()-shift2d.y());
    
    if(ClipTri.is_clockwise_oriented())
    {
      ClipTri.reverse_orientation();
      ClipTri2 = Triangle2(p3, p2, p1);
    } else {
      ClipTri2 = Triangle2(p1, p2, p3);
    }
  }

/*
  void _TriangleArrayIntersectionAfter(std::vector<std::vector<Point2>> &triArr, MQObject o, Vector3 *cliptri, transform3 &invMat, Point2 &shift2d, Vector3 &shiftFinal, int mqMatIdx, cv::Mat &matUVCalc, std::vector<Polygon2_K2> &retPolygon)
  {
    retPolygon.clear();
    Polygon2_K2 ClipTri;
    if(cliptri!=NULL)_MakeClipTri(cliptri, ClipTri);

    int numtri = triArr.size();
    for (int k=0;k<numtri;k++)
    {
      std::vector<Point2> &tri = triArr[k];
      retPolygon.push_back(_TriangleIntersectionAfter2(tri, o, ClipTri, invMat, shift2d, shiftFinal, mqMatIdx, matUVCalc));
    }
  }
  */

  bool _isDrawFace(cv::Point *points, cv::Mat &chkFace, int thresholdMennuki)
  {
    cv::Mat mask = cv::Mat::zeros(chkFace.rows, chkFace.cols, CV_8U);
    cv::fillConvexPoly(mask, points, 3, cv::Scalar(255));
    unsigned char v = cv::mean(chkFace, mask)[0];
    int nonzero = cv::countNonZero(mask);

    /*
    char b[151];
    sprintf(b, "v = %d\n", v);
    OutputDebugStringA(b);
    */
    /*
    const char* source_window = "chkFace";
    cv::namedWindow( source_window, cv::WINDOW_NORMAL );
    cv::Mat flipChkFace;
    cv::flip(chkFace, flipChkFace, -1);
    imshow( source_window, flipChkFace );
    const char* mask_window = "mask";
    cv::namedWindow( mask_window, cv::WINDOW_NORMAL );
    cv::Mat flipmask;
    cv::flip(mask, flipmask, -1);
    imshow( mask_window, flipmask );
    const char* mask2_window = "mask2";
    cv::namedWindow( mask2_window, cv::WINDOW_NORMAL );
    cv::Mat flipmask2;
    flipChkFace.copyTo(flipmask2, flipmask);
    imshow( mask2_window, flipmask2 );
    cv::waitKey(0);
    */
    return (v<thresholdMennuki && nonzero!=0) ? true:false;
  }
  bool isDrawFaceK2(Point2_K2 &p, Point2_K2 &p2, Point2_K2 &p3, cv::Mat &chkFace, int thresholdMennuki)
  {
    cv::Point points[3];
    points[0] = point2_K2ToCvPoint(p);
    points[1] = point2_K2ToCvPoint(p2);
    points[2] = point2_K2ToCvPoint(p3);
    //return true;
    return _isDrawFace(points, chkFace, thresholdMennuki);
  }
  bool isDrawFace(Point2 &p, Point2 &p2, Point2 &p3, cv::Mat &chkFace, int thresholdMennuki)
  {
    cv::Point points[3];
    points[0] = point2ToCvPoint(p);
    points[1] = point2ToCvPoint(p2);
    points[2] = point2ToCvPoint(p3);
    //return true;
    return _isDrawFace(points, chkFace, thresholdMennuki);
  }
  MyCoordinate CalcUVByAffineMat_K2(cv::Mat &matUVCalc, const Point2_K2 &p)
  {
    MyCoordinate ret;
    cv::Mat pMat(3, 1, CV_64F);
    pMat.at<double>(0,0) = CGAL::to_double(p.x());
    pMat.at<double>(1,0) = CGAL::to_double(p.y());
    pMat.at<double>(2,0) = 1.0;
    cv::Mat result = matUVCalc * pMat;
    return MyCoordinate(result.at<double>(0,0), result.at<double>(1,0));
  }
  MyCoordinate CalcUVByAffineMat(cv::Mat &matUVCalc, const Point2 &p)
  {
    MyCoordinate ret;
    cv::Mat pMat(3, 1, CV_64F);
    pMat.at<double>(0,0) = CGAL::to_double(p.x());
    pMat.at<double>(1,0) = CGAL::to_double(p.y());
    pMat.at<double>(2,0) = 1.0;
    cv::Mat result = matUVCalc * pMat;
    return MyCoordinate(result.at<double>(0,0), result.at<double>(1,0));
  }

};

#endif //_TAMARASTER2VECTOR_
