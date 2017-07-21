
#ifndef _TAMAMYPOLYGON_
#define _TAMAMYPOLYGON_



class MyPolygon
{
public:
  MyPolygon() { }
  MyPolygon(const MyPolygon &_p)
  {
    CopyFromMyPolygon(_p);
  }
  MyPolygon& operator=(const MyPolygon &_p)
  {
    CopyFromMyPolygon(_p);
    return *this;
  }
  void CopyFromMyPolygon(const MyPolygon &_p)
  {
    vert = _p.vert;
    coord = _p.coord;
    matid = _p.matid;
  }
  std::vector<Point3> vert;
  std::vector<MQCoordinate> coord;
  int matid;
};

#endif //_TAMAMYPOLYGON_
