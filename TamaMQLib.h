
#ifndef _TAMAMQLIB_
#define _TAMAMQLIB_

#include "MQPlugin.h"
#include "MQWidget.h"

void GetPointAndCoord(MQObject o, int fi, int numFV, MQPoint *pts, MQCoordinate *coord = NULL)
{
  //int numFV = o->GetFacePointCount(fi);
  std::vector<int> vidx2(numFV);
  int *vidx = &(*vidx2.begin());
  o->GetFacePointArray(fi, vidx);
  if(coord!=NULL)o->GetFaceCoordinateArray(fi, coord);
  if(pts!=NULL)
  {
    for(int i=0;i<numFV;i++)
    {
      pts[i] = o->GetVertex(vidx[i]);
    }
  }
}
void GetPoint(MQObject o, int fi, int numFV, std::vector<MQPoint> &pts)
{
  GetPointAndCoord(o, fi, numFV, &(*pts.begin()), NULL);
}
void GetPointAndCoord(MQObject o, int fi, int numFV, std::vector<MQPoint> &pts, std::vector<MQCoordinate> &coord)
{
  GetPointAndCoord(o, fi, numFV, &(*pts.begin()), &(*coord.begin()));
}
void GetPointAndCoord(MQObject o, int fi, int numFV, int *vidx, MQPoint *pts, MQCoordinate *coord = NULL)
{
  /*
  std::vector<int> vidx2;
  int *_vidx = vidx;
  if(_vidx==NULL)
  {
    int numFV = o->GetFacePointCount(fi);
    vidx2.resize(numFV);
    _vidx = &(*vidx2.begin())
  }
  */
  o->GetFacePointArray(fi, vidx);
  if(coord!=NULL)o->GetFaceCoordinateArray(fi, coord);
  if(pts!=NULL)
  {
    for(int i=0;i<numFV;i++)
    {
      pts[i] = o->GetVertex(vidx[i]);
    }
  }
}
void GetPoint(MQObject o, int fi, int numFV, std::vector<int> &vidx, std::vector<MQPoint> &pts)
{
  GetPointAndCoord(o, fi, numFV, &(*vidx.begin()), &(*pts.begin()), NULL);
}
void GetPointAndCoord(MQObject o, int fi, int numFV, std::vector<int> &vidx, std::vector<MQPoint> &pts, std::vector<MQCoordinate> &coord)
{
  GetPointAndCoord(o, fi, numFV, &(*vidx.begin()), &(*pts.begin()), &(*coord.begin()));
}
void GetPointAndCoord(MQObject o, int fi, int numFV, int *vidx, std::vector<MQPoint> &pts, std::vector<MQCoordinate> &coord)
{
  GetPointAndCoord(o, fi, numFV, vidx, &(*pts.begin()), &(*coord.begin()));
}

bool _Triangulate(MQDocument doc, int fi, MQObject o, bool bSelectNew = true, bool bDeleteOldFace = true)
{
  if(o==NULL)return false;
  
  int numFV = o->GetFacePointCount(fi);
  if(numFV<=3)return false;
  
  std::vector<int> vidxOld(numFV+1);
  std::vector<MQPoint> pts(numFV);
  std::vector<MQCoordinate> coord(numFV);
  GetPointAndCoord(o, fi, numFV, vidxOld, pts, coord);
  
  MQCoordinate newcoord[3];
  int idx[3];
  
  int matidx = o->GetFaceMaterial(fi);
  int numNewTri = numFV - 2;
  int triidx_size = numNewTri*3;
  std::vector<int> triidx(numNewTri*3);
  doc->Triangulate(&(*pts.begin()), numFV, &(*triidx.begin()), numNewTri*3);
  int oi = doc->GetObjectIndex(o);
  for(int i=0;i<numNewTri;i++)
  {
    int j = i*3;
    for(int k=0;k<3;k++)
    {
      int vii = triidx[j+k];
      idx[k] = vidxOld[vii];
      newcoord[k] = coord[vii];
    }
    int nfi = o->AddFace(3, idx);
    if(bSelectNew)doc->AddSelectFace(oi, nfi);
    o->SetFaceCoordinateArray(nfi, newcoord);
    o->SetFaceMaterial(nfi, matidx);
    for(int k=0;k<3;k++)o->SetFaceVertexColor(nfi, k, o->GetFaceVertexColor(fi, triidx[j+k]));
  }
  if(bDeleteOldFace)o->DeleteFace(fi);
  
  return true;
}
void TriangulateSelected(MQDocument doc, bool bSelectNew = true, bool bDeleteOldFace = true, bool bSkipLock = true, bool bSkipHidden = true)
{
  int numobj = doc->GetObjectCount();
  for(int oi=0;oi<numobj;oi++)
  {
    MQObject o = doc->GetObject(oi);
    if(o==NULL)continue;
    if(bSkipLock && o->GetLocking())continue;
    if(bSkipHidden && o->GetVisible()==0)continue;
    
    int numF = o->GetFaceCount();
    for(int fi=0;fi<numF;fi++)
    {
      if(doc->IsSelectFace(oi, fi)==FALSE)continue;
      
      _Triangulate(doc, fi, o, bSelectNew, bDeleteOldFace);
    }
  }
}
void TriangulateObj(MQDocument doc, MQObject o, bool bSelectNew = false, bool bDeleteOldFace = true, bool bSkipLock = false, bool bSkipHidden = false)
{
  if(o==NULL)return;
  if(bSkipLock && o->GetLocking())return;
  if(bSkipHidden && o->GetVisible()==0)return;
  
  int numF = o->GetFaceCount();
  for(int fi=0;fi<numF;fi++)
  {
    _Triangulate(doc, fi, o, bSelectNew, bDeleteOldFace);
  }
}
void TriangulateObj(MQDocument doc, int oi, bool bSelectNew = false, bool bDeleteOldFace = true, bool bSkipLock = false, bool bSkipHidden = false)
{
  MQObject o = doc->GetObject(oi);
  if(o!=NULL)TriangulateObj(doc, o, bSelectNew, bDeleteOldFace, bSkipLock, bSkipHidden);
}
void Triangulate1Poly(MQDocument doc, int oi, int fi, bool bSelectNew = true, bool bDeleteOldFace = true, bool bSkipLock = true, bool bSkipHidden = true)
{
  MQObject o = doc->GetObject(oi);
  if(o==NULL)return;
  if(bSkipLock && o->GetLocking())return;
  if(bSkipHidden && o->GetVisible()==0)return;
  
  int numF = o->GetFaceCount();
  if(fi>=numF)return;
  
  _Triangulate(doc, fi, o, bSelectNew, bDeleteOldFace);
}



MQObject _FindMQObjectByName(MQDocument doc, char *name, int *_oi = NULL)
{
  char tmp[_MAX_PATH*2+2];
  int numObj = doc->GetObjectCount();
  for(int oi=0;oi<numObj;oi++)
  {
    MQObject o = doc->GetObject(oi);
    if(o==NULL)continue;
    o->GetName(tmp, _MAX_PATH*2);
    if(strcmp(tmp, name)==0)
    {
      if(_oi!=NULL)*_oi = oi;
      return o;
    }
  }
  return NULL;
}

MQObject _FindChildMQObjectByName(MQDocument doc, char *name, MQObject oParent, int *_oi = NULL)
{
  char tmp[_MAX_PATH*2+2];
  int numObj = doc->GetChildObjectCount(oParent);
  for(int oi=0;oi<numObj;oi++)
  {
    MQObject o = doc->GetChildObject(oParent, oi);
    if(o==NULL)continue;
    o->GetName(tmp, _MAX_PATH*2);
    if(strcmp(tmp, name)==0)
    {
      if(_oi!=NULL)*_oi = oi;
      return o;
    }
  }
  return NULL;
}















#endif
