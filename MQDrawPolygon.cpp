

#if defined _WIN32 || defined __CYGWIN__
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <float.h>

#include <boost/filesystem.hpp>


#if defined _WIN32 || defined __CYGWIN__
#include <minmax.h>
#endif
#include "MQPlugin.h"
#include "MQWidget.h"

#include "My3DStruct.h"

#if defined _WIN32 || defined __CYGWIN__
//MQSDKのバージョンによってGetObject(atlgdi.hなど)が使えなくなる
//http://www.metaseq.net/bbs/metaseq/bbs.php?lang=jp&res=7044
#if MQPLUGIN_VERSION >= 0x0459
#ifndef GetObject
inline int GetObject(HGDIOBJ p1, int p2, LPVOID p3)
{
#ifdef UNICODE
  return GetObjectW(p1,p2,p3);
#else
  return GetObjectA(p1,p2,p3);
#endif
}
#endif
#endif
#include <atlbase.h>
#include <atlapp.h>
#include <atlmisc.h>
#endif

BOOL DrawPolygon(MQDocument doc);


#include <iostream>
#include <list>


#include "opencv2/highgui.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/imgproc.hpp"

#include "TamaMQLib.h"
#include "CacheTexInfo.h"

//#define MYOUTPUTDEBUG

#ifdef MYOUTPUTDEBUG
void MyOutputDebugStringA(const char *s) { OutputDebugStringA(s); }
#else
void MyOutputDebugStringA(const char *s) { }
#endif



#ifdef MYOUTPUTDEBUG
/*
void OutputDebugStringVector3(Vector3 v)
{
  char buf[1025];
  sprintf(buf, "Vector3: %f, %f, %f\n", float(CGAL::to_double(v.x())), float(CGAL::to_double(v.y())), float(CGAL::to_double(v.z())));
  OutputDebugStringA(buf);
}
void OutputDebugStringVector2(Vector2 v)
{
  char buf[1025];
  sprintf(buf, "Vector2: %f, %f\n", float(CGAL::to_double(v.x())), float(CGAL::to_double(v.y())));
  OutputDebugStringA(buf);
}
*/
void OutputDebugCVSize(cv::Size s)
{
  char buf[1025];
  sprintf(buf, "OutputDebugCVSize = (%d, %d)\n", s.width, s.height);
  OutputDebugStringA(buf);
}
#else
//inline void OutputDebugStringVector3(Vector3 v) {}
//inline void OutputDebugStringVector2(Vector2 v) {}
inline void OutputDebugCVSize(cv::Size s) {}
#endif


#include "DrawPolygonDlg.h"
//#include "Raster2Vector.h"
#include "MQTexManager.h"

//---------------------------------------------------------------------------
//  DllMain
//---------------------------------------------------------------------------
BOOL APIENTRY DllMain( HANDLE hModule, 
                       DWORD  ul_reason_for_call, 
                       LPVOID lpReserved
					 )
{
	//プラグインとしては特に必要な処理はないので、何もせずにTRUEを返す
    return TRUE;
}

//---------------------------------------------------------------------------
//  MQGetPlugInID
//    プラグインIDを返す。
//    この関数は起動時に呼び出される。
//---------------------------------------------------------------------------
MQPLUGIN_EXPORT void MQGetPlugInID(DWORD *Product, DWORD *ID)
{
	// プロダクト名(制作者名)とIDを、全部で64bitの値として返す
	// 値は他と重複しないようなランダムなもので良い
	*Product = 0xA8BEE201;
	*ID      = 0xCD9DA490;
}

//---------------------------------------------------------------------------
//  MQGetPlugInName
//    プラグイン名を返す。
//    この関数は[プラグインについて]表示時に呼び出される。
//---------------------------------------------------------------------------
MQPLUGIN_EXPORT const char *MQGetPlugInName(void)
{
	// プラグイン名
	return "MQDrawPolygon           Copyright(C) 2017, tamachan";
}

//---------------------------------------------------------------------------
//  MQGetPlugInType
//    プラグインのタイプを返す。
//    この関数は起動時に呼び出される。
//---------------------------------------------------------------------------
MQPLUGIN_EXPORT int MQGetPlugInType(void)
{
	// 選択部変形用プラグインである
	return MQPLUGIN_TYPE_SELECT;
}

//---------------------------------------------------------------------------
//  MQEnumString
//    ポップアップメニューに表示される文字列を返す。
//    この関数は起動時に呼び出される。
//---------------------------------------------------------------------------
MQPLUGIN_EXPORT const char *MQEnumString(int index)
{
	switch(index){
	case 0: return "DrawPolygon";
	}
	return NULL;
}

//---------------------------------------------------------------------------
//  MQModifySelect
//    メニューから選択されたときに呼び出される。
//---------------------------------------------------------------------------
MQPLUGIN_EXPORT BOOL MQModifySelect(int index, MQDocument doc)
{
  switch(index){
  case 0: return DrawPolygon(doc);
  }
  return FALSE;
}



bool GetDllDirA(char *path, int size)
{
  //char path[_MAX_PATH+16];
  HMODULE hModule = NULL;
  if(!GetModuleHandleEx(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,  (LPTSTR)&GetDllDirA, &hModule))return false;
  if(!GetModuleFileNameA(hModule, path, size))return false;
  PathRemoveFileSpecA(path);
  return true;
}
bool GetDllDirW(wchar_t *path, int size)
{
  //char path[_MAX_PATH+16];
  HMODULE hModule = NULL;
  if(!GetModuleHandleEx(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,  (LPTSTR)&GetDllDirW, &hModule))return false;
  if(!GetModuleFileNameW(hModule, path, size))return false;
  PathRemoveFileSpecW(path);
  return true;
}

std::string GetDrawPolygonPathA()
{
  char path[_MAX_PATH+16];
  path[0] = NULL;
  bool bRet = GetDllDirA(path, _MAX_PATH);
  if(!bRet)return "DrawPolygon.exe";
  std::string ret = path;
  return ret+"\\DrawPolygon.exe";
}

std::wstring GetDrawPolygonPathW()
{
  wchar_t path[_MAX_PATH+16];
  path[0] = NULL;
  bool bRet = GetDllDirW(path, _MAX_PATH);
  if(!bRet)return L"DrawPolygon.exe";
  std::wstring ret = path;
  return ret+L"\\DrawPolygon.exe";
}

std::string MyGetTempFilePathA()
{
  boost::filesystem::path path = boost::filesystem::temp_directory_path() / boost::filesystem::unique_path();
  return path.string();
}

std::wstring MyGetTempFilePathW()
{
  boost::filesystem::path path = boost::filesystem::temp_directory_path() / boost::filesystem::unique_path();
  return path.native();
}

MQPoint MyPointToMQPoint(MyPoint p)
{
  return MQPoint(p.x, p.y, p.z);
}

MQCoordinate MyCoordinateToMQCoordinate(MyCoordinate &_coord)
{
  return MQCoordinate(_coord.u, _coord.v);
}

void MyCoordinateToMQCoordinate(std::vector<MyCoordinate> &_coords, std::vector<MQCoordinate> &coords)
{
  int num = _coords.size();
  for(int i=0;i<num;i++)
  {
    coords.push_back(MyCoordinateToMQCoordinate(_coords[i]));
  }
}

void MyObjectToMQObject(MyObject o, MQObject mqoOut)
{
  int numF = o->GetFaceCount();
  for(int fi=0;fi<numF;fi++)
  {
    int numFV = o->GetFacePointCount(fi);
    if(numFV<3)continue;
    
    std::vector<MyPoint> vp;
    o->GetFaceMyPointArray(fi, vp);
    std::vector<int> vidx;
    for(int mpi=0;mpi<numFV;mpi++)
    {
      vidx.push_back(mqoOut->AddVertex(MyPointToMQPoint(vp[mpi])));
    }
    
    int newFaceIdx = mqoOut->AddFace(numFV, &(*vidx.begin()));
    int matId = o->GetFaceMaterial(fi);
    if(matId>=0)mqoOut->SetFaceMaterial(fi, matId);
    
    std::vector<MyCoordinate> _coords;
    o->GetFaceCoordinateArray(fi, _coords);
    std::vector<MQCoordinate> coords;
    MyCoordinateToMQCoordinate(_coords, coords);
    mqoOut->SetFaceCoordinateArray(fi, &(*coords.begin()));
  }
}

bool loadPoly(MQObject mqoOut, const char *path)
{
  _MyObject o;
  bool bRet = o.read(path);
  if(!bRet)return false;
  MyObjectToMQObject(&o, mqoOut);
  return true;
}
bool loadPoly(MQObject mqoOut, const wchar_t *path)
{
  _MyObject o;
  bool bRet = o.read(path);
  if(!bRet)return false;
  MyObjectToMQObject(&o, mqoOut);
  return true;
}

DWORD RunCmdA(std::string &cmd)
{
  STARTUPINFOA si = {0};
  si.cb = sizeof(si);
  PROCESS_INFORMATION pi = {0};
  
  if(!CreateProcessA(NULL, &cmd[0], NULL, NULL, FALSE, CREATE_NO_WINDOW, NULL, NULL, &si, &pi))
  {
    return -1;
  }

  DWORD result=-1;
  WaitForSingleObject(pi.hProcess, INFINITE);
  GetExitCodeProcess(pi.hProcess, &result);

  CloseHandle(pi.hProcess);
  CloseHandle(pi.hThread);
  return result;
}

DWORD RunCmdW(std::wstring &cmd)
{
  STARTUPINFOW si = {0};
  si.cb = sizeof(si);
  PROCESS_INFORMATION pi = {0};
  
  if(!CreateProcessW(NULL, &cmd[0], NULL, NULL, FALSE, CREATE_NO_WINDOW, NULL, NULL, &si, &pi))
  {
    return false;
  }

  DWORD result=-1;
  WaitForSingleObject(pi.hProcess, INFINITE);
  GetExitCodeProcess(pi.hProcess, &result);

  CloseHandle(pi.hProcess);
  CloseHandle(pi.hThread);
  return result;
}

std::string GetOptStrA(DrawPolygonDialog &dlg)
{
  int modeZScale = dlg.combo_zscale->GetCurrentIndex();
  double zscale = dlg.dblspin_zscale->GetPosition();
  int edgeproc = dlg.combo_edgeproc->GetCurrentIndex();
  int vectorconv = dlg.combo_vectorconv->GetCurrentIndex();
  double threshold = dlg.slider_threshold->GetPosition();
  int np = dlg.spin_np->GetPosition();
  bool bFacegen = dlg.check_facegen->GetChecked();
  bool bThresholdMennuki = dlg.check_thresholdMennuki->GetChecked();
  double thresholdMennuki = dlg.slider_thresholdMennuki->GetPosition();
  bool bOptimize = dlg.check_optimize->GetChecked();
  return " -a "+std::to_string((long long)modeZScale)+" -z "+std::to_string((long double)zscale)+" -e "+std::to_string((long long)edgeproc)+" -v "+std::to_string((long long)vectorconv)+" -b "+std::to_string((long double)threshold)+" -n "+std::to_string((long long)np)+(bFacegen?" -g":"")+(bThresholdMennuki?" -h":"")+" -c "+std::to_string((long double)thresholdMennuki)+(bOptimize?" -p":"")+" ";
}

std::wstring GetOptStrW(DrawPolygonDialog &dlg)
{
  int modeZScale = dlg.combo_zscale->GetCurrentIndex();
  double zscale = dlg.dblspin_zscale->GetPosition();
  int edgeproc = dlg.combo_edgeproc->GetCurrentIndex();
  int vectorconv = dlg.combo_vectorconv->GetCurrentIndex();
  double threshold = dlg.slider_threshold->GetPosition();
  int np = dlg.spin_np->GetPosition();
  bool bFacegen = dlg.check_facegen->GetChecked();
  bool bThresholdMennuki = dlg.check_thresholdMennuki->GetChecked();
  double thresholdMennuki = dlg.slider_thresholdMennuki->GetPosition();
  bool bOptimize = dlg.check_optimize->GetChecked();
  return L" -a "+std::to_wstring((long long)modeZScale)+L" -z "+std::to_wstring((long double)zscale)+L" -e "+std::to_wstring((long long)edgeproc)+L" -v "+std::to_wstring((long long)vectorconv)+L" -b "+std::to_wstring((long double)threshold)+L" -n "+std::to_wstring((long long)np)+(bFacegen?L" -g":L"")+(bThresholdMennuki?L" -h":L"")+L" -c "+std::to_wstring((long double)thresholdMennuki)+(bOptimize?L" -p":L"")+L" ";
}


bool ConvertFromClipboard(MQDocument doc, DrawPolygonDialog &dlg)
{
  int edgeproc = dlg.combo_edgeproc->GetCurrentIndex();
  //int vectorconv = dlg.combo_vectorconv->GetCurrentIndex();
  //double threshold = dlg.slider_threshold->GetPosition();
  //int np = dlg.spin_np->GetPosition();
  bool bFacegen = dlg.check_facegen->GetChecked();
  //int thresholdMennuki = dlg.slider_thresholdMennuki->GetPosition() * 255.0;
  bool bOptimize = dlg.check_optimize->GetChecked();
  bool bThresholdMennuki = dlg.check_thresholdMennuki->GetChecked();
  
  MQObject o = MQ_CreateObject();
  char objname[151];
  doc->GetUnusedObjectName(objname, 150, "draw");
  o->SetName(objname);
  
  std::wstring dppath = GetDrawPolygonPathW();
  std::wstring outpath = MyGetTempFilePathW();

  
  std::wstring cmd = L"\""+dppath+L"\" "+GetOptStrW(dlg)+L" --mode 0 --out \""+outpath+L"\"";

  DWORD result = RunCmdW(cmd);
  
/*  ::boost::process::environment env = ::boost::this_process::environment();
  ::boost::process::child c(dppath+" --mode 0 --out "+outpath, env, ::boost::process::windows::hide);
  c.wait();
  int result = c.exit_code();
  */
  if(result!=0)
  {
    OutputDebugStringA("DrawPolygon.exe failed!");
    _wremove(outpath.c_str());
    return false;
  }
  
  bool ret = loadPoly(o, outpath.c_str());
  
  _wremove(outpath.c_str());
  
  if(!ret)
  {
    OutputDebugStringA("loadPoly failed!");
    return false;
  }

  if(bOptimize)o->OptimizeVertex(0.0f, NULL);
  doc->AddObject(o);
  
  //const char* source_window = "Source";
  //cv::namedWindow( source_window, cv::WINDOW_NORMAL );
  //imshow( source_window, dstEdge );
  //cv::waitKey(0);
  
  return true;
}






MyPoint MQPointToMyPoint(MQPoint p)
{
  return MyPoint(p.x, p.y, p.z);
}
MyCoordinate MQCoordinateToMyCoordinate(MQCoordinate uv)
{
  return MyCoordinate(uv.u, uv.v);
}
void MQCoordinateToMyCoordinate(int num, MQCoordinate *uv, MyCoordinate *uv2)
{
  for(int i=0;i<num;i++)
  {
    uv2[i] = MQCoordinateToMyCoordinate(uv[i]);
  }
}

bool _OutputSelectedFaces(_MyObject &_oDst, std::vector<bool> &matUsed, MQDocument doc, std::vector<int> &matIds, MQObject ignoreParentObj = NULL)
{
  MyObject oDst = &_oDst;
  int newidx[3];
  MyCoordinate uv2[3];
  
  int vidx[3];
  MQCoordinate triuv[3];
  MQPoint pts[3];
  
  int numMaterials = doc->GetMaterialCount();
  matUsed.resize(numMaterials, false);
  
  int numobj = doc->GetObjectCount();
  for(int oi=0;oi<numobj;oi++)
  {
    MQObject o = doc->GetObject(oi);
    if(o==NULL || o->GetLocking() || o->GetVisible()==0)continue;
    
    if(ignoreParentObj!=NULL)
    {
      MQObject oParent = doc->GetParentObject(o);
      if(oParent==ignoreParentObj)continue;
    }
    
    int numV = o->GetVertexCount();
    int numF = o->GetFaceCount();
    for(int fi=0;fi<numF;fi++)
    {
      if(doc->IsSelectFace(oi, fi)==FALSE)continue;
      int numFV = o->GetFacePointCount(fi);
      if(numFV!=3)continue;
      
      GetPointAndCoord(o, fi, 3, vidx, pts, triuv);
      for(int i=0;i<3;i++)
      {
        newidx[i] = oDst->AddVertex(MQPointToMyPoint(pts[i]));
      }
      int newfi = oDst->AddFace(3, newidx);
      int mqMatIdx = o->GetFaceMaterial(fi);
      matUsed[mqMatIdx] = true;
      if(mqMatIdx>=0)oDst->SetFaceMaterial(newfi, mqMatIdx);
      MQCoordinateToMyCoordinate(3, triuv, uv2);
      oDst->SetFaceCoordinateArray(newfi, uv2);
    }
  }
  
  return true;
}
bool OutputSelectedFacesW(MQDocument doc, const wchar_t *outPath, std::vector<int> &matIds, MQObject ignoreParentObj = NULL)
{
  _MyObject _oDst;
  MyObject oDst = &_oDst;
  std::vector<bool> matUsed;
  
  bool bRet = _OutputSelectedFaces(_oDst, matUsed, doc, matIds, ignoreParentObj);
  if(!bRet)return false;
  
  bRet = oDst->write(outPath);
  if(!bRet)return false;
  
  int numMaterials = doc->GetMaterialCount();
  for(int i=0;i<numMaterials;i++)
  {
    if(matUsed[i])matIds.push_back(i);
  }
  return bRet;
}
bool OutputSelectedFacesA(MQDocument doc, const char *outPath, std::vector<int> &matIds, MQObject ignoreParentObj = NULL)
{
  _MyObject _oDst;
  MyObject oDst = &_oDst;
  std::vector<bool> matUsed;
  
  bool bRet = _OutputSelectedFaces(_oDst, matUsed, doc, matIds, ignoreParentObj);
  if(!bRet)return false;
  
  bRet = oDst->write(outPath);
  if(!bRet)return false;
  
  int numMaterials = doc->GetMaterialCount();
  for(int i=0;i<numMaterials;i++)
  {
    if(matUsed[i])matIds.push_back(i);
  }
  return bRet;
}


bool ConvertFromTextureSlow(MQDocument doc, DrawPolygonDialog &dlg)
{
  int edgeproc = dlg.combo_edgeproc->GetCurrentIndex();
  //int vectorconv = dlg.combo_vectorconv->GetCurrentIndex();
  //double threshold = dlg.slider_threshold->GetPosition();
  //int np = dlg.spin_np->GetPosition();
  bool bFacegen = dlg.check_facegen->GetChecked();
  //int thresholdMennuki = dlg.slider_thresholdMennuki->GetPosition() * 255.0;
  bool bOptimize = dlg.check_optimize->GetChecked();
  bool bThresholdMennuki = dlg.check_thresholdMennuki->GetChecked();
  
  TriangulateSelected(doc);
  
  //int numMaterial = doc->GetMaterialCount();
  //std::vector<int> texLoadStatus(numMaterial, -1);
  //std::vector<std::pair<std::string, cv::Mat>> texLoaded;
    
  MQObject oOut = MQ_CreateObject();
  char objname[151];
  doc->GetUnusedObjectName(objname, 150, "draw");
  oOut->SetName(objname);
  
  
  MQObject mqoCacheRoot = _FindMQObjectByName(doc, "MQDrawPolygonTexCache");
  
  std::wstring dppath = GetDrawPolygonPathW();
  
  std::wstring tripath = MyGetTempFilePathW();
  std::vector<int> matIds;
  bool bRet = OutputSelectedFacesW(doc, tripath.c_str(), matIds, mqoCacheRoot);
  if(!bRet)return false;
  
  std::wstring dpmrpath = MyGetTempFilePathW();
  bRet = MQTexManager::SaveDPMR(dpmrpath.c_str(), doc, matIds);
  if(!bRet)
  {
    _wremove(tripath.c_str());
    return false;
  }
  
  
  std::wstring outpath = MyGetTempFilePathW();
  
  
  std::wstring cmd = L"\""+dppath+L"\""+GetOptStrW(dlg)+L" --mode 2 --in \""+tripath+L"\" --texCache \""+dpmrpath+L"\" --out \""+outpath+L"\"";

  DWORD result = RunCmdW(cmd);
  if(result!=0)
  {
    OutputDebugStringA("DrawPolygon.exe failed!");
    MQDialog::MessageWarningBox(MQWindow::GetMainWindow(), L"低速変換に失敗", L"Error");
    _wremove(tripath.c_str());
    _wremove(dpmrpath.c_str());
    _wremove(outpath.c_str());
    return false;
  }
  
  bool ret = loadPoly(oOut, outpath.c_str());
  
  _wremove(tripath.c_str());
  _wremove(dpmrpath.c_str());
  _wremove(outpath.c_str());
  
  if(!ret)
  {
    OutputDebugStringA("loadPoly failed!");
    return false;
  }
  
  if(bOptimize)
  {
    MyOutputDebugStringA("OptimizeVertex\n");
    oOut->OptimizeVertex(0.0f, NULL);
  }
  MyOutputDebugStringA("AddObject\n");
  doc->AddObject(oOut);
  
  //const char* source_window = "Source";
  //cv::namedWindow( source_window, cv::WINDOW_NORMAL );
  //imshow( source_window, dstEdge );
  //cv::waitKey(0);
  MyOutputDebugStringA("End\n");
  
  return true;
}

bool ConvertFromTextureFast(MQDocument doc, DrawPolygonDialog &dlg)
{
  int edgeproc = dlg.combo_edgeproc->GetCurrentIndex();
  //int vectorconv = dlg.combo_vectorconv->GetCurrentIndex();
  //double threshold = dlg.slider_threshold->GetPosition();
  //int np = dlg.spin_np->GetPosition();
  bool bFacegen = dlg.check_facegen->GetChecked();
  //int thresholdMennuki = dlg.slider_thresholdMennuki->GetPosition() * 255.0;
  bool bOptimize = dlg.check_optimize->GetChecked();
  int modeZScale = dlg.combo_zscale->GetCurrentIndex();
  double zscale = dlg.dblspin_zscale->GetPosition();
  
  TriangulateSelected(doc);
  
  _MyObject oDst;
    
  MQObject oOut = MQ_CreateObject();
  char objname[151];
  doc->GetUnusedObjectName(objname, 150, "draw");
  oOut->SetName(objname);
  
  MQObject mqoCacheRoot = _FindMQObjectByName(doc, "MQDrawPolygonTexCache");
  
  std::string dppath = GetDrawPolygonPathA();
  
  std::string tripath = MyGetTempFilePathA();
  std::vector<int> matIds;
  bool bRet = OutputSelectedFacesA(doc, tripath.c_str(), matIds, mqoCacheRoot);
  if(!bRet)return false;
  
  std::string dpmvpath = MyGetTempFilePathA();
  bRet = MQTexManager::SaveDPMV(dpmvpath.c_str(), doc, dlg, matIds, bOptimize);
  if(!bRet)
  {
    remove(tripath.c_str());
    return false;
  }
  
  
  std::string outpath = MyGetTempFilePathA();
  
  
  std::string cmd = "\""+dppath+"\""+GetOptStrA(dlg)+" --mode 1 --in \""+tripath+"\" --texCache \""+dpmvpath+"\" --out \""+outpath+"\"";

  DWORD result = RunCmdA(cmd);
  if(result!=0)
  {
    OutputDebugStringA("DrawPolygon.exe failed!");
    MQDialog::MessageWarningBox(MQWindow::GetMainWindow(), L"ポリゴン高速変換に失敗", L"Error");
    remove(tripath.c_str());
    remove(dpmvpath.c_str());
    remove(outpath.c_str());
    return false;
  }
  
  bool ret = loadPoly(oOut, outpath.c_str());
  
  remove(tripath.c_str());
  remove(dpmvpath.c_str());
  remove(outpath.c_str());
  
  if(!ret)
  {
    OutputDebugStringA("loadPoly failed!");
    return false;
  }
  
  if(bOptimize)
  {
    MyOutputDebugStringA("OptimizeVertex\n");
    oOut->OptimizeVertex(0.0f, NULL);
  }
  MyOutputDebugStringA("AddObject\n");
  doc->AddObject(oOut);
  
  //const char* source_window = "Source";
  //cv::namedWindow( source_window, cv::WINDOW_NORMAL );
  //imshow( source_window, dstEdge );
  //cv::waitKey(0);
  MyOutputDebugStringA("End\n");
  
  return true;
}

BOOL DrawPolygon(MQDocument doc)
{
  MQWindow mainwin = MQWindow::GetMainWindow();
  DrawPolygonDialog dlg(mainwin);
  dlg.UpdateEnable(doc);
  if(dlg.Execute() != MQDialog::DIALOG_OK){
    return FALSE;
  }


  int srcType = dlg.combo_src->GetCurrentIndex();
  
  bool bRet = false;
  switch(srcType)
  {
  case 0:
  default:
    bRet = ConvertFromClipboard(doc, dlg);
    break;
  case 1:
    bRet = ConvertFromTextureFast(doc, dlg);
    break;
  case 2:
    bRet = ConvertFromTextureSlow(doc, dlg);
    break;
  }
  
  MQ_RefreshView(NULL);
  
  return bRet ? TRUE : FALSE;
}

