
#ifndef _TAMAMQTEXMANAGER_
#define _TAMAMQTEXMANAGER_

//#include "MyPolygon.h"
#include "TamaMQLib.h"



std::string GetDrawPolygonPathA();
std::string MyGetTempFilePathA();
DWORD RunCmdA(std::string &cmd);
bool loadPoly(MQObject mqoOut, const char *path);
std::string GetOptStrA(DrawPolygonDialog &dlg);

class MQTexManager
{
public:
  MQTexManager()
  {
  }
  //std::vector<int> index;
  //std::vector<CacheTexInfo> tex;
  //DrawPolygonOptions *dlg;
  
  
  /*
  bool LoadDPMV(const char *path)
  {
    FILE *fp = fopen(path, "rb");
    if(fp==NULL)return false;
    bool bRet = _LoadDPMV(fp);
    fclose(fp);
    return bRet;
  }
  
  bool _LoadDPMV(FILE *fp)
  {
    if(fgetc(fp)!='D')return false;
    if(fgetc(fp)!='P')return false;
    if(fgetc(fp)!='M')return false;
    if(fgetc(fp)!='V')return false;
    if(fgetc(fp)!=0x0)return false;
    
    unsigned int numMaterials = 0;
    if(fread(&numMaterials, 4, 1, fp)!=1)return false;
    
    //for(int i=0;i<numMaterials;i++)index.push_back(-2);//-1);
    index.resize(numMaterials, -2);//-1);
    
    for(int i=0;i<numMaterials;i++)
    {
      unsigned int w = 0;
      if(fread(&w, 4, 1, fp)!=1)return false;
      unsigned int h = 0;
      if(fread(&h, 4, 1, fp)!=1)return false;
      
      unsigned int numMatId = 0;
      if(fread(&numMatId, 4, 1, fp)!=1)return false;
      
      std::vector<int> matIDs;
      for(int n=0;n<numMatId;n++)
      {
        int id = 0;
        if(fread(&id, sizeof(int), 1, fp)!=1)return false;
        matIDs.push_back(id);
      }
      
      tex.push_back(CacheTexInfo());
      int newidx = tex.size()-1;
      CacheTexInfo &ci = tex[newidx];
      bool ret = ci.o.read(fp);
      if(!ret)return false;
      
      ci.w = w;
      ci.h = h;
      __UpdateMinMaxXY(ci);
      for(int n=0;n<numMatId;n++)
      {
        index[matIDs[n]] = newidx;
      }
    }
    return true;
  }
  
  
  bool LoadDPMR(const char *path)
  {
    FILE *fp = fopen(path, "rb");
    if(fp==NULL)return false;
    bool bRet = _LoadDPMR(fp);
    fclose(fp);
    return bRet;
  }
  
  bool _LoadDPMR(FILE *fp)
  {
    if(fgetc(fp)!='D')return false;
    if(fgetc(fp)!='P')return false;
    if(fgetc(fp)!='M')return false;
    if(fgetc(fp)!='R')return false;
    if(fgetc(fp)!=0x0)return false;
    
    unsigned int numMaterials = 0;
    if(fread(&numMaterials, 4, 1, fp)!=1)return false;
    
    //for(int i=0;i<numMaterials;i++)index.push_back(-2);//-1);
    index.resize(numMaterials, -2);//-1);
    
    for(int i=0;i<numMaterials;i++)
    {
      unsigned int numMatId = 0;
      if(fread(&numMatId, 4, 1, fp)!=1)return false;
      
      std::vector<int> matIDs;
      for(int n=0;n<numMatId;n++)
      {
        int id = 0;
        if(fread(&id, sizeof(int), 1, fp)!=1)return false;
        matIDs.push_back(id);
      }
      
      //filename
      unsigned int sizeFilename = 0;
      if(fread(&sizeFilename, 4, 1, fp)!=1)return false;
      
      char *filename = (char *)malloc(sizeFilename+17);
      if(filename==NULL)return false;
      
      if(fread(filename, sizeFilename, 1, fp)!=1)return false;
      memset(filename+sizeFilename, NULL, 16); //NULL
      
      cv::Mat mat = cv::imread(filename);
      if(mat.empty())return -2;
      
      free(filename);
      filename=NULL;
      
      tex.push_back(CacheTexInfo());
      int newidx = tex.size()-1;
      CacheTexInfo &ci = tex[newidx];
      
      ci.matRaster = mat;
      
      ci.w = mat.cols;
      ci.h = mat.rows;
      __UpdateMinMaxXY(ci);
      for(int n=0;n<numMatId;n++)
      {
        index[matIDs[n]] = newidx;
      }
    }
    return true;
  }
  */
  
  static bool MakeCacheNameSlow(const char *texFullPath, char *cachename, int *_w = NULL, int *_h = NULL)
  {
    strcpy(cachename, texFullPath);
    PathStripPathA(cachename);
    PathRemoveExtensionA(cachename);
    int len = strlen(cachename);
    char *pNull = cachename+len;
    cv::Mat mat = cv::imread(texFullPath);
    if(mat.empty())return false;
    int w = mat.cols, h = mat.rows;
    if(w<1 || h<1)return false;
    if(_w!=NULL)*_w = w;
    if(_h!=NULL)*_h = h;
    sprintf(pNull, "_%d_%d", w, h);
    return true;
  }
  static MQObject MakeCacheRoot(MQDocument doc)
  {
    MQObject mqoCacheRoot = MQ_CreateObject();
    if(mqoCacheRoot==NULL)return NULL;
    mqoCacheRoot->SetName("MQDrawPolygonTexCache");
    doc->AddObject(mqoCacheRoot);
    return mqoCacheRoot;
  }
  static bool _Rater2VectorByCmd(std::string &filename, MQObject oOut, DrawPolygonDialog &dlg, int matID = -1)
  {
    std::string dppath = GetDrawPolygonPathA();
    std::string outpath = MyGetTempFilePathA();
    long long matIDll = matID;
    std::string cmd = "\""+dppath+"\""+GetOptStrA(dlg)+" --mode 3 --matid "+std::to_string(matIDll)+" --in \""+filename+"\" --out \""+outpath+"\"";

    DWORD result = RunCmdA(cmd);
    if(result!=0)
    {
      OutputDebugStringA("DrawPolygon.exe failed!");
      remove(outpath.c_str());
      return false;
    }
    
    bool ret = loadPoly(oOut, outpath.c_str());
    
    remove(outpath.c_str());
    
    if(!ret)
    {
      OutputDebugStringA("loadPoly failed!");
      return false;
    }
    return true;
  }
  static MQObject Rater2VectorByCmd(const char *tmpname, std::string &filename, MQObject mqoCacheRoot, MQDocument doc, DrawPolygonDialog &dlg, bool bOptimize, int matID = -1)
  {
    MQObject o = MQ_CreateObject();
    if(o==NULL)return NULL;
    o->SetName(tmpname);
    
    bool bRet = _Rater2VectorByCmd(filename, o, dlg, matID);
    if(!bRet)return NULL;
    
    if(bOptimize)o->OptimizeVertex(0.0f, NULL);
    doc->InsertObject(o, mqoCacheRoot);
    o->SetDepth(1);
    return o;
  }
  static int maxBytes(std::vector<std::string> &filenames)
  {
    int numFile = filenames.size();
    int tmpnameSize = 1;
    for(int i=0;i<numFile;i++)
    {
      std::string &filename = filenames[i];
      int numFilenameSize = filename.size();
      if(numFilenameSize>tmpnameSize)tmpnameSize = numFilenameSize;
    }
    return tmpnameSize;
  }
  static bool MakeAllCache(MQDocument doc, DrawPolygonDialog &dlg, std::vector<std::string> &filenames, std::vector<std::vector<int> > &matidsPerImage, std::vector<MQObject> &mqoTexArr, std::vector<int> &wh, bool bOptimize)
  {
    MQObject mqoCacheRoot = _FindMQObjectByName(doc, "MQDrawPolygonTexCache");
    if(mqoCacheRoot==NULL)
    {
      mqoCacheRoot = MakeCacheRoot(doc);
      if(mqoCacheRoot==NULL)return false;
    }
    
    int tmpnameSize = maxBytes(filenames);
    char *tmpname = (char *)malloc(tmpnameSize+2);
    tmpname[tmpnameSize] = tmpname[tmpnameSize+1] = NULL;
    
    int numFile = filenames.size();
    for(int i=0;i<numFile;i++)
    {
      std::string &filename = filenames[i];
      int w=0, h=0;
      bool bRet = MakeCacheNameSlow(filename.c_str(), tmpname, &w, &h);
      if(!bRet)
      {
        free(tmpname);
        return false;
      }
      wh.push_back(w);
      wh.push_back(h);
      MQObject mqoTex = _FindChildMQObjectByName(doc, tmpname, mqoCacheRoot);
      if(mqoTex==NULL)
      {
        std::vector<int> &matIDs = matidsPerImage[i];
        int matID = -1;
        if(matIDs.size()>0)matID = matIDs[0];
        mqoTex = Rater2VectorByCmd(tmpname, filename, mqoCacheRoot, doc, dlg, bOptimize, matID);
        if(mqoTex==NULL)
        {
          free(tmpname);
          return false;
        }
      }
      mqoTexArr.push_back(mqoTex);
    }
    free(tmpname);
    return true;
  }
  
  static bool SaveDPMV(const char *path, MQDocument doc, DrawPolygonDialog &dlg, std::vector<int> &matIds, bool bOptimize)
  {
    FILE *fp = fopen(path, "wb");
    if(fp==NULL)return false;
    bool bRet = _SaveDPMV(fp, doc, dlg, matIds, bOptimize);
    fclose(fp);
    return bRet;
  }
  static void MakeUniqueFilenameArr(std::vector<int> &matIds, MQDocument doc, std::vector<std::string> &filenames, std::vector<std::vector<int> > &matidsPerImage)
  {
    int numMatIds = matIds.size();
    for(int i=0;i<numMatIds;i++)
    {
      int curmatid = matIds[i];
      MQMaterial mqmat = doc->GetMaterial(curmatid);
      if(mqmat==NULL)continue;
      
      
      char texPath[_MAX_PATH*2 + 2];
      mqmat->GetTextureName(texPath, _MAX_PATH*2);
      char texFullPath[_MAX_PATH*2 + 2];
      if(!doc->FindMappingFile(texFullPath, texPath, MQMAPPING_TEXTURE))continue;
      std::string f(texFullPath);
      int foundIdx = IsContainedString(f, filenames);
      if(foundIdx>=0)
      {
        matidsPerImage[foundIdx].push_back(curmatid);
      } else {
        filenames.push_back(f);
        std::vector<int> tmp;
        tmp.push_back(curmatid);
        matidsPerImage.push_back(tmp);
      }
    }
  }

  static MyPoint MQPointToMyPoint(MQPoint p)
  {
    return MyPoint(p.x, p.y, p.z);
  }
    
  static MyCoordinate MQCoordinateToMyCoordinate(MQCoordinate &_coord)
  {
    return MyCoordinate(_coord.u, _coord.v);
  }

  static void MQCoordinateToMyCoordinate(std::vector<MQCoordinate> &_coords, std::vector<MyCoordinate> &coords)
  {
    int num = _coords.size();
    for(int i=0;i<num;i++)
    {
      coords.push_back(MQCoordinateToMyCoordinate(_coords[i]));
    }
  }

  static void MQObjectToMyObject(MQObject mqo, MyObject oOut)
  {
    int numF = mqo->GetFaceCount();
    for(int fi=0;fi<numF;fi++)
    {
      int numFV = mqo->GetFacePointCount(fi);
      if(numFV<3)continue;
      
      std::vector<int> vp;
      vp.resize(numFV);
      mqo->GetFacePointArray(fi, &(*vp.begin()));
      std::vector<int> vidx;
      for(int mpi=0;mpi<numFV;mpi++)
      {
        MQPoint p = mqo->GetVertex(vp[mpi]);
        vidx.push_back(oOut->AddVertex(MQPointToMyPoint(p)));
      }
      
      int newFaceIdx = oOut->AddFace(numFV, &(*vidx.begin()));
      int matId = mqo->GetFaceMaterial(fi);
      if(matId>=0)oOut->SetFaceMaterial(fi, matId);
      
      std::vector<MQCoordinate> _coords;
      _coords.resize(numFV);
      mqo->GetFaceCoordinateArray(fi, &(*_coords.begin()));
      std::vector<MyCoordinate> coords;
      MQCoordinateToMyCoordinate(_coords, coords);
      oOut->SetFaceCoordinateArray(fi, &(*coords.begin()));
    }
  }
  
  static bool _SaveDPMV(FILE *fp, MQDocument doc, DrawPolygonDialog &dlg, std::vector<int> &matIds, bool bOptimize)
  {
    if(fputc('D', fp)==EOF)return false;
    if(fputc('P', fp)==EOF)return false;
    if(fputc('M', fp)==EOF)return false;
    if(fputc('V', fp)==EOF)return false;
    if(fputc(0x0, fp)==EOF)return false;
    
    std::vector<std::string> filenames;
    std::vector<std::vector<int> > matidsPerImage;
    std::vector<MQObject> mqoTexArr;
    std::vector<int> wh;
    
    MakeUniqueFilenameArr(matIds, doc, filenames, matidsPerImage);
    if(filenames.size()!=matidsPerImage.size())return false;
    
    bool bRet = MakeAllCache(doc, dlg, filenames, matidsPerImage, mqoTexArr, wh, bOptimize);
    if(!bRet)return false;
    if(filenames.size()!=mqoTexArr.size())return false;
    if(filenames.size()*2!=wh.size())return false;
    
    unsigned int numTextures = mqoTexArr.size();
    if(fwrite(&numTextures, 4, 1, fp)!=1)return false;
    
    for(int i=0;i<numTextures;i++)
    {
      unsigned int w = wh[i*2];
      if(fwrite(&w, 4, 1, fp)!=1)return false;
      unsigned int h = wh[i*2+1];
      if(fwrite(&h, 4, 1, fp)!=1)return false;
      
      std::vector<int> &matIDs = matidsPerImage[i];
      
      unsigned int numMatId = matIDs.size();
      if(fwrite(&numMatId, 4, 1, fp)!=1)return false;
      
      for(int n=0;n<numMatId;n++)
      {
        int id = matIDs[n];
        if(fwrite(&id, sizeof(int), 1, fp)!=1)return false;
      }
      
      MQObject mqoTex = mqoTexArr[i];
      TriangulateObj(doc, mqoTex);
      _MyObject _o;
      MQObjectToMyObject(mqoTex, &_o);
      bool ret = _o.write(fp);
      if(!ret)return false;
    }
    return true;
  }
  
  static bool SaveDPMR(const char *path, MQDocument doc, std::vector<int> &matIds)
  {
    FILE *fp = fopen(path, "wb");
    if(fp==NULL)return false;
    bool bRet = _SaveDPMR(fp, doc, matIds);
    fclose(fp);
    return bRet;
  }
  
  static bool SaveDPMR(const wchar_t *path, MQDocument doc, std::vector<int> &matIds)
  {
    FILE *fp = _wfopen(path, L"wb");
    if(fp==NULL)return false;
    bool bRet = _SaveDPMR(fp, doc, matIds);
    fclose(fp);
    return bRet;
  }
  
  static int IsContainedString(std::string &s, std::vector<std::string> &filenames)
  {
    int num = filenames.size();
    for(int i=0;i<num;i++)
    {
      if(s==filenames[i])return i;
    }
    return -1;
  }
  
  static bool _SaveDPMR(FILE *fp, MQDocument doc, std::vector<int> &matIds)
  {
    if(fputc('D', fp)==EOF)return false;
    if(fputc('P', fp)==EOF)return false;
    if(fputc('M', fp)==EOF)return false;
    if(fputc('R', fp)==EOF)return false;
    if(fputc(0x0, fp)==EOF)return false;
    
    std::vector<std::string> filenames;
    std::vector<std::vector<int>> matidsPerImage;
    
    int numMatIds = matIds.size();
    for(int i=0;i<numMatIds;i++)
    {
      int curmatid = matIds[i];
      MQMaterial mqmat = doc->GetMaterial(curmatid);
      if(mqmat==NULL)continue;
      
      
      char texPath[_MAX_PATH*2 + 2];
      mqmat->GetTextureName(texPath, _MAX_PATH*2);
      char texFullPath[_MAX_PATH*2 + 2];
      if(!doc->FindMappingFile(texFullPath, texPath, MQMAPPING_TEXTURE))return -2;
      std::string f(texFullPath);
      int foundIdx = IsContainedString(f, filenames);
      if(foundIdx>=0)
      {
        matidsPerImage[foundIdx].push_back(curmatid);
      } else {
        filenames.push_back(f);
        std::vector<int> tmp;
        tmp.push_back(curmatid);
        matidsPerImage.push_back(tmp);
      }
    }
    
    if(filenames.size()!=matidsPerImage.size())return false;
    
    unsigned int numTextures = matidsPerImage.size();
    if(fwrite(&numTextures, 4, 1, fp)!=1)return false;
    
    for(int i=0;i<numTextures;i++)
    {
      unsigned int numMatId = matidsPerImage[i].size();
      if(numMatId==0)continue;
      if(fwrite(&numMatId, 4, 1, fp)!=1)return false;
      
      std::vector<int> matIDs;
      for(int n=0;n<numMatId;n++)
      {
        int id = matidsPerImage[i][n];
        if(fwrite(&id, sizeof(int), 1, fp)!=1)return false;
      }
      
      //filename
      unsigned int sizeFilename = filenames[i].size();  /////////wstring‚Ìê‡‚Í‚±‚±‘‚«Š·‚¦‚é‚±‚Æ!!!!!!!!!!!!!
      if(fwrite(&sizeFilename, 4, 1, fp)!=1)return false;
      
      if(fwrite(&(filenames[i][0]), sizeFilename, 1, fp)!=1)return false;
    }
    return true;
  }
};

#endif //_TAMAMQTEXMANAGER_
