

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include "MQPlugin.h"
#include "MQWidget.h"

BOOL DrawPolygon(MQDocument doc);


#include <iostream>
#include <list>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Optimal_transportation_reconstruction_2.h>


//#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"

typedef CGAL::Simple_cartesian<double> K;

typedef K::FT FT;
typedef K::Ray_3 Ray;
typedef K::Line_3 Line;
typedef K::Point_3 Point3;
typedef K::Point_2 Point2;
typedef K::Triangle_3 Triangle;
typedef K::Direction_3 Direction;
typedef K::Vector_3 Vector3;
typedef K::Segment_2 Segment2;

typedef std::list<Triangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;

typedef std::pair<Point2, FT> PointMassPair;
typedef std::vector<PointMassPair> PointMassList;
typedef CGAL::First_of_pair_property_map <PointMassPair> Point_property_map;
typedef CGAL::Second_of_pair_property_map <PointMassPair> Mass_property_map;
typedef CGAL::Optimal_transportation_reconstruction_2<K, Point_property_map, Mass_property_map> Otr_2;
typedef CGAL::Optimal_transportation_reconstruction_2<K> Otr;


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

class DrawPolygonDialog : public MQDialog
{
public:
  DrawPolygonDialog(MQWindowBase& parent);
  ~DrawPolygonDialog();

  BOOL ComboChanged(MQWidgetBase *sender, MQDocument doc)
  {
    slider_threshold->SetEnabled( combo_vectorconv->GetCurrentIndex()==0 ? false:true);
    return FALSE;
  }

  MQComboBox *combo_edgeproc;
  MQComboBox *combo_vectorconv;
  MQSlider *slider_threshold;
  MQSpinBox *spin_np;
};

DrawPolygonDialog::DrawPolygonDialog(MQWindowBase& parent) : MQDialog(parent)
{
  SetTitle(L"DrawPolygon");

  MQFrame *mainFrame = CreateHorizontalFrame(this);

  MQFrame *paramFrame = CreateHorizontalFrame(mainFrame);
  paramFrame->SetMatrixColumn(2);
  
  CreateLabel(paramFrame, L"ライン抽出");
  combo_edgeproc = CreateComboBox(paramFrame);
  combo_edgeproc->AddItem(L"Canny (おすすめ)");
  combo_edgeproc->AddItem(L"Laplacian");
  combo_edgeproc->AddItem(L"Sobel");
  combo_edgeproc->AddItem(L"無し（中央線）");
  combo_edgeproc->SetCurrentIndex(0);
  
  CreateLabel(paramFrame, L"ベクトル化");
  combo_vectorconv = CreateComboBox(paramFrame);
  combo_vectorconv->AddItem(L"モノクロ値を重みとして使用");
  combo_vectorconv->AddItem(L"2値化（重み無し））");
  combo_vectorconv->SetCurrentIndex(0);
  combo_vectorconv->AddChangedEvent(this, &DrawPolygonDialog::ComboChanged);
  
  CreateLabel(paramFrame, L"2値化の判定値");
  slider_threshold = CreateSlider(paramFrame);
  slider_threshold->SetMin(0.01);
  slider_threshold->SetMax(0.99);
  slider_threshold->SetPosition(0.5);
  slider_threshold->SetEnabled(false);
  
  CreateLabel(paramFrame, L"出力頂点数");
  spin_np = CreateSpinBox(paramFrame);
  spin_np->SetMin(2);
  spin_np->SetMax(1000);
  spin_np->SetPosition(100);

  MQFrame *sideFrame = CreateVerticalFrame(mainFrame);

  MQButton *okbtn = CreateButton(sideFrame, L"OK");
  okbtn->SetDefault(true);
  okbtn->SetModalResult(MQDialog::DIALOG_OK);

  MQButton *cancelbtn = CreateButton(sideFrame, L"Cancel");
  cancelbtn->SetCancel(true);
  cancelbtn->SetModalResult(MQDialog::DIALOG_CANCEL);
}

DrawPolygonDialog::~DrawPolygonDialog()
{
}

#include <minmax.h>
#include <atlbase.h>
#include <atlapp.h>
#include <atlmisc.h>


BOOL GetClipboardBitmap(cv::Mat &mat)
{
  if(!::OpenClipboard(NULL))return FALSE;
  HBITMAP hBitmap = (HBITMAP)::GetClipboardData(CF_BITMAP);
  BITMAP bmp;
  if(hBitmap==NULL || !GetObject(hBitmap, sizeof(BITMAP), &bmp))
  {
    ::CloseClipboard();
    return FALSE;
  }
  BITMAPINFO bi = {0};
  BITMAPINFOHEADER *bh = &(bi.bmiHeader);
  bh->biSize = sizeof(BITMAPINFOHEADER);
  bh->biWidth = bmp.bmWidth;
  bh->biHeight = -bmp.bmHeight;
  bh->biPlanes = 1;
  bh->biBitCount = 32;
  bh->biCompression = BI_RGB;
  
  CDC dc;
  dc.CreateCompatibleDC();
  HBITMAP hBitmapOld = dc.SelectBitmap(hBitmap);
  mat.create(bmp.bmHeight, bmp.bmWidth, CV_8UC4);
  GetDIBits(dc, hBitmap, 0, bmp.bmHeight, mat.data, &bi, DIB_RGB_COLORS);
  
  dc.SelectBitmap(hBitmapOld);
  ::CloseClipboard();
  
  return TRUE;
}

void Edge(int edgeproc, cv::Mat &src, cv::Mat &srcf, cv::Mat &dst)
{
  cv::Mat tmp;
  switch(edgeproc)
  {
  case 2:
    //Sobel
    cv::Sobel(srcf, tmp, srcf.depth(),1,0,3);
    cvtColor(tmp, dst, CV_RGB2GRAY);
    break;
  case 0:
    cv::Canny(src, tmp, 100, 200);
    tmp.convertTo(dst, CV_32F, 1.0/255);
    break;
  case 1:
    cv::Laplacian(srcf, tmp, srcf.depth());
    cvtColor(tmp, dst, CV_RGB2GRAY);
    break;
  default:
    cvtColor(srcf, tmp, CV_RGB2GRAY);
    tmp.convertTo(tmp, CV_8U, 255.0);
    tmp = ~tmp;
    tmp.convertTo(dst, CV_32F, 1.0/255.0);
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
void ToPolygonThreshold(cv::Mat &src, std::vector<Point2>& points, float threshold)
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
      if(v>=threshold)
      {
        //float x2 = x;
        //x2 = (x2 / wf) * 100.0f;
        //float y2 = y;
        //y2 = (1.0f - (y2 / hf)) * 100.0f;
        pt = Point2(x, y);
        points.push_back(pt);
      }
      p++;
    }
  }
}

BOOL DrawPolygon(MQDocument doc)
{
  MQWindow mainwin = MQWindow::GetMainWindow();
  DrawPolygonDialog dlg(mainwin);
  if(dlg.Execute() != MQDialog::DIALOG_OK){
    return FALSE;
  }

  int edgeproc = dlg.combo_edgeproc->GetCurrentIndex();
  int vectorconv = dlg.combo_vectorconv->GetCurrentIndex();
  double threshold = dlg.slider_threshold->GetPosition();
  int np = dlg.spin_np->GetPosition();
  
  
  cv::Mat mat;
  if(GetClipboardBitmap(mat)==FALSE)return FALSE;

  cv::Mat matf;
  mat.convertTo(matf, CV_32F, 1.0/255);
  
  cv::Mat dstEdge;
  Edge(edgeproc, mat, matf, dstEdge);

  std::vector<Point2> otr2_v;
  std::vector<Segment2> otr2_line;
  
  Point_property_map point_pmap;
  Mass_property_map  mass_pmap;
  PointMassList pointsMass;
  std::vector<Point2> points;

  switch(vectorconv)
  {
  case 0:
    {
    ToPolygonMass(dstEdge, pointsMass);
    if(pointsMass.size()==0)return FALSE;
    Otr_2 otr2(pointsMass, point_pmap, mass_pmap);
    otr2.run_until(MIN(np, pointsMass.size()));
    otr2.list_output(std::back_inserter(otr2_v), std::back_inserter(otr2_line));
    }
    break;
  case 1:
    {
    ToPolygonThreshold(dstEdge, points, threshold);
    if(points.size()==0)return FALSE;
    Otr otr(points);
    otr.run_until(MIN(np, points.size()));
    otr.list_output(std::back_inserter(otr2_v), std::back_inserter(otr2_line));
    }
    break;
  default:
    return FALSE;
  }

  //printf("type = %d, depth = %d, channels = %d\n", dstEdge.type(), dstEdge.depth(), dstEdge.channels());

  MQObject o = MQ_CreateObject();
  char objname[151];
  doc->GetUnusedObjectName(objname, 150, "draw");
  o->SetName(objname);

  float wf = dstEdge.cols, hf = dstEdge.rows;

  int lnum = otr2_line.size();
  for(int i=0;i<lnum;i++)
  {
    Segment2 *l = &(otr2_line[i]);
    int idx[2];
    float x2 = l->source().x();
    //x2 = x2;
    float y2 = l->source().y();
    y2 = hf - y2;
    idx[0] = o->AddVertex(MQPoint(x2, y2, 0.0));
    x2 = l->target().x();
    //x2 = (x2 / wf) * 100.0f;
    y2 = l->target().y();
    y2 = hf - y2;
    idx[1] = o->AddVertex(MQPoint(x2, y2, 0.0));
    o->AddFace(2, idx);
  }

  o->OptimizeVertex(0.0f, NULL);
  doc->AddObject(o);
  
  //const char* source_window = "Source";
  //cv::namedWindow( source_window, cv::WINDOW_NORMAL );
  //imshow( source_window, dstEdge );
  //cv::waitKey(0);
  
  MQ_RefreshView(NULL);
  
  return TRUE;
}

