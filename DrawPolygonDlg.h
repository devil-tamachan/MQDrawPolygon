
#ifndef _TAMADRAWPOLYDLG_
#define _TAMADRAWPOLYDLG_

#include "MQSetting.h"

class MQDrawPolygonPlugin;

class DrawPolygonDialog : public MQDialog
{
public:
  DrawPolygonDialog(MQWindowBase& parent, MQDrawPolygonPlugin *plugin);
  ~DrawPolygonDialog();


  BOOL uiChanged(MQWidgetBase *sender, MQDocument doc)
  {
    UpdateEnable(doc);
    return FALSE;
  }
  
  void UpdateEnable(MQDocument doc)
  {
    //MakeObjList(doc);
    bool bCenterLine = combo_edgeproc->GetCurrentIndex()==3 ? false:true;
    bool bFaceGen = check_facegen->GetChecked();
    bool bThresholdMennuki = check_thresholdMennuki->GetChecked();
    slider_thresholdMennuki->SetEnabled(bCenterLine && bFaceGen && bThresholdMennuki);
    lbl_slider_thresholdMennuki->SetEnabled(slider_thresholdMennuki->GetEnabled());
    
    bool bConvModeFast = combo_src->GetCurrentIndex()==1;
    combo_zscale->SetEnabled(bConvModeFast);
    lbl_combo_zscale->SetEnabled(combo_zscale->GetEnabled());
    
    bool bZScaleModeInputValue = combo_zscale->GetCurrentIndex()==2;
    dblspin_zscale->SetEnabled(bZScaleModeInputValue);
    lbl_dblspin_zscale->SetEnabled(dblspin_zscale->GetEnabled());
    
    check_thresholdMennuki->SetEnabled(bCenterLine && bFaceGen);
    lbl_check_thresholdMennuki->SetEnabled(check_thresholdMennuki->GetEnabled());
    
    slider_threshold->SetEnabled( combo_vectorconv->GetCurrentIndex()==0 ? false:true);
    lbl_slider_threshold->SetEnabled(slider_threshold->GetEnabled());
    
    bool bClipboard = combo_src->GetCurrentIndex()==0 ? true : false;
    /*combo_filterobj->SetEnabled(!bClipboard);
    lbl_combo_filterobj->SetEnabled(combo_filterobj->GetEnabled());
    
    check_visibleObjOnly->SetEnabled(!bClipboard);
    lbl_check_visibleObjOnly->SetEnabled(check_visibleObjOnly->GetEnabled());
    */
  }
  BOOL cfgSettingsChanged(MQWidgetBase *sender, MQDocument doc)
  {
    MQSetting *config = OpenSetting();
    if(config==NULL)return FALSE;
    bool cfgAutoload = check_cfgAutoload->GetChecked();
    bool cfgSave = check_cfgSave->GetChecked();
    config->Save("cfgAutoload", cfgAutoload);
    config->Save("cfgSave", cfgSave);
    m_plugin->CloseSetting(config);
    return FALSE;
  }

/*
  void MakeObjList(MQDocument doc)
  {
    bool bVisibleObjOnly = check_visibleObjOnly->GetChecked();
    
    combo_filterobj->ClearItems();
    combo_filterobj->AddItem(L"すべてのオブジェクト");
    m_objIdx.clear();
    int numObj = doc->GetObjectCount();
    for(int i=0;i<numObj;i++)
    {
      MQObject o = doc->GetObject(i);
      if(o==NULL)continue;
      if(bVisibleObjOnly && o->GetVisible()==0)continue;
      m_objIdx.push_back(i);
      combo_filterobj->AddItem(o->GetNameW());
    }
  }

  int GetSelectObjIdx()
  {
    int i = combo_filterobj->GetCurrentIndex() - 1;
    if(i>=0 && i<m_objIdx.size())return m_objIdx[i];
    return -1;
  }
  */
  int GetOldIndex(int idxLatest, int SaveSize)
  {
    idxLatest--;
    if(idxLatest<0)idxLatest=SaveSize-1;
    if(idxLatest>=SaveSize)return -1;
    return idxLatest;
  }
  int GetNextIndex(int idxLatest, int SaveSize)
  {
    idxLatest++;
    if(idxLatest>=SaveSize)idxLatest=0;
    if(idxLatest<0)return -1;
    return idxLatest;
  }
  
  MQSetting* OpenSetting()
  {
    if(m_plugin==NULL)return NULL;
    MQSetting *config = m_plugin->OpenSetting();
    //if(config==NULL)return NULL;
    return config;
  }
  
  void SaveConfig()
  {
    int idxLatest, SaveSize;
    MQSetting *config = OpenSetting();
    if(config==NULL)return;
    
    config->Load("LatestIndex", idxLatest, -1);
    config->Load("SaveSize", SaveSize, 5);
    int newidx = GetNextIndex(idxLatest, SaveSize);
    _SaveConfig(config, newidx);
    
    m_plugin->CloseSetting(config);
  }
  void _SaveConfig(MQSetting *config, long long num = 0)
  {
    //if(m_plugin==NULL)return;
    //MQSetting *config = m_plugin->OpenSetting();
    //if(config==NULL)return;
    
    std::string snum = std::to_string(num);
    int itmp = 0;
    double dbltmp = 1.0;
    bool bltmp;
    
    itmp = combo_src->GetCurrentIndex();
    config->Save(("combo_src"+snum).c_str(), itmp);
    
    itmp = combo_zscale->GetCurrentIndex();
    config->Save(("combo_zscale"+snum).c_str(), itmp);
    
    dbltmp = dblspin_zscale->GetPosition();
    config->Save(("dblspin_zscale"+snum).c_str(), dbltmp);
    
    itmp = combo_edgeproc->GetCurrentIndex();
    config->Save(("combo_edgeproc"+snum).c_str(), itmp);
    
    itmp = combo_vectorconv->GetCurrentIndex();
    config->Save(("combo_vectorconv"+snum).c_str(), itmp);
    
    dbltmp = slider_threshold->GetPosition();
    config->Save(("slider_threshold"+snum).c_str(), dbltmp);
    
    itmp = spin_np->GetPosition();
    config->Save(("spin_np"+snum).c_str(), itmp);
    
    bltmp = check_facegen->GetChecked();
    config->Save(("check_facegen"+snum).c_str(), bltmp);
    
    bltmp = check_thresholdMennuki->GetChecked();
    config->Save(("check_thresholdMennuki"+snum).c_str(), bltmp);
    
    dbltmp = slider_thresholdMennuki->GetPosition();
    config->Save(("slider_thresholdMennuki"+snum).c_str(), dbltmp);
    
    bltmp = check_optimize->GetChecked();
    config->Save(("check_optimize"+snum).c_str(), bltmp);
    
    itmp = num;
    config->Save("LatestIndex", itmp);
    
    //m_plugin->CloseSetting(config);
  }
  
  /*
  void LoadConfig()
  {
    int idxLatest;
    
    if(m_plugin==NULL)return;
    MQSetting *config = m_plugin->OpenSetting();
    if(config==NULL)return;
    
    config->Load("LatestIndex", idxLatest, 0);
    
    _LoadConfig(config, idxLatest);
    
    m_plugin->CloseSetting(config);
  }*/
  bool LoadConfig(int before = 0)
  {
    int idxLatest, SaveSize;
    if(before<0)return false;
    
    MQSetting *config = OpenSetting();
    if(config==NULL)return false;
    
    config->Load("LatestIndex", idxLatest, 0);
    if(before>0)
    {
      config->Load("SaveSize", SaveSize, 5);
      if(before>SaveSize-1)return false;
      for(int i=0;i<before;i++)
      {
        idxLatest = GetOldIndex(idxLatest, SaveSize);
      }
    }
    
    bool bRet = _LoadConfig(config, idxLatest);
    
    m_plugin->CloseSetting(config);
    return bRet;
  }
  bool _LoadConfig(MQSetting *config, long long num = 0)
  {
    //if(m_plugin==NULL)return;
    //MQSetting *config = m_plugin->OpenSetting();
    //if(config==NULL)return;
    
    std::string snum = std::to_string(num);
    int itmp = 0;
    double dbltmp = 1.0;
    bool bltmp;
    //MQSetting::Load(const char *name, int& value, int default_value=0);
    
    config->Load(("combo_src"+snum).c_str(), itmp, 0);
    combo_src->SetCurrentIndex(itmp);
    
    config->Load(("combo_zscale"+snum).c_str(), itmp, 1);
    combo_zscale->SetCurrentIndex(itmp);
    
    config->Load(("dblspin_zscale"+snum).c_str(), dbltmp, 1.0);
    dblspin_zscale->SetPosition(dbltmp);
    
    config->Load(("combo_edgeproc"+snum).c_str(), itmp, 0);
    combo_edgeproc->SetCurrentIndex(itmp);
    
    config->Load(("combo_vectorconv"+snum).c_str(), itmp, 0);
    combo_vectorconv->SetCurrentIndex(itmp);
    
    config->Load(("slider_threshold"+snum).c_str(), dbltmp, 0.5);
    slider_threshold->SetPosition(dbltmp);
    
    config->Load(("spin_np"+snum).c_str(), itmp, 100);
    spin_np->SetPosition(itmp);
    
    config->Load(("check_facegen"+snum).c_str(), bltmp, true);
    check_facegen->SetChecked(bltmp);
    
    config->Load(("check_thresholdMennuki"+snum).c_str(), bltmp, true);
    check_thresholdMennuki->SetChecked(bltmp);
    
    config->Load(("slider_thresholdMennuki"+snum).c_str(), dbltmp, 0.5);
    slider_thresholdMennuki->SetPosition(dbltmp);
    
    config->Load(("check_optimize"+snum).c_str(), bltmp, true);
    check_optimize->SetChecked(bltmp);
    
    //m_plugin->CloseSetting(config);
    return true;
  }
  
  BOOL OnRecent0(MQWidgetBase *sender, MQDocument doc)
  {
    if(!LoadConfig(0))InitConfig();
    return TRUE;
  }
  BOOL OnRecent1(MQWidgetBase *sender, MQDocument doc)
  {
    if(!LoadConfig(1))InitConfig();
    return TRUE;
  }
  BOOL OnRecent2(MQWidgetBase *sender, MQDocument doc)
  {
    if(!LoadConfig(2))InitConfig();
    return TRUE;
  }
  BOOL OnRecent3(MQWidgetBase *sender, MQDocument doc)
  {
    if(!LoadConfig(3))InitConfig();
    return TRUE;
  }
  BOOL OnRecent4(MQWidgetBase *sender, MQDocument doc)
  {
    if(!LoadConfig(4))InitConfig();
    return TRUE;
  }
  void InitConfig()
  {
    combo_src->SetCurrentIndex(0);
    combo_zscale->SetCurrentIndex(1);
    dblspin_zscale->SetPosition(1.0);
    combo_edgeproc->SetCurrentIndex(0);
    combo_vectorconv->SetCurrentIndex(0);
    slider_threshold->SetPosition(0.5);
    spin_np->SetPosition(100);
    check_facegen->SetChecked(true);
    check_thresholdMennuki->SetChecked(true);
    slider_thresholdMennuki->SetPosition(0.5);
    check_optimize->SetChecked(true);
  }
  BOOL OnCfgInit(MQWidgetBase *sender, MQDocument doc)
  {
    InitConfig();
    return TRUE;
  }
  
  BOOL OnOK(MQWidgetBase *sender, MQDocument doc)
  {
    if(!check_cfgSave->GetChecked())return FALSE;
    SaveConfig();
    return FALSE;
  }
  
  MQDrawPolygonPlugin *m_plugin;

  MQComboBox *combo_src;
  MQComboBox *combo_filterobj;
  MQLabel *lbl_combo_filterobj;
  MQCheckBox *check_visibleObjOnly;
  MQLabel *lbl_check_visibleObjOnly;
  MQComboBox *combo_zscale;
  MQLabel *lbl_combo_zscale;
  MQDoubleSpinBox *dblspin_zscale;
  MQLabel *lbl_dblspin_zscale;
  MQComboBox *combo_edgeproc;
  MQComboBox *combo_vectorconv;
  MQSlider *slider_threshold;
  MQLabel *lbl_slider_threshold;
  MQSpinBox *spin_np;
  MQCheckBox *check_facegen;
  MQSlider *slider_thresholdMennuki;
  MQLabel *lbl_slider_thresholdMennuki;
  MQCheckBox *check_thresholdMennuki;
  MQLabel *lbl_check_thresholdMennuki;
  MQCheckBox *check_optimize;
  std::vector<int> m_objIdx;
  
  
  MQCheckBox *check_cfgAutoload;
  MQCheckBox *check_cfgSave;
};

DrawPolygonDialog::DrawPolygonDialog(MQWindowBase& parent, MQDrawPolygonPlugin *plugin) : MQDialog(parent)
{
  m_plugin = plugin;
  SetTitle(L"DrawPolygon");

  MQFrame *mainFrame = CreateHorizontalFrame(this);

  MQFrame *paramFrame = CreateHorizontalFrame(mainFrame);
  paramFrame->SetMatrixColumn(2);
  
  MQSetting *config = OpenSetting();
  if(config==NULL)return;
  bool cfgAutoload, cfgSave;
  config->Load("cfgAutoload", cfgAutoload, true);
  config->Load("cfgSave", cfgSave, true);
  m_plugin->CloseSetting(config);
  
  CreateLabel(paramFrame, L"入力画像");
  combo_src = CreateComboBox(paramFrame);
  combo_src->AddItem(L"コピペ画像");
  combo_src->AddItem(L"面テクスチャ/高速");
  combo_src->AddItem(L"面テクスチャ/低速");
  combo_src->SetCurrentIndex(0);
  combo_src->AddChangedEvent(this, &DrawPolygonDialog::uiChanged);
  
  /*
  lbl_combo_filterobj = CreateLabel(paramFrame, L"入力オブジェクト");
  combo_filterobj = CreateComboBox(paramFrame);
  combo_filterobj->SetCurrentIndex(0);

  lbl_check_visibleObjOnly = CreateLabel(paramFrame, L"\"入力オブジェクト\"から非表示オブジェを除く");
  check_visibleObjOnly = CreateCheckBox(paramFrame);
  check_visibleObjOnly->SetChecked(true);
  check_visibleObjOnly->AddChangedEvent(this, &DrawPolygonDialog::uiChanged);
  */
  
  lbl_combo_zscale = CreateLabel(paramFrame, L"立体テクスチャのZ伸縮");
  combo_zscale = CreateComboBox(paramFrame);
  combo_zscale->AddItem(L"伸縮しない");
  combo_zscale->AddItem(L"UV面積比でZ伸縮");
  combo_zscale->AddItem(L"固定値");
  combo_zscale->SetCurrentIndex(1);
  combo_zscale->AddChangedEvent(this, &DrawPolygonDialog::uiChanged);
  
  lbl_dblspin_zscale = CreateLabel(paramFrame, L"立体テクスチャのZ伸縮率");
  dblspin_zscale = CreateDoubleSpinBox(paramFrame);
  dblspin_zscale->SetMin(0.0000001);
  dblspin_zscale->SetPosition(1.0);
  
  CreateLabel(paramFrame, L"ライン抽出");
  combo_edgeproc = CreateComboBox(paramFrame);
  combo_edgeproc->AddItem(L"Canny (おすすめ)");
  combo_edgeproc->AddItem(L"Laplacian");
  combo_edgeproc->AddItem(L"Sobel");
  combo_edgeproc->AddItem(L"無し（中央線）");
  combo_edgeproc->SetCurrentIndex(0);
  combo_edgeproc->AddChangedEvent(this, &DrawPolygonDialog::uiChanged);
  
  CreateLabel(paramFrame, L"ベクトル化");
  combo_vectorconv = CreateComboBox(paramFrame);
  combo_vectorconv->AddItem(L"モノクロ値を重みとして使用");
  combo_vectorconv->AddItem(L"2値化（重み無し））");
  combo_vectorconv->SetCurrentIndex(0);
  combo_vectorconv->AddChangedEvent(this, &DrawPolygonDialog::uiChanged);
  
  lbl_slider_threshold = CreateLabel(paramFrame, L"2値化の判定値");
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

  CreateLabel(paramFrame, L"面を生成");
  check_facegen = CreateCheckBox(paramFrame);
  check_facegen->SetChecked(true);
  check_facegen->AddChangedEvent(this, &DrawPolygonDialog::uiChanged);

  lbl_check_thresholdMennuki = CreateLabel(paramFrame, L"面抜きする");
  check_thresholdMennuki = CreateCheckBox(paramFrame);
  check_thresholdMennuki->SetChecked(true);
  check_thresholdMennuki->AddChangedEvent(this, &DrawPolygonDialog::uiChanged);
  
  
  lbl_slider_thresholdMennuki = CreateLabel(paramFrame, L"面抜き判定値 (多く残る→)");
  slider_thresholdMennuki = CreateSlider(paramFrame);
  slider_thresholdMennuki->SetMin(0.01);
  slider_thresholdMennuki->SetMax(0.99);
  slider_thresholdMennuki->SetPosition(0.5);
  slider_thresholdMennuki->SetEnabled(true);

  CreateLabel(paramFrame, L"近接する頂点を接合");
  check_optimize = CreateCheckBox(paramFrame);
  check_optimize->SetChecked(true);

  MQFrame *sideFrame = CreateVerticalFrame(mainFrame);

  MQButton *okbtn = CreateButton(sideFrame, L"OK");
  okbtn->SetDefault(true);
  okbtn->SetModalResult(MQDialog::DIALOG_OK);
  okbtn->AddClickEvent(this, &DrawPolygonDialog::OnOK, true);

  MQButton *cancelbtn = CreateButton(sideFrame, L"Cancel");
  cancelbtn->SetCancel(true);
  cancelbtn->SetModalResult(MQDialog::DIALOG_CANCEL);

  CreateLabel(sideFrame, L"起動時に前回設定");
  check_cfgAutoload = CreateCheckBox(sideFrame);
  check_cfgAutoload->SetChecked(cfgAutoload);
  check_cfgAutoload->AddChangedEvent(this, &DrawPolygonDialog::cfgSettingsChanged);

  CreateLabel(sideFrame, L"設定を保存");
  check_cfgSave = CreateCheckBox(sideFrame);
  check_cfgSave->SetChecked(cfgSave);
  check_cfgSave->AddChangedEvent(this, &DrawPolygonDialog::cfgSettingsChanged);

  MQButton *recent0btn = CreateButton(sideFrame, L"前回設定");
  recent0btn->AddClickEvent(this, &DrawPolygonDialog::OnRecent0);

  MQButton *recent1btn = CreateButton(sideFrame, L"前前回設定");
  recent1btn->AddClickEvent(this, &DrawPolygonDialog::OnRecent1);

  MQButton *recent2btn = CreateButton(sideFrame, L"3つ前の設定");
  recent2btn->AddClickEvent(this, &DrawPolygonDialog::OnRecent2);

  MQButton *recent3btn = CreateButton(sideFrame, L"4つ前の設定");
  recent3btn->AddClickEvent(this, &DrawPolygonDialog::OnRecent3);

  MQButton *recent4btn = CreateButton(sideFrame, L"5つ前の設定");
  recent4btn->AddClickEvent(this, &DrawPolygonDialog::OnRecent4);

  MQButton *cfgInitBtn = CreateButton(sideFrame, L"初期設定");
  cfgInitBtn->AddClickEvent(this, &DrawPolygonDialog::OnCfgInit);
  
  if(cfgAutoload)LoadConfig();
}

DrawPolygonDialog::~DrawPolygonDialog()
{
}

#endif //_TAMADRAWPOLYDLG_
