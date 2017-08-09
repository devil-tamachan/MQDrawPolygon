
#ifndef _TAMADRAWPOLYDLG_
#define _TAMADRAWPOLYDLG_



class DrawPolygonDialog : public MQDialog
{
public:
  DrawPolygonDialog(MQWindowBase& parent);
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
};

DrawPolygonDialog::DrawPolygonDialog(MQWindowBase& parent) : MQDialog(parent)
{
  SetTitle(L"DrawPolygon");

  MQFrame *mainFrame = CreateHorizontalFrame(this);

  MQFrame *paramFrame = CreateHorizontalFrame(mainFrame);
  paramFrame->SetMatrixColumn(2);
  
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

  MQButton *cancelbtn = CreateButton(sideFrame, L"Cancel");
  cancelbtn->SetCancel(true);
  cancelbtn->SetModalResult(MQDialog::DIALOG_CANCEL);
}

DrawPolygonDialog::~DrawPolygonDialog()
{
}

#endif //_TAMADRAWPOLYDLG_
