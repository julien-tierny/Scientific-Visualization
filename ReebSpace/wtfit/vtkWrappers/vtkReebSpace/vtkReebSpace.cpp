#include                  <vtkReebSpace.h>

vtkStandardNewMacro(vtkReebSpace)

vtkReebSpace::vtkReebSpace(){

  // init
  SetNumberOfOutputPorts(4);
  VaryingConnectivity = false;
  VaryingValues = false;
  perturbated_ = false;
  
  ZeroSheetId = true;
  ZeroSheetType = true;
  ZeroSheetValue = true;
  ZeroSheetVertexId = true;
  
  OneSheetId = true;
  OneSheetType = true;
  OneSheetValue = true;
  OneSheetVertexId = true;
  OneSheetEdgeId = true;
  
  TwoSheets = true;
  TwoSheetCaseId = true;
  TwoSheetEdgeId = true;
  TwoSheetEdgeType = true;
  TwoSheetId = true;
  TwoSheetParameterization = true;
  TwoSheetTetId = true;
  TwoSheetValue = true;
  
  ThreeSheetTetNumber = true;
  ThreeSheetVertexNumber = true;
  ThreeSheetExpansion = true;
  ThreeSheetDomainVolume = true;
  ThreeSheetRangeArea = true;
  ThreeSheetHyperVolume = true;
  
  SimplificationThreshold = 0;
  SimplificationCriterion = 1;
  
  UseAllCores = false;
}

vtkReebSpace::~vtkReebSpace(){

}


// transmit abort signals -- to copy paste in other wrappers
bool vtkReebSpace::needsToAbort(){
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int vtkReebSpace::updateProgress(const float &progress){

  {
    stringstream msg;
    msg << "[vtkReebSpace] " << progress*100 
      << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }
  
  UpdateProgress(progress);
  return 0;
}

template<class dataTypeU, class dataTypeV> int vtkReebSpace::baseCall(
  vtkUnstructuredGrid *input, vtkDataArray *uField, vtkDataArray *vField){

  Timer t;
  bool preProcess = false;
  long long int *tetList = (long long int *) input->GetCells()->GetPointer();
  
  reebSpace_.setWrapper(this);
  reebSpace_.setTetList(tetList);
  
  reebSpace_.setVertexNumber(input->GetNumberOfPoints());
  reebSpace_.setTetNumber(input->GetNumberOfCells());
  // NOTE: assuming VTK always use float for points
  reebSpace_.setPointSet((float *) input->GetPoints()->GetVoidPointer(0));
  
  // connectivity pre-processing
  if((edgeList_.empty())
    ||(VaryingConnectivity)){
    
    vtkOneSkeleton oneSkeleton;
    oneSkeleton.setThreadNumber(threadNumber_);
    oneSkeleton.buildEdgeList(input, edgeList_, true);
    preProcess = true;
  }
  reebSpace_.setEdgeList(&edgeList_);
  
  if((vertexEdgeList_.empty())
    ||(VaryingConnectivity)){
    
    vtkZeroSkeleton zeroSkeleton;
    zeroSkeleton.setThreadNumber(threadNumber_);
    zeroSkeleton.buildVertexEdges(input, 
      vertexEdgeList_,
      &edgeList_, true);
    preProcess = true;
  }
  reebSpace_.setVertexEdgeList(&vertexEdgeList_);
  
  if((edgeFans_.empty())
    ||(edgeFanLinkEdgeList_.empty())
    ||(VaryingConnectivity)){
    
    vtkOneSkeleton oneSkeleton;
    oneSkeleton.setThreadNumber(threadNumber_);
    
    vertexStars_.clear();
    oneSkeleton.buildEdgeStars(
      input, edgeStars_, &edgeList_, &vertexStars_, true);
    
    reebSpace_.setVertexStars(&vertexStars_);
    reebSpace_.setEdgeStars(&edgeStars_);
    reebSpace_.connectivityPreprocessing<dataTypeU, dataTypeV>(edgeStars_,
      edgeFanLinkEdgeList_, edgeFans_, sosOffsets_);
    preProcess = true;
  }
  
  if((tetNeighbors_.empty())||(VaryingConnectivity)){
    
    // this guy is an actual vtk filter
    vtkSmartPointer<vtkThreeSkeleton> threeSkeleton 
      = vtkSmartPointer<vtkThreeSkeleton>::New();
      
    threeSkeleton->setThreadNumber(threadNumber_);
    threeSkeleton->buildCellNeighbors(input, tetNeighbors_, 
      &vertexStars_, NULL, true);
    reebSpace_.setTetNeighbors(&tetNeighbors_);
    preProcess = true;
  }
  
  reebSpace_.setEdgeFans(&edgeFans_);
  reebSpace_.setEdgeFanLinkEdgeList(&edgeFanLinkEdgeList_);
  reebSpace_.setSosOffsets(&sosOffsets_);
  reebSpace_.setExpand3Sheets(ThreeSheetExpansion);
  
    
  if(preProcess){
    stringstream msg;
    msg << "[vtkReebSpace] VTK connectivity pre-processed." << endl;
    msg << "[vtkReebSpace] First field: '"
      << uField->GetName() << "'" << endl;
    msg << "[vtkReebSpace] Second field: '"
      << vField->GetName() << "'" << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  reebSpace_.setInputField(uField->GetVoidPointer(0),
    vField->GetVoidPointer(0));
  
  // simulation of simplicity
  if((!perturbated_)||(VaryingValues)){
    
    switch(uField->GetDataType()){
    
      case VTK_DOUBLE:
        switch(vField->GetDataType()){
          case VTK_DOUBLE:
            reebSpace_.perturbate(
              (dataTypeU) pow(10, -DBL_DIG+1), 
              (dataTypeV) pow(10, -DBL_DIG+1));
            break;
            
          case VTK_FLOAT:
            reebSpace_.perturbate(
              (dataTypeU) pow(10, -DBL_DIG+1), 
              (dataTypeV) pow(10, -FLT_DIG+1));
            break;
            
          default:
            // same epsilon for char, unsigned char, int, etc.
            reebSpace_.perturbate(
              (dataTypeU) pow(10, -DBL_DIG+1), 
              (dataTypeV) 2);
            break;
        }
        break;

      case VTK_FLOAT:
        switch(vField->GetDataType()){
          case VTK_DOUBLE:
            reebSpace_.perturbate(
              (dataTypeU) pow(10, -FLT_DIG+1), 
              (dataTypeV) pow(10, -DBL_DIG+1));
            break;
            
          case VTK_FLOAT:
            reebSpace_.perturbate(
              (dataTypeU) pow(10, -FLT_DIG+1), 
              (dataTypeV) pow(10, -FLT_DIG+1));
            break;
            
          default:
            // same epsilon for char, unsigned char, int, etc.
            reebSpace_.perturbate(
              (dataTypeU) pow(10, -FLT_DIG+1), 
              (dataTypeV) 2);
            break;
        }
        break;
        
      default:
        // same epsilon for char, unsigned char, int, etc.
        switch(vField->GetDataType()){
          
          case VTK_DOUBLE:
            reebSpace_.perturbate(
              (dataTypeU) 2, 
              (dataTypeV) pow(10, -DBL_DIG+1));
            break;
            
          case VTK_FLOAT:
            reebSpace_.perturbate(
              (dataTypeU) 2, 
              (dataTypeV) pow(10, -FLT_DIG+1));
            break;
            
          default:
            reebSpace_.perturbate(
              (dataTypeU) 2, 
              (dataTypeV) 2);
            break;
        }
        break;
    }
    
    perturbated_ = true;
    preProcess = true;
  }
  
  // go!
  if((reebSpace_.empty())||(VaryingConnectivity)||(VaryingValues)){
    {
      stringstream msg;
      msg << "[vtkReebSpace] Starting computation..." << endl;
      
      dMsg(cout, msg.str(), timeMsg);
    }
    reebSpace_.execute<dataTypeU, dataTypeV>();
  }
  
  if(SimplificationThreshold > 0){
    reebSpace_.simplify<dataTypeU, dataTypeV>(
      SimplificationThreshold, 
      (ReebSpace::SimplificationCriterion) SimplificationCriterion);
  }
  
  return 0;
}

int vtkReebSpace::doIt(vtkUnstructuredGrid *input, 
  vtkUnstructuredGrid *sheet0, vtkUnstructuredGrid *sheet1, 
  vtkUnstructuredGrid *sheet2, vtkUnstructuredGrid *sheet3){

  // check data components
  vtkDataArray *uComponent = NULL, *vComponent = NULL;
  
  if(Ucomponent.length()){
    uComponent = input->GetPointData()->GetArray(Ucomponent.data());
  }
  else{
    // default
    uComponent = input->GetPointData()->GetArray(0);
  }
  if(!uComponent)
    return -1;
  
  if(Vcomponent.length()){
    vComponent = input->GetPointData()->GetArray(Vcomponent.data());
  }
  else{
    // default
    vComponent = input->GetPointData()->GetArray(0);
  }
  if(!vComponent)
    return -2;
  
  // set the Reeb space functor
  switch(uComponent->GetDataType()){
    
    case VTK_CHAR:
      switch(vComponent->GetDataType()){
        vtkTemplateMacro((
          {
            baseCall<char, VTK_TT>(input, uComponent, vComponent);
          }
        ));
      }
      break;
      
    case VTK_DOUBLE:
      switch(vComponent->GetDataType()){
        vtkTemplateMacro((
          {
            baseCall<double, VTK_TT>(input, uComponent, vComponent);
          }
        ));
      }
      break;
      
    case VTK_FLOAT:
      switch(vComponent->GetDataType()){
        vtkTemplateMacro((
          {
            baseCall<float, VTK_TT>(input, uComponent, vComponent);
          }
        ));
      }
      break;
      
    case VTK_INT:
      switch(vComponent->GetDataType()){
        vtkTemplateMacro((
          {
            baseCall<int, VTK_TT>(input, uComponent, vComponent);
          }
        ));
      }
      break;
      
    case VTK_UNSIGNED_CHAR:
      switch(vComponent->GetDataType()){
        vtkTemplateMacro((
          {
            baseCall<unsigned char, VTK_TT>(input, uComponent, vComponent);
          }
        ));
      }
      break;
      
    case VTK_UNSIGNED_SHORT:
      switch(vComponent->GetDataType()){
        vtkTemplateMacro((
          {
            baseCall<unsigned short, VTK_TT>(input, uComponent, vComponent);
          }
        ));
      }
      break;
      
    default:
      {
        stringstream msg;
        msg << "[vtkReebSpace] Unsupported U-component data type :( ["
          << uComponent->GetDataType() << "]" << endl;
        dMsg(cerr, msg.str(), 1);
      }
      break;
  }
  
  // prepare the output
  {
    stringstream msg;
    msg << "[vtkReebSpace] Preparing the VTK-output..." << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
 
  // 0-sheets - 
  // Optional additional fields:
  // PointData; u, v, vertexId, type, sheetId
  const vector<int> *sheet0segmentation = reebSpace_.get0sheetSegmentation();
  int vertexNumber = 0;
  for(int i = 0; i < (int) sheet0segmentation->size(); i++){
    int sheet0Id = (*sheet0segmentation)[i];
    if(sheet0Id != -1){
      const ReebSpace::Sheet0 *sheet = reebSpace_.get0sheet(sheet0Id);
      if(!sheet->pruned_){
        vertexNumber++;
      }
    }
  }
  vtkSmartPointer<vtkPoints> sheet0Points = 
    vtkSmartPointer<vtkPoints>::New();
  
  if(!sheet0->GetPoints())
    sheet0->SetPoints(sheet0Points);
    
  sheet0->GetPoints()->SetNumberOfPoints(vertexNumber);
  
  vtkSmartPointer<vtkDoubleArray> vertexScalarsU = 
    vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> vertexScalarsV = 
    vtkSmartPointer<vtkDoubleArray>::New();
  if(ZeroSheetValue){
    vertexScalarsU->SetNumberOfTuples(vertexNumber);
    vertexScalarsU->SetName(uComponent->GetName());
    vertexScalarsV->SetNumberOfTuples(vertexNumber);
    vertexScalarsV->SetName(vComponent->GetName());
  }
  else{
    sheet0->GetPointData()->RemoveArray(uComponent->GetName());
    sheet0->GetPointData()->RemoveArray(vComponent->GetName());
  }
  
  vtkSmartPointer<vtkIntArray> vertexIds = 
    vtkSmartPointer<vtkIntArray>::New();
  if(ZeroSheetVertexId){
    vertexIds->SetNumberOfTuples(vertexNumber);
    vertexIds->SetName("VertexIds");
  }
  else{
    sheet0->GetPointData()->RemoveArray("VertexIds");
  }
  
  vtkSmartPointer<vtkCharArray> vertexTypes = 
    vtkSmartPointer<vtkCharArray>::New();
  if(ZeroSheetType){
    vertexTypes->SetNumberOfTuples(vertexNumber);
    vertexTypes->SetName("SheetType");
  }
  else{
    sheet0->GetPointData()->RemoveArray("SheetType");
  }
  
  vtkSmartPointer<vtkIntArray> vertexSheetId =
    vtkSmartPointer<vtkIntArray>::New();
  if(ZeroSheetId){
    vertexSheetId->SetNumberOfTuples(vertexNumber);
    vertexSheetId->SetName("0-SheetId");
  }
  else{
    sheet0->GetPointData()->RemoveArray("0-SheetId");
  }
  
  vertexNumber = 0;
  double *p = NULL;
  for(int i = 0; i < (int) sheet0segmentation->size(); i++){
    int sheet0Id = (*sheet0segmentation)[i];
    if(sheet0Id != -1){
      
      const ReebSpace::Sheet0 *sheet = reebSpace_.get0sheet(sheet0Id);
      
      if(!sheet->pruned_){
        p = input->GetPoint(i);
        
        sheet0->GetPoints()->SetPoint(vertexNumber, p);
        
        if(ZeroSheetId){
          vertexSheetId->SetTuple1(vertexNumber, 
            (*sheet0segmentation)[i]);
        }
        if(ZeroSheetVertexId){
          vertexIds->SetTuple1(vertexNumber, i);
        }
        if(ZeroSheetValue){
          double u, v;
          
          uComponent->GetTuple(i, &u);
          vComponent->GetTuple(i, &v);
          
          vertexScalarsU->SetTuple1(vertexNumber, u);
          vertexScalarsV->SetTuple1(vertexNumber, v);
        }
        if(ZeroSheetType){
          const ReebSpace::Sheet0* sheet = 
            reebSpace_.get0sheet((*sheet0segmentation)[i]);
          vertexTypes->SetTuple1(vertexNumber, sheet->type_);
        }
      
        vertexNumber++;
      }
    }
  }
  if(ZeroSheetId)
    sheet0->GetPointData()->AddArray(vertexSheetId);
  if(ZeroSheetVertexId)
    sheet0->GetPointData()->AddArray(vertexIds);
  if(ZeroSheetValue){
    sheet0->GetPointData()->AddArray(vertexScalarsU);
    sheet0->GetPointData()->AddArray(vertexScalarsV);
  }
  if(ZeroSheetType)
    sheet0->GetPointData()->AddArray(vertexTypes);
    
  
  // 1-sheets
  // Optional additional fields:
  // PointData: u, v, vertexId, 
  // CellData: edgeId, type, sheetId
  const vector<int> *sheet1segmentation = reebSpace_.get1sheetSegmentation();
  
  vtkSmartPointer<vtkPoints> sheet1Points = 
    vtkSmartPointer<vtkPoints>::New();
  sheet1->SetPoints(sheet1Points);
  
  vtkSmartPointer<vtkCellArray> sheet1Edges =
    vtkSmartPointer<vtkCellArray>::New();
  
  vtkSmartPointer<vtkDoubleArray> edgeScalarsU = 
    vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> edgeScalarsV = 
    vtkSmartPointer<vtkDoubleArray>::New();
    
  if(OneSheetValue){
    edgeScalarsU->SetName(uComponent->GetName());
    edgeScalarsV->SetName(vComponent->GetName());
  }
  else{
    sheet1->GetPointData()->RemoveArray(uComponent->GetName());
    sheet1->GetPointData()->RemoveArray(vComponent->GetName());
  }
 
  vtkSmartPointer<vtkIntArray> edgeVertexIds = 
    vtkSmartPointer<vtkIntArray>::New();
  if(OneSheetVertexId){
    edgeVertexIds->SetName("VertexIds");
  }
  else{
    sheet1->GetPointData()->RemoveArray("VertexIds");
  }
  
  vtkSmartPointer<vtkIntArray> edgeType = 
    vtkSmartPointer<vtkIntArray>::New();
  if(OneSheetType){
    edgeType->SetName("EdgeType");
  }
  else{
    sheet1->GetCellData()->RemoveArray("EdgeType");
  }
  
  vtkSmartPointer<vtkIntArray> edgeIds = 
    vtkSmartPointer<vtkIntArray>::New();
  if(OneSheetEdgeId){
    edgeIds->SetName("EdgeIds");
  }
  else{
    sheet1->GetCellData()->RemoveArray("EdgeIds");
  }
  
  vtkSmartPointer<vtkIntArray> edgeSheetIds = 
    vtkSmartPointer<vtkIntArray>::New();
  if(OneSheetId){
    edgeSheetIds->SetName("1-SheetId");
  }
  else{
    sheet1->GetCellData()->RemoveArray("1-SheetId");
  }
  
  vertexNumber = 0;
  double p0[3], p1[3];
  vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
  idList->SetNumberOfIds(2);
  const vector<int> *edgeTypes = reebSpace_.getEdgeTypes();
  
  for(int i = 0; i < (int) sheet1segmentation->size(); i++){
    
    int sheet1Id = (*sheet1segmentation)[i];
    
    if(sheet1Id != -1){
      
      const ReebSpace::Sheet1 *sheet = reebSpace_.get1sheet(sheet1Id);
      
      if(!sheet->pruned_){
      
        input->GetPoint(edgeList_[i].first, p0);
        input->GetPoint(edgeList_[i].second, p1);
        
        sheet1->GetPoints()->InsertNextPoint(p0);
        sheet1->GetPoints()->InsertNextPoint(p1);
        
        if(OneSheetValue){
          double u, v;
          uComponent->GetTuple(edgeList_[i].first, &u);
          vComponent->GetTuple(edgeList_[i].first, &v);
          
          edgeScalarsU->InsertNextTuple1(u);
          edgeScalarsV->InsertNextTuple1(v);
          
          uComponent->GetTuple(edgeList_[i].second, &u);
          vComponent->GetTuple(edgeList_[i].second, &v);
          
          edgeScalarsU->InsertNextTuple1(u);
          edgeScalarsV->InsertNextTuple1(v);
        }
        
        if(OneSheetVertexId){
          edgeVertexIds->InsertNextTuple1(edgeList_[i].first);
          edgeVertexIds->InsertNextTuple1(edgeList_[i].second);
        }
      
        idList->SetId(0, vertexNumber);
        idList->SetId(1, vertexNumber + 1);
        sheet1Edges->InsertNextCell(idList);
        vertexNumber += 2;
        
        if(OneSheetEdgeId){
          edgeIds->InsertNextTuple1(i);
        }
        if(OneSheetType){
          edgeType->InsertNextTuple1((*edgeTypes)[i]);
        }
        if(OneSheetId){
          edgeSheetIds->InsertNextTuple1((*sheet1segmentation)[i]);
        }
      }
    }
  }
  sheet1->SetCells(VTK_LINE, sheet1Edges);
  
  if(OneSheetValue){
    sheet1->GetPointData()->AddArray(edgeScalarsU);
    sheet1->GetPointData()->AddArray(edgeScalarsV);
  }
  if(OneSheetVertexId){
    sheet1->GetPointData()->AddArray(edgeVertexIds);
  }
  if(OneSheetId){
    sheet1->GetCellData()->AddArray(edgeSheetIds);
  }
  if(OneSheetEdgeId){
    sheet1->GetCellData()->AddArray(edgeIds);
  }
  if(OneSheetType){
    sheet1->GetCellData()->AddArray(edgeType);
  }
  
  // 2-sheets
  // optional fields:
  // pointdata: twoSheetValues, twoSheetParameterization
  if(TwoSheets){
    const vector<FiberSurface::Vertex> *vertexList = 
      reebSpace_.getFiberSurfaceVertices();
      
    vtkSmartPointer<vtkPoints> sheet2Points = vtkSmartPointer<vtkPoints>::New();
    sheet2Points->SetNumberOfPoints(vertexList->size());
    sheet2->SetPoints(sheet2Points);
      
    vtkSmartPointer<vtkDoubleArray> triangleScalarsU = 
      vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> triangleScalarsV =
      vtkSmartPointer<vtkDoubleArray>::New();
      
    if(TwoSheetValue){
      triangleScalarsU->SetName(uComponent->GetName());
      triangleScalarsU->SetNumberOfTuples(vertexList->size());
      triangleScalarsV->SetName(vComponent->GetName());
      triangleScalarsV->SetNumberOfTuples(vertexList->size());
    }
    else{
      sheet2->GetPointData()->RemoveArray(uComponent->GetName());
      sheet2->GetPointData()->RemoveArray(vComponent->GetName());
    }
    
    vtkSmartPointer<vtkDoubleArray> triangleParameterization = 
      vtkSmartPointer<vtkDoubleArray>::New();
    if(TwoSheetParameterization){
      triangleParameterization->SetName("EdgeParameterization");
      triangleParameterization->SetNumberOfTuples(vertexList->size());
    }
    else{
      sheet2->GetPointData()->RemoveArray("EdgeParameterization");
    }
    
    int sheet2TriangleNumber = 0;
    for(int i = 0; i < reebSpace_.getNumberOf2sheets(); i++){
      const ReebSpace::Sheet2 *sheet = reebSpace_.get2sheet(i);
      
      if(!sheet->pruned_){
        for(int j = 0; j  < (int) sheet->triangleList_.size(); j++){
          sheet2TriangleNumber += sheet->triangleList_[j].size();
        }
      }
    }
    
    vtkSmartPointer<vtkCellArray> sheet2Triangles = 
    vtkSmartPointer<vtkCellArray>::New();

    // celldata: twoSheetId, twoSheetEdgeId, twoSheetTetId
    vtkSmartPointer<vtkIntArray> triangleSheetIds = 
      vtkSmartPointer<vtkIntArray>::New();
    if(TwoSheetId){
      triangleSheetIds->SetName("2-SheetId");
      triangleSheetIds->SetNumberOfTuples(sheet2TriangleNumber);
    }
    else{
      sheet2->GetCellData()->RemoveArray("2-SheetId");
    }
    
    vtkSmartPointer<vtkIntArray> triangleEdgeIds = 
      vtkSmartPointer<vtkIntArray>::New();
    if(TwoSheetEdgeId){
      triangleEdgeIds->SetName("EdgeIds");
      triangleEdgeIds->SetNumberOfTuples(sheet2TriangleNumber);
    }
    else{
      sheet2->GetCellData()->RemoveArray("EdgeIds");
    }
    
    vtkSmartPointer<vtkIntArray> triangleEdgeType = 
      vtkSmartPointer<vtkIntArray>::New();
    if(TwoSheetEdgeType){
      triangleEdgeType->SetName("EdgeType");
      triangleEdgeType->SetNumberOfTuples(sheet2TriangleNumber);
    }
    else{
      sheet2->GetCellData()->RemoveArray("EdgeType");
    }
    
    vtkSmartPointer<vtkIntArray> triangleTetIds = 
      vtkSmartPointer<vtkIntArray>::New();
    if(TwoSheetTetId){
      triangleTetIds->SetName("TetIds");
      triangleTetIds->SetNumberOfTuples(sheet2TriangleNumber);
    }
    else{
      sheet2->GetCellData()->RemoveArray("TetIds");
    }
    
    vtkSmartPointer<vtkIntArray> triangleCaseIds = 
      vtkSmartPointer<vtkIntArray>::New();
    if(TwoSheetCaseId){
      triangleCaseIds->SetName("CaseIds");
      triangleCaseIds->SetNumberOfTuples(sheet2TriangleNumber);
    }
    else{
      sheet2->GetCellData()->RemoveArray("CaseIds");
    }
    
    for(int i = 0; i < (int) vertexList->size(); i++){
      sheet2->GetPoints()->SetPoint(i, 
        (*vertexList)[i].p_[0], (*vertexList)[i].p_[1], (*vertexList)[i].p_[2]);
      
      if(TwoSheetValue){
        triangleScalarsU->SetTuple1(i, (*vertexList)[i].uv_.first);
        triangleScalarsV->SetTuple1(i, (*vertexList)[i].uv_.second);
      }
      if(TwoSheetParameterization){
        triangleParameterization->SetTuple1(i, (*vertexList)[i].t_);
      }
    }
    if(TwoSheetValue){
      sheet2->GetPointData()->AddArray(triangleScalarsU);
      sheet2->GetPointData()->AddArray(triangleScalarsV);
    }
    if(TwoSheetParameterization){
      sheet2->GetPointData()->AddArray(triangleParameterization);
    }

    int triangleNumber = 0;
    idList->SetNumberOfIds(3);
    for(int i = 0; i < reebSpace_.getNumberOf2sheets(); i++){
      const ReebSpace::Sheet2 *sheet = reebSpace_.get2sheet(i);
     
      if(!sheet->pruned_){
        for(int j = 0; j < (int) sheet->triangleList_.size(); j++){
          
          for(int k = 0; k < (int) sheet->triangleList_[j].size(); k++){
            
            for(int l = 0; l < 3; l++){
              idList->SetId(l, sheet->triangleList_[j][k].vertexIds_[l]);
            }
            
            sheet2Triangles->InsertNextCell(idList);
            
            if(TwoSheetId){
              triangleSheetIds->SetTuple1(triangleNumber, i);
            }
            
            if(TwoSheetEdgeId){
              const ReebSpace::Sheet1 *sheet1 = 
                reebSpace_.get1sheet(sheet->sheet1Id_);
              triangleEdgeIds->SetTuple1(triangleNumber,
                sheet1->edgeList_[j]);
            }
            
            if(TwoSheetEdgeType){
              int polygonEdgeId = sheet->triangleList_[j][k].polygonEdgeId_;
              int edgeId = reebSpace_.getJacobi2Edge(polygonEdgeId);
              triangleEdgeType->SetTuple1(triangleNumber, (*edgeTypes)[edgeId]);
            }
            
            if(TwoSheetTetId){
              triangleTetIds->SetTuple1(triangleNumber, 
                sheet->triangleList_[j][k].tetId_);
            }
            if(TwoSheetCaseId){
              triangleCaseIds->SetTuple1(triangleNumber,
                sheet->triangleList_[j][k].caseId_);
            }
            
            triangleNumber++;
          }
        }
      }
    }
    sheet2->SetCells(VTK_TRIANGLE, sheet2Triangles);
    
    if(TwoSheetId){
      sheet2->GetCellData()->AddArray(triangleSheetIds);
    }
    if(TwoSheetEdgeId){
      sheet2->GetCellData()->AddArray(triangleEdgeIds);
    }
    if(TwoSheetEdgeType){
      sheet2->GetCellData()->AddArray(triangleEdgeType);
    }
    if(TwoSheetTetId){
      sheet2->GetCellData()->AddArray(triangleTetIds);
    }
    if(TwoSheetCaseId){
      sheet2->GetCellData()->AddArray(triangleCaseIds);
    }
  }
 
  // now take care of the 3 sheets
//   vector<float> *triangulationPoints 
//     = reebSpace_.getSheetTriangulationPoints();
//   vector<long long int> *triangulationCells
//     = reebSpace_.getSheetTriangulationCells();
//     
//   vtkSmartPointer<vtkPoints> sheet3Points = vtkSmartPointer<vtkPoints>::New();
//   vtkSmartPointer<vtkFloatArray> pointData = 
//     vtkSmartPointer<vtkFloatArray>::New();
//   pointData->SetNumberOfComponents(3);
//   pointData->SetVoidArray(
//     triangulationPoints->data(), triangulationPoints->size(), 1);
//   sheet3Points->SetData(pointData);
//   sheet3->SetPoints(sheet3Points);
//   
//   vtkSmartPointer<vtkCellArray> sheet3Cells 
//     = vtkSmartPointer<vtkCellArray>::New();
//   vtkSmartPointer<vtkIdTypeArray> idArray 
//     = vtkSmartPointer<vtkIdTypeArray>::New();
//   idArray->SetVoidArray(
//     triangulationCells->data(), triangulationCells->size(), 1);
//   sheet3Cells->SetCells(triangulationCells->size()/5, idArray);
//   sheet3->SetCells(VTK_TETRA, sheet3Cells);

  // now take care of the 3 sheets
  sheet3->ShallowCopy(input);
  const vector<int> *vertex3sheets = reebSpace_.get3sheetVertexSegmentation();
  
  vtkSmartPointer<vtkIntArray> vertexNumberField
    = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> tetNumberField
    = vtkSmartPointer<vtkIntArray>::New();
  
  if(ThreeSheetTetNumber){
    tetNumberField->SetNumberOfTuples(input->GetNumberOfPoints());
    tetNumberField->SetName("3-SheetTetNumber");
    for(int i = 0; i < input->GetNumberOfPoints(); i++){
      const ReebSpace::Sheet3 *sheet3
        = reebSpace_.get3sheet((*vertex3sheets)[i]);
      if((sheet3)&&(!sheet3->pruned_))
        tetNumberField->SetTuple1(i, sheet3->tetList_.size());
      else
        tetNumberField->SetTuple1(i, 0);
    }
    sheet3->GetPointData()->AddArray(tetNumberField);
  }
  else{
    sheet3->GetPointData()->RemoveArray("3-SheetTetNumber");
  }
  
  if(ThreeSheetVertexNumber){
    vertexNumberField->SetNumberOfTuples(input->GetNumberOfPoints());
    vertexNumberField->SetName("3-SheetVertexNumber");
    for(int i = 0; i < input->GetNumberOfPoints(); i++){
      const ReebSpace::Sheet3 *sheet3
        = reebSpace_.get3sheet((*vertex3sheets)[i]);
      if((sheet3)&&(!sheet3->pruned_))
        vertexNumberField->SetTuple1(i, sheet3->vertexList_.size());
      else
        vertexNumberField->SetTuple1(i, 0);
    }
    sheet3->GetPointData()->AddArray(vertexNumberField);
  }
  else{
    sheet3->GetPointData()->RemoveArray("3-SheetTetNumber");
  }
  
  vtkSmartPointer<vtkDoubleArray> domainVolume = 
    vtkSmartPointer<vtkDoubleArray>::New();
  if(ThreeSheetDomainVolume){
    domainVolume->SetNumberOfTuples(input->GetNumberOfPoints());
    domainVolume->SetName("3-SheetDomainVolume");
    
    for(int i = 0; i < input->GetNumberOfPoints(); i++){
      const ReebSpace::Sheet3 *sheet3 = 
        reebSpace_.get3sheet((*vertex3sheets)[i]);
      if((sheet3)&&(!sheet3->pruned_)){
        domainVolume->SetTuple1(i, sheet3->domainVolume_);
      }
      else{
        domainVolume->SetTuple1(i, 0);
      }
    }
    
    sheet3->GetPointData()->AddArray(domainVolume);
  }
  else{
    sheet3->GetPointData()->RemoveArray("3-SheetDomainVolume");
  }
  
  vtkSmartPointer<vtkDoubleArray> rangeArea = 
    vtkSmartPointer<vtkDoubleArray>::New();
  if(ThreeSheetDomainVolume){
    rangeArea->SetNumberOfTuples(input->GetNumberOfPoints());
    rangeArea->SetName("3-SheetRangeArea");
    
    for(int i = 0; i < input->GetNumberOfPoints(); i++){
      const ReebSpace::Sheet3 *sheet3 = 
        reebSpace_.get3sheet((*vertex3sheets)[i]);
      if((sheet3)&&(!sheet3->pruned_)){
        rangeArea->SetTuple1(i, sheet3->rangeArea_);
      }
      else{
        rangeArea->SetTuple1(i, 0);
      }
    }
    
    sheet3->GetPointData()->AddArray(rangeArea);
  }
  else{
    sheet3->GetPointData()->RemoveArray("3-SheetRangeArea");
  }
  
  vtkSmartPointer<vtkDoubleArray> hyperVolume = 
    vtkSmartPointer<vtkDoubleArray>::New();
  if(ThreeSheetDomainVolume){
    hyperVolume->SetNumberOfTuples(input->GetNumberOfPoints());
    hyperVolume->SetName("3-SheetHyperVolume");
    
    for(int i = 0; i < input->GetNumberOfPoints(); i++){
      const ReebSpace::Sheet3 *sheet3 = 
        reebSpace_.get3sheet((*vertex3sheets)[i]);
      if((sheet3)&&(!sheet3->pruned_)){
        hyperVolume->SetTuple1(i, sheet3->hyperVolume_);
      }
      else{
        hyperVolume->SetTuple1(i, 0);
      }
    }
    
    sheet3->GetPointData()->AddArray(hyperVolume);
  }
  else{
    sheet3->GetPointData()->RemoveArray("3-SheetHyperVolume");
  }
  
  vtkSmartPointer<vtkIntArray> vertexSegmentation 
    = vtkSmartPointer<vtkIntArray>::New();
  vertexSegmentation->SetName("3-SheetId");
  vertexSegmentation->SetNumberOfTuples(input->GetNumberOfPoints());
  for(int i = 0; i < input->GetNumberOfPoints(); i++){
    const ReebSpace::Sheet3 *sheet = reebSpace_.get3sheet((*vertex3sheets)[i]);
    if(sheet){
      vertexSegmentation->SetTuple1(i, sheet->simplificationId_);
    }
    else{
      vertexSegmentation->SetTuple1(i, (*vertex3sheets)[i]);
    }
  }
  sheet3->GetPointData()->AddArray(vertexSegmentation);
  
  const vector<int> *tet3sheets = reebSpace_.get3sheetTetSegmentation();
  vtkSmartPointer<vtkIntArray> tetSegmentation 
    = vtkSmartPointer<vtkIntArray>::New();
  tetSegmentation->SetName("3-SheetId");
  tetSegmentation->SetNumberOfTuples(input->GetNumberOfCells());
  for(int i = 0; i < input->GetNumberOfCells(); i++){
    const ReebSpace::Sheet3 *sheet = reebSpace_.get3sheet((*tet3sheets)[i]);
    if(sheet){
      tetSegmentation->SetTuple1(i, sheet->simplificationId_);
    }
    else{
      tetSegmentation->SetTuple1(i, (*tet3sheets)[i]);
    }
  }
  sheet3->GetCellData()->AddArray(tetSegmentation);
  
  return 0;
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int vtkReebSpace::RequestData(vtkInformation *request, 
  vtkInformationVector **inputVector, vtkInformationVector *outputVector){

  Memory m;
  
  // here the vtkDataSet type should be changed to whatever type you consider.
  vtkUnstructuredGrid *input = vtkUnstructuredGrid::GetData(inputVector[0]);
  vtkUnstructuredGrid *sheet0 = 
    vtkUnstructuredGrid::SafeDownCast(
      outputVector->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT()));
  vtkUnstructuredGrid *sheet1 = 
    vtkUnstructuredGrid::SafeDownCast(
      outputVector->GetInformationObject(1)->Get(vtkDataObject::DATA_OBJECT()));
  vtkUnstructuredGrid *sheet2 = 
    vtkUnstructuredGrid::SafeDownCast(
      outputVector->GetInformationObject(2)->Get(vtkDataObject::DATA_OBJECT()));
  vtkUnstructuredGrid *sheet3 = 
    vtkUnstructuredGrid::SafeDownCast(
      outputVector->GetInformationObject(3)->Get(vtkDataObject::DATA_OBJECT()));
    
  doIt(input, sheet0, sheet1, sheet2, sheet3);
  
  {
    stringstream msg;
    msg << "[vtkReebSpace] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 1;
}