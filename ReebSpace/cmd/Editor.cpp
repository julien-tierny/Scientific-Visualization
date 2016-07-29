#include                  <Editor.h>

Editor::Editor(){

  input_ = NULL;
  lastObject_ = true;
  debugLevel_ = 0;

  reebSpace_ = vtkReebSpace::New();
}

Editor::~Editor(){

  // delete the mesh reader and the associated data at once
  if(reader_) 
    reader_->Delete();
  if(reebSpace_)
    reebSpace_->Delete();
}

int Editor::execute(){
  
  reebSpace_->Update();
  
  return 0;
}

int Editor::init(int &argc, char **argv){

  double threshold;
  int criterion;
  string inputFilePath;
  CommandLineParser parser;
 
  parser.setStringArgument("i", &inputFilePath,
    "Path to the input VTU file");
  
  parser.setDoubleArgument("s", &threshold,
    "Simplification threshold, in [0, 1] (default: 0)", true);
  parser.setIntArgument("c", &criterion,
    "Simplification criterion, in [0, 2] (default: 1)", true);
  
  parser.setIntArgument("u", &uId_, 
    "Identifier of the u-component field (default: 0)", true);
  parser.setIntArgument("v", &vId_, 
    "Identifier of the v-component field (default: 1)", true);
  {
    stringstream msg;
    msg << "Thread number (default: " << OsCall::getNumberOfCores() << ")";
    parser.setIntArgument("t", &threadNumber_, msg.str(), true);
  }
  
  // now parse the command line
  parser.parse(argc, argv);
   
  // set default values
  if(globalDebugLevel_ == -INT_MAX){
    globalDebugLevel_ = infoMsg;
  }
  debugLevel_ = globalDebugLevel_;
  if(threadNumber_ == -INT_MAX){
    threadNumber_ = OsCall::getNumberOfCores();
  }
  if(uId_ == -INT_MAX){
    uId_ = 0;
  }
  if(vId_ == -INT_MAX){
    vId_ = 1;
  }
  if(criterion == -INT_MAX){
    criterion = 1;
  }
  if(threshold == -DBL_MAX){
    threshold = 0;
  }
  
  // now load the data to the editor
  int ret = loadData(inputFilePath);
  if(ret != 0)
    return ret;
  
  reebSpace_->SetThreadNumber(threadNumber_);
  reebSpace_->SetdebugLevel_(debugLevel_);
  reebSpace_->SetSimplificationCriterion(criterion);
  reebSpace_->SetSimplificationThreshold(threshold);
  
  return 0;
}

int Editor::loadData(const string &fileName){

  // create a reader object
  reader_ = vtkXMLUnstructuredGridReader::New();
  reader_->SetFileName(fileName.data());  
  
  // handle debug messages
  {
    stringstream msg;
    msg << "[Editor] Reading input mesh..." << endl;
    // choose where to display this message (cout, cerr, a file)
    // choose the priority of this message (1, nearly always displayed, 
    // higher values mean lower priorities)
    dMsg(cout, msg.str(), 1);
  }
  
  reader_->Update();
  input_ = reader_->GetOutput();

  // set the input
  reebSpace_->SetInputData(input_);
  if(!input_->GetPointData()){
    stringstream msg;
    msg << "[Editor] No scalar field attached to the input tet mesh!" << endl;
    msg << "[Editor] Exiting..." << endl;
    dMsg(cerr, msg.str(), 0);
    return -1;
  }
  
  vtkDataArray *component = input_->GetPointData()->GetArray(uId_);
  if(!component){
    stringstream msg;
    msg << "[Editor] Could not find scalar field with Id `" << uId_ << "'!"
      << endl;
    msg << "[Editor] Exiting..." << endl;
    dMsg(cerr, msg.str(), 0);
    return -2;
  }
  reebSpace_->SetUcomponent(component->GetName());
  
  component = input_->GetPointData()->GetArray(vId_);
  if(!component){
    stringstream msg;
    msg << "[Editor] Could not find scalar field with Id `" << vId_ << "'!"
      << endl;
    msg << "[Editor] Exiting..." << endl;
    dMsg(cerr, msg.str(), 0);
    return -3;
  }
  reebSpace_->SetVcomponent(component->GetName());
  
  {
    stringstream msg;
    msg << "[Editor]   done! (read " 
      << input_->GetNumberOfPoints()
      << " vertices, "
      << input_->GetNumberOfCells() 
      << " cells)" << endl;
    dMsg(cout, msg.str(), 1);
  }

  return 0;
}

int Editor::saveData(const string &fileName) const{
  
  for(int i = 0; i < 4; i++){
    vtkXMLUnstructuredGridWriter *writer = vtkXMLUnstructuredGridWriter::New();
    stringstream filePath;
    filePath << fileName << "_" << i << ".vtu";
    writer->SetFileName(filePath.str().data());
    writer->SetInputData(reebSpace_->GetOutput(i));
    writer->Write();
    writer->Delete();
  }
  
  return 0;
}
