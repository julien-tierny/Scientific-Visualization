/*
 * file:                  Editor.cpp
 * description:           Data-structures and processing.
 * author:                Your Name Here <Your Email Address Here>.
 * date:                  The Date Here.
 */

#include                  <Editor.h>

#include <vtksys/SystemTools.hxx>

Editor::Editor(){

  grid_ = NULL;
  lastObject_ = true;

  debug_ = 1;
  contourTree_ = vtkContourForests::New();
}

Editor::~Editor(){

  // delete the mesh reader and the associated data at once
  contourTree_->Delete();
}

int Editor::execute(){

  contourTree_->setDebugLevel(debug_);
  contourTree_->SetdebugLevel_(debug_);
  contourTree_->SetThreadNumber(core_);

  contourTree_->SetInputData(grid_);
  contourTree_->SetFieldId(fieldId_);
  contourTree_->ShowSegmentation(true);
  contourTree_->SetLessPartition(lessPartitions_);
  contourTree_->SetSimplificationType(method_);
  contourTree_->SetSimplificationThreshold(threshold_);
  contourTree_->SetTreeType(treeType_);
  contourTree_->Update();

  grid_->ShallowCopy(contourTree_->GetOutput());

  return 0;
}

int Editor::init(int &argc, char **argv){

  CommandLineParser parser;

  // specify argument "-g" : The path of the grid
  parser.setStringArgument("g", &inputFilePath_, "Path to the input 3D grid");
  //parser.setIntArgument("d", &debug_, "Debug Level", true);
  parser.setIntArgument("n", &core_, "Thread number", true);
  parser.setIntArgument("f", &fieldId_, "Field identifier", true);
  parser.setDoubleArgument("s", &threshold_, "Simplification threshold between 0 and 1", true);
  parser.setIntArgument("m", &method_, "Simplification method : 0 persist, 1 Vertices ...", true);
  parser.setIntArgument("t", &treeType_, "type of tree : 2 is CT", true);
  parser.setOption("l", &lessPartitions_, "Use 2 time less partitions than nb threads");

  // now parse the command line
  parser.parse(argc, argv);

  // set default values
  debug_ = wtfit::globalDebugLevel_;
  if(debug_ == -INT_MAX){
    debug_ = 1;
  }

  if(fieldId_ == -INT_MAX){
    fieldId_ = 0;
  }

  if(core_ == -INT_MAX){
    core_ = 1;
  }

  if(threshold_ == -DBL_MAX){
    threshold_ = 0;
  }

  if (treeType_ == -INT_MAX) {
     treeType_ = 2;
  }

  lessPartitions_ &= core_ != 1;

  // now load the data to the editor
  loadData();

  return 0;
}

int Editor::loadData(){

  // create a reader object
  std::string extension =
          vtksys::SystemTools::GetFilenameLastExtension(inputFilePath_);

  if (extension == ".vtu"){
    grid_ = ReadAnXMLFile<vtkXMLUnstructuredGridReader> (inputFilePath_.c_str());
  } else if (extension == ".vti"){
    grid_ = ReadAnXMLFile<vtkXMLImageDataReader> (inputFilePath_.c_str());
  } else {
    cerr << "Bad format, need vtu" << endl;
    return -1;
  }

  {
    stringstream msg;
    msg << "[Editor]   done! (read "
      << grid_->GetNumberOfPoints()
      << " vertices, "
      << grid_->GetNumberOfCells()
      << " cells)" << endl;
    dMsg(cout, msg.str(), 1);
  }

  return 0;
}

int Editor::saveData(const string &fileName) const{

  vtkXMLUnstructuredGridWriter *writer = vtkXMLUnstructuredGridWriter::New();

  writer->SetFileName(fileName.data());
  writer->SetInputData(
          (vtkUnstructuredGrid*)contourTree_->GetOutput());
  writer->Write();

  writer->Delete();

  return 0;
}
