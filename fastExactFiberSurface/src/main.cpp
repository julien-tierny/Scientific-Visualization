/* 
 * file:                  main.cpp
 * description:           Command-line program to compute fiber surfaces.
 * author(s):             Pavol Klacansky <pavol@klacansky.com>.
 *                        Julien Tierny <julien.tierny@lip6.fr>
 * date:                  March 2015.
 */

#include                  <CommandLineParser.h>

#include                  <vtkFiberSurfaceFilter.h>
#include                  <vtkSmartPointer.h>
#include                  <vtkXMLPolyDataWriter.h>
#include                  <vtkXMLUnstructuredGridReader.h>

int main(int argc, char **argv){
  
  // 1) parse command line options and arguments
  CommandLineParser       parser;
  
  bool                    manifoldOutput, fiberTexture, edgeSegmentation;
  int                     algorithm,
                          dataF0, dataF1,
                          polygonF0, polygonF1, 
                          implementation,
                          threadNumber,
                          threadingStrategy,
                          bvhMinimumTetNumber;
  double                  octreeMinimalRangeAreaRatio,
                          automaticThreadingStrategyBalance;
  string                  dataSetPath, polygonPath;
  
  parser.setStringArgument("i", &dataSetPath, "Path to VTU input data-set");
  parser.setStringArgument("p", &polygonPath, "Path to VTU polygon");
  
  parser.setIntArgument("iF0", &dataF0, 
    "Scalar field Id for first component (data-set, default: 0)", true);
  parser.setIntArgument("iF1", &dataF1, 
    "Scalar field Id for second component (data-set, default: 1)", true);
  parser.setIntArgument("pF0", &polygonF0, 
    "Scalar field Id for first component (polygon, default: 0)", true);
  parser.setIntArgument("pF1", &polygonF1, 
    "Scalar field Id for second component (polygon, default: 1)", true);
  
  parser.setIntArgument("I", &implementation, 
    "Implementation: [0, Regular] [1, Octree] [2, BVH] (default: 0)", true);
  
  parser.setIntArgument("t", &threadNumber,
    "Thread number (default: 1)", true);
  
  {
    stringstream msg;
    msg << "Threading strategy";
    msg << " [0, Automatic] [1, Domain tetrahedra] [2, Polygon Edges]";
    msg << " (default: 1)";
    parser.setIntArgument("s", &threadingStrategy, msg.str(), true);
  }
  
  parser.setIntArgument("a", &algorithm,
  "Extraction algorithm: [0, normal] [1, grey] [2, grey minimal] (default: 2)",
  true);
  
  parser.setRealArgument("o", &octreeMinimalRangeAreaRatio,
    "Octree minimum range area ratio (default: 0)", true);
  
  parser.setIntArgument("b", &bvhMinimumTetNumber,
    "BVH minimum tet number (default: 1)", true);
  
  parser.setRealArgument("B", &automaticThreadingStrategyBalance,
    "Automatic threading strategy balance (default: 1)", true);
  
  parser.setOption("m", &manifoldOutput,
    "Manifold output (default: off)");
  
  parser.setOption("f", &fiberTexture,
    "Compute fiber texture coordinates (default: off)");
  
  parser.setOption("e", &edgeSegmentation,
    "Compute polygon edge segmentation (default: off)");
  
  parser.parse(argc, argv);
  
  // default values
  if(dataF0 == -INT_MAX){
    dataF0 = 0;
  }
  if(dataF1 == -INT_MAX){
    dataF1 = 1;
  }
  
  if(polygonF0 == -INT_MAX){
    polygonF0 = 0;
  }
  if(polygonF1 == -INT_MAX){
    polygonF1 = 1;
  }
  
  if(implementation == -INT_MAX){
    implementation = 0;
  }
  
  if(algorithm == -INT_MAX){
    algorithm = 2;
  }
 
  if(threadNumber <= 0){
    threadNumber = 1;
  }

  if(threadingStrategy == -INT_MAX){
    threadingStrategy = 1;
  }
  
  if(octreeMinimalRangeAreaRatio == -REAL_MAX){
    octreeMinimalRangeAreaRatio = 0;
  }
  
  if(bvhMinimumTetNumber < 1){
    bvhMinimumTetNumber = 1;
  }
  if(bvhMinimumTetNumber > 16){
    bvhMinimumTetNumber = 16;
  }
  
  if(automaticThreadingStrategyBalance < 0){
    automaticThreadingStrategyBalance = 1;
  }
  
  
  // 2) reading
  cout << "[fiberSurface] Reading input data-set..." << endl;
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader0 = 
    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader0->SetFileName(dataSetPath.data());
  reader0->Update();
  cout << "[fiberSurface]   "
    << reader0->GetOutput()->GetNumberOfCells() 
    << " tet(s) read." << endl; 
    
  cout << "[fiberSurface] Reading input polygon..." << endl;
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader1 = 
    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader1->SetFileName(polygonPath.data());
  reader1->Update();
  cout << "[fiberSurface]   "
    << reader1->GetOutput()->GetNumberOfCells() 
    << " segment(s) read." << endl; 
  
  vtkSmartPointer<vtkFiberSurfaceFilter> fiberSurface = 
    vtkSmartPointer<vtkFiberSurfaceFilter>::New();
    
  // set the input
  fiberSurface->SetInputData(0, reader0->GetOutput());
  fiberSurface->SetInputData(1, reader1->GetOutput());
  
  
  // 3) setting options
  
  // field selection
  fiberSurface->SetDataUComponent(
    reader0->GetOutput()->GetPointData()->GetArray(dataF0)->GetName());
  fiberSurface->SetDataVComponent(
    reader0->GetOutput()->GetPointData()->GetArray(dataF1)->GetName());
  cout << "[fiberSurface] Data fields: '"
    << reader0->GetOutput()->GetPointData()->GetArray(dataF0)->GetName()
    << "' and '"
    << reader0->GetOutput()->GetPointData()->GetArray(dataF1)->GetName()
    << "'" << endl;
  
  fiberSurface->SetPolygonUComponent(
    reader1->GetOutput()->GetPointData()->GetArray(polygonF0)->GetName());
  fiberSurface->SetPolygonVComponent(
    reader1->GetOutput()->GetPointData()->GetArray(polygonF1)->GetName());
  cout << "[fiberSurface] Polygon fields: '"
    << reader1->GetOutput()->GetPointData()->GetArray(polygonF0)->GetName()
    << "' and '"
    << reader1->GetOutput()->GetPointData()->GetArray(polygonF1)->GetName()
    << "'" << endl;
  
  fiberSurface->SetImplementation(implementation);
  cout << "[fiberSurface] Implementation: " << implementation << endl;
  fiberSurface->SetAlgorithm(algorithm);
  cout << "[fiberSurface] Algorithm: " << algorithm << endl;
  fiberSurface->SetThreadNumber(threadNumber);
  cout << "[fiberSurface] Thread number: " << threadNumber << endl;
  fiberSurface->SetThreadStrategy(threadingStrategy);
  cout << "[fiberSurface] Threading strategy: " << threadingStrategy << endl;
  fiberSurface->SetOctreeRangeRatio(octreeMinimalRangeAreaRatio);
  cout << "[fiberSurface] Octree minimum range area ratio: "
    << octreeMinimalRangeAreaRatio << endl;
  fiberSurface->SetBvhMinimumTetNumber(bvhMinimumTetNumber);
  cout << "[fiberSurface] BVH minimum tet number: "
    << bvhMinimumTetNumber << endl;
  fiberSurface->SetThreadBalance(automaticThreadingStrategyBalance);
  cout << "[fiberSurface] Automatic threading strategy balance: "
    << automaticThreadingStrategyBalance << endl;
  fiberSurface->SetManifold(manifoldOutput);
  cout << "[fiberSurface] Manifold output: " << manifoldOutput << endl;
  fiberSurface->SetVisibleFibers(fiberTexture);
  cout << "[fiberSurface] Fiber texture coordinates: " << fiberTexture << endl;
  fiberSurface->SetPolygonSegmentation(edgeSegmentation);
  cout << "[fiberSurface] Polygon edge segmentation: " 
    << edgeSegmentation << endl;
  cout << "[fiberSurface] Let's go!" << endl;
  cout << "[fiberSurface]" << endl;
  
  
  // 4) Computing
  fiberSurface->Update();
  
  
  // 5) Saving
  cout << "[fiberSurface] Saving output to output.vtp..." << endl;
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = 
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetInputData(fiberSurface->GetOutput());
  writer->SetFileName("output.vtp");
  writer->Write();
  cout << "[fiberSurface]   done!" << endl;
  
  return 0;
}
