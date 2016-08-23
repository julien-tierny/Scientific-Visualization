#include                  <vtkTriangulation.h>

vtkTriangulation::vtkTriangulation(){

}

int vtkTriangulation::setInputData(vtkDataSet* dataSet){

  if(dataSet->GetDataObjectType() == VTK_UNSTRUCTURED_GRID){
    setInputPoints(dataSet->GetNumberOfPoints(),
      (float *)
      ((vtkUnstructuredGrid *) dataSet)->GetPoints()->GetVoidPointer(0));
    setInputCells(dataSet->GetNumberOfCells(),
      ((vtkUnstructuredGrid *) dataSet)->GetCells()->GetPointer());
  }
  else if(dataSet->GetDataObjectType() == VTK_POLY_DATA){
    setInputPoints(dataSet->GetNumberOfPoints(),
      (float *)
      ((vtkPolyData *) dataSet)->GetPoints()->GetVoidPointer(0));
    setInputCells(dataSet->GetNumberOfCells(),
      ((vtkPolyData *) dataSet)->GetPolys()->GetPointer());
  }
  else if(dataSet->GetDataObjectType() == VTK_IMAGE_DATA){
    vtkImageData *imageData = (vtkImageData *) dataSet;

    double origin[3];
    imageData->GetOrigin(origin);
    
    double spacing[3];
    imageData->GetSpacing(spacing);
    
    int dimensions[3];
    imageData->GetDimensions(dimensions);
    
    setInputGrid(
      origin[0], origin[1], origin[2],
      spacing[0], spacing[1], spacing[2],
      dimensions[0], dimensions[1], dimensions[2]);
  }
  else{
    stringstream msg;
    msg << "[vtkTriangulation] Unsupported input VTK class (ref=" 
      << dataSet->GetDataObjectType() << ")" << endl;
    msg << "[vtkTriangulation] Leaving an empty triangulation..." << endl;
    dMsg(cerr, msg.str(), Debug::fatalMsg);
  }
  
  return 0;
}

