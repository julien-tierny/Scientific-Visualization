/// \date February 2016.
///
/// \brief Data-structures and processing.

#ifndef EDITOR_H
#define EDITOR_H

// base code includes
#include                  <CommandLineParser.h>

// vtk wrappers
#include                  <vtkReebSpace.h>

// VTK includes
#include                  <vtkUnstructuredGrid.h>
#include                  <vtkXMLUnstructuredGridReader.h>
#include                  <vtkXMLUnstructuredGridWriter.h>

class Editor : public Debug{

  public:
      
    Editor();
      
    ~Editor();
    
    int execute();
    
    vtkDataSet* getData() const{ return input_;};

    int init(int &argc, char **argv);
    
    int saveData(const string &fileName) const;
   
    
  protected:

    int                           uId_, vId_;
    vtkUnstructuredGrid           *input_;
    vtkReebSpace                  *reebSpace_;
    vtkXMLUnstructuredGridReader  *reader_;
    
    int loadData(const string &fileName);
    
};

#endif // EDITOR_H
