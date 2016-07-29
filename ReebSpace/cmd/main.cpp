/// \date February 2016.
///
/// \brief Command-line program that computes the Reeb space segmentation of 
/// a piecewise linear bivariate scalar field defined on a piecewise linear 
/// 3-manifold.

// include the local headers
#include                  <Editor.h>

int main(int argc, char **argv) {

  // init editor
  Editor editor;
  int ret = editor.init(argc, argv);
 
  if(ret != 0)
    return ret;
  
  // execute data processing
  ret = editor.execute();
  
  if(ret != 0)
    return ret;
  
  // save the output
  editor.saveData("output");
  
  return 0;
}
