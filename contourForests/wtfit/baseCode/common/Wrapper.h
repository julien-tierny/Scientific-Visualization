/// \ingroup baseCode
/// \class wtfit::Wrapper
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date October 2014.
/// 
/// \brief Wrapper class to wrap wtfit code.

#ifndef                 _WRAPPER_H
#define                 _WRAPPER_H

#include                <Debug.h>

namespace wtfit{

  class Wrapper : public Debug{
    
    public:
      
      Wrapper(){
        processingProgress_ = 0;
      };
      
      ~Wrapper(){};
     
      virtual bool needsToAbort() = 0;
     
      virtual int updateProgress(const float &progress) = 0;
      
    protected:
      
      float             processingProgress_;
      
  };
}

#endif
