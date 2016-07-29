/// \ingroup baseCode
/// \class wtfit::JacobiSet 
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date June 2015.
///
/// \brief %JacobiSet processing package.
///
/// %JacobiSet is a processing package that computes the Jacobi Set of a 
/// bivariate scalar function.
/// \param dataTypeU Data type of the input first component field (char, float, 
/// etc.).
/// \param dataTypeV Data type of the input second component field (char, float,
/// etc.)
/// \sa vtkJacobiSet

#ifndef _JACOBISET_H
#define _JACOBISET_H

// base code includes
#include                  <ZeroSkeleton.h>
#include                  <OneSkeleton.h>
#include                  <ScalarFieldCriticalPoints.h>
#include                  <Wrapper.h>


namespace wtfit{
  
  template <class dataTypeU, class dataTypeV> class JacobiSet : public Debug{

    public:
        
      JacobiSet();
      
      ~JacobiSet();

      int connectivityPreprocessing(const vector<vector<int> > &edgeStarList,
        vector<vector<pair<int, int> > > &edgeFanLinkEdgeLists,
        vector<vector<long long int> > &edgeFans,
        vector<int> &sosOffsets) const;
      
      int execute(vector<pair<int, char> > &jacobiSet);
    
      int perturbate(const dataTypeU &uEpsilon = pow(10, -DBL_DIG),
        const dataTypeV &vEpsilon = pow(10, -DBL_DIG)) const;
      
      int setEdgeFans(const vector<vector<long long int> > *edgeFans){
        edgeFans_ = edgeFans;
        return 0;
      }
      
      int setEdgeFanLinkEdgeList(
        const vector<vector<pair<int, int> > > *edgeFanLinkEdgeLists){
        edgeFanLinkEdgeLists_ = edgeFanLinkEdgeLists;
        return 0;
      }
      
      int setEdgeList(const vector<pair<int, int> > *edgeList){
        edgeList_ = edgeList;
        return 0;
      }
      
      int setInputField(const void *uField, const void *vField){
        
        uField_ = uField;
        vField_ = vField;
        return 0;
      }
      
      int setSosOffsets(vector<int> *sosOffsets){
        sosOffsets_ = sosOffsets;
        return 0;
      }

      // NOTE: here it's not clear how vtk builds vtkIdType 
      // to check on bigger data-sets
      int setTetList(const long long int *tetList){
        tetList_ = tetList;
        return 0;
      }
      
      int setVertexNumber(const int &vertexNumber){
        vertexNumber_ = vertexNumber;
        return 0;
      }
      
    protected:
    
      int                   vertexNumber_;
      const long long int   *tetList_;
      const void            *uField_, *vField_;
      const vector<pair<int, int> > *edgeList_;
      // for each edge, one skeleton of its triangle fan
      const vector<vector<pair<int, int> > > *edgeFanLinkEdgeLists_;
      // for each edge, the one skeleton of its triangle fan
      const vector<vector<long long int> > *edgeFans_;
      vector<int>           *sosOffsets_;
  };
}

// if the package is not a template, comment the following line
#include                  <JacobiSet.cpp>

#endif // JACOBISET_H
