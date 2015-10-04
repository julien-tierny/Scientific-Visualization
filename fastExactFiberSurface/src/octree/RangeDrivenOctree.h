/* 
 * file:                  RangeDrivenOctree.h
 * description:           RangeDrivenOctree processing package.
 * author:                Julien Tierny <julien.tierny@lip6.fr>.
 * date:                  March 2015.
 */

/// \ingroup baseCode
/// \class wtfit::RangeDrivenOctree 
/// \brief %RangeDrivenOctree processing package.
///
/// \param dataType Data type of the input scalar field (char, float, 
/// etc.).

#ifndef _RANGE_DRIVEN_OCTREE_H
#define _RANGE_DRIVEN_OCTREE_H

// base code includes
#include                  <Wrapper.h>


namespace wtfit{
  
  template <class dataType0, class dataType1> class RangeDrivenOctree 
    : public Debug{

    public:
        
      RangeDrivenOctree();
      
      ~RangeDrivenOctree();

      int build();
      
      inline int getTet2NodeMap(vector<int> &map, 
        const bool &forSegmentation = false) const;
      
      inline int rangeSegmentQuery(const double *p0, const double *p1, 
        vector<int> &cellList) const;
      
      inline void setCellList(const long long int *cellList){
        cellList_ = cellList;
      }
      
      inline void setCellNumber(const int &cellNumber){
        cellNumber_ = cellNumber;
      }

      inline void setLeafMinimumCellNumber(const int &number){
        leafMinimumCellNumber_ = number;
      }
     
      inline void setLeafMinimumDomainVolumeRatio(const float &ratio){
        leafMinimumRangeAreaRatio_ = ratio;
      }
     
      inline void setLeafMinimumRangeAreaRatio(const float &ratio){
        leafMinimumRangeAreaRatio_ = ratio;
      }
     
      inline void setPointList(const float *pointList){
        pointList_ = pointList;
      }
      
      inline void setRange(const dataType0 *u, const dataType1 *v){
        u_ = u;
        v_ = v;
      }
      
      inline void setVertexNumber(const int &vertexNumber){
        vertexNumber_ = vertexNumber;
      }
      
      int stats(ostream &stream);
      
      int statNode(const int &nodeId, ostream &stream);
      
    protected:
   
      class OctreeNode{
       
        public:
          pair<pair<dataType0, dataType0>, pair<dataType1, dataType1> >
                          rangeBox_;
          vector<int>     cellList_;
          vector<int>     childList_;
          vector<pair<float, float> > 
                          domainBox_;
      };
      
      int buildNode(const vector<int> &cellList, 
        const vector<pair<float, float> > &domainBox,
        const pair<pair<dataType0, dataType0>, pair<dataType1, dataType1> >
          &rangeBox, int &nodeId);
      
      inline int rangeSegmentQuery(const double *p0, const double *p1,
        const int &nodeId,
        vector<int> &cellList) const;
      
      inline bool segmentIntersection(const double *p0, const double *p1,
        const double *q0, const double *q1) const;
      
      const dataType0     *u_;
      const dataType1     *v_;
      const float         *pointList_;
      const long long int *cellList_;
      float               domainVolume_,
                          leafMinimumDomainVolumeRatio_, 
                          leafMinimumRangeAreaRatio_, 
                          rangeArea_;
      int                 cellNumber_, vertexNumber_, 
                          leafMinimumCellNumber_, rootId_;
      mutable int         queryResultNumber_;
      vector<OctreeNode>  nodeList_;
      vector<vector<pair<float, float> > > 
                          cellDomainBox_;
      vector<pair<pair<dataType0, dataType0>, pair<dataType1, dataType1> > >
                          cellRangeBox_;
  };
}

// if the package is not a template, comment the following line
#include                  <RangeDrivenOctree.cpp>

#endif // _RANGE_DRIVEN_OCTREE_H
