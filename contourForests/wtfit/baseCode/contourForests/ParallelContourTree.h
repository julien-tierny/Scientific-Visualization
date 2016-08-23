/*
 * file:                  ContourTree.h
 * description:           ContourTree processing package.
 * author:                Gueunet Charles
 * date:                  aout 2015
 */

///\ingroup baseCode
///\class wtfit::ParallelContourTree
///\brief %ParallelContourTree processing package.
///
///%ParallelContourTree compute the contour tree for a given input
/// scalar in parallel
///\param dataType Data type of the input scalar field (char, float,
/// etc.).
///\sa vtkContourTree

#ifndef _PARALLELCONTOURTREE_H
#define _PARALLELCONTOURTREE_H

#include <ContourTree.h>

#include <typeinfo>

namespace wtfit
{
   // Classes Interface
   //         Partition
   //         ParallelContourTree

   class Interface
   {
     private:
      // same size : number of conn. components.
      idVertex seed_;

      // Overlap
      vector<idVertex> lowerOverlap_;
      vector<idVertex> upperOverlap_;

     public:
      Interface(const idVertex &seed);

      // Getter & Setter
      // {
      inline const idVertex &getSeed(void) const
      {
         return seed_;
      }

      inline vector<idVertex> &getUpper(void)
      {
         return upperOverlap_;
      }

      inline vector<idVertex> &getLower(void)
      {
         return lowerOverlap_;
      }
      inline const idVertex getNbUpper(void) const
      {
         return upperOverlap_.size();
      }

      inline const idVertex getNbLower(void) const
      {
         return lowerOverlap_.size();
      }

      inline void setSeed(const idVertex &local_seed)
      {
         seed_ = local_seed;
      }

      inline void addUpper(const idVertex &lb)
      {
         upperOverlap_.emplace_back(lb);
      }

      inline void addLower(const idVertex &lb)
      {
         lowerOverlap_.emplace_back(lb);
      }

      inline void upReserve(const idVertex &u)
      {
          upperOverlap_.reserve(u);
      }

      inline void loReserve(const idVertex &l)
      {
        lowerOverlap_.reserve(l);
      }

      inline void appendUpper(const vector<idVertex> &vertices)
      {
         upperOverlap_.insert(upperOverlap_.end(), vertices.cbegin(), vertices.cend());
      }

      inline void appendLower(const vector<idVertex> &vertices)
      {
         lowerOverlap_.insert(lowerOverlap_.end(), vertices.cbegin(), vertices.cend());
      }

      // }
   };

   class ParallelContourTree : public ContourTree
   {
     private:
      // Parallel info
      numThread   nbThread_;
      idInterface nbInterfaces_;
      idPartition nbPartitions_;
      int partitionNum_;
      bool lessPartition_;

      vector<Interface>   vect_interfaces_;
      vector<ContourTree> vect_ct_;

      // vector<vector<Node> *> vect_corrNodes_;
      // vector<vector<SuperArc> *> vect_corrSuperArc_;

     public:
      ParallelContourTree(const numThread &nbThread);

      ParallelContourTree();

      virtual ~ParallelContourTree();

      // Getters & Setters
      // {
      inline void setNbThread(const unsigned short nbThread)
      {
         if (nbThread) {
            nbThread_ = nbThread;
         } else {
            nbThread_ = OsCall::getNumberOfCores();
         }

#ifndef withKamikaze
         if (nbThread_ == 0) {
            nbThread_     = 1;
            nbInterfaces_ = 0;
            nbPartitions_ = 1;
            stringstream debug;
            debug << "[ParallelContourTree] won't use 0 thread, set to 1" << endl;
            err(debug.str(), fatalMsg);
         }
#endif

         nbPartitions_ = nbThread_;
         nbInterfaces_ = nbPartitions_ - 1;

         vect_interfaces_.reserve(nbInterfaces_);
         vect_ct_.resize(nbPartitions_);
      }

      inline void setPartitionNum(int p)
      {
          partitionNum_ = p;
      }

      inline void setLessPartition(bool l)
      {
          lessPartition_ = l;
      }

      // }

      // Init
      // {
      void initInterfaces();

      void initOverlap();

      template <typename scalarType>
      void initLocalCT(decltype(nbPartitions_) i);

      void flush(void);

      //}

      // Process
      // {

      template <typename scalarType>
      int build(const bool ct, const bool segment, const double threshold);

      template <typename scalarType>
      int parallelBuild(vector<vector<ExtendedUnionFind *>> &baseUF_JT,
                        vector<vector<ExtendedUnionFind *>> &baseUF_ST, const double threshold);

      void stitch(void);
      void stitchTree(const char tree);

      // replace distributed tree by a global one, will be removed
      void unify();
      void unifyCT();
      void unifyMT();
      // }

      // Print
      // {
      void printDebug(DebugTimer &timer, const string &str);

      void printVectCT();
      // }
   };

#include <ParallelContourTreeTemplate.h>
}

#endif
