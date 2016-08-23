/**
 * @file ParallelContourTreeTemplate.h
 * @brief Template function for parallel contour tree processing
 * @author Gueunet Charles
 * @version 1
 * @date 2016-03-30
 */

#ifndef PARALLELCONTOURTREETEMPLATE_H
#define PARALLELCONTOURTREETEMPLATE_H

#include <ParallelContourTree.h>

// TODO
// Template wrapper

/// ------------------- ParallelCT

// Init
// {

template <typename scalarType>
void ParallelContourTree::initLocalCT(decltype(nbPartitions_) i)
{
   vect_ct_[i].setDebugLevel(debugLevel_);

   vect_ct_[i].setTriangulation(mesh_);
   vect_ct_[i].setPartition(i);

   vect_ct_[i].setVertexScalars<scalarType>((scalarType *)scalars_);
   vect_ct_[i].setVertexSoSoffsets(soSOffsets_);
   vect_ct_[i].setSorted(sortedVertices_);
   vect_ct_[i].setMirror(mirrorOffsets_);

   vect_ct_[i].setSegmentation(segmentation_);
   vect_ct_[i].setComputeContourTree(computeContourTree_);
   vect_ct_[i].setSimplificationMethod(simplifyMethod_);

   // take care of vect2tree if not set before
   vect_ct_[i].flush();
}

// }

// Process
// {

template <typename scalarType>
int ParallelContourTree::build(const bool ct, const bool segment, const double threshold)
{
#ifdef withOpenMP
    if(lessPartition_)
        omp_set_num_threads(nbPartitions_*2);
    else
        omp_set_num_threads(nbPartitions_);

// Get number of proc <- find number of socket for NUMA ??
// omp_get_num_procs();

// std::thread::hardware_concurrency();
#endif

   if (debugLevel_ > 3) {
       cout << " nb partition : " << static_cast<unsigned>(nbPartitions_);
       if(lessPartition_){
         cout << " with independant merge trees ";
       }
       cout << " nb vertices " << mesh_->getNumberOfVertices() << endl;
       cout << endl;
   }

   DebugTimer timerAllocSeq;

   setSegmentation(segment);
   setComputeContourTree(ct);

   // reserve interface partition et cts
   flush();
   vector<idVertex> seeds;

   printDebug(timerAllocSeq, "Allocations Sequentielles        ");

   DebugTimer timerTOTAL;

   DebugTimer timerSort;
   sortInput<scalarType>();
   printDebug(timerSort, "Sort                             ");

   initInterfaces();

   // {
   // cross edges to initialize overlap bounds
   //DebugTimer timerSubdivide;
   //seeds.reserve(nbInterfaces_);
   //for (Interface &interface : vect_interfaces_) {
      //seeds.push_back(interface.getSeed());
   //}
   //mesh_->subdivide(seeds, sortedVertices_, mirrorOffsets_);

   //for (idInterface i = 0; i < nbInterfaces_; i++) {
      //vect_interfaces_[i].setUpper(mirrorOffsets_[seeds[i]]);
      //vect_interfaces_[i].setLower(mirrorOffsets_[seeds[i]]);
   //}

   //printDebug(timerSubdivide, "Subdivize mesh                   ");
   // }

   DebugTimer timerInitOverlap;
   initOverlap();
   if(debugLevel_ > 2){
        for (unsigned i = 0; i < nbInterfaces_; i++) {
            cout << "interface : " << i;
            cout << " seed : " << vect_interfaces_[i].getSeed();
            cout << endl;
        }
   }
   printDebug(timerInitOverlap, "Initialize Overlap               ");

   DebugTimer                          timerAllocPara;
   vector<vector<ExtendedUnionFind *>> vect_baseUF_JT(nbPartitions_), vect_baseUF_ST(nbPartitions_);
   const auto nbVert = mesh_->getNumberOfVertices();

#pragma omp parallel for num_threads(nbPartitions_) schedule(static)
   for (decltype(nbPartitions_) tree = 0; tree < nbPartitions_; ++tree) {
      // CT
      initLocalCT<scalarType>(tree);
      vect_ct_[tree].initDataMT<scalarType>();

      // JT
      vect_ct_[tree].jt_->vect_nodes_.reserve((nbVert / nbPartitions_) / 10);
      vect_ct_[tree].jt_->vect_superArcs_.reserve((nbVert / nbPartitions_) / 10);

      // ST
      vect_ct_[tree].st_->vect_nodes_.reserve((nbVert / nbPartitions_) / 10);
      vect_ct_[tree].st_->vect_superArcs_.reserve((nbVert / nbPartitions_) / 10);

      // Big reserve
      vect_baseUF_JT[tree].resize(nbVert);
      vect_baseUF_ST[tree].resize(nbVert);
   }
   printDebug(timerAllocPara, "Alloc parallel                   ");

   DebugTimer timerbuild;
   parallelBuild<scalarType>(vect_baseUF_JT, vect_baseUF_ST, threshold);
   printDebug(timerbuild, "ParallelBuild                    ");

   DebugTimer timerZip;
   if (computeContourTree_ && nbPartitions_ > 1 && partitionNum_ == -1) {
      stitch();
      for (unsigned p = 0; p < nbPartitions_; ++p) {
         vect_ct_[p].parallelInitNodeValence(nbPartitions_);
      }
   }
   printDebug(timerZip, "Stitch                           ");

   if (debugLevel_ >= 4) {
      printVectCT();
   }


   const unsigned nbThreadReal = (lessPartition_)?nbPartitions_*2:nbPartitions_;

   // Unify to create a normal tree for wrapper / test
   DebugTimer timerUnify;
   if (ct) {
       if(partitionNum_ != -1){
           if(partitionNum_ > nbInterfaces_){
              shallowCopy(&vect_ct_[nbInterfaces_]);
           } else {
              shallowCopy(&vect_ct_[partitionNum_]);
           }
       } else if (nbPartitions_ == 1) {
          shallowCopy(&vect_ct_[0]);
       } else {
          unifyCT();
          // for global simlify
          parallelInitNodeValence(nbThreadReal);
       }

   } else {
       if(partitionNum_ != -1){
           if(partitionNum_ > nbInterfaces_){
              jt_ = vect_ct_[nbInterfaces_].jt_->clone();
              st_ = vect_ct_[nbInterfaces_].st_->clone();
           } else {
              jt_ = vect_ct_[partitionNum_].jt_->clone();
              st_ = vect_ct_[partitionNum_].st_->clone();
           }
       } else {
          if (nbThread_ > 1) {
             cout << "Only partition 1 shown !" << endl;
          }
          jt_ = vect_ct_[0].jt_->clone();
          st_ = vect_ct_[0].st_->clone();
       }
       jt_->parallelInitNodeValence(nbThreadReal);
       st_->parallelInitNodeValence(nbThreadReal);
   }

   printDebug(timerUnify, "Create Contour tree              ");

   if(ct && partitionNum_ == -1 && nbPartitions_ > 1 && threshold){
      DebugTimer timerGlobalSimplify;
      idEdge simplifed = globalSimplify<scalarType>(-1, nullVertex, threshold);
      if(debugLevel_ >=1){
         cout << "Global Simplification                                       " << timerGlobalSimplify.getElapsedTime();
         cout << " ( " << simplifed << " pairs merged )" << endl;
      }
   }

   printDebug(timerTOTAL, "TOTAL                            ");

   //verifyTree();

    if (debugLevel_ >= 5)  {
        if(ct)
           printTree2();
        else {
           cout << "JT :" << endl;
           jt_->printTree2();
           cout << "ST :" << endl;
           st_->printTree2();
        }
    }

   // reclaim memory
   {
      for (decltype(nbPartitions_) tree = 0; tree < nbPartitions_; ++tree) {
         vect_ct_[tree].jt_->vect_nodes_.shrink_to_fit();
         vect_ct_[tree].jt_->vect_superArcs_.shrink_to_fit();
         vect_ct_[tree].st_->vect_nodes_.shrink_to_fit();
         vect_ct_[tree].st_->vect_superArcs_.shrink_to_fit();
      }
   }

   if(ct){
      updateSegmentation(true);
   } else {
      jt_->updateSegmentation();
      st_->updateSegmentation();
   }

   cout << "Contour Tree computed " << endl;

   return 0;
}

template <typename scalarType>
int ParallelContourTree::parallelBuild(vector<vector<ExtendedUnionFind *>> &vect_baseUF_JT,
                                       vector<vector<ExtendedUnionFind *>> &vect_baseUF_ST,
                                       const double                         threshold)
{
   const auto &nbVert = mesh_->getNumberOfVertices();
   vector<float> timeSimplify(nbPartitions_, 0);
   vector<float> speedProcess(nbPartitions_*2, 0);
   idEdge nbPairMerged = 0;

#ifdef withOpenMP
   omp_set_nested(1);
#endif

//cout << "NO PARALLEL DEBUG MODE" << endl;
#pragma omp parallel for num_threads(nbPartitions_) schedule(static)
   for (int i = 0; i < nbPartitions_; ++i) {
      DebugTimer timerMergeTree;

       //if(i<nbInterfaces_)
       //cout << vect_interfaces_[i].getSeed()
      //<< " , pos : " << mirrorOffsets_[vect_interfaces_[i].getSeed()] << endl;

      const idVertex &startJT =
          (i == 0) ? 0 : mirrorOffsets_[vect_interfaces_[i - 1].getSeed()];
      const idVertex &startST =
          (i == nbInterfaces_) ? nbVert - 1 : mirrorOffsets_[vect_interfaces_[i].getSeed()]-1;

      const idVertex &endJT =
          (i == nbInterfaces_) ? nbVert : mirrorOffsets_[vect_interfaces_[i].getSeed()];
      const idVertex &endST =
          (i == 0) ? -1 : mirrorOffsets_[vect_interfaces_[i - 1].getSeed()] - 1;

      const idVertex &posSeed0 = (i == 0) ? -1 : mirrorOffsets_[vect_interfaces_[i - 1].getSeed()];
      const idVertex &posSeed1 =
          (i == nbInterfaces_) ? nullVertex : mirrorOffsets_[vect_interfaces_[i].getSeed()];

      const vector<idVertex> &lowerOverlap =
          (i == 0) ? vector<idVertex>() : vect_interfaces_[i - 1].getLower();

      const vector<idVertex> &upperOverlap =
          (i == nbInterfaces_) ? vector<idVertex>() : vect_interfaces_[i].getUpper();

      const idVertex partitionSize =
          abs(endJT - startJT) + lowerOverlap.size() + upperOverlap.size();

      if (lessPartition_) {
#pragma omp parallel sections num_threads(2)
         {
#pragma omp section
            {
                DebugTimer timerSimplify;
                DebugTimer timerBuild;
               vect_ct_[i].getJoinTree()->build(vect_baseUF_JT[i], lowerOverlap, upperOverlap,
                                                startJT, endJT, posSeed0, posSeed1);
               speedProcess[i] = partitionSize / timerBuild.getElapsedTime();

               timerSimplify.reStart();
               const idEdge tmpMerge =
                   vect_ct_[i].getJoinTree()->localSimplify<scalarType>(posSeed0, posSeed1, threshold);
#pragma omp atomic update
               timeSimplify[i] += timerSimplify.getElapsedTime();
#pragma omp atomic update
               nbPairMerged += tmpMerge;

            }

#pragma omp section
            {
               DebugTimer timerSimplify;
               DebugTimer timerBuild;
               vect_ct_[i].getSplitTree()->build(vect_baseUF_ST[i], upperOverlap, lowerOverlap,
                                                 startST, endST, posSeed0, posSeed1);
               speedProcess[nbPartitions_ + i] = partitionSize / timerBuild.getElapsedTime();

               timerSimplify.reStart();
               const idEdge tmpMerge =
                   vect_ct_[i].getSplitTree()->localSimplify<scalarType>(posSeed0, posSeed1, threshold);
#pragma omp atomic update
               timeSimplify[i] += timerSimplify.getElapsedTime();
#pragma omp atomic update
               nbPairMerged += tmpMerge;
            }
         }
      }

      else {
         DebugTimer timerSimplify;
         DebugTimer timerBuild;
         vect_ct_[i].getJoinTree()->build(vect_baseUF_JT[i], lowerOverlap, upperOverlap, startJT,
                                          endJT, posSeed0, posSeed1);
         speedProcess[i] = partitionSize / timerBuild.getElapsedTime();

         //if(nbPartitions_ > 1){
             timerSimplify.reStart();
             nbPairMerged +=
                 vect_ct_[i].getJoinTree()->localSimplify<scalarType>(posSeed0, posSeed1, threshold);
             timeSimplify[i] += timerSimplify.getElapsedTime();
         //}

         vect_ct_[i].getSplitTree()->build(vect_baseUF_ST[i], upperOverlap, lowerOverlap, startST,
                                           endST, posSeed0, posSeed1);
         speedProcess[nbPartitions_ + i] = partitionSize / timerBuild.getElapsedTime();

         //if(nbPartitions_ > 1){
             timerSimplify.reStart();
             nbPairMerged +=
                 vect_ct_[i].getSplitTree()->localSimplify<scalarType>(posSeed0, posSeed1, threshold);
             timeSimplify[i] += timerSimplify.getElapsedTime();
         //}
      }

      {
         stringstream mt;
         mt << "[ParallelBuild] Merge Tree " << static_cast<unsigned>(i)
            << " constructed in : " << timerMergeTree.getElapsedTime() << endl;
         dMsg(cout, mt.str(), infoMsg);
      }

      if (segmentation_ && (threshold || !computeContourTree_)) {
         DebugTimer timerUpdateSegm;
         vect_ct_[i].getJoinTree()->updateSegmentation();
         vect_ct_[i].getSplitTree()->updateSegmentation();

         if (debugLevel_ >= 3) {
            cout << "Local MT : updated in " << timerUpdateSegm.getElapsedTime() << endl;
         }
      }

      if (0 || computeContourTree_) {
         DebugTimer timerCombine;

         auto jt = vect_ct_[i].getJoinTree();
         auto st = vect_ct_[i].getSplitTree();

         // Copy missing nodes of a tree to the other one
         // Maintain this traversal order for good insertion
         for (int t = 0; t < st->getNumberOfNodes(); ++t) {
            if (!st->getNode(t)->isHidden()) {
               // cout << "insert in jt : " << st->getNode(t)->getVertexId() << endl;
               jt->insertNode(st->getNode(t), segmentation_);
            }
         }
         // and vice versa
         for (int t = 0; t < jt->getNumberOfNodes(); ++t) {
            if (!jt->getNode(t)->isHidden()) {
               // cout << "insert in st : " << jt->getNode(t)->getVertexId() << endl;
               st->insertNode(jt->getNode(t), segmentation_);
            }
         }

         if (debugLevel_ >= 6) {
            cout << "Local JT :" << endl;
            vect_ct_[i].getJoinTree()->printTree2();
            cout << "Local ST :" << endl;
            vect_ct_[i].getSplitTree()->printTree2();
            cout << "combine" << endl;
         }

         if (0 || computeContourTree_) {
            // DebugTimer timerNoise;

            // here Arcs & Nodes of the merge tree are destroyed
            vect_ct_[i].combine(posSeed0, posSeed1);

            vect_ct_[i].updateSegmentation(true);

            //vect_ct_[i].initNodeValence();

            if(debugLevel_ >2){
                printDebug(timerCombine, "Trees combined   in    ");
            }

            if (debugLevel_ >= 4) {
               vect_ct_[i].printTree2();
            }

            //if (nbPartitions_ > 1) {
               //DebugTimer   timerSimplify;
               //const idEdge tmpMerge =
                   //vect_ct_[i].globalSimplify<scalarType>(posSeed0, posSeed1, threshold);
//#pragma omp atomic update
               //timeSimplify[i] += timerSimplify.getElapsedTime();
//#pragma omp atomic update
               //nbPairMerged += tmpMerge;

               //vect_ct_[i].updateSegmentation(true);

               //vect_ct_[i].initNodeValence();
            //}
         }
      } else {
         // cout << "JT crossing below : " << endl;
         // for(auto & arc : vect_ct_[i].getJoinTree()->vect_arcsCrossingBelow_){
         // SuperArc * a = vect_ct_[i].getJoinTree()->getSuperArc(arc);
         // cout << vect_ct_[i].getJoinTree()->getNode(a->getDownNodeId())->getVertexId();
         // cout << " - ";
         // cout << vect_ct_[i].getJoinTree()->getNode(a->getUpNodeId())->getVertexId();
         // cout << " :: " << a->getOverlapBelow();
         // cout << endl;
         //}

         // cout << "JT crossing above : " << endl;
         // for(auto & arc : vect_ct_[i].getJoinTree()->vect_arcsCrossingAbove_){
         // SuperArc * a = vect_ct_[i].getJoinTree()->getSuperArc(arc);
         // cout << vect_ct_[i].getJoinTree()->getNode(a->getDownNodeId())->getVertexId();
         // cout << " - ";
         // cout << vect_ct_[i].getJoinTree()->getNode(a->getUpNodeId())->getVertexId();
         // cout << " :: " << a->getOverlapAbove();
         // cout << endl;
         //}

         // cout << "ST crossing below: " << endl;
         // for(auto & arc : vect_ct_[i].getSplitTree()->vect_arcsCrossingBelow_){
         // SuperArc * a = vect_ct_[i].getSplitTree()->getSuperArc(arc);
         // cout << vect_ct_[i].getSplitTree()->getNode(a->getDownNodeId())->getVertexId();
         // cout << " - ";
         // cout << vect_ct_[i].getSplitTree()->getNode(a->getUpNodeId())->getVertexId();
         // cout << " :: " << a->getOverlapBelow();
         // cout << endl;
         //}

         // cout << "ST crossing above: " << endl;
         // for(auto & arc : vect_ct_[i].getSplitTree()->vect_arcsCrossingAbove_){
         // SuperArc * a = vect_ct_[i].getSplitTree()->getSuperArc(arc);
         // cout << vect_ct_[i].getSplitTree()->getNode(a->getDownNodeId())->getVertexId();
         // cout << " - ";
         // cout << vect_ct_[i].getSplitTree()->getNode(a->getUpNodeId())->getVertexId();
         // cout << " :: " << a->getOverlapAbove();
         // cout << endl;
         //}

         if (debugLevel_ >= 5) {
            cout << "Local JT :" << endl;
            vect_ct_[i].getJoinTree()->printTree2();
            cout << "Local ST :" << endl;
            vect_ct_[i].getSplitTree()->printTree2();
            cout << "combine" << endl;
         }
      }
   }

   if (debugLevel_ >= 1 && threshold && nbPartitions_ > 1) {
      auto  maxSimplifIt = max_element(timeSimplify.cbegin(), timeSimplify.cend());
      float maxSimplif   = *maxSimplifIt;
      cout << "Local simplification maximum time :                         " << maxSimplif;
      cout << " ( " << nbPairMerged << " pairs merged )" << endl;
   }

   if(debugLevel_ > 1) {
      auto maxProcSpeed = max_element(speedProcess.cbegin(), speedProcess.cend());
      auto minProcSpeed = min_element(speedProcess.cbegin(), speedProcess.cend());
      cout << "process speed : ";
      cout << " min is " << *minProcSpeed << " vert/sec";
      cout << " max is " << *maxProcSpeed << " vert/sec";
      cout << endl;
   }

   return 0;
}

//}

#endif /* end of include guard: PARALLELCONTOURTREETEMPLATE_H */
