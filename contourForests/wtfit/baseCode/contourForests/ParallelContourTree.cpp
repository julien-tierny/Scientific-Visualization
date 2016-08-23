
/*
 * file:                  ContourTree.cpp
 * description:           ContourTree processing package.
 * author:                Gueunet Charles
 * date:                  Aout 2015
 */

#include <ParallelContourTree.h>

Interface::Interface(const idVertex &seed) : seed_(seed)
{
}

/// ------------------------- ParallelContourTree

ParallelContourTree::ParallelContourTree(const numThread &nbThread)
    : ContourTree(),
      nbThread_(nbThread),
      nbInterfaces_(nbThread - 1),
      nbPartitions_(nbThread),
      partitionNum_(-1),
      lessPartition_(false)
{
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

   vect_interfaces_.reserve(nbThread - 1);
   // vect_corrNodes_.reserve(nbThread);
   // vect_corrSuperArc_.reserve(nbThread);
   vect_ct_.resize(nbThread);
}

ParallelContourTree::ParallelContourTree() : ParallelContourTree(1)
{
}

ParallelContourTree::~ParallelContourTree()
{
   // Shared info
   if (destroyVectorSortedVertices_ && sortedVertices_.size()) {
      sortedVertices_.clear();
   }

   mirrorOffsets_.clear();

   if (destroyVectorSoS_ && soSOffsets_) {
      delete[] soSOffsets_;
      soSOffsets_ = nullptr;
   }

   if (destroyVect2Tree_ && vect_vert2tree_) {
      delete vect_vert2tree_;
      vect_vert2tree_ = nullptr;
   }
   // Personal info
   // Problem com from the destuction of the parallelCT, calling ~CT and so delete this->jt_ that
   // have been set at vect_ct_[0].jt_ => unify problem
}

// Init
// {
void ParallelContourTree::initInterfaces()
{
   // We have nbThread_ partition of the same size through all vertices
   size_t partitionSize = mesh_->getNumberOfVertices() / nbPartitions_;

   // ------------------
   // Seeds
   // ------------------
   // {

   // We initiate interface with their seed (isovalue) and their adjacent partition
   //  and each partition with it size and bounds.
   for (decltype(nbInterfaces_) i = 0; i < nbInterfaces_; ++i) {
      // interfaces have their first vertex of the sorted array as seed
      vect_interfaces_.emplace_back(sortedVertices_[partitionSize * (i + 1)]);
   }

   // }
   // ------------------
   // Print Debug
   // ------------------
   // {

   if (debugLevel_ >= 4) {
      stringstream partition;
      partition << "seeds :";
      for (const auto &i : vect_interfaces_) {
         partition << i.getSeed() << " ";
      }
      partition << endl;
      dMsg(cout, partition.str(), 3);
   }

   // }
}

void ParallelContourTree::initOverlap()
{
   const idEdge nbEdges = mesh_->getNumberOfEdges();

   // if we choose to have less partition, we still want to use all thread for overlap init.
   const unsigned nbThreadOverlap = (lessPartition_)?nbPartitions_*2:nbPartitions_;

   // ------------------
   // Parallel find border vertices
   // ------------------
   // {

   vector<vector<vector<idVertex>>> lowers(nbThreadOverlap);
   vector<vector<vector<idVertex>>> uppers(nbThreadOverlap);

   for (unsigned p = 0; p < nbThreadOverlap; p++) {
        lowers[p].resize(nbInterfaces_);
        uppers[p].resize(nbInterfaces_);
   }

#pragma omp parallel for num_threads(nbThreadOverlap) schedule(static)
   for (idEdge e = 0; e < nbEdges; e++) {

#ifdef withOpenMP
       unsigned part = omp_get_thread_num();
#else
       unsigned part = 0;
#endif

       vector<vector<idVertex>> &localUppers = uppers[part];
       vector<vector<idVertex>> &localLowers = lowers[part];

       idVertex v0, v1;
       mesh_->getEdgeVertex(e,0,v0);
       mesh_->getEdgeVertex(e,1,v1);

       for (unsigned i = 0; i < nbInterfaces_; i++) {
           const bool side0 = isEqHigher(v0, vect_interfaces_[i].getSeed());
           const bool side1 = isEqHigher(v1, vect_interfaces_[i].getSeed());

           if(side0 != side1){
              // edge cross this interface, add both extrema in it
              if(side0){
                 // The seed is already in the partition, we do not want to have it twice
                 //if (v0 != vect_interfaces_[i].getSeed()) {
                    localUppers[i].emplace_back(v0);
                 //}
                 localLowers[i].emplace_back(v1);
              } else {
                 //if (v1 != vect_interfaces_[i].getSeed()) {
                    localUppers[i].emplace_back(v1);
                 //}
                 localLowers[i].emplace_back(v0);
              }
           }
       }
   }

   // }
   // --------------------------
   // Insert in interfaces
   // --------------------------
   // {

   // reserve
   vector<idVertex> sizeReserveUp(nbInterfaces_, 0);
   vector<idVertex> sizeReserveLo(nbInterfaces_, 0);
   for (unsigned p = 0; p < nbThreadOverlap; p++) {
       for (unsigned i = 0; i < nbInterfaces_; i++) {
           sizeReserveUp[i] += uppers[p][i].size();
           sizeReserveLo[i] += lowers[p][i].size();
       }
   }

   for (unsigned i = 0; i < nbInterfaces_; i++) {
       vect_interfaces_[i].upReserve(sizeReserveUp[i]);
       vect_interfaces_[i].loReserve(sizeReserveLo[i]);
   }

   // append
   for (unsigned p = 0; p < nbThreadOverlap; p++) {
       for (unsigned i = 0; i < nbInterfaces_; i++) {
           vect_interfaces_[i].appendUpper(uppers[p][i]);
           vect_interfaces_[i].appendLower(lowers[p][i]);
       }
   }

   // }
   // -----------------
   // Sort the overlap
   // ----------------
   // {

   // As we don't care about noise in this part, do we need to sort ?

   auto isLowerComp = [&](const idVertex &a, const idVertex &b){
        return isLower(a,b);
   };

#pragma omp parallel for num_threads(nbInterfaces_) schedule(static)
   for (unsigned i = 0; i < nbInterfaces_; i++) {
      vector<idVertex> &upOverlap = vect_interfaces_[i].getUpper();
      vector<idVertex> &loOverlap = vect_interfaces_[i].getLower();

      sort(upOverlap.begin(), upOverlap.end(), isLowerComp);
      auto upLast = unique(upOverlap.begin(), upOverlap.end());
      upOverlap.erase(upLast, upOverlap.end());
      sort(loOverlap.begin(), loOverlap.end(), isLowerComp);
      auto loLast = unique(loOverlap.begin(), loOverlap.end());
      loOverlap.erase(loLast, loOverlap.end());
   }

   // }
   // -----------
   // Debug print
   // -----------
   // {

   //for (unsigned i = 0; i < nbInterfaces_; i++) {
      //cout << "interface : " << i << endl;

      //cout << "upper" << endl;
      //for (const idVertex &v : vect_interfaces_[i].getUpper()) {
         //cout << v << ", ";
      //}

      //cout << endl << "lower" << endl;
      //for (const idVertex &v : vect_interfaces_[i].getLower()) {
         //cout << v << ", ";
      //}

      //cout << endl;
   //}

   // }

}

void ParallelContourTree::flush(void)
{
   MergeTree::flush();
   vect_interfaces_.reserve(nbInterfaces_);
   vect_ct_.resize(nbPartitions_);
}
// }

// Process
// {

void ParallelContourTree::stitch()
{
   // We need goods informations here befor starting
   // Get the arc/NODE correspoding to the seed + crossingEdges

   //cout << " Stitch : " << endl;
   if (computeContourTree_) {
      //cout << "ct" << endl;
      stitchTree(2);
   } else {
      //cout << "st" << endl;
      stitchTree(0);
      //cout << "jt" << endl;
      stitchTree(1);
   }
}

void ParallelContourTree::stitchTree(const char treetype)
{

    const bool DEBUG = false;
    // CAUTION, printArc may segfault because arc may have nodes frome several partitions

   // For the instance assume ct :
   //    insert node with true (jt / ct ok)
   auto getTree = [&, treetype](const idPartition &i) -> MergeTree * {
      if (treetype == 0) {
         return vect_ct_[i].getJoinTree();
      }
      if (treetype == 1) {
         return vect_ct_[i].getSplitTree();
      }
      return &vect_ct_[i];
   };

   // the last partition can't stitch above
   for (int i = 0; i < nbPartitions_-1; i++) {
      MergeTree *tree = getTree(i);
      // we stitch with above tree
      for (const auto &arc : tree->vect_arcsCrossingAbove_) {
         SuperArc *crossing = tree->getSuperArc(arc);

         // ignore useless arc (should not happend)
         if (!crossing->isVisible()) {
            if (DEBUG) {
               cout << "ignore hidden arc : " << tree->printArc(arc) << endl;
               cout << endl;
            }
            continue;
         }

         if (DEBUG) {
            cout << "crossing above ->" << tree->printArc(arc) << endl;
            cout << " (part: " << i << ")" << endl;
         }

         const pair<idVertex, bool> &seed = make_pair(vect_interfaces_[i].getSeed(), false);
         if (isLower(tree->getNode(crossing->getUpNodeId())->getVertexId(), seed.first)) {
            if (DEBUG) {
               cout << "up is below seed !! " << endl;
            }
            continue;
         }

         // info about the other side
         // the arc cross above, that mean we are not on the last partition
         const idPartition &otherPartition = i + 1;
         const idVertex stitchVert = tree->cutArcAboveSeed(arc, seed);

         // Curent tree stitch node
         if (!tree->isCorrespondingNode(stitchVert)) {
            auto newNodeId = tree->makeNode(stitchVert);
            tree->updateCorrespondingArc(stitchVert, arc);

            if (DEBUG) {
               idSuperArc sa = tree->getCorrespondingSuperArcId(stitchVert);
               cout << "insert node in current in : ";
               cout << tree->printArc(sa);
               cout << endl;
            }

            tree->insertNode(tree->getNode(newNodeId), false, true);
            tree->getNode(newNodeId)->setUpValence(1);
            tree->getNode(newNodeId)->setDownValence(1);
         }
         const idNode curTreeStitchNodeId = tree->getCorrespondingNode(stitchVert);

         if(DEBUG){
            cout << "stitch vertex is : " << stitchVert << endl;
         }

         // Above tree stitch node
         idSuperArc toHideAbove = nullSuperArc;
         if (getTree(otherPartition)->isCorrespondingArc(stitchVert)) {
            if (DEBUG) {
               idSuperArc arcOtherTree =
                   getTree(otherPartition)->getCorrespondingSuperArcId(stitchVert);
               cout << "rev. insert node in: ";
               cout << getTree(otherPartition)->printArc(arcOtherTree);
               cout << endl;
            }

            // insert stitchpoint
            toHideAbove =
                getTree(otherPartition)
                    ->reverseInsertNode(tree->getNode(curTreeStitchNodeId), segmentation_, true);
            tree->getNode(curTreeStitchNodeId)->setUpValence(1);
            tree->getNode(curTreeStitchNodeId)->setDownValence(1);
         }
         const idNode otherTreeStitchNodeId = getTree(otherPartition)->getCorrespondingNode(stitchVert);

         if(DEBUG){
            cout << "verify vertex : " << otherTreeStitchNodeId << " id : "
                 << getTree(otherPartition)->getNode(otherTreeStitchNodeId)->getVertexId() << endl;

            cout << "degree : " << static_cast<unsigned>(getTree(otherPartition)
                                                             ->getNode(otherTreeStitchNodeId)
                                                             ->getNumberOfDownSuperArcs());
            cout << " - " << static_cast<unsigned>(getTree(otherPartition)
                                                       ->getNode(otherTreeStitchNodeId)
                                                       ->getNumberOfUpSuperArcs());
            cout << " hidden : "
                 << static_cast<unsigned>(
                        getTree(otherPartition)->getNode(otherTreeStitchNodeId)->isHidden())
                 << endl;

            if (getTree(otherPartition)->getNode(otherTreeStitchNodeId)->getVertexId() != stitchVert) {
                cout << "stitching problem : " << endl;
                cout << getTree(otherPartition)->getNode(otherTreeStitchNodeId)->getVertexId()
                    << " should be " << stitchVert << endl;
            }
         }

         if (getTree(otherPartition)->getNode(otherTreeStitchNodeId)->isHidden()) {
            cout << "------> is hidden !" << endl;
            continue;
         }

         //if(toHideHere != nullSuperArc){
             //tree->hideArc(toHideHere);
             //if (DEBUG) {
                //cout << " hide in current : " << tree->printArc(toHideHere) << endl;
             //}
         //}

         idSuperArc leading = nullSuperArc;

         // Hide in tree above, arc leading to stitch
         if(toHideAbove != nullSuperArc){
             getTree(otherPartition)->hideArc(toHideAbove);
             leading = 0;
             if (DEBUG) {
                cout << " hide in above : " << getTree(otherPartition)->printArc(toHideAbove)
                     << endl;
             }
         } else if(crossing->getDownCT() == i){

            // hide and disconnect arc leading to : current -> can be node or segmentation
            const idVertex below =
                tree->cutArcBelowSeed(arc, seed, getTree(otherPartition)->vect_vert2tree_);
            leading = getTree(otherPartition)->hideAndClearLeadingTo(otherTreeStitchNodeId, below);

            if (DEBUG) {
               cout << " hide in above leading to at vert : " << static_cast<unsigned>(below);
               if(leading != nullSuperArc){
                  cout << " arc : " << getTree(otherPartition)->printArc(leading) << endl;
               } else {
                  cout << " no arc found " << endl;
               }
            }
         }

         // avoid duplicate if : /\      |
         if (tree->alreadyExtLinked(curTreeStitchNodeId, otherPartition, otherTreeStitchNodeId)) {
            if (DEBUG) {
               cout << "This arc already exist" << endl << endl;
            }
            continue;
         }

         // For gaussian 2 threads : hide all noise
         tree->getNode(curTreeStitchNodeId)->clearUpSuperArcs();

         // Insert

         // Arc from top to bottom
         // TODO TODO Might overlap above !! not necessarly false
         // If true, add in the vector !!
         bool stillAbove = false;

         if(i < nbPartitions_-2){
           const idVertex nextSeed = vect_interfaces_[otherPartition].getSeed();
           stillAbove = isLower(nextSeed, stitchVert);
         }
         getTree(otherPartition)
             ->vect_superArcs_.emplace_back(curTreeStitchNodeId, otherTreeStitchNodeId, true,
                                            stillAbove, i, otherPartition);
         getTree(otherPartition)
             ->getNode(otherTreeStitchNodeId)
             ->addDownSuperArcId(getTree(otherPartition)->vect_superArcs_.size() - 1);
         if(stillAbove){
            const idSuperArc & lastArc = getTree(otherPartition)->vect_superArcs_.size() - 1;
            getTree(otherPartition)->addCrossingAbove(lastArc);

            if (DEBUG) {
               cout << "add new crossing above arc" << endl;
            }
         }

         // VERY UGLY FIX
         if (leading == nullSuperArc) {
            getTree(otherPartition)->getSuperArc(
                    getTree(otherPartition)->getNumberOfSuperArcs() -1)->hide();
         }

         // Arc bottom to top
         tree->vect_superArcs_.emplace_back(curTreeStitchNodeId, otherTreeStitchNodeId, false, true,
                                            i, otherPartition);
         tree->getNode(curTreeStitchNodeId)->addUpSuperArcId(tree->vect_superArcs_.size() - 1);

         if (DEBUG) {
             cout << "added in " << static_cast<unsigned>(otherPartition) << " and current ";
             cout << "arc "
                  << getTree(otherPartition)->getNode(otherTreeStitchNodeId)->getVertexId();
             cout << " - " << tree->getNode(curTreeStitchNodeId)->getVertexId();
             cout << " stitch with above done" << endl;
             cout << endl;
         }
      }  // for each arc of this tree

   }  // for each partition

   if (DEBUG) {
     printVectCT();
   }
}

void ParallelContourTree::unify()
{
}

void ParallelContourTree::unifyCT()
{

   //int treetype = 2;
   //// For the instance assume ct :
   ////    insert node with true (jt / ct ok)
   //auto getTree = [&, treetype](const idPartition &i) -> MergeTree * {
      //if (treetype == 0) {
         //return vect_ct_[i].getSplitTree();
      //}
      //if (treetype == 1) {
         //return vect_ct_[i].getJoinTree();
      //}
      //return &vect_ct_[i];
   //};
   //

   MergeTree tmpTree;
   tmpTree.vect_vert2tree_ = new vector<idCorresp>(mesh_->getNumberOfVertices(), nullCorresp);
   tmpTree.vect_nodes_.reserve(vect_vert2tree_->size() / 10);
   tmpTree.vect_superArcs_.reserve(vect_vert2tree_->size() / 10);
   tmpTree.mesh_ = mesh_;

   // partition, node in partion, is a leaf
   queue<tuple<idInterface, idNode,bool>> queue_LeavesNodes;
   vector<unsigned char> vect_alreadySeen(mesh_->getNumberOfVertices(), 0);
   //const idVertex & maxVert = mesh_->getNumberOfOriginalVertices();

   const bool DEBUG = false;

   set<idVertex> addedVerts;

   for (idPartition partition = 0; partition < nbPartitions_; ++partition) {
      for (auto &l : vect_ct_[partition].vect_leaves_) {
         if (!vect_ct_[partition].getNode(l)->isHidden()) {
            idVertex leafVert = vect_ct_[partition].getNode(l)->getVertexId();

            // if not in partition
            if (partition != 0 && isLower(leafVert, vect_interfaces_[partition - 1].getSeed()))
               continue;
            if (partition != nbInterfaces_ &&
                isHigher(leafVert, vect_interfaces_[partition].getSeed()))
               continue;

            // if have stitched, will be intergrated in an arc
            if (vect_ct_[partition].getNode(l)->getNumberOfUpSuperArcs() &&
                vect_ct_[partition].getNode(l)->getNumberOfDownSuperArcs())
               continue;

            // if alone
            if (!vect_ct_[partition].getNumberOfVisibleArcs(l))
               continue;

            if(addedVerts.find(leafVert) == addedVerts.end()){
               queue_LeavesNodes.emplace(partition, l, true);

               if (DEBUG) {
                  cout << "will see : currentTree : " << static_cast<unsigned>(partition)
                       << " leave "
                       << static_cast<unsigned>(vect_ct_[partition].getNode(l)->getVertexId());
                  cout << " state hidden : "
                       << static_cast<unsigned>(vect_ct_[partition].getNode(l)->isHidden()) << endl;
               }

               ++vect_alreadySeen[leafVert];
               addedVerts.emplace(leafVert);
            }
         }
      }
   }

   addedVerts.clear();

   idPartition currentTree, retainTree;
   idNode      currentNode, retainNode;
   bool isLeaf;

   idNode oldUpNodeId;
   idNode newNodeId;

   idSuperArc oldSuperArcId;
   idSuperArc newSuperArcId;

   idVertex closingVertex;

   unsigned char upArc;

   pair<idVertex, bool> **appendVertList = nullptr;
   idVertex *appendSize = nullptr;
   unsigned  nbAppend;

   if (segmentation_) {
      appendVertList = new pair<idVertex, bool> *[nbInterfaces_ * 2];
      appendSize     = new idVertex[nbInterfaces_ * 2];
   }

   while (!queue_LeavesNodes.empty()) {
      tie(currentTree, currentNode, isLeaf) = queue_LeavesNodes.front();

      if (DEBUG) {
         cout << endl;
         cout << " process : currentTree : " << static_cast<unsigned>(currentTree) << " node "
              << static_cast<unsigned>(vect_ct_[currentTree].getNode(currentNode)->getVertexId())
              << endl;
      }

      // We are on the down node of a superArc
      // make traversal on currentNode
      // oldNodeId -> currentNode in currentTree

      newNodeId = tmpTree.makeNode(vect_ct_[currentTree].getNode(currentNode));
      if(isLeaf) tmpTree.vect_leaves_.emplace_back(newNodeId);

      if (DEBUG) {
         cout << "tree " << static_cast<unsigned>(currentTree) << " node "
              << static_cast<unsigned>(vect_ct_[currentTree].getNode(currentNode)->getVertexId());
         cout << " have : "
              << static_cast<unsigned>(
                     vect_ct_[currentTree].getNode(currentNode)->getNumberOfUpSuperArcs())
              << " up arcs";
         cout << " and "
              << static_cast<unsigned>(
                     vect_ct_[currentTree].getNode(currentNode)->getNumberOfDownSuperArcs())
              << " down arcs " << endl;
      }

      // openArc for each old upArc of this node
      retainTree       = currentTree;
      retainNode       = currentNode;
      idSuperArc nbArc = vect_ct_[retainTree].getNode(retainNode)->getNumberOfUpSuperArcs();

      for (upArc = 0; upArc < nbArc; ++upArc) {
         // Traverse regular node of this branch, update vect_vert2tree with currentArc
         currentTree   = retainTree;
         currentNode   = retainNode;
         oldSuperArcId = vect_ct_[currentTree].getNode(currentNode)->getUpSuperArcId(upArc);
         oldUpNodeId   = vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getUpNodeId();

         // ignore already processed arcs
         if (!vect_ct_[currentTree].getSuperArc(oldSuperArcId)->isVisible() ){
             //vect_ct_[currentTree]
                     //.getNode(vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getUpNodeId())
                     //->getVertexId() > mesh_->getNumberOfOriginalVertices()) {
             if (DEBUG) {
                idNode      up   = vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getUpNodeId();
                idPartition upct = vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getUpCT();

                idNode down = vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getDownNodeId();
                idPartition downct = vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getDownCT();

                cout << "ignore : " << currentTree << " arc "
                     << vect_ct_[downct].getNode(down)->getVertexId() << " - "
                     << vect_ct_[upct].getNode(up)->getVertexId() << endl;
            }
            continue;
         }

         if (DEBUG) {
            idNode      up   = vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getUpNodeId();
            idPartition upct = vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getUpCT();

            idNode      down   = vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getDownNodeId();
            idPartition downct = vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getDownCT();

            cout << "start : " << static_cast<unsigned>(currentTree) << " arc "
                 << vect_ct_[downct].getNode(down)->getVertexId() << " - "
                 << vect_ct_[upct].getNode(up)->getVertexId() << endl;
         }

         // there is no more overlapping when unifying
         newSuperArcId = tmpTree.openSuperArc(newNodeId, false, false);
         if (segmentation_) {
            tmpTree.getSuperArc(newSuperArcId)
                ->setVertList(vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getVertList());
            tmpTree.getSuperArc(newSuperArcId)
                ->setVertSize(vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getVertSize());
            // cout << "initial vert size : " << static_cast<unsigned>(
            // vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getVertSize() ) << endl;
         }

         nbAppend = 0;
         vect_ct_[currentTree].getSuperArc(oldSuperArcId)->hide();
         currentTree = vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getUpCT();

         if (DEBUG) {
            cout << "  arc in "
                 << static_cast<unsigned>(vect_ct_[currentTree].getNode(oldUpNodeId)->getVertexId())
                 << endl;
         }

         // traversal through regular nodes (and through partitions if needed)
         while (vect_ct_[currentTree].getNode(oldUpNodeId)->getDownValence() <= 1 &&
                vect_ct_[currentTree].getNode(oldUpNodeId)->getUpValence() == 1) {
            // update vect2tree
            tmpTree.updateCorrespondingArc(
                vect_ct_[currentTree].getNode(oldUpNodeId)->getVertexId(), newSuperArcId);
            // become the upNode get info in SuperArc!!
            oldSuperArcId = vect_ct_[currentTree].getNode(oldUpNodeId)->getUpSuperArcId(0);
            oldUpNodeId   = vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getUpNodeId();
            if (segmentation_) {
               appendVertList[nbAppend] =
                   vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getVertList();
               appendSize[nbAppend] =
                   vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getVertSize();
            }
            ++nbAppend;

            if (DEBUG) {
               idNode      up   = vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getUpNodeId();
               idPartition upct = vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getUpCT();

               idNode      down = vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getDownNodeId();
               idPartition downct = vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getDownCT();

               cout << "hide in " << static_cast<unsigned>(currentTree) << " arc "
                    << vect_ct_[downct].getNode(down)->getVertexId() << " - "
                    << vect_ct_[upct].getNode(up)->getVertexId() << " nbVert="
                    << static_cast<unsigned>(
                           vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getVertSize())
                    << endl;
            }

            vect_ct_[currentTree].getSuperArc(oldSuperArcId)->hide();
            currentTree = vect_ct_[currentTree].getSuperArc(oldSuperArcId)->getUpCT();

            if (DEBUG) {
               cout << "    traversal : currentTree : " << static_cast<unsigned>(currentTree)
                    << " node " << static_cast<unsigned>(
                                       vect_ct_[currentTree].getNode(oldUpNodeId)->getVertexId())
                    << endl;
            }
         }

         if (segmentation_) {
            tmpTree.getSuperArc(newSuperArcId)
                ->appendVertLists(appendVertList, appendSize, nbAppend);

            // cout << "size : " << static_cast<unsigned>(
            // tmpTree.getSuperArc(newSuperArcId)->getVertSize() ) << endl;
         }

         // closeOpenedArc with makeNode to retrieve node if already exist
         tmpTree.closeSuperArc(newSuperArcId,
                               tmpTree.makeNode(vect_ct_[currentTree].getNode(oldUpNodeId)), false,
                               false);
         closingVertex = vect_ct_[currentTree].getNode(oldUpNodeId)->getVertexId();
         idVertex tmpOpenedVert = tmpTree.getNode(tmpTree.getSuperArc(newSuperArcId)->getDownNodeId())->getVertexId();
         // hide useless arc
         if (closingVertex == tmpOpenedVert) {
            if (DEBUG) {
               cout << "unify remove arc : " << tmpTree.printArc(newSuperArcId) << endl;
            }
            tmpTree.hideArc(newSuperArcId);
         }

         // push closing node in queue
         ++vect_alreadySeen[closingVertex];

         if (DEBUG) {
            cout << "create :";
            cout << tmpTree.printArc(newSuperArcId) << " id : " << newSuperArcId;
            cout << endl;

            cout << " seen " << static_cast<unsigned>(vect_alreadySeen[closingVertex]) << endl;
            cout << " nb down arcs "
                 << static_cast<unsigned>(
                        vect_ct_[currentTree].getNode(oldUpNodeId)->getDownValence())
                 << endl;

            cout << "node " << vect_ct_[currentTree].getNode(oldUpNodeId)->getVertexId() << " seen "
                 << static_cast<unsigned>(vect_alreadySeen[closingVertex]) << " on "
                 << static_cast<unsigned>(
                        vect_ct_[currentTree].getNode(oldUpNodeId)->getDownValence())
                 << " tree " << static_cast<unsigned>(currentTree)
                 << endl;
         }

         if (vect_alreadySeen[closingVertex] >=
                 vect_ct_[currentTree].getNode(oldUpNodeId)->getDownValence()){

               queue_LeavesNodes.emplace(currentTree, oldUpNodeId, false);

               if (DEBUG) {
                  cout << "-- push : currentTree : " << static_cast<unsigned>(currentTree)
                       << " node " << static_cast<unsigned>(
                                          vect_ct_[currentTree].getNode(oldUpNodeId)->getVertexId())
                       << endl;
               }

         }
      }  // end for upArc

      queue_LeavesNodes.pop();
      if (DEBUG) {
          cout << endl;
      }
   }

   tmpTree.vect_superArcs_.shrink_to_fit();
   tmpTree.vect_nodes_.shrink_to_fit();

   // Do swaps vector
   vect_superArcs_ = tmpTree.vect_superArcs_;
   vect_nodes_     = tmpTree.vect_nodes_;
   vect_vert2tree_->swap(*tmpTree.vect_vert2tree_);
   vect_leaves_.swap(tmpTree.vect_leaves_);
}

void ParallelContourTree::unifyMT()
{
   cout << "UNIFY MT TODO" << endl;
}

// }

// Print
// {
void ParallelContourTree::printDebug(DebugTimer &timer, const string &str)
{
   stringstream msg;
   msg << "[ParallelContourTree] " << str << " " << static_cast<unsigned>(partition_) << " : "
       << timer.getElapsedTime() << endl;
   dMsg(cout, msg.str(), timeMsg);
}

void ParallelContourTree::printVectCT()
{
   int arcCTUp, arcCTDown;

   for (size_t nb = 0; nb < vect_ct_.size(); ++nb) {
      cout << "CT " << nb << endl;
      cout << "Nodes" << endl;

      for (const auto &n : vect_ct_[nb].getNodes()) {
         if (!n.isHidden()) {
            cout << "Node  " << n.getVertexId();

            if (n.isHidden())
               cout << " X ";

            cout << endl;
            cout << "  arc up : ";

            for (int i = 0; i < n.getNumberOfUpSuperArcs(); ++i) {
               cout << n.getUpSuperArcId(i) << " ";
            }

            cout << endl << " arc down : ";

            for (int i = 0; i < n.getNumberOfDownSuperArcs(); ++i) {
               cout << n.getDownSuperArcId(i) << " ";
            }

            cout << endl;
         }
      }

      cout << "Arcs" << endl;

      for (const auto &sa : vect_ct_[nb].getSuperArc()) {
         if (!sa.isHidden()) {
            arcCTDown = sa.getDownCT();
            arcCTUp   = sa.getUpCT();

            if (sa.getDownNodeId() == nullNodes) {
               cout << "||";

            } else {
               cout << static_cast<unsigned>(arcCTDown) << ":";
               cout << vect_ct_[arcCTDown].getNode(sa.getDownNodeId())->getVertexId();
            }

            if (sa.isHidden())
               cout << " <X> ";
            else if(!sa.isVisible())
               cout << " <-> ";
            else
               cout << " <>> ";

            if (sa.getUpNodeId() == nullNodes) {
               cout << "||";

            } else {
               cout << static_cast<unsigned>(arcCTUp) << ":";
               cout << vect_ct_[arcCTUp].getNode(sa.getUpNodeId())->getVertexId();
            }

            cout << endl;
         }
      }

      if (1) {
         cout << "Leaves" << endl;

         for (const auto &l : vect_ct_[nb].vect_leaves_)
            cout << " " << l;

         cout << endl;

         cout << "Roots" << endl;

         for (const auto &r : vect_ct_[nb].vect_roots_)
            cout << " " << r;

         cout << endl;
      }
   }
}
// }

