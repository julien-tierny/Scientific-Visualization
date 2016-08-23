/**
 * @file ContourTreeTemplate.h
 * @brief Template function for Contour Tree processing
 * @author Gueunet Charles
 * @version 1
 * @date 2016-03-30
 */

#ifndef CONTOURTREETEMPLATE_H
#define CONTOURTREETEMPLATE_H

#include <ContourTree.h>

/// ---------------------------- MT

// Init
// {

template <typename scalarType>
void MergeTree::sortInput(void)
{
   auto nbVertices = mesh_->getNumberOfVertices();

   if (!sortedVertices_.size()) {
      auto indirect_sort = [&](const size_t &a, const size_t &b) {
         return isLower<scalarType>(a, b);
      };

      sortedVertices_.resize(nbVertices, 0);

      iota(sortedVertices_.begin(), sortedVertices_.end(), 0);

#ifdef withOpenMP
      __gnu_parallel::sort(sortedVertices_.begin(), sortedVertices_.end(), indirect_sort);
#else
      sort(sortedVertices_.begin(), sortedVertices_.end(), indirect_sort);
#endif
      destroyVectorSortedVertices_ = true;
   }

   if (!mirrorOffsets_.size()) {
      DebugTimer timerMirror;
      mirrorOffsets_.resize(nbVertices);

#pragma omp parallel for
      for (int i = 0; i < nbVertices; i++) {
         mirrorOffsets_[sortedVertices_[i]] = i;
      }

      // stringstream msgMir;
      // msgMir << "Mirror    (include in sort) took           " << timerMirror.getElapsedTime()
      //<< endl;
      // dMsg(cout, msgMir.str(), timeMsg);
   }
}
// }

// Process
// {

// Simplify

template <typename scalarType>
idEdge MergeTree::localSimplify(const idVertex &posSeed0, const idVertex &posSeed1, const double threshold)
{
    const bool DEBUG = false;

    // -----------------
    // Persistance pairs
    // -----------------
    // {

    vector<tuple<idVertex, idVertex, scalarType, bool>> pairs;
    computePersistencePairs<scalarType>(pairs);

    if (DEBUG) {
       cout << "pairs : ( threshold : " << threshold << "  )" << endl;
       for (const auto &p : pairs) {
          const idVertex &thisOriginVert = get<0>(p);
          const idVertex &thisEndVert    = get<1>(p);
          const double &  thisPersist    = get<2>(p);
          const bool      thisNeedUp     = get<3>(p);
          cout << thisOriginVert << " - " << thisEndVert << " : " << thisPersist;
          cout << " ( " << thisNeedUp << " )" << endl;
       }
       cout << endl;
    }

    // }
    // --------------
    // Simplify
    // --------------
    // {

    return simplifyTree<scalarType>(posSeed0, posSeed1, threshold, pairs);

    // }
}

template <typename scalarType>
idEdge MergeTree::globalSimplify(const idVertex posSeed0, const idVertex posSeed1, const double threshold)
{

    // if null threshold, leave
    if (!threshold) {
        return 0;
    }

    //---------------------
    // Sort Nodes
    //---------------------
    //{

    auto isLowerComp = [&](const idNode &n1, const idNode &n2) {
        return isLower(getNode(n1)->getVertexId(), getNode(n2)->getVertexId());
    };

    const auto nbNode = vect_nodes_.size();

    vector<idNode> sortedNodes(nbNode);
    iota(sortedNodes.begin(), sortedNodes.end(), 0);
// Sort nodes by vertex scalar
//{
#ifdef withOpenMP
    __gnu_parallel::sort(sortedNodes.begin(), sortedNodes.end(), isLowerComp);
#else
    sort(sortedNodes.begin(), sortedNodes.end(), isLowerComp);
#endif
//}

    //}
    //---------------------
    // Make pairs
    //---------------------
    //{

    // origin, end, persistance, needToGoUp
    vector<tuple<idVertex, idVertex, scalarType,bool>> pairsJT;
    vector<tuple<idVertex, idVertex, scalarType,bool>> pairsST;

    computePersistencePairsMT<scalarType>(sortedNodes, pairsJT, pairsST);

    //}
    //---------------------
    // Fusionne & Sort pairs
    //---------------------
    //{

        auto pairComp = [](const tuple<idVertex, idVertex, scalarType, bool> &a,
                           const tuple<idVertex, idVertex, scalarType, bool> &b) {
            // sort by persistence
           return get<2>(a) < get<2>(b);
        };

        vector<tuple<idVertex, idVertex, scalarType,bool>> sortedPairs;
        size_t sizePJT = pairsJT.size(), sizePST = pairsST.size();

        sortedPairs.reserve(sizePJT + sizePST);

        sortedPairs.insert(sortedPairs.end(), pairsJT.begin(), pairsJT.end());
        sortedPairs.insert(sortedPairs.end(), pairsST.begin(), pairsST.end());

        // Sort pairs by persistence
        //{

#ifdef withOpenMP
   __gnu_parallel::sort(sortedPairs.begin(), sortedPairs.end(), pairComp);
#else
   sort(sortedPairs.begin(), sortedPairs.end(), pairComp);
#endif

    auto last = unique(sortedPairs.begin(), sortedPairs.end());
    sortedPairs.erase(last, sortedPairs.end());

    // Debug
     //if(!partition_){
        //cout << "threshold : " << threshold << endl;
        //for (const auto &p : sortedPairs) {
            //cout << "pair : " << get<0>(p) << " - " << get<1>(p);
            //cout << " persist : " << get<2>(p) << " / " << threshold;
            //cout << " from jt : " << get<3>(p) << endl;
        //}
     //}

        //}

    //}
    //---------------------
    // Traverse pairs and merge on the tree
    //---------------------
    //{

        // identify subtrees and merge them in recept'arcs
        return simplifyTree<scalarType>(posSeed0, posSeed1, threshold, sortedPairs);
    //}
}

template <typename scalarType>
idEdge MergeTree::simplifyTree( const idVertex &posSeed0, const idVertex &posSeed1,
    const double threshold, const vector<tuple<idVertex, idVertex, scalarType, bool>> &sortedPairs)
{
   const auto nbNode = vect_nodes_.size();
   const auto nbArcs = vect_superArcs_.size();
   // Retain the relation between merge coming from st, jt
   // also retain info about what we keep
   vector<ExtendedUnionFind *> subtreeUF(nbNode, nullptr);

   // nb arc seen below / above this node
   vector<pair<short, short>> valenceOffset(nbNode, make_pair(0,0));
   idEdge nbPairMerged = 0;

   const bool DEBUG = false;

   if (DEBUG) {
      cout << "Imapct simplify on tree btwn : " << posSeed0 << " and " << posSeed1 << endl;
   }

   //----------
   // Make subtrees
   //-----------
   //{

   queue<tuple<idNode,bool>> node2see;

   // Add the origin of all pairs that need to merge
   for (const auto & pp : sortedPairs) {
      if (get<2>(pp) < threshold ) {
         const idVertex &thisOriginVert = get<0>(pp);
         const idVertex &thisEndVert    = get<1>(pp);
         const idNode &  thisOriginId   = getCorrespondingNode(thisOriginVert);
         //const idNode &  thisEndId      = getCorrespondingNode(thisEndVert);

         if (mirrorOffsets_[thisOriginVert] <= posSeed0 ||
             mirrorOffsets_[thisOriginVert] >= posSeed1 ||
             mirrorOffsets_[thisEndVert] <= posSeed0 || mirrorOffsets_[thisEndVert] >= posSeed1) {
            continue;
         }

         node2see.emplace(thisOriginId, get<3>(pp));
         subtreeUF[thisOriginId] = new ExtendedUnionFind(0);
         ++nbPairMerged;
         if(DEBUG){
            cout << "willSee " << thisOriginVert << " pos " << mirrorOffsets_[thisOriginVert] << endl;
         }
      } else break;
   }

   //---
   // In UF :
   //  Origin is the size of the segmentation
   //  Data is negative : -idNodeRoot-1
   //  Data is positive : Receptacle Arc id
   //--
   // Use the queue to mark Arc that will be merged and UF to identify
   // subtree. When a node have only one way out : enqueue it to continue travresall
   while (!node2see.empty()) {
       idNode curNodeId;
       bool needToGoUp; // identify up/down traversall

       tie(curNodeId, needToGoUp) = node2see.front(); // should have only one arc valid to take
       node2see.pop();

       // Here we take the only available arc :
       idSuperArc mergingArcId;
       idNode parentNodeId;
       // continue traversall
       if(needToGoUp){
          mergingArcId = newUpArc(curNodeId, subtreeUF);
          parentNodeId = getSuperArc(mergingArcId)->getUpNodeId();
          ++valenceOffset[curNodeId].second;
          ++valenceOffset[parentNodeId].first;
       } else {
          mergingArcId = newDownArc(curNodeId, subtreeUF);
          parentNodeId = getSuperArc(mergingArcId)->getDownNodeId();
          ++valenceOffset[curNodeId].first;
          ++valenceOffset[parentNodeId].second;
       }

       markThisArc(subtreeUF, curNodeId, mergingArcId, parentNodeId);

       // if we have processed all but one arc of this node, we nee to continue traversall
       // throug it
       if (valenceOffset[parentNodeId].first + valenceOffset[parentNodeId].second + 1 ==
           getNode(parentNodeId)->getValence()) {
           // only one way out, is it up ?
           node2see.emplace(parentNodeId, valenceOffset[parentNodeId].second + 1 ==
                                              getNode(parentNodeId)->getUpValence());
           if(DEBUG){
               cout << " add to see " << getNode(parentNodeId)->getVertexId() << endl;
           }
       }
   } // end while node2see

   // for each node valenceOffset is the number of arc attached to this node that will merge

   // Debug print
   if (DEBUG) {
      cout << "node subtrees before creating receptarc " << endl;
      for (idNode nid = 0; nid < nbNode; nid++) {
         if (subtreeUF[nid]) {
            cout << "node " << getNode(nid)->getVertexId() << " is in subtree rooted :";
            const idNode &root = -subtreeUF[nid]->find()->getData() - 1;
            cout << getNode(root)->getVertexId();
            const idVertex &segmSize = subtreeUF[nid]->find()->getOrigin();
            cout << " with segmentation of " << segmSize << endl;
         }
      }
   }

   //}
   //----------
   // Create the recept'arcs
   //-----------
   //{

   // Add the origin of all pairs that need to merge
   for (const auto & pp : sortedPairs) {
      if (get<2>(pp) < threshold ) {
         const idVertex &thisOriginVert = get<0>(pp);
         const idVertex &thisEndVert    = get<1>(pp);
         const idNode &  thisOriginId   = getCorrespondingNode(thisOriginVert);
         //const idNode &  thisEndId      = getCorrespondingNode(thisEndVert);

         if (mirrorOffsets_[thisOriginVert] <= posSeed0 ||
             mirrorOffsets_[thisOriginVert] >= posSeed1 ||
             mirrorOffsets_[thisEndVert] <= posSeed0 || mirrorOffsets_[thisEndVert] >= posSeed1) {
            continue;
         }

         if (subtreeUF[thisOriginId]->find()->getData() < 0) {
            // create receptarc
            const idNode &subtreeRoot = -subtreeUF[thisOriginId]->find()->getData() - 1;
            // The id of the next arc to be created : NOT PARALLEL
            const idSuperArc receptArcId = vect_superArcs_.size();
            // down , up, segmentation size
            // create the receptacle arc and merge arc not in sub-tree in it
            const tuple<idNode, idNode, idVertex> &receptArc =
                createReceptArc(subtreeRoot, receptArcId, subtreeUF, valenceOffset);

            // make superArc and do the makeAlloc on it
            const bool overlapB = mirrorOffsets_[getNode(get<0>(receptArc))->getVertexId()] < posSeed0;
            const bool overlapA = mirrorOffsets_[getNode(get<1>(receptArc))->getVertexId()] >= posSeed1;
            const idSuperArc na = makeSuperArc(get<0>(receptArc), get<1>(receptArc), overlapB, overlapA, nullptr, -1);

            if(overlapB){
                vect_arcsCrossingBelow_.emplace_back(na);
            }

            if(overlapA){
                vect_arcsCrossingAbove_.emplace_back(na);
            }

            subtreeUF[thisOriginId]->find()->setData(receptArcId);
            getSuperArc(receptArcId)->makeAllocGlobal(get<2>(receptArc));

            if (DEBUG) {
               cout << "create arc : " << printArc(receptArcId)
                    << " with segm : " << get<2>(receptArc) << endl;
            }
         }
      } else break;
   }

   //}
   //----------
   // Merge these subtrees
   //-----------
   //{

   // The merge is done after because we don't want to forget arcs parallel to the recept'arc
   // but we cannot make the difference before the former are created

   // nbArcs is before the insertion of receptarcs so they will no be crossed here
   for (idSuperArc arc = 0; arc < nbArcs; arc++) {

       if(getSuperArc(arc)->isMerged()){
          const idSuperArc &receptacleArcId = getSuperArc(arc)->getReplacantArcId();
          // take care of connectivity
          mergeArc(arc, receptacleArcId);
          if(DEBUG){
            cout << " parallel merge in " << printArc(receptacleArcId) << endl;
            cout << " arc " << printArc(arc) << " size : " << getSuperArc(arc)->getVertSize()
                 << endl;
          }
          getSuperArc(receptacleArcId)
              ->addSegmentationGlobal(getSuperArc(arc)->getVertList(),
                                      getSuperArc(arc)->getVertSize());
       } else {
          const idNode &downNode = getSuperArc(arc)->getDownNodeId();
          const idNode &upNode   = getSuperArc(arc)->getUpNodeId();

          if (!(subtreeUF[downNode] && subtreeUF[upNode]))
             continue;

          if (subtreeUF[downNode] && subtreeUF[upNode] &&
              subtreeUF[downNode]->find() != subtreeUF[upNode]->find()) {
             if (DEBUG) {
                cout << "Arc between 2 degenerate with mergin " << printArc(arc) << endl;
                cout << "below recept : " << printArc(subtreeUF[downNode]->find()->getData());
                cout << endl;
                cout << "Above recept : " << printArc(subtreeUF[upNode]->find()->getData()) << endl;
                cout << endl;
             }

             continue;
          }

          ExtendedUnionFind *curUF =
              (subtreeUF[upNode]) ? subtreeUF[upNode]->find() : subtreeUF[downNode]->find();

          const idSuperArc &receptacleArcId = curUF->getData();

          if(DEBUG){
             cout << "merge in " << printArc(receptacleArcId) << endl;
             cout << "  arc " << printArc(arc) << " size : " << getSuperArc(arc)->getVertSize();
             cout << endl;
          }

          getSuperArc(arc)->merge(receptacleArcId);
          if (getSuperArc(arc)->getVertSize()) {
             getSuperArc(receptacleArcId)
                 ->addSegmentationGlobal(getSuperArc(arc)->getVertList(),
                                         getSuperArc(arc)->getVertSize());

             getSuperArc(receptacleArcId)->addSegmentationGlobal(getNode(downNode)->getVertexId());
             getSuperArc(receptacleArcId)->addSegmentationGlobal(getNode(upNode)->getVertexId());
          }

          // Tree topology
          getNode(downNode)->removeUpSuperArc(arc);
          getNode(downNode)->decUpValence();
          getNode(upNode)->removeDownSuperArc(arc);
          getNode(upNode)->decDownValence();

          if(!getNode(downNode)->getNumberOfUpSuperArcs())
              hideNode(downNode);
          if(!getNode(upNode)->getNumberOfDownSuperArcs())
              hideNode(upNode);

       }
   } // end for arcs

   //}

        return nbPairMerged;
}

// Persistence

template <typename scalarType>
int MergeTree::computePersistencePairs(vector<pair<pair<int, int>, double>> &pairs)
{
#ifndef withKamikaze
   if (!vect_superArcs_.size()) {
      return -1;
   }
#endif

   pairs.reserve(vect_leaves_.size());

   for (const idNode &leave : vect_leaves_) {
      Node *   curNode  = getNode(leave);
      idVertex curVert  = curNode->getVertexId();
      idVertex termVert = getNode(curNode->getTerminaison())->getVertexId();

      addPair<scalarType>(pairs,curVert,termVert);
   }

   auto pair_sort = [](const pair<pair<int, int>, double> &a,
                        const pair<pair<int, int>, double> &b) { return a.second < b.second; };

   sort(pairs.begin(), pairs.end(), pair_sort);

   return 0;
}

template <typename scalarType>
int MergeTree::computePersistencePairs(vector<tuple<idVertex, idVertex, scalarType, bool>> &pairs)
{
#ifndef withKamikaze
   if (!vect_superArcs_.size()) {
      return -1;
   }
#endif

   pairs.reserve(vect_leaves_.size());

   for (const idNode &leave : vect_leaves_) {
      Node *   curNode  = getNode(leave);
      idVertex curVert  = curNode->getVertexId();
      idVertex termVert = getNode(curNode->getTerminaison())->getVertexId();

      addPair<scalarType>(pairs, curVert, termVert, true);
   }

   auto pair_sort = [](const tuple<idVertex,idVertex,scalarType,bool> &a,
                        const tuple<idVertex,idVertex,scalarType,bool> &b) { return get<2>(a) < get<2>(b); };

   sort(pairs.begin(), pairs.end(), pair_sort);

   return 0;
}

template <typename scalarType>
void MergeTree::computePersistencePairsMT(const vector<idNode> &sortedNodes,
                                        vector<tuple<idVertex, idVertex, scalarType,bool>> &pairsJT,
                                        vector<tuple<idVertex, idVertex, scalarType,bool>> &pairsST)
{
    const auto nbNode = vect_nodes_.size();

    vector<ExtendedUnionFind *> vect_JoinUF(nbNode, nullptr);
    vector<ExtendedUnionFind *> vect_SplitUF(nbNode, nullptr);

#pragma omp parallel sections num_threads(2)
    {
#pragma omp section
       {
          pairsJT.reserve(vect_leaves_.size());
          // For the biggest pair of the component
          map<idVertex, idVertex> pendingMinMax;

          for (auto it = sortedNodes.cbegin(); it != sortedNodes.cend(); ++it) {
             const idNode &  n = *it;
             const idVertex &v = getNode(n)->getVertexId();

             const auto &nbUp = getNode(n)->getNumberOfUpSuperArcs();
             const auto &nbDown = getNode(n)->getNumberOfDownSuperArcs();

             if (nbDown == 0) {
                // leaf
                vect_JoinUF[n] = new ExtendedUnionFind(getNode(n)->getVertexId());
                vect_JoinUF[n]->setOrigin(v);
                //cout << " jt origin : " << v << endl;
             } else {
                // first descendant
                const idSuperArc   firstSaId = getNode(n)->getDownSuperArcId(0);
                const SuperArc *   firstSA   = getSuperArc(firstSaId);
                const idNode &firstChildNodeId = firstSA->getDownNodeId();

                ExtendedUnionFind *merge     = vect_JoinUF[firstChildNodeId]->find();
                idVertex           further   = merge->getOrigin();
                unsigned           furtherI  = 0;

                // Find the most persistant way
                for (unsigned ni = 1; ni < nbDown; ++ni) {
                   const idSuperArc curSaId = getNode(n)->getDownSuperArcId(ni);
                   const SuperArc * curSA   = getSuperArc(curSaId);
                   const idNode &   neigh   = curSA->getDownNodeId();

                   // fix
                   if (neigh == n)
                      continue;

                   ExtendedUnionFind *neighUF = vect_JoinUF[neigh]->find();

                   if (isLower(neighUF->getOrigin(), further)) {
                      further = neighUF->getOrigin();
                      furtherI = ni;
                   }
                }

                if (nbDown > 1) {
                   // close finish pair and make union
                   for (unsigned ni = 0; ni < nbDown; ++ni) {
                      const idSuperArc curSaId = getNode(n)->getDownSuperArcId(ni);
                      const SuperArc * curSA   = getSuperArc(curSaId);
                      const idNode &   neigh   = curSA->getDownNodeId();

                      // fix
                      if (neigh == n)
                         continue;

                      ExtendedUnionFind *neighUF = vect_JoinUF[neigh]->find();

                      if (ni != furtherI) {  // keep the more persitent pair
                         addPair<scalarType>(pairsJT, neighUF->getOrigin(), v, true);
                         pendingMinMax.erase(neighUF->getOrigin());

                         // cout << " jt make pair : " << neighUF->getOrigin() << " - " << v <<
                         // endl;
                      }

                      ExtendedUnionFind::makeUnion(merge, neighUF)->setOrigin(further);
                   }
                }

                merge->find()->setOrigin(further);
                vect_JoinUF[n] = merge->find();

                if(!nbUp){
                   // potential close of the component
                   //cout << "pending for " << further << " is " << v << endl;
                   pendingMinMax[further] = v;
                }
             }
          } // end for each node

          // Add the pending biggest pair of each component
          for (const auto & pair_vert : pendingMinMax) {
              if (isCorrespondingNode(pair_vert.second)) {
                 //cout << " add : " << pair_vert.first << " - " << pair_vert.second << endl;
                 addPair<scalarType>(pairsJT, pair_vert.first, pair_vert.second, true);
              }
          }
       } // end para section

#pragma omp section
       {
          pairsST.reserve(vect_leaves_.size());
          // For the biggest pair of the component
          map<idVertex, idVertex> pendingMinMax;

          for (auto it = sortedNodes.crbegin(); it != sortedNodes.crend(); ++it) {
             const idNode &  n = *it;
             const idVertex &v = getNode(n)->getVertexId();

             const auto &nbUp = getNode(n)->getNumberOfUpSuperArcs();
             const auto &nbDown = getNode(n)->getNumberOfDownSuperArcs();

             if (nbUp  == 0) {
                // leaf
                vect_SplitUF[n] = new ExtendedUnionFind(getNode(n)->getVertexId());
                vect_SplitUF[n]->setOrigin(v);
                //cout << " st origin : " << v << endl;
             } else {
                // first descendant
                const idSuperArc   firstSaId = getNode(n)->getUpSuperArcId(0);
                const SuperArc *   firstSA   = getSuperArc(firstSaId);
                const idNode &firstChildNodeId = firstSA->getUpNodeId();

                ExtendedUnionFind *merge     = vect_SplitUF[firstChildNodeId]->find();
                idVertex           further   = merge->getOrigin();
                unsigned           furtherI  = 0;

                for (unsigned ni = 1; ni < nbUp; ++ni) {
                    // find the more persistant way
                   const idSuperArc curSaId = getNode(n)->getUpSuperArcId(ni);
                   const SuperArc * curSA   = getSuperArc(curSaId);
                   // Ignore hidden / simplified arc
                   if(!curSA->isVisible()) continue;
                   const idNode &   neigh   = curSA->getUpNodeId();

                   // fix
                   if (neigh == n)
                      continue;

                   //cout << "visit neighbor : " << ni << " which is " << getNode(neigh)->getVertexId() << endl;

                   ExtendedUnionFind *neighUF = vect_SplitUF[neigh]->find();

                   if (isHigher(neighUF->getOrigin(), further)) {
                      further = neighUF->getOrigin();
                      furtherI = ni;
                   }
                }

                if (nbUp > 1) {
                   // close finsh pair and make union
                   for (unsigned ni = 0; ni < nbUp; ++ni) {
                      const idSuperArc curSaId = getNode(n)->getUpSuperArcId(ni);
                      const SuperArc * curSA   = getSuperArc(curSaId);
                      const idNode &   neigh   = curSA->getUpNodeId();

                      // fix
                      if (neigh == n)
                         continue;

                      ExtendedUnionFind *neighUF = vect_SplitUF[neigh]->find();

                      if(ni != furtherI){
                         addPair<scalarType>(pairsST, neighUF->getOrigin(), v, false);

                         pendingMinMax.erase(neighUF->getOrigin());

                         //cout << " st make pair : " << neighUF->getOrigin() << " - " << v
                              //<< " for neighbor " << getNode(neigh)->getVertexId() << endl;
                      }

                      ExtendedUnionFind::makeUnion(merge, neighUF)->setOrigin(further);
                      // Re-visit after merge lead to add the most persistant pair....
                   }
                }
                merge->find()->setOrigin(further);
                vect_SplitUF[n] = merge->find();

                if (!nbDown) {
                   pendingMinMax[further] = v;
                }

             } // end nbUp == 0 else
          } // end for each node

          // Add the pending biggest pair of each component
          for (const auto & pair_vert : pendingMinMax) {
             addPair<scalarType>(pairsST, pair_vert.first, pair_vert.second, false);
          }
       } // end para section
    } // end para
}


template <typename scalarType>
int MergeTree::computePersistencePlot(vector<pair<double, int>> &plot,
                                      vector<pair<pair<int, int>, double>> *persistencePairs)
{
   bool isLocal = false;
   if (!persistencePairs) {
      persistencePairs = new vector<pair<pair<int, int>, double>>;
      isLocal          = true;
   }

   if (!persistencePairs->size())
      computePersistencePairs<scalarType>(*persistencePairs);
   unsigned int nbElmnt = persistencePairs->size();
   plot.resize(nbElmnt);

   for (unsigned int i = 0; i < nbElmnt; ++i) {
      plot[i].first  = max((*persistencePairs)[i].second, pow(10, -REAL_SIGNIFICANT_DIGITS));
      plot[i].second = persistencePairs->size() - i;
   }

   if (isLocal)
      delete persistencePairs;
   return 0;
}

template <typename scalarType>
int MergeTree::computePersistenceDiagram(vector<pair<double, double>> &diagram,
                                         vector<pair<pair<int, int>, double>> *pairs)
{
   bool isLocal = false;
   if (!pairs) {
      pairs   = new vector<pair<pair<int, int>, double>>;
      isLocal = true;
   }

   if (!pairs->size())
      computePersistencePairs<scalarType>(*pairs);

   diagram.resize(pairs->size());

   for (int i = 0; i < (int)pairs->size(); i++) {
      if (isJT) {
         diagram[i].first  = ((scalarType *)scalars_)[(*pairs)[i].first.first];
         diagram[i].second = ((scalarType *)scalars_)[(*pairs)[i].first.second];
      } else {
         diagram[i].second = ((scalarType *)scalars_)[(*pairs)[i].first.first];
         diagram[i].first  = ((scalarType *)scalars_)[(*pairs)[i].first.second];
      }
   }

   auto pair_sort2 = [&](const pair<double, double> &a, const pair<double, double> &b) {
      return a.first < b.first;
   };

   sort(diagram.begin(), diagram.end(), pair_sort2);

   if (isLocal)
      delete pairs;

   return 0;
}

// }
// Tools
// {

template <typename scalarType>
void MergeTree::addPair(vector<tuple<idVertex, idVertex, scalarType, bool>> &pairs,
                        const idVertex &orig, const idVertex &term, const bool goUp)
{
   if (simplifyMethod_ == SimplifMethod::Persist) {
      pairs.emplace_back(
          orig, term, fabs(((scalarType *)scalars_)[orig] - ((scalarType *)scalars_)[term]), goUp);
   } else if (simplifyMethod_ == SimplifMethod::Span) {
      float coordOrig[3], coordTerm[3], span;
      mesh_->getVertexPoint(orig, coordOrig[0], coordOrig[1], coordOrig[2]);
      mesh_->getVertexPoint(term, coordTerm[0], coordTerm[1], coordTerm[2]);
      span = Geometry::distance(coordOrig, coordTerm);
      pairs.emplace_back(orig, term, span, goUp);
   }
}

template <typename scalarType>
void MergeTree::addPair(vector<pair<pair<int, int>, double>> &pairs,
                        const idVertex &orig, const idVertex &term)
{
   if (simplifyMethod_ == SimplifMethod::Persist) {
      pairs.emplace_back(reorderEdgeRel(make_pair(orig, term)),
                         fabs(((scalarType *)scalars_)[orig] - ((scalarType *)scalars_)[term]));
   } else if (simplifyMethod_ == SimplifMethod::Span) {
      float coordOrig[3], coordTerm[3], span;
      mesh_->getVertexPoint(orig, coordOrig[0], coordOrig[1], coordOrig[2]);
      mesh_->getVertexPoint(term, coordTerm[0], coordTerm[1], coordTerm[2]);
      span = Geometry::distance(coordOrig, coordTerm);
      pairs.emplace_back(reorderEdgeRel(make_pair(orig, term)), span);
   }

}

// }
/// ---------------------------- CT

// Process
// {

template <typename scalarType>
void ContourTree::initDataMT(void)
{
   jt_->setTriangulation(mesh_);
   st_->setTriangulation(mesh_);

   jt_->setPartition(partition_);
   st_->setPartition(partition_);

   jt_->setVertexScalars<scalarType>((scalarType *)scalars_);
   st_->setVertexScalars<scalarType>((scalarType *)scalars_);

   initSoS();
   jt_->setVertexSoSoffsets(soSOffsets_);
   st_->setVertexSoSoffsets(soSOffsets_);

   // initVert2Tree();
   // jt_->initVert2Tree();
   // st_->initVert2Tree();

   // parallel ct already do the sort
   // sortInput<scalarType>();
   jt_->setSorted(sortedVertices_);
   st_->setSorted(sortedVertices_);

   jt_->setMirror(mirrorOffsets_);
   st_->setMirror(mirrorOffsets_);

   jt_->setSegmentation(segmentation_);
   st_->setSegmentation(segmentation_);

   jt_->setComputeContourTree(computeContourTree_);
   st_->setComputeContourTree(computeContourTree_);

   jt_->setSimplificationMethod(simplifyMethod_);
   st_->setSimplificationMethod(simplifyMethod_);
}

// Persistance

template <typename scalarType>
int ContourTree::computePersistencePairs(vector<pair<pair<int, int>, double>> *pairs,
                                         vector<pair<pair<int, int>, double>> *mergePairs,
                                         vector<pair<pair<int, int>, double>> *splitPairs)
{
   if (pairs->size())
      return 0;

   bool isLocalJoin = false, isLocalSplit = false;

   if (!mergePairs) {
      mergePairs  = new vector<pair<pair<int, int>, double>>;
      isLocalJoin = true;
   }

   if (!splitPairs) {
      splitPairs   = new vector<pair<pair<int, int>, double>>;
      isLocalSplit = true;
   }

   if (!mergePairs->size())
      jt_->computePersistencePairs<scalarType>(*mergePairs);

   if (!splitPairs->size())
      jt_->computePersistencePairs<scalarType>(*splitPairs);

   int mpSize = mergePairs->size();
   int msSize = splitPairs->size();

   pairs->resize(mpSize + msSize);

   for (int i = 0; i < mpSize; i++) {
      (*pairs)[i] = (*mergePairs)[i];
   }

   for (int i = 0; i < msSize; i++) {
      (*pairs)[mpSize + i] = (*splitPairs)[i];
   }

   auto pair_sort = [&](const pair<pair<int, int>, double> &a,
                        const pair<pair<int, int>, double> &b) { return a.second < b.second; };

   sort(pairs->begin(), pairs->end(), pair_sort);

   //for (auto &p : *pairs) {
      //cout << p.first.first << "-" << p.first.second << "   ";
   //}
   //cout << endl;

   if (isLocalJoin)
      delete mergePairs;

   if (isLocalSplit)
      delete splitPairs;

   return 0;
}

template <typename scalarType>
int ContourTree::computePersistencePlot(vector<pair<double, int>> &plot,
                                        vector<pair<pair<int, int>, double>> *mergePairs,
                                        vector<pair<pair<int, int>, double>> *splitPairs,
                                        vector<pair<pair<int, int>, double>> *pairs)
{
   bool isLocal = false;
   if (!pairs) {
      pairs   = new vector<pair<pair<int, int>, double>>;
      isLocal = true;
   }

   if (!pairs->size()) {
      computePersistencePairs<scalarType>(pairs, mergePairs, splitPairs);
   }

   plot.resize(pairs->size());

   for (int i = 0; i < (int)plot.size(); i++) {
      plot[i].first  = max((*pairs)[i].second, pow(10, -REAL_SIGNIFICANT_DIGITS));
      plot[i].second = pairs->size() - i;
   }

   if (isLocal)
      delete pairs;

   return 0;
}

template <typename scalarType>
int ContourTree::computePersistenceDiagram(vector<pair<double, double>> &diagram,
                                           vector<pair<pair<int, int>, double>> *mergePairs,
                                           vector<pair<pair<int, int>, double>> *splitPairs,
                                           vector<pair<pair<int, int>, double>> *pairs)
{
   bool isLocal = false;
   if (!pairs) {
      pairs   = new vector<pair<pair<int, int>, double>>;
      isLocal = true;
   }

   if (!pairs->size())
      computePersistencePairs<scalarType>(pairs, mergePairs, splitPairs);

   diagram.resize(pairs->size());

   for (int i = 0; i < (int)pairs->size(); i++) {
      diagram[i].first  = ((scalarType *)scalars_)[(*pairs)[i].first.first];
      diagram[i].second = ((scalarType *)scalars_)[(*pairs)[i].first.second];
   }

   auto pair_sort2 = [&](const pair<double, double> &a, const pair<double, double> &b) {
      return a.first < b.first;
   };

   sort(diagram.begin(), diagram.end(), pair_sort2);

   if (isLocal)
      delete pairs;

   return 0;
}

// }

#endif /* end of include guard: CONTOURTREETEMPLATE_H */

