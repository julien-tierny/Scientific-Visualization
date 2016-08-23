/*
 * file: ContourTree.cpp
 * description: ContourTree processing package.
 * author: Gueunet Charles
 * date: Juin 2015
 */

#include <ContourTree.h>
#include <string>
#include<iterator>

// ------------------ Merge Tree

// Constructors & destructors
MergeTree::MergeTree(const unsigned char part) : MergeTree(false, part)
{
}

MergeTree::MergeTree(const bool t, const unsigned char part)
    : partition_{part},
      scalars_{nullptr},
      soSOffsets_{nullptr},
      sortedVertices_{},
      mirrorOffsets_{},
      isJT{t},
      vect_superArcs_{},
      vect_nodes_{},
      vect_leaves_{},
      vect_roots_{},
      vect_vert2tree_{nullptr},
      debugLevel_{3},
      simplifyMethod_{SimplifMethod::Persist}
{
}

MergeTree::~MergeTree()
{
   // The parallel contour Tree is in charge of the destruction of shared data !
   // if (destroyVectorSortedVertices_ && sortedVertices_) {
   // delete[] sortedVertices_;
   // sortedVertices_ = nullptr;
   //}

   // if (destroyVectorSoS_ && soSOffsets_) {
   // delete[] soSOffsets_;
   // soSOffsets_ = nullptr;
   //}

   // if (destroyVect2Tree_ && vect_vert2tree_) {
   // delete vect_vert2tree_;
   // vect_vert2tree_ = nullptr;
   //}
}

// Init
// {

void MergeTree::flush(void)
{
   initSoS();
   initVert2Tree();

   vect_superArcs_.clear();
   vect_nodes_.clear();
   vect_leaves_.clear();
}

void MergeTree::initSoS(void)
{
   if (!soSOffsets_) {
      soSOffsets_             = new idVertex[mesh_->getNumberOfVertices()];
      idVertex *soSOffsetsEnd = soSOffsets_ + mesh_->getNumberOfVertices();
      iota(soSOffsets_, soSOffsetsEnd, 0);
      destroyVectorSoS_ = true;
   }
}

void MergeTree::initVert2Tree(void)
{
   if (!vect_vert2tree_) {
      vect_vert2tree_ = new vector<idCorresp>(mesh_->getNumberOfVertices(), nullCorresp);
   } else {
      cout << "flush vert2tree : recompute" << endl;
      vect_vert2tree_->resize(mesh_->getNumberOfVertices());
      fill(vect_vert2tree_->begin(), vect_vert2tree_->end(), nullCorresp);
   }
}

//}

// Process
// {


int MergeTree::build(vector<ExtendedUnionFind *> &vect_baseUF,
                     const vector<idVertex> &overlapBefore, const vector<idVertex> &overlapAfter,
                     idVertex start, idVertex end, const idVertex &posSeed0,
                     const idVertex &posSeed1)
{
   // idea, work on the neighbohood instead of working on the node itsef.
   // Need lower / higher star construction.
   //
   // at this time, ST have no root except in Adjacency list. These root are isolated vertices.
   //
   // clear and reset tree data (this step should take almost no time)
   flush();

   DebugTimer begin;

   // -----------------
   // Find boundaries
   // -----------------
   // {

   idVertex sortedNode;
   const idVertex step = (isJT)?1:-1;

   // main
   const idVertex mainStart = start;
   const idVertex mainEnd   = end;

   // overlap before
   const idVertex beforeStart = (isJT) ? 0 : overlapBefore.size() - 1;
   const idVertex beforeEnd   = (isJT) ? overlapBefore.size() : -1;

   // overlap after
   const idVertex afterStart = (isJT) ? 0 : overlapAfter.size() - 1;
   const idVertex afterEnd   = (isJT) ? overlapAfter.size() : -1;

   // print debug
   if(debugLevel_ >= 3){
#pragma omp critical
      {
         cout << "partition : " << static_cast<unsigned>(partition_);
         cout << ", isJT : "    << isJT;
         cout << ",  size : ";
         cout << "before  : " << abs(beforeEnd - beforeStart);
         cout << ", main : "  << abs(mainEnd - mainStart);
         cout << ", after : " << abs(afterEnd - afterStart);
         cout << ", Total : ";
         cout << abs(beforeEnd - beforeStart) + abs(mainEnd - mainStart) +
                     abs(afterEnd - afterStart);
         cout << endl;
      }
   }

   // }
   // --------------
   // Overlap Before
   // --------------
   // {

   // for each vertex of our triangulation
   for (sortedNode = beforeStart; sortedNode != beforeEnd; sortedNode += step) {
      const idVertex currentVertex = overlapBefore[sortedNode];
      const bool overlapB = isJT;
      const bool overlapA = !isJT;
      processVertex(currentVertex, vect_baseUF, overlapB, overlapA, begin);
   }  // foreach node

   // }
   // ---------------
   // Partition
   // ---------------
   // {

   // for each vertex of our triangulation
   for (sortedNode = mainStart; sortedNode != mainEnd; sortedNode += step) {
      const idVertex currentVertex = sortedVertices_[sortedNode];
      processVertex(currentVertex, vect_baseUF, false, false, begin);
   }  // foreach node

   // }
   // ---------------
   // Overlap After
   // ---------------
   // {

   // for each vertex of our triangulation
   for (sortedNode = afterStart; sortedNode != afterEnd; sortedNode += step) {
      const idVertex currentVertex = overlapAfter[sortedNode];
      const bool overlapB = !isJT;
      const bool overlapA = isJT;
      processVertex(currentVertex, vect_baseUF, overlapB, overlapA, begin);
   }  // foreach node

   // }
   // ---------------
   // Close root arcs
   // ---------------
   // {

   // DebugTimer timeClose;

   // Closing step for openedSuperArc.
   // Not aesthetic but efficient
   // More efficient that using nullity of the neighborhood to detect extrema of the opponent tree
   idNode     rootNode;
   idVertex   corrVertex, origin;
   idSuperArc tmp_sa;

   // It can't be more connected component that leaves so test for each leaves (even virtual
   // extrema)
   for (const auto &l : vect_leaves_) {
      corrVertex = getNode(l)->getVertexId();

      // case of an isolated point TODO remove here and below
      if (!mesh_->getVertexNeighborNumber(corrVertex)) {
         tmp_sa = getNode(l)->getUpSuperArcId(0);
      } else {
         tmp_sa = (idSuperArc)((vect_baseUF[corrVertex])->find()->getData());
         origin = (idSuperArc)((vect_baseUF[corrVertex])->find()->getOrigin());
      }

      if (vect_superArcs_[tmp_sa].getUpNodeId() == nullNodes) {
         rootNode         = makeNode(vect_superArcs_[tmp_sa].getLastVisited(), origin);
         Node *originNode = getNode(getNode(getSuperArc(tmp_sa)->getDownNodeId())->getOrigin());
         originNode->setTerminaison(rootNode);

         const bool overlapB = mirrorOffsets_[getNode(rootNode)->getVertexId()] <= posSeed0;
         const bool overlapA = mirrorOffsets_[getNode(rootNode)->getVertexId()] >= posSeed1;

         closeSuperArc(tmp_sa, rootNode, overlapB, overlapA);

         // cout << getNode(rootNode)->getVertexId() << " - "
         //<< getNode(vect_superArcs_[tmp_sa].getDownNodeId())->getVertexId() << endl;

         // in the case we have 1 vertex domain,
         // hide the close SuperArc wich is point1 <>> point1
         if (getSuperArc(tmp_sa)->getDownNodeId() == getSuperArc(tmp_sa)->getUpNodeId()) {
            hideArc(tmp_sa);
         }

         vect_roots_.emplace_back(rootNode);
      }
   }

   // }
   // -----------
   // Timer print
   // ------------
   // {

   //if (debugLevel_ >= infoMsg) {
      //stringstream msgClose;

      //if (isJT)
         //msgClose << "[MergeTree] Join Tree ";
      //else
         //msgClose << "[MergeTree] Split Tree ";

      //msgClose << static_cast<unsigned>(partition_) << " closed in : " << timeClose.getElapsedTime()
               //<< endl;
      //dMsg(cout, msgClose.str(), infoMsg);
   //}

   if (debugLevel_ >= infoMsg) {
      stringstream msg;

      if (isJT)
         msg << "[MergeTree] Join  Tree ";
      else
         msg << "[MergeTree] Split Tree ";

      msg << static_cast<unsigned>(partition_) << " ";
      msg << "computed   in        " << begin.getElapsedTime();
      msg << "              \t( nb arcs : " << vect_superArcs_.size() << " )" << endl;
      dMsg(cout, msg.str(), infoMsg);
   }

   // }

   return 0;
}

void MergeTree::processVertex(const idVertex &             currentVertex,
                              vector<ExtendedUnionFind *> &vect_baseUF, const bool overlapB,
                              const bool overlapA, DebugTimer &begin)
{
   vector<ExtendedUnionFind *> vect_neighUF;
   ExtendedUnionFind *         seed = nullptr, *tmpseed;

   u_char    neighSize;
   const int neighborNumber = mesh_->getVertexNeighborNumber(currentVertex);

   idSuperArc currentArc;
   idNode     closingNode, currentNode;
   idVertex   neighbor;

   // Check UF in neighborhood
   for (int n = 0; n < neighborNumber; ++n) {
      mesh_->getVertexNeighbor(currentVertex, n, neighbor);
      // if the vertex is out: consider it null
      tmpseed = vect_baseUF[neighbor];
      // unvisited vertex, we continue.
      if (tmpseed == nullptr) {
         continue;
      }

      tmpseed = tmpseed->find();

      // get all different UF in neighborhood
      if (find(vect_neighUF.cbegin(), vect_neighUF.cend(), tmpseed) == vect_neighUF.end()) {
         vect_neighUF.emplace_back(tmpseed);
         seed = tmpseed;
      }
   }

   (neighSize = vect_neighUF.size());

   // idVertex test = 1;
   // if (currentVertex == test)
   // cout << test << " : " << vect_neighUF.size() << " " << vect_interfaceUF.size() << endl;
   // Make output
   if (!neighSize) {
      // we are on a real extrema we have to create a new UNION FIND and a branch
      // a real extrema can't be a virtual extrema

      seed = new ExtendedUnionFind(currentVertex);
      // When creating an extrema we create a pair ending on this node.
      currentNode = makeNode(currentVertex);
      getNode(currentNode)->setOrigin(currentNode);
      currentArc = openSuperArc(currentNode, overlapB, overlapA);
      // if(overlap && partition_ == 1) cout << currentVertex << endl;
      vect_leaves_.emplace_back(currentNode);

      if (debugLevel_ >= advancedInfoMsg) {
         stringstream msg;
         msg << "[MergeTree] ";

         if (isJT) {
            msg << "min node id:" << currentVertex;
         } else {
            msg << "max node id:" << currentVertex;
         }

         msg << " at : " << begin.getElapsedTime() << endl;
         dMsg(cout, msg.str(), advancedInfoMsg);
      }
   } else if (neighSize > 1) {
      // Is a saddle if have more than one UF in neighborhood
      // Or is linked to an interface (virtual extrema -> leaves or
      // interface cc representant : regular node )

      // Merge operation => close all arriving SuperArc (ignoring UF from interface)
      // And open a new one
      closingNode = makeNode(currentVertex);
      currentArc  = openSuperArc(closingNode, overlapB, overlapA);

      // TODO why use find() for neighboring UF ??
      idVertex farOrigin = vect_neighUF[0]->find()->getOrigin();

      // close each SuperArc finishing here
      for (auto *neigh : vect_neighUF) {
         closeSuperArc((idSuperArc)neigh->find()->getData(), closingNode, overlapB, overlapA);
         // persistance pair closing here.
         // For the one who will continue, it will be overide later
         getNode(getCorrespondingNode(neigh->find()->getOrigin()))->setTerminaison(closingNode);

         // cout << getNode(getCorrespondingNode(neigh->find()->getOrigin()))->getVertexId()
         //<< " terminate on " << getNode(closingNode)->getVertexId() << endl;

         if ((isJT && isLower(neigh->find()->getOrigin(), farOrigin)) ||
             (!isJT && isHigher(neigh->find()->getOrigin(), farOrigin))) {
            // here we keep the continuing the most persitant pair.
            // It means a pair end when a parent have another origin thant the current leaf (or
            // is the root)
            // It might be not intuitive but it is more convenient for degenerate cases
            farOrigin = neigh->find()->getOrigin();
            // cout << "find origin  " << farOrigin << " for " << currentVertex << "   " << isJT
            //<< endl;
         }
      }

      // Union correspond to the merge
      seed = ExtendedUnionFind::makeUnion(vect_neighUF);
      seed->setOrigin(farOrigin);
      getNode(closingNode)->setOrigin(getCorrespondingNode(farOrigin));

      // cout << "  " << getNode(closingNode)->getVertexId() << " have origin at "
      //<< getNode(getCorrespondingNode(farOrigin))->getVertexId() << endl;

      if (debugLevel_ >= advancedInfoMsg) {
         stringstream msg;

         msg << "[MergeTree] ";
         msg << "saddle node id:" << currentVertex;
         dMsg(cout, msg.str(), advancedInfoMsg);
      }

   } else {
      // regular node
      currentArc = (idSuperArc)seed->find()->getData();
      updateCorrespondingArc(currentVertex, currentArc);
   }
   // common
   seed->setData((ufDataType)currentArc);
   getSuperArc(currentArc)->setLastVisited(currentVertex, segmentation_);
   vect_baseUF[currentVertex] = seed;
}

// update lately

void MergeTree::updateSegmentation(const bool ct)
{
   auto compL = [&](const pair<idVertex, bool> &a, const pair<idVertex, bool> &b) {
      return isLower(a.first, b.first);
   };

   auto compH = [&](const pair<idVertex, bool> &a, const pair<idVertex, bool> &b) {
      return isHigher(a.first, b.first);
   };

   const idSuperArc nbArc = vect_superArcs_.size();
   if (isJT || ct) {
      for (idSuperArc sa = 0; sa < nbArc; sa++) {
         SuperArc *superArc     = getSuperArc(sa);
         if(!superArc->isVisible()) continue;

         auto *segmentation = superArc->getVertList();
         auto &segmSize     = superArc->getVertSize();

         //fix
         if(!segmentation)
             continue;

         sort(segmentation, segmentation + segmSize, compL);
         for (idVertex i = 0; i < segmSize; i++) {
            const idVertex &vert = segmentation[i].first;
            if (!segmentation[i].second) {
               updateCorrespondingArc(vert, sa);
            }
         }
      }
   } else {
      for (idSuperArc sa = 0; sa < nbArc; sa++) {
         SuperArc *superArc     = getSuperArc(sa);
         if(!superArc->isVisible()) continue;

         auto *    segmentation = superArc->getVertList();
         auto &    segmSize     = superArc->getVertSize();

         //fix
         if(!segmentation)
             continue;

         sort(segmentation, segmentation + segmSize, compH);
         for (idVertex i = 0; i < segmSize; i++) {
            const idVertex &vert = segmentation[i].first;
            updateCorrespondingArc(vert, sa);
         }
      }
   }

   const idNode nbNode = vect_nodes_.size();
   for (idNode n = 0; n < nbNode; n++) {
       if(!getNode(n)->isHidden()){
          updateCorrespondingNode(getNode(n)->getVertexId(), n);
       }
   }
}

void MergeTree::parallelUpdateSegmentation(const bool ct)
{
   auto compL = [&](const pair<idVertex, bool> &a, const pair<idVertex, bool> &b) {
      return isLower(a.first, b.first);
   };

   auto compH = [&](const pair<idVertex, bool> &a, const pair<idVertex, bool> &b) {
      return isHigher(a.first, b.first);
   };

   const idSuperArc nbArc = vect_superArcs_.size();
   if (isJT || ct) {
#pragma omp parallel for
      for (idSuperArc sa = 0; sa < nbArc; sa++) {
         SuperArc *superArc     = getSuperArc(sa);
         if(!superArc->isVisible()) continue;

         auto *segmentation = superArc->getVertList();
         auto &segmSize     = superArc->getVertSize();

         //fix
         if(!segmentation)
             continue;

         sort(segmentation, segmentation + segmSize, compL);
         for (idVertex i = 0; i < segmSize; i++) {
            const idVertex &vert = segmentation[i].first;
            if (!segmentation[i].second) {
               updateCorrespondingArc(vert, sa);
            }
         }
      }
   } else {
#pragma omp parallel for
      for (idSuperArc sa = 0; sa < nbArc; sa++) {
         SuperArc *superArc     = getSuperArc(sa);
         if(!superArc->isVisible()) continue;

         auto *    segmentation = superArc->getVertList();
         auto &    segmSize     = superArc->getVertSize();

         //fix
         if(!segmentation)
             continue;

         sort(segmentation, segmentation + segmSize, compH);
         for (idVertex i = 0; i < segmSize; i++) {
            const idVertex &vert = segmentation[i].first;
            updateCorrespondingArc(vert, sa);
         }
      }
   }

   const idNode nbNode = vect_nodes_.size();
   for (idNode n = 0; n < nbNode; n++) {
       if(!getNode(n)->isHidden()){
          updateCorrespondingNode(getNode(n)->getVertexId(), n);
       }
   }
}

void MergeTree::parallelInitNodeValence(const int nbThreadValence)
{
   //cout << "SENTINEL : Parallel Init Node Valence " << endl;
   const auto &nbNodes = vect_nodes_.size();

#pragma omp parallel for num_threads(nbThreadValence)
   for (idNode n = 0; n < nbNodes; n++) {
      short downVal = 0, upVal = 0;
      Node *node = getNode(n);

      const auto nbDown = node->getNumberOfDownSuperArcs();
      for (unsigned i = 0; i < nbDown; i++) {
         const SuperArc *sa = getSuperArc(node->getDownSuperArcId(i));
         if (sa->isVisible()) {
            ++downVal;
         }
      }

      const auto nbUp = node->getNumberOfUpSuperArcs();
      for (unsigned i = 0; i < nbUp; i++) {
         const SuperArc *sa = getSuperArc(node->getUpSuperArcId(i));
         if (sa->isVisible()) {
            ++upVal;
         }
      }

      node->setDownValence(downVal);
      node->setUpValence(upVal);
      //if(getNode(n)->getVertexId() == 2){
        //cout << "toto : " << static_cast<unsigned>(partition_) << " : " << downVal << endl << endl;
        //for (unsigned i = 0; i < nbDown; i++) {
           //cout << printArc ( node->getDownSuperArcId(i) ) << " id : " << static_cast<unsigned>(node->getDownSuperArcId(i));
           //cout << " " << static_cast<unsigned>(getSuperArc(node->getDownSuperArcId(i))->getDownCT());
           //cout << " - " << static_cast<unsigned>(getSuperArc(node->getDownSuperArcId(i))->getUpCT());
           //cout << endl;
        //}
      //}
   }
}


void MergeTree::initNodeValence()
{
   cout << "SENTINEL : Init Node Valence " << endl;
   // This function is called before simplification, nothing is hidden yet
   const auto &nbNodes = vect_nodes_.size();

   for (idNode n = 0; n < nbNodes; n++) {
      Node *node = getNode(n);

      const auto nbDown = node->getNumberOfDownSuperArcs();
      const auto nbUp = node->getNumberOfUpSuperArcs();

      node->setDownValence(nbDown);
      node->setUpValence(nbUp);
   }
}


// }
// Arcs and node manipulations
// {

// SuperArcs
idSuperArc MergeTree::openSuperArc(const idNode &downNodeId, const bool overlapB,
                                   const bool overlapA)
{
#ifndef withKamikaze
   if (downNodeId < 0 || (size_t)downNodeId >= vect_nodes_.size()) {
      cout << "[Merge Tree] openSuperArc on a inexisting node !" << endl;
      return -2;
   }
#endif

   idSuperArc newSuperArcId = (idSuperArc)vect_superArcs_.size();
   vect_superArcs_.emplace_back(downNodeId, nullNodes, overlapB, overlapA, partition_, partition_);
   vect_nodes_[downNodeId].addUpSuperArcId(newSuperArcId);
   vect_nodes_[downNodeId].incUpValence();

   return newSuperArcId;
}

idSuperArc MergeTree::makeSuperArc(const idNode &downNodeId, const idNode &upNodeId,
                                   const bool overlapB, const bool overlapA,
                                   pair<idVertex, bool> *vertexList, int vertexSize)
{
   idSuperArc newSuperArcId = (idSuperArc)vect_superArcs_.size();

   if (downNodeId != upNodeId)
      vect_superArcs_.emplace_back(downNodeId, upNodeId, overlapB, overlapA, partition_,
                                   partition_);
   else
      vect_superArcs_.emplace_back(downNodeId, upNodeId, overlapB, overlapA, partition_, partition_,
                                   ComponentState::HIDDEN);

   vect_superArcs_[newSuperArcId].setVertList(vertexList);
   vect_superArcs_[newSuperArcId].setVertSize(vertexSize);

   vect_nodes_[downNodeId].addUpSuperArcId(newSuperArcId);
   vect_nodes_[downNodeId].incUpValence();
   vect_nodes_[upNodeId].addDownSuperArcId(newSuperArcId);
   vect_nodes_[upNodeId].incDownValence();

   return newSuperArcId;
}

void MergeTree::closeSuperArc(const idSuperArc &superArcId, const idNode &upNodeId,
                              const bool overlapB, const bool overlapA)
{
#ifndef withKamikaze

   if (superArcId < 0 || (size_t)superArcId >= vect_superArcs_.size()) {
      cout << "[Merge Tree] closeSuperArc on a inexisting arc !" << endl;
      return;
   }

   if (upNodeId < 0 || (size_t)upNodeId >= vect_nodes_.size()) {
      cout << "[Merge Tree] closeOpenedArc on a inexisting node !" << endl;
      return;
   }

#endif
   vect_superArcs_[superArcId].setUpNodeId(upNodeId);
   // TODO why do we need to re-set last-visited ? (maybe Saddle, check it out)
   vect_superArcs_[superArcId].setLastVisited(getNode(upNodeId)->getVertexId(), false);
   vect_nodes_[upNodeId].addDownSuperArcId(superArcId);
   vect_nodes_[upNodeId].incDownValence();

   vect_superArcs_[superArcId].setOverlapBelow(vect_superArcs_[superArcId].getOverlapBelow() !=
                                               overlapB);
   if (vect_superArcs_[superArcId].getOverlapBelow()) {
      vect_arcsCrossingBelow_.emplace_back(superArcId);
   }

   vect_superArcs_[superArcId].setOverlapAbove(vect_superArcs_[superArcId].getOverlapAbove() !=
                                               overlapA);
   if (vect_superArcs_[superArcId].getOverlapAbove()) {
      vect_arcsCrossingAbove_.emplace_back(superArcId);
   }
}

const idVertex MergeTree::cutArcAboveSeed(const idSuperArc &arc, const pair<idVertex, bool> &seed)
{
   auto isLowerComp = [&](const pair<idVertex, bool> &a, const pair<idVertex, bool> &b) {
      return isLower(a.first, b.first);
   };

   SuperArc *crossing = getSuperArc(arc);
   idVertex  stitchVert;
   // get stitching node
   const auto &vertList = crossing->getVertList();
   const auto &vertSize = crossing->getVertSize();

   if (vertSize) {
      auto posVert = lower_bound(vertList, vertList + vertSize, seed, isLowerComp);
      if (posVert == vertList + vertSize) {
         stitchVert = getNode(crossing->getUpNodeId())->getVertexId();
      } else {
         while (posVert < vertList + vertSize && posVert->second) {
            ++posVert;
         }

         if (posVert == vertList + vertSize) {
            stitchVert = getNode(crossing->getUpNodeId())->getVertexId();
         } else {
             // update segmentation
            crossing->setVertSize(posVert - vertList);
            stitchVert = posVert->first;
            // cout << "new size is " << crossing->getVertSize() << endl;
         }
      }
   } else {
      stitchVert = getNode(crossing->getUpNodeId())->getVertexId();
   }

   return stitchVert;
}

const idVertex MergeTree::cutArcBelowSeed(const idSuperArc &arc, const pair<idVertex, bool> &seed, const vector<idCorresp>* vert2treeOther)
{
   auto isLowerComp = [&](const pair<idVertex, bool> &a, const pair<idVertex, bool> &b) {
      return isLower(a.first, b.first);
   };

   SuperArc *crossing = getSuperArc(arc);
   idVertex  stitchVert;
   // get stitching node
   const auto &vertList = crossing->getVertList();
   const auto &vertSize = crossing->getVertSize();

   if (vertSize) {
      auto posVert = lower_bound(vertList, vertList + vertSize, seed, isLowerComp);
      if (posVert == vertList) {
         stitchVert = getNode(crossing->getDownNodeId())->getVertexId();
      } else {
         --posVert;  // we want below and not hidden (TODO remove nullCorresp ?)
         while (posVert > vertList &&
                (posVert->second || (*vert2treeOther)[posVert->first] == nullCorresp)) {
            --posVert;
          }

          if (posVert == 0) {
             stitchVert = getNode(crossing->getDownNodeId())->getVertexId();
         } else {
             // we only find, do not touch the segmentation
            stitchVert = posVert->first;
         }
      }
   } else {
      stitchVert = getNode(crossing->getDownNodeId())->getVertexId();
   }

   return stitchVert;
}

void MergeTree::removeHiddenDownArcs(const idNode &n)
{

    Node * node = getNode(n);

    const auto nbDown = node->getNumberOfDownSuperArcs();
    for (unsigned i = 0; i < nbDown; i++) {
       const SuperArc *sa = getSuperArc(node->getDownSuperArcId(i));
       if(!sa->isVisible()){
        node->removeDownSuperArcId(i);
       }
    }

}

// test

unsigned MergeTree::getNumberOfVisibleArcs(const idNode &n)
{

    Node * node = getNode(n);

    const auto nbDown = node->getNumberOfDownSuperArcs();
    const auto nbUp = node->getNumberOfUpSuperArcs();
    unsigned res = 0;
    for (unsigned i = 0; i < nbDown; i++) {
       const SuperArc *sa = getSuperArc(node->getDownSuperArcId(i));
       if(sa->isVisible()){
           ++res;
       }
    }
    for (unsigned i = 0; i < nbUp; i++) {
       const SuperArc *sa = getSuperArc(node->getUpSuperArcId(i));
       if(sa->isVisible()){
           ++res;
       }
    }

    return res;
}

unsigned MergeTree::getNumberOfUnmergedDownArcs(const idNode &n)
{

    Node * node = getNode(n);

    const auto nbDown = node->getNumberOfDownSuperArcs();
    unsigned res = 0;
    for (unsigned i = 0; i < nbDown; i++) {
       const SuperArc *sa = getSuperArc(node->getDownSuperArcId(i));
       if(!sa->isMerged()){
           ++res;
       }
    }

    //cout << "partition : " << static_cast<unsigned>(partition_);
    //cout << " vertex " << node->getVertexId() << " have " << res << " real downs";
    //cout << " on " << static_cast<unsigned>(nbDown) << endl;

    return res;
}

const bool MergeTree::alreadyExtLinked(const idNode &node, const idPartition &tree, const idNode &treeNode)
{
    Node* n = getNode(node);
    const auto nbUp = n->getNumberOfUpSuperArcs();


    for (unsigned a = 0; a < nbUp; a++) {
        const idSuperArc sa = n->getUpSuperArcId(a);

        if (getSuperArc(sa)->getUpCT() == tree && getSuperArc(sa)->getUpNodeId() == treeNode) {
           return true;
        }
    }

    return false;
}

// state

void MergeTree::hideArc(const idSuperArc &sa)
{
   vect_superArcs_[sa].hide();
   vect_nodes_[vect_superArcs_[sa].getUpNodeId()].removeDownSuperArc(sa);
   vect_nodes_[vect_superArcs_[sa].getUpNodeId()].decDownValence();
   vect_nodes_[vect_superArcs_[sa].getDownNodeId()].removeUpSuperArc(sa);
   vect_nodes_[vect_superArcs_[sa].getDownNodeId()].decUpValence();
}

void MergeTree::mergeArc(const idSuperArc &sa, const idSuperArc &recept,
                         const bool changeConnectivity)
{
   vect_superArcs_[sa].merge(recept);

   if (changeConnectivity) {
      vect_nodes_[vect_superArcs_[sa].getUpNodeId()].removeDownSuperArc(sa);
      vect_nodes_[vect_superArcs_[sa].getUpNodeId()].decDownValence();
      vect_nodes_[vect_superArcs_[sa].getDownNodeId()].removeUpSuperArc(sa);
      vect_nodes_[vect_superArcs_[sa].getDownNodeId()].decUpValence();
   }
}

void MergeTree::hideNode(const idNode &node)
{
    vect_nodes_[node].hide();
}

// Nodes

idNode MergeTree::makeNode(const idVertex &vertexId, const idVertex &term)
{
#ifndef withKamikaze
   if (vertexId < 0 || vertexId >= mesh_->getNumberOfVertices()) {
      cout << "[Merge Tree] make node, wrong vertex :" << vertexId << " on "
           << mesh_->getNumberOfVertices() << endl;
      return -1;
   }
#endif

   if (isCorrespondingNode(vertexId)) {
      return getCorrespondingNode(vertexId);
   }

   int size_base = (int)vect_nodes_.size();
   vect_nodes_.emplace_back(vertexId, term);
   updateCorrespondingNode(vertexId, size_base);

   return size_base;
}

idNode MergeTree::makeNode(const Node *const n, const idVertex &term)
{
   const idVertex &vertexId = n->getVertexId();

   if (isCorrespondingNode(vertexId)) {
      // Node already exist
      return getCorrespondingNode(vertexId);
   }

   int size_base = (int)vect_nodes_.size();
   vect_nodes_.emplace_back(vertexId, term);
   updateCorrespondingNode(vertexId, size_base);

   return size_base;
}

void MergeTree::delNode(const idNode &node, const pair<idVertex, bool> *markVertices,
                        const int &nbMark)
{
   Node *mainNode = getNode(node);

   if (mainNode->getNumberOfUpSuperArcs() == 0) {
      // Root
      if (mainNode->getNumberOfDownSuperArcs() != 1) {
         cout << endl << "[MergeTree]:delNode won't delete ";
         cout << mainNode->getVertexId() << " (root) with ";
         cout << static_cast<unsigned>(mainNode->getNumberOfDownSuperArcs()) << " down ";
         cout << static_cast<unsigned>(mainNode->getNumberOfUpSuperArcs()) << " up ";
         cout << " partition : " << static_cast<unsigned>(partition_) << endl;
         return;
      }

      idSuperArc downArc  = mainNode->getDownSuperArcId(0);
      Node *     downNode = getNode(vect_superArcs_[downArc].getDownNodeId());
      downNode->removeUpSuperArc(downArc);
      mainNode->clearDownSuperArcs();
      // USELESS ??? TODO remove ?? HIDE
      vect_superArcs_[downArc].hide();
      // vect_superArcs_[downArc].setUpNodeId(0);
      // vect_superArcs_[downArc].setDownNodeId(0);

   } else if (mainNode->getNumberOfDownSuperArcs() < 2) {
      idSuperArc upArc  = mainNode->getUpSuperArcId(0);
      idNode     upId   = vect_superArcs_[upArc].getUpNodeId();
      Node *     upNode = getNode(upId);

      // Have a parent
      upNode->removeDownSuperArc(upArc);
      // vect_superArcs_[upArc].setUpNodeId(0);
      // vect_superArcs_[upArc].setDownNodeId(0);
      vect_superArcs_[upArc].hide();
      mainNode->clearUpSuperArcs();
      // Have child(s)
      // Should be 0 or 1, verify
      if (mainNode->getNumberOfDownSuperArcs() == 1) {
         const idSuperArc &downArc = mainNode->getDownSuperArcId(0);

         // IN CASE OF SEGM
         if (markVertices != nullptr) {
            // if contiguous vert list
            if ((vect_superArcs_[downArc].getVertList() + vect_superArcs_[downArc].getVertSize()) ==
                vect_superArcs_[upArc].getVertList()) {
               // down arc <- union of the two list
               vect_superArcs_[downArc].setVertSize(vect_superArcs_[downArc].getVertSize() +
                                                    vect_superArcs_[upArc].getVertSize());
               // mark removed ones (passed as parameter of this function)
               // the markVertices is sorted in reverse order as it come form the other tree
               idVertex acc = -1;
               for (idVertex i = nbMark - 1; i >= 0; --i) {
                  if (!vect_superArcs_[downArc].getVertSize())
                     break;
                  while (vect_superArcs_[downArc].getRegularNodeId(++acc) !=
                         markVertices[i].first) {
                     if (acc == vect_superArcs_[downArc].getVertSize())
                        break;
                  }
                  if (acc == vect_superArcs_[downArc].getVertSize())
                     break;

                  vect_superArcs_[downArc].setMasqued(acc);
               }

            } else {
               // manually concatene segmentation of both arc
               // and put it in the bottom arc.
               const auto &upSize   = vect_superArcs_[upArc].getVertSize();
               const auto &downSize = vect_superArcs_[downArc].getVertSize();

               const auto *upSegm   = vect_superArcs_[upArc].getVertList();
               const auto *downSegm = vect_superArcs_[downArc].getVertList();

               pair<idVertex, bool>* newSegmentation = new pair<idVertex, bool>[upSize+downSize];

               for (idVertex i = 0; i < downSize; i++) {
                  newSegmentation[i] = downSegm[i];
               }

               for (idVertex i = 0; i < upSize; i++) {
                  newSegmentation[i + downSize] = upSegm[i];
               }

               vect_superArcs_[downArc].setVertList(newSegmentation);
               vect_superArcs_[downArc].setVertSize(downSize+upSize);
            }
         }

         vect_superArcs_[downArc].setUpNodeId(upId);
         vect_superArcs_[upArc].setDownNodeId(getSuperArc(downArc)->getDownNodeId());
         upNode->addDownSuperArcId(downArc);
      }

      mainNode->clearDownSuperArcs();
   }
}

// Normal insert : existing arc stay below inserted (JT example)
//  *   - <- upNodeId
//  | \ |   <- newSA
//  |   * <- newNodeId
//  |   |   <- currentSA
//  - - -
const idSuperArc MergeTree::insertNode(Node *node, const bool segment, const bool isCT)
{
   // already present
   if (isCorrespondingNode(node->getVertexId())) {
      Node *myNode = getNode(getCorrespondingNode(node->getVertexId()));
      // If it has been hidden / replaced we need to re-make it
      if (myNode->isHidden()) {
         SuperArc *sa = getSuperArc(myNode->getUpSuperArcId(0));
         const idSuperArc &correspondingArcId = (sa->getReplacantArcId() == nullSuperArc)
                                                     ? myNode->getUpSuperArcId(0)
                                                     : sa->getReplacantArcId();
         updateCorrespondingArc(myNode->getVertexId(), correspondingArcId);
         //cout << "update in insert : " << myNode->getVertexId() << " is on " << correspondingArcId << endl;
      } else
         return nullSuperArc;
   }

   idNode     upNodeId, newNodeId;
   idSuperArc currentSA, newSA;
   idVertex   origin;

   // Create new node
   currentSA = getCorrespondingSuperArcId(node->getVertexId());
   upNodeId  = vect_superArcs_[currentSA].getUpNodeId();
   origin    = vect_nodes_[vect_superArcs_[currentSA].getDownNodeId()].getOrigin();
   newNodeId = makeNode(node, origin);

   // Connectivity
   // TODO use makeSuperArc for the newArc
   // Insert only node inside the partition : created arc don t cross
   newSA = openSuperArc(newNodeId, false, false);

   vect_superArcs_[newSA].setUpNodeId(upNodeId);
   vect_nodes_[upNodeId].removeDownSuperArc(currentSA);
   vect_nodes_[upNodeId].addDownSuperArcId(newSA);

   vect_superArcs_[currentSA].setUpNodeId(newNodeId);
   vect_nodes_[newNodeId].addDownSuperArcId(currentSA);

   if (segment) {
      // cut the vertex list at the node position and
      // give each arc its part.
      SuperArc *tmpSA = getSuperArc(currentSA);
      pair<idVertex, bool> *newNodePosPtr =
          (isJT || isCT)
              ? lower_bound(tmpSA->getVertList(), tmpSA->getVertList() + tmpSA->getVertSize(),
                            make_pair(node->getVertexId(), false),
                            [&](const pair<idVertex, bool> &a, const pair<idVertex, bool> &b) {
                               return isLower(a.first, b.first);
                            })
              : lower_bound(tmpSA->getVertList(), tmpSA->getVertList() + tmpSA->getVertSize(),
                            make_pair(node->getVertexId(), false),
                            [&](const pair<idVertex, bool> &a, const pair<idVertex, bool> &b) {
                               return isHigher(a.first, b.first);
                            });

      unsigned newNodePos = newNodePosPtr - tmpSA->getVertList();

      getSuperArc(newSA)->setVertList(newNodePosPtr);
      getSuperArc(newSA)->setVertSize(tmpSA->getVertSize() - newNodePos);

      tmpSA->setVertSize(newNodePos);
   }

   return newSA;
}

// Reverse insert : existing arc stay above inserted (JT example)
//  *   - <- upNodeId
//  |   |   <- currentSa
//  |   * <- newNodeId
//  | / |   <- newSA
//  - - -
const idSuperArc MergeTree::reverseInsertNode(Node *node, const bool segment, const bool isCT)
{
   // already present
   if (isCorrespondingNode(node->getVertexId())) {
      Node *myNode = getNode(getCorrespondingNode(node->getVertexId()));
      // If it has been hidden / replaced we need to re-make it
      if (myNode->isHidden()) {
         cout << "reverse insert don t  deal with hidden" << endl;
      } else
         return nullSuperArc;
   }

   idNode     downNodeId, newNodeId;
   idSuperArc currentSA, newSA;
   idVertex   origin;

   // Create new node
   currentSA  = getCorrespondingSuperArcId(node->getVertexId());
   downNodeId = vect_superArcs_[currentSA].getDownNodeId();
   origin    = vect_nodes_[vect_superArcs_[currentSA].getDownNodeId()].getOrigin();
   newNodeId = makeNode(node, origin);

   // Connectivity
   // TODO use makeSuperArc for the newArc
   // Insert only node inside the partition : created arc don t cross
   newSA = openSuperArc(downNodeId, false, false);

   vect_superArcs_[newSA].setUpNodeId(newNodeId);
   vect_nodes_[downNodeId].removeUpSuperArc(currentSA);
   vect_nodes_[downNodeId].addUpSuperArcId(newSA);

   vect_superArcs_[currentSA].setDownNodeId(newNodeId);
   vect_nodes_[newNodeId].addUpSuperArcId(currentSA);

   vect_nodes_[newNodeId].addDownSuperArcId(newSA);

   if (segment) {
      // cut the vertex list at the node position and
      // give each arc its part.
      SuperArc *tmpSA = getSuperArc(currentSA);
      pair<idVertex, bool> *newNodePosPtr =
          (isJT || isCT)
              ? lower_bound(tmpSA->getVertList(), tmpSA->getVertList() + tmpSA->getVertSize(),
                            make_pair(node->getVertexId(), false),
                            [&](const pair<idVertex, bool> &a, const pair<idVertex, bool> &b) {
                               return isLower(a.first, b.first);
                            })
              : lower_bound(tmpSA->getVertList(), tmpSA->getVertList() + tmpSA->getVertSize(),
                            make_pair(node->getVertexId(), false),
                            [&](const pair<idVertex, bool> &a, const pair<idVertex, bool> &b) {
                               return isHigher(a.first, b.first);
                            });

      unsigned newNodePos = newNodePosPtr - tmpSA->getVertList();

      getSuperArc(newSA)->setVertList(tmpSA->getVertList());
      getSuperArc(newSA)->setVertSize(newNodePos);

      tmpSA->setVertList(newNodePosPtr);
      tmpSA->setVertSize(tmpSA->getVertSize() - newNodePos);
   }

   return newSA;
}


// traverse

inline Node *MergeTree::getDownNode(const SuperArc *a)
{
   return &(vect_nodes_[a->getDownNodeId()]);
}

inline Node *MergeTree::getUpNode(const SuperArc *a)
{
   return &(vect_nodes_[a->getUpNodeId()]);
}

idNode MergeTree::getParent(const idNode &n)
{
   return getSuperArc(getNode(n)->getUpSuperArcId(0))->getUpNodeId();
}

// Here the return of the vector use the move constructor
const vector<idNode> MergeTree::getNodeNeighbors(const idNode &n)
{
   Node *node   = getNode(n);
   auto  nbUp   = node->getNumberOfUpSuperArcs();
   auto  nbDown = node->getNumberOfDownSuperArcs();

   vector<idNode> res(nbDown + nbUp);

   unsigned currentPos = 0;

   // Nodes below
   for (int i = 0; i < nbDown; i++) {
      const idSuperArc &corArc  = node->getDownSuperArcId(i);
      const idNode &    corNode = getSuperArc(corArc)->getDownNodeId();
      res[currentPos++]         = corNode;
   }

   // Nodes above
   for (int i = 0; i < nbUp; i++) {
      const idSuperArc &corArc  = node->getUpSuperArcId(i);
      const idNode &    corNode = getSuperArc(corArc)->getUpNodeId();
      res[currentPos++]         = corNode;
   }

   return move(res);
}

const vector<idNode> MergeTree::getNodeUpNeighbors(const idNode &n)
{
   Node *node   = getNode(n);
   auto  nbUp   = node->getNumberOfUpSuperArcs();

   vector<idNode> res(nbUp);

   // Nodes above
   for (int i = 0; i < nbUp; i++) {
      const idSuperArc &corArc  = node->getUpSuperArcId(i);
      const idNode &    corNode = getSuperArc(corArc)->getUpNodeId();
      res[i]                    = corNode;
   }

   return move(res);
}

const vector<idNode> MergeTree::getNodeDownNeighbors(const idNode &n)
{
   Node *node   = getNode(n);
   auto  nbDown   = node->getNumberOfDownSuperArcs();

   vector<idNode> res(nbDown);

   // Nodes above
   for (int i = 0; i < nbDown; i++) {
      const idSuperArc &corArc  = node->getDownSuperArcId(i);
      const idNode &    corNode = getSuperArc(corArc)->getDownNodeId();
      res[i]                    = corNode;
   }

   return move(res);
}

// hide / clear

void MergeTree::hideAndClearArcsAbove(const idNode &baseNode)
{
   const idSuperArc nbArc = getNode(baseNode)->getNumberOfUpSuperArcs();
   for (unsigned i = 0; i < nbArc; ++i) {
      const idSuperArc &upsaid  = getNode(baseNode)->getUpSuperArcId(i);
      SuperArc *curArc  = getSuperArc(upsaid);
      if(curArc->getUpCT() != partition_) continue;
      // Node *    topNode = getNode(curArc->getUpNodeId());

      hideArc(upsaid);
      hideNode(curArc->getUpNodeId());

      //cout << " above hide vertex " << topNode->getVertexId() << " in "
           //<< static_cast<unsigned>(partition_) << endl;
      //cout << " above hide arc  " << printArc(getNode(baseNode)->getUpSuperArcId(i)) << endl;
   }

   getNode(baseNode)->clearUpSuperArcs();
}

void MergeTree::hideAndClearArcsBelow(const idNode &baseNode, const idVertex &seed)
{
   const idSuperArc nbArc = getNode(baseNode)->getNumberOfDownSuperArcs();
   for (unsigned i = 0; i < nbArc; ++i) {
      const idSuperArc &downsaid = getNode(baseNode)->getDownSuperArcId(i);
      SuperArc *curArc   = getSuperArc(downsaid);
      if(curArc->getDownCT() != partition_) continue;
      Node *    downNode = getNode(curArc->getDownNodeId());

      if (isLower(downNode->getVertexId(), seed) ){
          //downNode->getVertexId() >= mesh_->getNumberOfOriginalVertices())
          hideArc(downsaid);
          hideNode(curArc->getDownNodeId());

         //cout << " below hide vertex " << downNode->getVertexId() << " in "
              //<< static_cast<unsigned>(partition_) << endl;
         //cout << " below hide arc  " << printArc(getNode(baseNode)->getDownSuperArcId(i)) << endl;
      }
   }

   removeHiddenDownArcs(baseNode);
}

const idSuperArc MergeTree::hideAndClearLeadingTo(const idNode &baseNode, const idVertex &v)
{
    if(isCorrespondingNode(v)){
        const auto nbDown = getNode(baseNode)->getNumberOfDownSuperArcs();
        for (unsigned aid = 0; aid < nbDown; aid++) {
           const idSuperArc a        = getNode(baseNode)->getDownSuperArcId(aid);
           if(getSuperArc(a)->getDownCT() != partition_ || !getSuperArc(a)->isVisible()){
              // keep external ?
              continue;
           }
           const idNode     downNode = getSuperArc(a)->getDownNodeId();

           if (getNode(downNode)->getVertexId() == v) {
              hideArc(a);
              return a;
           }
        }
    } else {
       if (isCorrespondingArc(v)) {
          const idSuperArc a = getCorrespondingSuperArcId(v);
          hideArc(a);
          return a;
       }
    }

    return nullSuperArc;
}

// }

// Operators : find, print & clone
// {

// Print
void MergeTree::printTree2()
{
#ifdef withOpenMP
#pragma omp critical
#endif
   {
      cout << "Partition : " << (int)partition_ << endl;
      cout << "Nodes" << endl;

      for (const auto &n : getNodes()) {
         if (!n.isHidden()) {
            cout << "Node " << n.getVertexId();

            if (n.isHidden())
               cout << " X ";

            cout << endl;
            cout << " arc up (" << n.getUpValence() << ") : ";

            for (int i = 0; i < n.getNumberOfUpSuperArcs(); ++i) {
               if (getSuperArc(n.getUpSuperArcId(i))->isVisible()) {
                  cout << "+";
               } else {
                  cout << "-";
               }
               cout << n.getUpSuperArcId(i) << " ";
            }

            cout << endl << " arc down (" << n.getDownValence() << ") : ";

            for (int i = 0; i < n.getNumberOfDownSuperArcs(); ++i) {
               if (getSuperArc(n.getDownSuperArcId(i))->isVisible()) {
                  cout << "+";
               } else {
                  cout << "-";
               }
               cout << n.getDownSuperArcId(i) << " ";
            }

            cout << endl;
         }
      }

      cout << "Arcs" << endl;

      for (idSuperArc said = 0; said < getNumberOfSuperArcs(); ++said) {
          const SuperArc &sa = vect_superArcs_[said];
         if (sa.isVisible()) {
            if (sa.getDownNodeId() == nullNodes) {
               cout << "||";

            } else {
               cout << getNode(sa.getDownNodeId())->getVertexId();
            }

            if (sa.isHidden() || sa.isMerged())
               cout << " <X> ";
            else
               cout << " <<" << said << ">> ";

            if (sa.getUpNodeId() == nullNodes) {
               cout << "||";

            } else {
               cout << getNode(sa.getUpNodeId())->getVertexId();
            }

            cout << endl;
         }
      }

      cout << "Leaves" << endl;

      for (const auto &l : vect_leaves_)
         cout << " " << vect_nodes_[l].getVertexId();

      cout << endl;

      cout << "Roots" << endl;

      for (const auto &r : vect_roots_)
         cout << " " << vect_nodes_[r].getVertexId();

      cout << endl;
   }
}

// Clone
MergeTree *MergeTree::clone() const
{
   MergeTree *newMT = new MergeTree(isJT, partition_);
   // shared / const data
   newMT->partition_      = partition_;
   newMT->debugLevel_     = debugLevel_;
   newMT->mesh_           = mesh_;
   newMT->scalars_        = scalars_;
   newMT->soSOffsets_     = soSOffsets_;
   newMT->sortedVertices_ = sortedVertices_;
   newMT->mirrorOffsets_  = mirrorOffsets_;
   // unique data
   newMT->vect_superArcs_         = vect_superArcs_;
   newMT->vect_nodes_             = vect_nodes_;
   newMT->vect_leaves_            = vect_leaves_;
   newMT->vect_roots_             = vect_roots_;
   newMT->vect_arcsCrossingBelow_ = vect_arcsCrossingBelow_;
   newMT->vect_arcsCrossingAbove_ = vect_arcsCrossingAbove_;
   // -> is a memory leak (vert2tree destructed by paraCT as a shared data, but new one here)
   newMT->vect_vert2tree_ = new vector<idCorresp>(vect_vert2tree_->size());
   copy(vect_vert2tree_->cbegin(), vect_vert2tree_->cend(), newMT->vect_vert2tree_->begin());

   return newMT;
}

void MergeTree::clone(const MergeTree *mt)
{
   // we already have common data
   vect_superArcs_         = mt->vect_superArcs_;
   vect_nodes_             = mt->vect_nodes_;
   vect_leaves_            = mt->vect_leaves_;
   vect_roots_             = mt->vect_roots_;
   vect_arcsCrossingBelow_ = mt->vect_arcsCrossingBelow_;
   vect_arcsCrossingAbove_ = mt->vect_arcsCrossingAbove_;
   vect_vert2tree_         = new vector<idCorresp>(vect_vert2tree_->size());
   copy(mt->vect_vert2tree_->cbegin(), mt->vect_vert2tree_->cend(), vect_vert2tree_->begin());
}

void MergeTree::shallowCopy(const MergeTree *mt)
{
   vect_superArcs_         = mt->vect_superArcs_;
   vect_nodes_             = mt->vect_nodes_;
   vect_leaves_            = mt->vect_leaves_;
   vect_roots_             = mt->vect_roots_;
   vect_vert2tree_         = mt->vect_vert2tree_;
   vect_arcsCrossingBelow_ = mt->vect_arcsCrossingBelow_;
   vect_arcsCrossingAbove_ = mt->vect_arcsCrossingAbove_;
}

// }
// Simplification
// {

// Hide the basenode and its upArc, merging the segmentation with masterArc
void MergeTree::hideAndMerge(const idSuperArc &mergingArcId, const idSuperArc &receptacleArcId,
                             const bool preserveNode)
{
   SuperArc *mergingArc    = getSuperArc(mergingArcId);
   SuperArc *receptacleArc = getSuperArc(receptacleArcId);

   // merge     HERE WE LOST THE SORTED SEGM
   if(mergingArc->getVertSize() != -1){
      receptacleArc->appendSegmentation(mergingArc->getSegmentation());
   }

   if (!preserveNode) {
      hideNode(mergingArc->getDownNodeId());
   }

   mergeArc(mergingArcId, receptacleArcId);
}

void MergeTree::markThisArc(vector<ExtendedUnionFind *> &ufArray, const idNode &curNodeId,
                            const idSuperArc &mergingArcId, const idNode &parentNodeId)
{
    // size of this subtree segmentation + segmentation of this arc
    const auto &curSegmenSize =
        getSuperArc(mergingArcId)->getVertSize() + ufArray[curNodeId]->find()->getOrigin() + 2;
    // +2 for the merging nodes

    // UF propagation
    if (ufArray[parentNodeId] == nullptr) {
       // Parent have never been seen : recopy UF
       ufArray[parentNodeId] = ufArray[curNodeId]->find();
       ufArray[parentNodeId]->find()->setOrigin(curSegmenSize);
       // cout << "will merge " << getNode(curNodeId)->getVertexId() << endl;
   } else {
      // The parent have already been visited : merge UF and segmentation
      const auto &oldSegmentationSize = ufArray[parentNodeId]->find()->getOrigin();
      ExtendedUnionFind::makeUnion(ufArray[curNodeId]->find(), ufArray[parentNodeId]->find())
          ->setOrigin(oldSegmentationSize + curSegmenSize);
      //cout << "Union on " << getNode(parentNodeId)->getVertexId();
      //cout << " from " << getNode(curNodeId)->getVertexId() << endl;
   }

   // The last parentNode is the root of the subtree
   //cout << "for " << getNode(curNodeId)->getVertexId() << " set root " << getNode(parentNodeId)->getVertexId() << endl;
   ufArray[parentNodeId]->find()->setData(-((ufDataType)parentNodeId) - 1);
}

const idSuperArc MergeTree::newUpArc(const idNode &curNodeId, vector<ExtendedUnionFind *> &ufArray)
{

    idSuperArc keepArc = nullSuperArc;
    const auto nbUp = getNode(curNodeId)->getNumberOfUpSuperArcs();
    for (int d = 0; d < nbUp; d++) {
        const idSuperArc &curArc = getNode(curNodeId)->getUpSuperArcId(d);
        if (!getSuperArc(curArc)->isVisible())
           continue;

        keepArc = curArc;

        const idNode &newUp = getSuperArc(curArc)->getUpNodeId();
        if (!ufArray[newUp] || ufArray[curNodeId]->find() != ufArray[newUp]->find()) {
           return curArc;
        }
    }

    //cout << "node " << getNode(curNodeId)->getVertexId() << " have no up arc to take" << endl;
    return keepArc;
}

const idSuperArc MergeTree::newDownArc(const idNode &               curNodeId,
                                       vector<ExtendedUnionFind *> &ufArray)
{

    idSuperArc keepArc = nullSuperArc;
    const auto nbDown = getNode(curNodeId)->getNumberOfDownSuperArcs();
    for (int d = 0; d < nbDown; d++) {
        const idSuperArc &curArc = getNode(curNodeId)->getDownSuperArcId(d);
        if (!getSuperArc(curArc)->isVisible())
           continue;

        keepArc = curArc;

        const idNode &newDown = getSuperArc(curArc)->getDownNodeId();
        if (!ufArray[newDown] || ufArray[curNodeId]->find() != ufArray[newDown]->find()) {
           return curArc;
        }
    }

    //cout << "node " << getNode(curNodeId)->getVertexId() << " have no down arc to take" << endl;
    return keepArc;
}

const tuple<idNode, idNode, idVertex> MergeTree::createReceptArc(
    const idNode &root, const idSuperArc &receptacleArcId, vector<ExtendedUnionFind *> &ufArray,
    const vector<pair<short, short>> &valenceOffsets)
{

   const bool DEBUG = false;

   ExtendedUnionFind *ufRoot   = ufArray[root]->find();
   idNode        downNode = root;
   idNode        upNode   = root;

   if (DEBUG) {
      cout << " create receptarc for root : " << printNode(root) << endl;
      cout << " custom valence : " << valenceOffsets[root].first;
      cout << " + " << valenceOffsets[root].second << endl;
   }

   // descend in the tree until valence is not 2
   idVertex segmentationSize = ufRoot->find()->getOrigin();
   //cout << "init size " << segmentationSize << endl;

    // We need a valence of 2 (we don't want to cross a futur saddle
    // But we want to avoid up && down = root

   while (getNode(downNode)->getUpValence() - valenceOffsets[downNode].second == 1 &&
          getNode(downNode)->getDownValence() - valenceOffsets[downNode].first == 1) {
      // take the down node not leading to the current subtree
      // if have an UF, merge with current subtree
      // (else init it?)
      const idSuperArc &downArc = newDownArc(downNode, ufArray);

      // deal with arc segmentation
      segmentationSize += getSuperArc(downArc)->getVertSize()+2;
      // + 2 for merging nodes

      downNode = getSuperArc(downArc)->getDownNodeId();
      const idNode tmpUp = getSuperArc(downArc)->getUpNodeId();

      if (DEBUG) {
         cout << "change down to " << getNode(downNode)->getVertexId() << endl;
         cout << " new segmentation : " << segmentationSize << endl;
      }

      // UF
      if (ufArray[downNode]) {
         segmentationSize += ufArray[downNode]->find()->getOrigin();
         //ExtendedUnionFind::makeUnion(ufArray[downNode], ufRoot);
         if (ufArray[downNode]->find()->getData() < 0) {
            ufArray[downNode]->find()->setData(receptacleArcId);
         }
       } else {
         //ufArray[downNode] = ufRoot->find();
       }
      mergeArc(downArc, receptacleArcId, false);
      hideNode(tmpUp);
   }

   if (DEBUG) {
      cout << " continue receptarc for root : " << printNode(root) << endl;
      cout << " custom valence : " << valenceOffsets[root].first;
      cout << " + " << valenceOffsets[root].second << endl;
   }

    // for a node to be regular, it must have a down valence = 1 but
    // but we process the down SO ADAPT HERE
   while (getNode(upNode)->getUpValence() - valenceOffsets[upNode].second == 1 &&
          getNode(upNode)->getDownValence() - valenceOffsets[upNode].first == 1) {

      const idSuperArc &upArc = newUpArc(upNode, ufArray);

      segmentationSize += getSuperArc(upArc)->getVertSize()+2;

      upNode = getSuperArc(upArc)->getUpNodeId();
      const idNode tmpDown = getSuperArc(upArc)->getDownNodeId();

      if (DEBUG) {
         cout << "change up to " << getNode(upNode)->getVertexId() << endl;
         cout << " new segmentation : " << segmentationSize << endl;
      }

      if (ufArray[upNode]) {
         segmentationSize += ufArray[upNode]->find()->getOrigin();
         //ExtendedUnionFind::makeUnion(ufArray[upNode], ufRoot);
         if (ufArray[upNode]->find()->getData() < 0) {
            ufArray[upNode]->find()->setData(receptacleArcId);
         }
      } else {
         //ufArray[upNode] = ufRoot->find();
      }
      mergeArc(upArc, receptacleArcId, false);
      hideNode(tmpDown);
   }

   // if upNode == downNode, take one none merging arc randomly
   // (this case is possbile if several degenerate node are following)
   if (upNode == downNode) {
      // several degen. nodes adjacent
      // Prefer down for JT / ST
      idSuperArc tmpDown = newDownArc(downNode, ufArray);
      idSuperArc tmpUp   = newUpArc(upNode, ufArray);

      if (tmpDown == nullSuperArc) {
         upNode = getSuperArc(tmpUp)->getUpNodeId();
         if (ufArray[upNode]) {
            segmentationSize += ufArray[upNode]->find()->getOrigin();
            //ExtendedUnionFind::makeUnion(ufArray[downNode], ufRoot);
         } else {
            ufArray[upNode] = ufRoot->find();
         }
         getSuperArc(tmpUp)->merge(receptacleArcId);
      } else {
         downNode = getSuperArc(tmpDown)->getDownNodeId();
         if (ufArray[downNode]) {
            segmentationSize += ufArray[downNode]->find()->getOrigin();
            //ExtendedUnionFind::makeUnion(ufArray[upNode], ufRoot);
         } else {
            ufArray[downNode] = ufRoot->find();
         }
         getSuperArc(tmpDown)->merge(receptacleArcId);
      }

      if(tmpDown != nullSuperArc)
         segmentationSize += getSuperArc(tmpDown)->getVertSize() + 2;
      if(tmpUp != nullSuperArc)
         segmentationSize += getSuperArc(tmpUp)->getVertSize() + 2;

      //cout << " special : new segmentation : " << segmentationSize << endl;
   }

   return make_tuple(downNode, upNode, segmentationSize);
}

// }
// ---------------------- Debug
// {

bool MergeTree::verifyTree(void)
{
   bool res = true;

#pragma omp critical
   {
      const unsigned &nbArcs  = getNumberOfSuperArcs();
      const unsigned &nbNodes = getNumberOfNodes();

      cout << "Verify Tree : " << endl;
      cout << "nbNode initial : " << nbNodes << endl;
      cout << "nbArcs initial : " << nbArcs << endl;

      unsigned int nbArcsVisibles  = 0;
      unsigned int nbNodesVisibles = 0;

      // for each visible arc, verify he is in the node

      for (unsigned aid = 0; aid < nbArcs; aid++) {
         const SuperArc &arc = vect_superArcs_[aid];
         if (arc.isVisible()) {
            ++nbArcsVisibles;
            const idNode &up  = arc.getUpNodeId();
            const idNode &down  = arc.getDownNodeId();
            if(up == nullSuperArc || down == nullSuperArc){
                res = false;
                cout << "[Verif]: arc id : " << aid << "have a null boundary :";
                cout << " down :" << down << " up:" << up << endl;
            } else {
               bool isIn = false;

               // Arc is present in its upNode
               const Node &nup = vect_nodes_[up];
               const auto &upNbDown = nup.getNumberOfDownSuperArcs();
               for (unsigned int d = 0; d < upNbDown; d++) {
                  if (nup.getDownSuperArcId(d) == aid) {
                     isIn = true;
                     break;
                  }
               }
               if (!isIn) {
                  res = false;
                  cout << "[Verif]: arc " << printArc(aid) << " is not known by its up node :";
                  cout << vect_nodes_[up].getVertexId() << endl;
               }

               isIn = false;

               // Arc is present in its upNode
               const Node &ndown = vect_nodes_[arc.getDownNodeId()];
               const auto &upNbUp = ndown.getNumberOfUpSuperArcs();
               for (unsigned int u = 0; u < upNbUp; u++) {
                  if (ndown.getUpSuperArcId(u) == aid) {
                     isIn = true;
                     break;
                  }
               }
               if (!isIn) {
                  res = false;
                  cout << "[Verif]: arc " << printArc(aid) << " is not known by its down node :";
                  cout << vect_nodes_[down].getVertexId() << endl;
               }
            }
         }
      }

      // for each node, verify she is in the arc

      for (unsigned nid = 0; nid < nbNodes; nid++) {
          const Node &node = vect_nodes_[nid];
          if(!node.isHidden()){
            ++nbNodesVisibles;

            // Verify up arcs
            const unsigned &nbup = node.getNumberOfUpSuperArcs();
            for (unsigned ua = 0; ua < nbup; ua++) {
               const SuperArc &arc         = vect_superArcs_[node.getUpSuperArcId(ua)];
               const idNode    arcDownNode = arc.getDownNodeId();
               if(arcDownNode != nid || !arc.isVisible()){
                   res = false;
                   const idNode upnode = arc.getUpNodeId();
                   const idNode downnode = arc.getDownNodeId();
                   if(upnode == nullSuperArc || downnode == nullSuperArc){
                      cout << "[Verif]: arc id : " << node.getUpSuperArcId(ua);
                      cout << "have a null boundary :";
                      cout << " down :" << downnode << " up:" << upnode << endl;
                   } else {
                      cout << "[Verif] Node " << node.getVertexId() << " id : " << nid;
                      cout << " Problem with up arc : " << printArc(node.getUpSuperArcId(ua))
                           << endl;
                   }
               }
            }

            // Verify down arcs
            const unsigned &nbdown = node.getNumberOfDownSuperArcs();
            for (unsigned da = 0; da < nbdown; da++) {
               const SuperArc &arc         = vect_superArcs_[node.getDownSuperArcId(da)];
               const idNode    arcUpNode = arc.getUpNodeId();
               if(arcUpNode != nid || ! arc.isVisible()){
                   res = false;
                   const idNode upnode = arc.getUpNodeId();
                   const idNode downnode = arc.getDownNodeId();
                   if(upnode == nullSuperArc || downnode == nullSuperArc){
                      cout << "[Verif]: arc id : " << node.getDownSuperArcId(da);
                      cout << "have a null boundary :";
                      cout << " down :" << downnode << " up:" << upnode << endl;
                   } else {
                      cout << "[Verif] Node " << node.getVertexId() << " id : " << nid;
                      cout << " Problem with down arc : " << printArc(node.getUpSuperArcId(da))
                           << endl;
                   }
               }
            }

          }
      }

      // verify segmentation information
      const auto& nbVert = mesh_->getNumberOfVertices();
      vector<bool> segmSeen(nbVert, false);

      for (unsigned aid = 0; aid < nbArcs; aid++) {
          SuperArc &arc = vect_superArcs_[aid];

          if (!arc.isVisible())
             continue;

          const unsigned segmSize = arc.getVertSize();
          const pair<idVertex, bool> *segmVect = arc.getVertList();

          if (segmSize && !segmVect) {
             res = false;
             cout << "[Verif] Inconsistant segmentation for arc : ";
             cout << printArc(aid);
             cout << " have size of " << segmSize;
             cout << " and a null list" << endl;
         }

         for (unsigned v = 0; v < segmSize; v++) {
             if(!segmVect[v].second) {
                segmSeen.at(segmVect[v].first) = true;
             }
         }
      }

      for(const Node & node : vect_nodes_) {
        if(node.isHidden())
            continue;

        segmSeen.at(node.getVertexId()) = true;
      }

      cout << "Segm missing : ";
      for (int v = 0; v < nbVert; v++) {
          if(!segmSeen[v]) {
             res = false;
             cout << v << ", ";
          }
      }
      cout << endl;

      cout << "Nb visible Node : " << nbNodesVisibles << endl;
      cout << "Nb visible Arcs : " << nbArcsVisibles << endl;
   }
   return res;
}

// }

// ------------------ Contour Tree

ContourTree::ContourTree() : MergeTree(false), jt_(new MergeTree(true)), st_(new MergeTree(false))
{
}

ContourTree::~ContourTree()
{
   if (jt_) {
      delete jt_;
      jt_ = nullptr;
   }
   if (st_) {
      delete st_;
      st_ = nullptr;
   }
}

// Process
// {

int ContourTree::combine(const idVertex &seed0, const idVertex &seed1)
{
   deque<pair<bool, idNode>> queue_growingNodes;
   pair<bool, idNode>        head;

   MergeTree *xt = nullptr, *yt = nullptr;
   idNode     correspondingNodeId, parentId, node1, node2;
   Node *     currentNode, *parentNode;

   const bool DEBUG = false;

   // If a tree add only one leaf, the other tree is the result
   idVertex nbAddedleavesST = 0,
            nbAddedleavesJT = 0;

   // Add leves to growing nodes
   // We insert non hidden nodes, only those of the interface or those
   // just beyond, linked to a crossing arc
   for (const auto &nId : st_->getLeaves()) {
      if (!st_->getNode(nId)->isHidden() && st_->getNode(nId)->getNumberOfSuperArcs()) {
         queue_growingNodes.emplace_back(false, nId);
         //cout << "add : id " << nId << " node : " << st_->getNode(nId)->getVertexId() << endl;
         ++nbAddedleavesST;
      }
   }

   if (nbAddedleavesST == 1) {
      queue_growingNodes.clear();
   }

   for (auto &nId : jt_->getLeaves()) {
      if (!jt_->getNode(nId)->isHidden() && jt_->getNode(nId)->getNumberOfSuperArcs()) {
         queue_growingNodes.emplace_back(true, nId);
         ++nbAddedleavesST;
      }
   }

   if (nbAddedleavesJT == 1) {
      queue_growingNodes.pop_back();
   }

   if(DEBUG) {
      cout << "queue_growingNodes : " << queue_growingNodes.size() << endl;
   }

   if (nbAddedleavesST == 1 && nbAddedleavesJT == 1) {
       // ultra simplistic case where both tree a filliform
       shallowCopy(jt_);
       return 0;
   }

   // Warning, have a reserve here, can't make it at the begnining, need build output
   vect_leaves_.reserve(jt_->getLeaves().size() + st_->getLeaves().size());

   if (queue_growingNodes.empty()) {
      cout << "[ContourTree::combine ] Nothing to combine" << endl;
   }

   // seed : to keep crossing edges;
   const idVertex &s0 = (seed0 == -1) ? nullVertex : sortedVertices_[seed0];
   const idVertex &s1 =
       (seed1 > mesh_->getNumberOfVertices()) ? nullVertex : sortedVertices_[seed1];

   while (!queue_growingNodes.empty()) {
      // i <- Get(Q)
      head = queue_growingNodes.front();

      if (head.first) {
         // node come frome jt
         xt = jt_;
         yt = st_;
         if (DEBUG) {
            cout << endl << "JT :";
         }
      } else {
         // node come from st
         xt = st_;
         yt = jt_;
         if (DEBUG) {
            cout << endl << "ST :";
         }
      }

      currentNode = xt->getNode(head.second);

      if (DEBUG) {
         cout << "node : " << currentNode->getVertexId() << endl;
      }

      correspondingNodeId = yt->getCorrespondingNode(currentNode->getVertexId());

      if (isCorrespondingNode(currentNode->getVertexId())) {
         // already a node in the tree
         node1 = getCorrespondingNode(currentNode->getVertexId());
      } else {
         // create a new node
         node1 = makeNode(currentNode);

         if (!currentNode->getNumberOfDownSuperArcs())
            vect_leaves_.emplace_back(node1);
         else if (!currentNode->getNumberOfUpSuperArcs())
            vect_leaves_.emplace_back(node1);
      }

      // "choose a non-root leaf that is not a split in ST" so we ignore such nodes
      if (currentNode->getNumberOfUpSuperArcs() == 0 ||
          yt->getNode(correspondingNodeId)->getNumberOfDownSuperArcs() > 1) {

         queue_growingNodes.pop_front();

         if (yt->getNode(correspondingNodeId)->getNumberOfDownSuperArcs() > 1) {
            if (DEBUG) {
               cout << "re-enqueue and";
            }

            queue_growingNodes.emplace_back(head.first, head.second);
         }

         if (DEBUG) {
            cout << " ignore" << endl;
         }

         continue;
      }

      // j <- GetAdj(XT, i)
      idSuperArc curUpArc = currentNode->getUpSuperArcId(0);
      if(xt->getSuperArc(curUpArc)->isMerged()) curUpArc = xt->getSuperArc(curUpArc)->getReplacantArcId();
      parentId = xt->getSuperArc(curUpArc)->getUpNodeId();
      if (parentId == nullNodes) {
         // security : if not closed arc, close it here
         // can append when 1 vertex is isolated
         parentId =
             xt->makeNode(xt->getSuperArc(currentNode->getUpSuperArcId(0))->getLastVisited());
         // single not si not crossing anything
         xt->closeSuperArc(currentNode->getUpSuperArcId(0), parentId, false, false);
      }

      parentNode = xt->getNode(parentId);

      // cout << " parent node :" << parentNode->getVertexId() << endl;

      // HERE parent is null ...
      if (isCorrespondingNode(parentNode->getVertexId())) {
         // already a node in the tree
         node2 = getCorrespondingNode(parentNode->getVertexId());
      } else {
         // create a new node
         node2 = makeNode(parentNode);
         if (!parentNode->getNumberOfUpSuperArcs())
            vect_leaves_.emplace_back(node2);
      }

      // AddArc(CT, ij)
      pair<idVertex, bool> *arcVertList = nullptr;
      idVertex arcVertSize = 0;
      {
         // Retrieve segmentation info
         if (segmentation_) {
            arcVertList = xt->getSuperArc(currentNode->getUpSuperArcId(0))->getVertList();
            arcVertSize = xt->getSuperArc(currentNode->getUpSuperArcId(0))->getVertSize();
         }

         bool overlapB = false, overlapA = false;

         // If the created arc cross tha above or below interface, keep this info
         if (s0 != nullVertex) {
            overlapB = (isEqHigher(currentNode->getVertexId(), s0) !=
                        isEqHigher(parentNode->getVertexId(), s0));
         }
         if (s1 != nullVertex) {
            overlapA = (isEqHigher(currentNode->getVertexId(), s1) !=
                        isEqHigher(parentNode->getVertexId(), s1));
         }

         idSuperArc createdArc;
         // create the arc
         if (isLower(currentNode->getVertexId(), parentNode->getVertexId())) {  // take care of the order
            createdArc = makeSuperArc(node1, node2, overlapB, overlapA, arcVertList, arcVertSize);
         } else {
            createdArc = makeSuperArc(node2, node1, overlapB, overlapA, arcVertList, arcVertSize);
         }

         if (overlapB) {
            vect_arcsCrossingBelow_.emplace_back(createdArc);
         }
         if (overlapA)
            vect_arcsCrossingAbove_.emplace_back(createdArc);

         int nbv = 0;
         // update segmentation info
         if (segmentation_) {
            // pour chaque vertex de l'arc fraichement ajoute : regard vect_vert2tree
            // si deja set : masque le vertex, sinon set to current arc
            for (idVertex vert = 0; vert < arcVertSize; vert++) {
               const idVertex &v = arcVertList[vert].first;
               if (isCorrespondingNull(v)) {
                  updateCorrespondingArc(v, createdArc);
                  ++nbv;
               } else {
                  getSuperArc(createdArc)->setMasqued(vert);
               }
            }
         }

         if (DEBUG) {
            cout << " arc added : (segm: " << nbv << ") ";
            cout << printArc(createdArc) << endl;
         }
      }

      // DelNode(XT, i)
      {
         if (DEBUG) {
         cout << " delete xt (" << isJT << ") node :" << xt->getNode(head.second)->getVertexId()
              << endl;
         }

         xt->delNode(head.second);
      }

      // DelNode(YT, i)
      {
         if (yt->getNode(correspondingNodeId)->getNumberOfDownSuperArcs() < 2) {
            if (DEBUG) {
               cout << " delete yt (" << head.first << ") node :";
               cout << yt->getNode(correspondingNodeId)->getVertexId();
               cout << " have : ";
               cout << static_cast<unsigned>(
                   yt->getNode(correspondingNodeId)->getNumberOfDownSuperArcs());
               cout << " down";
               cout << " and : " << static_cast<unsigned>(
                                        yt->getNode(correspondingNodeId)->getNumberOfUpSuperArcs())
                    << " up" << endl;
            }

            yt->delNode(correspondingNodeId, arcVertList, arcVertSize);
         }
      }

      if (parentNode->getNumberOfDownSuperArcs() == 0 && parentNode->getNumberOfUpSuperArcs()) {
         queue_growingNodes.emplace_back(head.first, parentId);

         if(DEBUG){
             cout << "will see : " << parentNode->getVertexId() << endl;
         }
      }

      queue_growingNodes.pop_front();
   }

   return 0;
}

//}

// -------------------- Print arc and node

ostream &wtfit::operator<<(ostream &o, SuperArc const &a)
{
   o << a.getDownNodeId() << " <>> " << a.getUpNodeId();
   return o;
}

ostream &wtfit::operator<<(ostream &o, Node const &n)
{
   o << n.getNumberOfDownSuperArcs() << " .-. " << n.getNumberOfUpSuperArcs();
   return o;
}

// ----------------------
// NOTE : Can we make this class parent of the wrapper ? (public => only build and access the tree)
