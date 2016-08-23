/*
 * file: ContourTree.h
 * description: ContourTree processing package.
 * author: Gueunet Charles
 * date: Juin 2015
 */

/// \ingroup baseCode
/// \class wtfit::ContourTree
/// \brief %ContourTree processing package.
///
/// %ContourTree compute sequentially the contour tree for a given input
/// scalar
/// \param dataType Data type of the input scalar field (char, float,
/// etc.).
/// \sa vtkContourTree

#ifndef _CONTOURTREE_H
#define _CONTOURTREE_H

// base code includes
#include <DataTypes.h>
#include <ExtendedUF.h>
#include <Geometry.h>
#include <Triangulation.h>
#include <Wrapper.h>

#include <queue>
#include <set>
#include <tuple>
#include <vector>

#include <memory>
#include <parallel/algorithm>

namespace wtfit
{
   // Classes : SuperArc
   //           Node
   //           MergeTree
   //           ContourTree

   class SuperArc
   {
     private:
      // Extrema
      idNode downNodeId_, upNodeId_;
      // SuperArc are in charge of interCT communication
      idPartition downCT_, upCT_, replacantCT_;
      // Before stitching, an arc can cross the interface below the partition.
      // It can also cross the interface above.
      // We keep these information
      bool overlapBelow_, overlapAbove_;
      // Keep th last vertex seen by this arc
      // After the build a a merge tree, a close step is
      // done, using this field to close each root arc
      idVertex lastVisited_;
      // Stat of this arc, if replaced...
      ComponentState state_;
      // ... use this field to know by wich other arc.
      // Caution, we do not want chained replacant !
      idSuperArc replacantId_;

      // Regular nodes in this arc
      // Vector for initialisation only (mergetree::build & simplify)
      // Use vertList and sizeVertList_ to acces Arc's vertices after combine
      vector<pair<idVertex, bool>> vertices_;
      // initialized with vertices_.data() and vertices_.size()
      // these two variable allow to split the segmentation when
      // inserting node in the arc without moving memory
      // |------------| Arc
      // [............] Size N
      //
      // |---*--------| Arc with new node
      // [..][........] Size N1 + N2 = N
      // So we split static array with a2 = a1 + N1
      //
      // It works because vertices are sorted by construction in these arrays

      // Array version, also Retain masqued regular nodes :
      // when a node is removed from an arc of a tree,
      // we mark as masqed vertices that are in the arc we added in CT
      // and also in the staying arc to avoid duplicate.
      // _/!\_ used at combine time => create this array using sizeVertList_
      pair<idVertex, bool> *vertList_;
      // The size of *vertList_
      idVertex sizeVertList_;
#ifndef withKamikaze
      // add a size verification for global simplify step
      idVertex allocSgm_=-1;
#endif

     public:
      // CONSTRUCT
      // -----------------
      // {
      SuperArc(const idNode &d, const idNode &u, const bool overB, const bool overA,
               const unsigned char &ctd = 0, const unsigned char &ctu = 0, const size_t &resv = 0ul,
               const ComponentState &state = ComponentState::VISIBLE)
          : downNodeId_(d),
            upNodeId_(u),
            downCT_(ctd),
            upCT_(ctu),
            overlapBelow_(overB),
            overlapAbove_(overA),
            lastVisited_(nullVertex),
            state_(state),
            replacantId_(nullSuperArc),
            vertList_(NULL),
            sizeVertList_(-1)
      {
         vertices_.reserve(resv);
      }

      // }
      // ------------------
      // ACCESSOR
      // --------------------
      // {
      //
      // node
      // .................................{

      inline const idNode &getUpNodeId(void) const
      {
         return upNodeId_;
      }

      inline const idNode &getDownNodeId(void) const
      {
         return downNodeId_;
      }

      inline void setUpNodeId(const idNode &upId)
      {
         upNodeId_ = upId;
      }

      inline void setDownNodeId(const idNode &downId)
      {
         downNodeId_ = downId;
      }

      // }
      // tree
      // .................................{
      //
      // idPartition u_char so lighter than a ref
      inline const idPartition getDownCT(void) const
      {
         return downCT_;
      }

      inline const idPartition getUpCT(void) const
      {
         return upCT_;
      }

      inline void setUpCT(decltype(upCT_) &ct)
      {
         upCT_ = ct;
      }
      // }
      // overlap
      // .................................{

      inline const bool getOverlapAbove(void) const
      {
         return overlapAbove_;
      }

      inline const bool getOverlapBelow(void) const
      {
         return overlapBelow_;
      }

      inline const bool isCrossing(void) const
      {
          return overlapBelow_ || overlapAbove_;
      }

      inline void setOverlapAbove(const bool local_overlapAbove)
      {
         overlapAbove_ = local_overlapAbove;
      }

      inline void setOverlapBelow(const bool local_overlapBelow)
      {
         overlapBelow_ = local_overlapBelow;
      }

      // }
      // last vertex seen
      // .................................{

      inline const idVertex &getLastVisited(void) const
      {
         return lastVisited_;
      }

      inline void setLastVisited(const idVertex &vertId, const bool segment)
      {
         lastVisited_ = vertId;
         if (segment) {
            vertices_.emplace_back(vertId, false);
         }
      }

      // }
      // state
      // .................................{

      inline bool isHidden(void) const
      {
         return state_ == ComponentState::HIDDEN;
      }

      inline bool isPruned(void) const
      {
         return state_ == ComponentState::MERGED;
      }

      inline bool isMerged(void) const
      {
         return state_ == ComponentState::MERGED;
      }

      inline bool isVisible(void) const
      {
         return state_ == ComponentState::VISIBLE;
      }

      inline void hide(void)
      {
         state_ = ComponentState::HIDDEN;
      }

      inline void merge(const idSuperArc &arc, const idPartition ct = 255)
      {
         replacantCT_ = (ct == 255) ? upCT_ : ct;
         replacantId_ = arc;
         state_       = ComponentState::MERGED;
      }

      // }
      // replacant arc/tree (merge)
      // .................................{

      inline const idSuperArc &getReplacantArcId(void) const
      {
         return replacantId_;
      }

      inline const idPartition getReplacantCT(void) const
      {
         return replacantCT_;
      }

      // }
      // regular nodes (segmentation)
      // .................................{

      inline int getNumberOfRegularNodes(void)
      {
         return getVertSize();
      }

      inline const idVertex &getRegularNodeId(const int &idx)
      {
         return getVertList()[idx].first;
      }

      inline const bool isMasqued(const idVertex &v) const
      {
         return vertList_[v].second;
      }

      // The vector

      inline const idVertex getSegmentationSize(void) const
      {
         return vertices_.size();
      }

      // not const for sort in simplify
      inline vector<pair<idVertex, bool>> &getSegmentation(void)
      {
         return vertices_;
      }

      // The array

      inline pair<idVertex, bool> *getVertList()
      {
         if (sizeVertList_ == -1) {
            vertList_     = vertices_.data();
            sizeVertList_ = vertices_.size();
         }
         return vertList_;
      }

      inline const idVertex &getVertSize()
      {
         if (sizeVertList_ == -1) {
            vertList_     = vertices_.data();
            sizeVertList_ = vertices_.size();
         }
         return sizeVertList_;
      }

      inline void setMasqued(const unsigned &v)
      {
         vertList_[v].second = true;
      }

      inline void setVertList(pair<idVertex, bool> *vl)
      {
         vertList_ = vl;
      }

      inline void setVertSize(const int &s)
      {
         sizeVertList_ = s;
      }

      // append regular nodes :

      // From array : alloc should be already done
      inline void appendVertLists(pair<idVertex, bool> **vertList, idVertex *vertSize,
                                  const unsigned &nb)
      {
         // size local
         size_t newSize = sizeVertList_;
         // size added
         for (unsigned i = 0; i < nb; ++i) {
            newSize += vertSize[i];
         }

         // alloc
         pair<idVertex, bool> *tmpVert = new pair<idVertex, bool>[newSize];
         size_t pos = 0;

         // values local
         for (idVertex i = 0; i < sizeVertList_; ++i) {
            tmpVert[pos++] = vertList_[i];
         }

         // values added
         for (unsigned i = 0; i < nb; ++i) {
            for (idVertex j = 0; j < vertSize[i]; ++j) {
               tmpVert[pos++] = vertList[i][j];
            }
         }

         // set values
         vertList_     = tmpVert;
         sizeVertList_ = newSize;
      }

      // From array : alloc should be already done (chck boundary no kamikaze only)
      inline int addSegmentationGlobal(const pair<idVertex, bool> *arr, const idVertex &size)
      {

#ifndef withKamikaze
          //cout << "size " << sizeVertList_ << " add " << size << " on " << allocSgm_ << endl;
         if (sizeVertList_ + size >= allocSgm_) {
            cout << "SEGMENTATION SIZE PROBLEM :" << endl;
            cout << "alloc : " << allocSgm_ << endl;
            cout << "size : " << sizeVertList_ << endl;
            cout << "add : " << size << endl;
            // gdb
            return 1;
         }
#endif

         for (idVertex v = 0; v < size; v++) {
             if(!arr[v].second) {
                vertList_[sizeVertList_++] = arr[v];
             }
         }

         return 0;
      }

      // from vector
      inline void appendSegmentation(const vector<pair<idVertex, bool>> &other)
      {
         vertices_.insert(vertices_.end(), other.begin(), other.end());
      }

      // from vector with a move
      inline void setVertices(vector<pair<idVertex, bool>>::iterator &a,
                              vector<pair<idVertex, bool>>::iterator  b)
      {
         vertices_.insert(vertices_.end(), make_move_iterator(a), make_move_iterator(b));
      }

      // add one vertex, alloc should be already done
      inline void addSegmentationGlobal(const idVertex &v){
          vertList_[sizeVertList_++] = make_pair(v, false);
      }

      // }

      // }
      // -----------
      // RESERVE
      // ----------
      // {
      // Trick
      // During simplification we use sizeVertList to keep the number of vertices that will merge in
      // materArc
      // Doing so we ensure we have only one reserve on the vector
      // At the end, we set sizeVertList_ back to minus one for normal utilisation.

      inline void addFuturReserve(const idVertex &nb)
      {
         sizeVertList_ += nb;
         // We have an offset of -1 due to the initial value of sizeVertList_
      }

      inline void makeAllocGlobal(const idVertex& size)
      {
         vertList_     = new pair<idVertex, bool>[size];
         sizeVertList_ = 0;
#ifndef withKamikaze
         allocSgm_     = size;
#endif
      }

      // Real reserve
      inline void makeReserve(void)
      {
         if (sizeVertList_ != -1) {
            // We have an offset of -1 due to the initial value of sizeVertList_
            vertices_.reserve(sizeVertList_ + 1);
            sizeVertList_ = -1;
         }
      }

      //}
    };

   class Node
   {
      friend class MergeTree;

     private:
      // mesh vertex where this node is
      idVertex vertexId_;
      // For leaves, linkedNode is the saddle ending the persistance pair
      // For saddle, linked is the leaf starting the persistance pair in which they are
      idVertex linkedNode_;
      // link with superArc above and below
      vector<idSuperArc> vect_downSuperArcList_, vect_upSuperArcList_;
      // Won't be displayed if hidden
      bool hidden_;
      // valence down / up
      tuple<short,short> valence_;

     public:
      // -----------------
      // CONSTRUCTOR
      // -----------------
      // {

      Node(const idVertex &id, const idVertex &linked)
          : vertexId_(id), linkedNode_(linked), hidden_(false), valence_(0,0)
      {
      }

      // }
      // -----------------
      // ACCESSOR
      // ------------------
      // {

      // Vertex id
      // ........................{

      inline idVertex getVertexId() const
      {
         return vertexId_;
      }

      inline void setVertexId(const idVertex &vertexId)
      {
         vertexId_ = vertexId;
      }

      // }
      // Linked node
      // ........................{

      inline const idVertex &getOrigin(void) const
      {
         return linkedNode_;
      }

      inline const idVertex &getTerminaison(void) const
      {
         return linkedNode_;
      }

      inline void setOrigin(const idVertex &linked)
      {
         linkedNode_ = linked;
      }

      inline void setTerminaison(const idVertex &linked)
      {
         linkedNode_ = linked;
      }

      // }
      // vector arcs
      // ............................{

      inline unsigned char getNumberOfDownSuperArcs() const
      {
         return (unsigned char)vect_downSuperArcList_.size();
      }


      inline unsigned char getNumberOfUpSuperArcs() const
      {
         return (unsigned char)vect_upSuperArcList_.size();
      }


      inline unsigned char getNumberOfSuperArcs() const
      {
         return (unsigned char)(vect_upSuperArcList_.size() + vect_downSuperArcList_.size());
      }

      inline idSuperArc getDownSuperArcId(const unsigned char &neighborId) const
      {
#ifndef withKamikaze
         if ((neighborId < 0) || ((size_t)neighborId >= vect_downSuperArcList_.size())) {
            cout << "[Merge Tree:Node] get down on bad neighbor !";
            cout << endl;
            return 0;
         }
#endif
         return vect_downSuperArcList_[neighborId];
      };

      inline idSuperArc getUpSuperArcId(const unsigned char &neighborId) const
      {
#ifndef withKamikaze
         if ((unsigned)neighborId >= vect_upSuperArcList_.size()) {
            cerr << "[MergeTree:Node] No SuperArc to access " << static_cast<unsigned>(neighborId);
            cerr << endl;
         }
#endif
         if (vect_upSuperArcList_.size() == 0) {
            return nullSuperArc;
         }
         return vect_upSuperArcList_[neighborId];
      }

      inline void addDownSuperArcId(const idSuperArc &downSuperArcId)
      {
         vect_downSuperArcList_.emplace_back(downSuperArcId);
      }

      inline void addUpSuperArcId(const idSuperArc &upSuperArcId)
      {
         vect_upSuperArcList_.emplace_back(upSuperArcId);
      }

      inline unsigned clearDownSuperArcs(void)
      {
         unsigned s = vect_downSuperArcList_.size();
         vect_downSuperArcList_.clear();
         return s;
      }

      inline unsigned clearUpSuperArcs(void)
      {
         unsigned s = vect_upSuperArcList_.size();
         vect_upSuperArcList_.clear();
         return s;
      }

      // remove the i^th arc
      inline void removeDownSuperArcId(unsigned i)
      {
          vect_downSuperArcList_[i] = vect_downSuperArcList_.back();
          vect_downSuperArcList_.pop_back();

          decDownValence();
      }

      // Find and remove the arc
      inline void removeDownSuperArc(const idSuperArc &idSa)
      {
         for (unsigned char i = 0; i < vect_downSuperArcList_.size(); ++i) {
            if (vect_downSuperArcList_[i] == idSa) {
               vect_downSuperArcList_[i] = vect_downSuperArcList_.back();
               vect_downSuperArcList_.pop_back();

               decDownValence();
               return;
            }
         }
      }

      // Find and remove the arc (better perf for young added arc)
      inline void removeDownSuperArcFromLast(const idSuperArc &idSa)
      {
         for (unsigned char i = vect_downSuperArcList_.size() - 1; i >= 0; --i) {
            if (vect_downSuperArcList_[i] == idSa) {
               vect_downSuperArcList_[i] = vect_downSuperArcList_.back();
               vect_downSuperArcList_.pop_back();

               decDownValence();
               return;
            }
         }
      }


      // Find and remove the arc
      inline void removeUpSuperArc(const idSuperArc &idSa)
      {
         for (unsigned char i = 0; i < vect_upSuperArcList_.size(); ++i) {
            if (vect_upSuperArcList_[i] == idSa) {
               vect_upSuperArcList_[i] = vect_upSuperArcList_.back();
               vect_upSuperArcList_.pop_back();

               decUpValence();
               return;
            }
         }
      }

      // Find and remove the arc (better perf for young added arc)
      inline void removeUpSuperArcFromLast(const idSuperArc &idSa)
      {
         for (unsigned char i = vect_upSuperArcList_.size(); i >= 0; --i) {
            if (vect_upSuperArcList_[i] == idSa) {
               vect_upSuperArcList_[i] = vect_upSuperArcList_.back();
               vect_upSuperArcList_.pop_back();

               decUpValence();
               return;
            }
         }
      }

      // }
      // hidden node
      // ...........................................{

      inline bool isHidden() const
      {
         return hidden_;
      }

      inline bool isVisible() const
      {
        return !hidden_;
      }

      inline void hide()
      {
         hidden_ = true;
      }

      inline void setHidden(const bool local_hidden)
      {
         hidden_ = local_hidden;
      }


      // }
      // Valence
      // .......................................... {

      inline const short getUpValence(void) const
      {
        return get<1>(valence_);
      }

      inline const short getDownValence(void) const
      {
        return get<0>(valence_);
      }

      inline const short getValence(void) const
      {
        return get<0>(valence_) + get<1>(valence_);
      }

      inline void setUpValence(const short v)
      {
        get<1>(valence_) = v;
      }

      inline void setDownValence(const short v)
      {
        get<0>(valence_) = v;
      }

      inline void incUpValence(void)
      {
        ++get<1>(valence_);
      }

      inline void incDownValence(void)
      {
        ++get<0>(valence_);
      }

      inline void decUpValence(void)
      {
        --get<1>(valence_);
      }

      inline void decDownValence(void)
      {
        --get<0>(valence_);
      }

      // }

      // }
   };

   class MergeTree : virtual public Debug
   {
      friend class ParallelContourTree;
      friend class ContourTree;

     protected:
      // MESH DATA ---------------------------------------- if change : update clone()

      Triangulation *mesh_;
      unsigned char  partition_;

      void *    scalars_;
      idVertex *soSOffsets_;

      vector<idVertex> sortedVertices_, mirrorOffsets_;

      // TREE DATA -----------------------------------------

      const bool       isJT;

      vector<SuperArc> vect_superArcs_;
      vector<Node>     vect_nodes_;
      vector<idNode>   vect_leaves_,  vect_roots_;

      // on arc can be in both
      vector<idSuperArc> vect_arcsCrossingBelow_, vect_arcsCrossingAbove_;

      vector<idCorresp> *vect_vert2tree_;

      // parameters
      bool segmentation_;
      bool computeContourTree_;
      int  debugLevel_;
      SimplifMethod simplifyMethod_;

     public:

      // CONSTRUCT
      // -----------
      // {

      MergeTree(decltype(partition_) part = 0);

      MergeTree(const bool t, const unsigned char part = 0);

      virtual ~MergeTree();

      //}
      // -------------
      // ACCESSOR
      // ------------
      //{

      // mesh
      // .....................{

      inline void setTriangulation( Triangulation * const m)
      {
          mesh_ = m;
      }

      //}
      // partition
      // .....................{

      inline const unsigned char getPartition(void) const
      {
         return partition_;
      }

      inline void setPartition(decltype(partition_) part)
      {
         partition_ = part;
      }

      // }
      // scalar
      // .....................{

      template <typename scalarType>
      inline const scalarType &getValue(const idVertex &idNode) const
      {
         return (((scalarType *)scalars_))[idNode];
      }

      template <typename scalarType>
      inline void setVertexScalars(scalarType *scalars)
      {
         scalars_ = (void *)scalars;
      }


      // }
      // offset
      // .....................{

      inline void setVertexSoSoffsets(decltype(soSOffsets_) const offsets)
      {
         destroyVectorSoS_ = false;
         soSOffsets_       = offsets;
      }

      // }
      // sorted / mirror
      // .....................{

      inline void setSorted(decltype(sortedVertices_) const sorted)
      {
         destroyVectorSortedVertices_ = false;
         sortedVertices_              = sorted;
      }

      inline void setMirror(decltype(mirrorOffsets_) const mir)
      {
         mirrorOffsets_ = mir;
      }

      // }
      // arcs
      // .....................{

      inline const idSuperArc getNumberOfSuperArcs(void) const
      {
         return (int)vect_superArcs_.size();
      }

      inline const idSuperArc getNumberOfVisibleArcs(void) const
      {
         idSuperArc visibleArc = 0;
         for (const SuperArc &arc : vect_superArcs_) {
            if (arc.isVisible())
               ++visibleArc;
         }
         return visibleArc;
      }

      inline const vector<SuperArc> &getSuperArc(void) const
      {
         return vect_superArcs_;
      }

      inline SuperArc *getSuperArc(const idSuperArc &i)
      {
#ifndef withKamikaze
         if ((size_t)i >= vect_superArcs_.size()) {
            cout << "[Merge Tree] get superArc on bad id :" << i;
            cout << " / " << vect_superArcs_.size() << endl;
            return nullptr;
         }
#endif
         return &(vect_superArcs_[i]);
      }

      inline idVertex getNumberOfVisibleRegularNode(const idSuperArc &sa)
      {
         idVertex  res   = 0;
         SuperArc *a     = getSuperArc(sa);
         const auto nbReg = a->getNumberOfRegularNodes();
         for (idVertex v = 0; v < nbReg; v++) {
            if (!a->isMasqued(v))
               ++res;
         }

         return res;
      }

      inline void addCrossingAbove(const idSuperArc &sa)
      {
         vect_arcsCrossingAbove_.emplace_back(sa);
      }

      // }
      // nodes
      // .....................{

      inline const int getNumberOfNodes(void) const
      {
         return (int)vect_nodes_.size();
      }

      inline const vector<Node>& getNodes(void) const
      {
         return vect_nodes_;
      }

      inline Node *getNode(const idNode &nodeId)
      {
         // NOTE needed for jt / st detection
         if (/*nodeId < 0 || */ nodeId >= (idNode)vect_nodes_.size()) {
            return nullptr;
         }
         return &(vect_nodes_[nodeId]);
      }

      // }
      // leaves / root
      // .....................{

      inline const idVertex getNumberOfLeaves(void) const
      {
         return vect_leaves_.size();
      }

      inline const vector<idNode>& getLeaves(void) const
      {
         return vect_leaves_;
      }

      inline const idNode &getLeave(const int &id) const
      {
#ifndef withKamikaze
         if ((id < 0) || (size_t)id > (vect_leaves_.size())) {
            stringstream msg;
            msg << "[MergTree] getLeaves out of bounds : " << id << endl;
            err(msg.str(), fatalMsg);
            return vect_leaves_[0];
         }
#endif
         return vect_leaves_[id];
      }

      // }
      // vert2tree
      // .....................{

      inline void setVert2Tree(decltype(vect_vert2tree_) const vect2tree)
      {
         destroyVect2Tree_ = false;
         vect_vert2tree_   = vect2tree;
      }

      // + Section spectial VERT 2 TREE after ACCESSOR
      // }
      // segmentation
      // .....................{

      inline void setSegmentation(const bool local_segmentation)
      {
         segmentation_ = local_segmentation;
      }

      // }
      // compute ct
      // .....................{

      inline void setComputeContourTree(const bool local_computeContourTree)
      {
         computeContourTree_ = local_computeContourTree;
      }

      // }
      // debuglevel
      // .....................{

      inline void setDebugLevel(decltype((debugLevel_)) d)
      {
         Debug::setDebugLevel(d);
         debugLevel_ = d;
      }

      // }
      // Simplification method
      // .....................{

      inline void setSimplificationMethod(int m) {
         switch (m) {
            case 0:
               simplifyMethod_ = SimplifMethod::Persist;
               break;
            case 1:
               simplifyMethod_ = SimplifMethod::Span;
               break;
            case 2:
               simplifyMethod_ = SimplifMethod::NbVert;
               break;
            case 3:
               simplifyMethod_ = SimplifMethod::NbArc;
               break;
            default:
               simplifyMethod_ = SimplifMethod::Persist;
               break;
         }
      }

      // }

      // }
      // --------------------
      // VERT 2 TREE Special functions
      // --------------------
      //{

      // test vertex correpondance
      // ...........................{

      inline const bool isCorrespondingArc(const int &val) const
      {
         return !isCorrespondingNull(val) && (*vect_vert2tree_)[val] >= 0;
      }

      inline const bool isCorrespondingNode(const int &val) const
      {
         return (*vect_vert2tree_)[val] < 0;
      }

      inline const bool isCorrespondingNull(const int &val) const
      {
         return (*vect_vert2tree_)[val] == nullCorresp;
      }

      //}
      // Get vertex info
      // ...........................{

      inline const idNode getCorrespondingDownNode(const idVertex &vertexId)
      {
         idCorresp val = (*vect_vert2tree_)[vertexId];
         if (val < 0) {
            return idNodeToCorr(val);
         }
         return vect_superArcs_[val].getDownNodeId();
      }

      inline const idNode getCorrespondingUpNode(const idVertex &vertexId)
      {
         idCorresp val = (*vect_vert2tree_)[vertexId];
         if (val < 0) {
            return idNodeToCorr(val);
         }
         return vect_superArcs_[val].getUpNodeId();
      }

      inline const idNode getCorrespondingNode(const int &val) const
      {
#ifndef withKamikaze
         if (!isCorrespondingNode(val)) {
            stringstream debug;
            debug << "[MergeTree] : getCorrespondingNode, ";
            debug << "Vertex :" << val << " is not a node";
            debug << (*vect_vert2tree_)[val] << endl;
            err(debug.str(), fatalMsg);
            //gdb
         }
#endif
         return -(*vect_vert2tree_)[val] - 1;
      }

      inline const idSuperArc getCorrespondingSuperArcId(const int &val) const
      {
#ifndef withKamikaze
         if (!isCorrespondingArc(val)) {
            stringstream debug;
            debug << "[MergeTree] : getCorrespondingSuperArcId, ";
            debug << "Vertex :" << val << " is not on an arc" << endl;
            err(debug.str(), fatalMsg);
         }
#endif
         return (*vect_vert2tree_)[val];
      }

      // }
      // Get vertex correponding object
      // ................................{

      inline SuperArc *vertex2SuperArc(const idVertex &vert)
      {
         return &(vect_superArcs_[getCorrespondingSuperArcId(vert)]);
      }

      inline Node *vertex2Node(const idVertex &vert)
      {
         return &(vect_nodes_[getCorrespondingNode(vert)]);
      }

      // }
      // Update vertex info
      // ................................{

      inline void updateCorrespondingArc(const idVertex &arc, const idSuperArc &val)
      {
         (*vect_vert2tree_)[arc] = val;
      }

      inline void updateCorrespondingNode(const idVertex &vert, const int &val)
      {
         (*vect_vert2tree_)[vert] = idNodeToCorr(val);
      }

      inline const idCorresp idNodeToCorr(const idNode &id) const
      {
         // transform idNode to special value for the array : -idNode -1
         return -(idCorresp)(id + 1);
      }


      // }

      // }
      // --------------------
      // Init
      // --------------------
      // {
      void flush(void);

      /// \brief init Simulation of Simplicity datastructure if not set
      void initSoS(void);

      /// \brief ensure the size of vert2tree
      void initVert2Tree(void);

      /// \brief if sortedVertices_ is null, define and fill it
      template <typename scalarType>
      void sortInput(void);

      //}
      // -------------------
      // Process
      // -------------------
      //{

      // build
      // ..........................{

      // Merge tree processing of a vertex during build
      void processVertex(const idVertex &vertex, vector<ExtendedUnionFind *> &vect_baseUF,
                         const bool overlapB, const bool overlapA, DebugTimer &begin);

      /// \brief Compute the merge tree using Carr's algorithm
      int build(vector<ExtendedUnionFind *> &vect_baseUF, const vector<idVertex> &overlapBefore,
                const vector<idVertex> &overlapAfter, idVertex start, idVertex end,
                const idVertex &posSeed0, const idVertex &posSeed1);

      // }
      // Simplify
      // ...........................{

      // BFS simplification for local CT
      template <typename scalarType>
      idEdge localSimplify(const idVertex &podSeed0, const idVertex &podSeed1,
                           const double threshold);

      // BFS simpliciation for global CT
      template <typename scalarType>
      idEdge globalSimplify(const idVertex posSeed0, const idVertex posSeed1, const double threshold);

      // Having sorted pairs, simplify the current tree
      // in accordance with threashol, between the two seeds.
      template <typename scalarType>
      idEdge simplifyTree(const idVertex &posSeed0, const idVertex &posSeed1,
                          const double threshold,
                          const vector<tuple<idVertex, idVertex, scalarType, bool>> &sortedPairs);

      // add this arc in the subtree which is in the parentNode
      void markThisArc(vector<ExtendedUnionFind *> &ufArray, const idNode &curNodeId,
                       const idSuperArc &mergingArcId, const idNode &parentNodeId);
      // }
      // PersistencePairs
      // ...........................{

      template <typename scalarType>
      int computePersistencePairs(vector<pair<pair<int, int>, double>> &pairs);

      template <typename scalarType>
      int computePersistencePairs(vector<tuple<idVertex, idVertex, scalarType, bool>> &pairs);

      // Construct abstract JT / ST on a CT and fill pairs in accordance.
      // used for global simplification
      template <typename scalarType>
      void computePersistencePairsMT(const vector<idNode> &sortedNodes,
                                   vector<tuple<idVertex, idVertex, scalarType,bool>> &pairsJT,
                                   vector<tuple<idVertex, idVertex, scalarType,bool>> &pairsST);

      template <typename scalarType>
      int computePersistencePlot(vector<pair<double, int>> &plot,
                                 vector<pair<pair<int, int>, double>> *persistencePairs);

      template <typename scalarType>
      int computePersistenceDiagram(vector<pair<double, double>> &diagram,
                                    vector<pair<pair<int, int>, double>> *pairs);


      // }

      // }
      // --------------------------------
      // Arcs and node manipulations
      // --------------------------------
      // {

      // SuperArcs
      // .......................{

      idSuperArc openSuperArc(const idNode &downNodeId, const bool overlapB, const bool overlapA);

      idSuperArc makeSuperArc(const idNode &downNodeId, const idNode &upNodeId, const bool overlapB,
                              const bool overlapA, pair<idVertex, bool> *vertexList = nullptr,
                              int vertexSize = -1);

      void closeSuperArc(const idSuperArc &superArcId, const idNode &upNodeId, const bool overlapB,
                         const bool overlapA);

      void hideArc(const idSuperArc &sa);

      void mergeArc(const idSuperArc &sa, const idSuperArc &recept, const bool changeConnectivity = true );

      const idVertex cutArcAboveSeed(const idSuperArc &arc, const pair<idVertex, bool> &seed);

      const idVertex cutArcBelowSeed(const idSuperArc &arc, const pair<idVertex, bool> &seed,const vector<idCorresp>* vert2treeOther);

      // is there an external arc linkind node with treeNode in tree
      const bool alreadyExtLinked(const idNode &node, const idPartition &tree,
                                  const idNode &treeNode);

      // TODO Remove that

      void removeHiddenDownArcs(const idNode &n);

      unsigned getNumberOfVisibleArcs(const idNode &n);

      unsigned getNumberOfUnmergedDownArcs(const idNode &n);

      // }
      // Nodes
      // ...........................{

      idNode makeNode(const idVertex &vertexId, const idVertex &linked = nullVertex);

      idNode makeNode(const Node *const n, const idVertex &linked = nullVertex);

      const idSuperArc insertNode(Node *node, const bool segment, const bool isCT = false);

      const idSuperArc reverseInsertNode(Node *node, const bool segment, const bool isCT = false);

      inline Node *getDownNode(const SuperArc *a);

      inline Node *getUpNode(const SuperArc *a);

      idNode getParent(const idNode &n);

      void delNode(const idNode &node, const pair<idVertex, bool> *mv = nullptr,
                   const int &nbm = 0);

      void hideNode(const idNode &node);

      // For persistance pair on CT
      // these function allow to make a JT / ST od the CT
      const vector<idNode> getNodeNeighbors(const idNode &node);

      const vector<idNode> getNodeUpNeighbors(const idNode &n);

      const vector<idNode> getNodeDownNeighbors(const idNode &n);

      // Remove part not in partition

      void hideAndClearArcsAbove(const idNode &baseNode);

      void hideAndClearArcsBelow(const idNode &baseNode, const idVertex &seed);

      const idSuperArc hideAndClearLeadingTo(const idNode &baseNode, const idVertex &v);

      // }
      // Update informations
      // ...........................{

      void updateSegmentation(const bool ct = false);

      void parallelUpdateSegmentation(const bool ct = false);


      // will disapear
      void initNodeValence();

      // will disapear
      void parallelInitNodeValence(const int nbThreadValence);

      // }

      // }
      // ---------------------------
      // Operators : print & clone
      // ---------------------------
      // {

      // Print
      void printTree2(void);

      string printArc(const idSuperArc &a)
      {
         const SuperArc *sa = getSuperArc(a);
         stringstream    res;
         res << getNode(sa->getDownNodeId())->getVertexId() << " - ";
         res << getNode(sa->getUpNodeId())->getVertexId();
         res << " (" << sa->isVisible() << ")";
         return res.str();
      }

      string printNode(const idNode &n)
      {
         const Node *node = getNode(n);
         stringstream res;
         res << node->getVertexId() << " : ";
         res << "( " << node->isHidden() << " ) ";
         res << " valence : " << node->getDownValence() << " + ";
         res << node->getUpValence();
         return res.str();
      }

      // Clone
      MergeTree *clone() const;

      void clone(const MergeTree *mt);

      // Share even the vert2tree
      void shallowCopy(const MergeTree *mt);

      //}

     protected:
      // ------------------
      // Comparisons
      // -----------------
      // {

      // Strict

      inline const bool isLower(const idVertex &a, const idVertex &b) const
      {
         return mirrorOffsets_[a] < mirrorOffsets_[b];
      }

      inline const bool isHigher(const idVertex &a, const idVertex &b) const
      {
         return mirrorOffsets_[a] > mirrorOffsets_[b];
      }

      // Large

      inline const bool isEqLower(const idVertex &a, const idVertex &b) const
      {
         return mirrorOffsets_[a] <= mirrorOffsets_[b];
      }

      inline const bool isEqHigher(const idVertex &a, const idVertex &b) const
      {
         return mirrorOffsets_[a] >= mirrorOffsets_[b];
      }

      //}

     private:

      // ------------------
      // Comparisons
      // -----------------
      // {
      // Compare using the scalar array : only for sort step

      template <typename scalarType>
      inline const bool isLower(const idVertex &a, const idVertex &b) const
      {
         return ((scalarType *)scalars_)[a] < ((scalarType *)scalars_)[b] ||
                (((scalarType *)scalars_)[a] == ((scalarType *)scalars_)[b] &&
                 soSOffsets_[a] < soSOffsets_[b]);
      }

      template <typename scalarType>
      inline const bool isHigher(const idVertex &a, const idVertex &b) const
      {
         return ((scalarType *)scalars_)[a] > ((scalarType *)scalars_)[b] ||
                (((scalarType *)scalars_)[a] == ((scalarType *)scalars_)[b] &&
                 soSOffsets_[a] > soSOffsets_[b]);
      }

      template <typename scalarType>
      inline const bool isEqLower(const idVertex &a, const idVertex &b) const
      {
         return ((scalarType *)scalars_)[a] < ((scalarType *)scalars_)[b] ||
                (((scalarType *)scalars_)[a] == ((scalarType *)scalars_)[b] &&
                 soSOffsets_[a] <= soSOffsets_[b]);
      }

      template <typename scalarType>
      inline const bool isEqHigher(const idVertex &a, const idVertex &b) const
      {
         return ((scalarType *)scalars_)[a] > ((scalarType *)scalars_)[b] ||
                (((scalarType *)scalars_)[a] == ((scalarType *)scalars_)[b] &&
                 soSOffsets_[a] >= soSOffsets_[b]);
      }

      // }
      // ----------------
      // Simplification
      // ----------------
      // {

      // preserve = do no hide it.
      void hideAndMerge(const idSuperArc &mergingArcId, const idSuperArc &receptacleArcId,
                        const bool preserveDownNode = false);

      // Use BFS from root to find down and up of the receptarc (maintaining segmentation information)
      const tuple<idNode, idNode, idVertex> createReceptArc(
          const idNode &root, const idSuperArc &receptArcId, vector<ExtendedUnionFind *> &arrayUF,
          const vector<pair<short, short>> &valenceOffsets);

      // during this BFS nodes should have only one arc up/down : find it :
      const idSuperArc newUpArc(const idNode &curNodeId, vector<ExtendedUnionFind *> &ufArray);

      const idSuperArc newDownArc(const idNode &curNodeId, vector<ExtendedUnionFind *> &ufArray);

      // }
      // --------------
      // Tool
      // --------------
      // {
      // create a pair with relative order : child vertex first

      inline pair<idVertex, idVertex> reorderEdgeRel(const pair<idVertex, idVertex> &vert)
      {
         if (isJT) {
            if (isLower(vert.first, vert.second))
               return vert;

            return make_pair(vert.second, vert.first);
         }  // else

         if (isHigher(vert.first, vert.second)) {
            return vert;
         }

         return make_pair(vert.second, vert.first);
      }

      bool verifyTree(void);

      // for simplification
      template <typename scalarType>
      void addPair(vector<tuple<idVertex, idVertex, scalarType, bool>> &pairs, const idVertex &orig,
                   const idVertex &term, const bool goUp);

      template <typename scalarType>
      void addPair(vector<pair<pair<int, int>, double>> &pairs, const idVertex &orig,
                   const idVertex &term);

      // }

      // Know what to destroy
      bool destroyVectorSortedVertices_ = true;
      bool destroyVectorSoS_            = true;
      bool destroyVect2Tree_            = true;
   };

   class ContourTree : public MergeTree
   {
      friend class ParallelContourTree;

     protected:
      MergeTree *jt_, *st_;


     public:

      // -----------------
      // Constructors
      // -----------------
      // {

      ContourTree();
      virtual ~ContourTree();

      // }
      // -----------------
      // ACCESSOR
      // -----------------
      // {

      inline void setTriangulation(decltype(mesh_) m)
      {
         MergeTree::setTriangulation(m);
         jt_->setTriangulation(m);
         st_->setTriangulation(m);
      }

      inline auto getJoinTree(void) const -> decltype(jt_) const
      {
         return jt_;
      }

      inline auto getSplitTree(void) const -> decltype(st_) const
      {
         return st_;
      }

      inline void setPartition(decltype((partition_)) part)
      {
         partition_ = part;
         jt_->setPartition(part);
         st_->setPartition(part);
      }

      inline void setDebugLevel(int &d)
      {
         MergeTree::setDebugLevel(d);
         jt_->setDebugLevel(d);
         st_->setDebugLevel(d);
      }

      // }
      // -----------------
      // PROCESS
      // -----------------
      // {

      /// \brief Combine tree with Natarajan's algorithm
      int combine(const idVertex &seed0, const idVertex &seed1);

      // Persistenc pairs

      template <typename scalarType>
      int computePersistencePairs(vector<pair<pair<int, int>, double>> *pairs,
                                  vector<pair<pair<int, int>, double>> *mergePairs,
                                  vector<pair<pair<int, int>, double>> *splitPairs);

      template <typename scalarType>
      int computePersistencePlot(vector<pair<double, int>> &plot,
                                 vector<pair<pair<int, int>, double>> *mergePairs,
                                 vector<pair<pair<int, int>, double>> *splitPairs,
                                 vector<pair<pair<int, int>, double>> *pairs);

      template <typename scalarType>
      int computePersistenceDiagram(vector<pair<double, double>> &diagram,
                                    vector<pair<pair<int, int>, double>> *mergePairs,
                                    vector<pair<pair<int, int>, double>> *splitPairs,
                                    vector<pair<pair<int, int>, double>> *pairs);
      // }

     private:

      // -----------------
      // PROCESS
      // -----------------
      // {

      /// \brief initialize data of the Merge Trees jt & st
      template <typename scalarType>
      void initDataMT(void);

      // }
   };

   ostream &operator<<(ostream &o, Node const &n);
   ostream &operator<<(ostream &o, SuperArc const &a);

#include <ContourTreeTemplate.h>
}

#endif  // CONTOURTREE_H
