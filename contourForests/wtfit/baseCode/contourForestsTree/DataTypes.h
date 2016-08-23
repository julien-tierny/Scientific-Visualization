/**
 * @file DataTypes.h
 * @brief Additional types for ContourTree.
 * @author Gueunet Charles
 * @version 1
 * @date 2016-03-30
 */

#ifndef DATATYPES_H
#define DATATYPES_H

#include <tuple>
#include <limits>


namespace wtfit
{
   /// \brief SuperArc index in vect_superArcs_
   using idSuperArc = long unsigned int;
   /// \brief Node index in vect_nodes_
   using idNode = unsigned int;
   /// \brief Vertex index in scalars_
   using idVertex = int;
   /// \brief Edge index in vect_edgeList_
   using idEdge = int;
   /// \brief Cell index in vect_cellList_
   using idCell = int;

   /// \brief type used to recover Node/Arc in vert2tree SIGNED ONLY
   // Warning, in long long int the max super arc is -1, might not be able to deal with
   // too large data
   using idCorresp = long long int;

   /// \brief type use to store threads related numbers
   using numThread = unsigned char;
   /// \brief index of the interface/partition in vect_interfaces_
   using idInterface = numThread;
   using idPartition = numThread;
   /// \brief index of a connected component of interface (max 65 535 CC)
   using idCC = unsigned long int;
   /// \brief short type of an edge
   using EdgeType = std::pair<idVertex, idVertex>;
   /// \brief vtkIdCell
   using vtkIdCell = long long int;

   /// \brief type stored by UnionFind, TODO make it unsigned
   using ufDataType = long int;

   // Global null

   // QUESTION impact on performance using max (0 would be faster alloacted)
   static const idSuperArc     nullSuperArc       = std::numeric_limits<idSuperArc>::max();
   static const idNode         nullNodes          = std::numeric_limits<idNode>::max();
   static const idVertex       nullVertex         = std::numeric_limits<idVertex>::max();
   static const idCorresp      nullCorresp        = std::numeric_limits<idCorresp>::max();
   static const idInterface    nullInterface      = std::numeric_limits<idInterface>::max();
   static const idPartition    nullPartition      = std::numeric_limits<idPartition>::max();
   static const ufDataType     nullUfData         = std::numeric_limits<ufDataType>::max();
   static constexpr ufDataType specialUfData      = std::numeric_limits<ufDataType>::max() - 1;

   enum ComponentState : char { VISIBLE, HIDDEN, PRUNED, MERGED };

   enum SimplifMethod : char { Persist=0, Span,  NbVert, NbArc };
}

#endif /* end of include guard: DATATYPES_H */
