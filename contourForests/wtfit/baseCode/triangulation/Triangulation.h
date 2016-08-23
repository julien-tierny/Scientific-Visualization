/// \ingroup baseCode
/// \class wtfit::Triangulation
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date January 2016.
///
/// \brief Triangulation is a class that provides time and memory efficient 
/// traversal methods on triangulations of piecewise linear manifolds. It 
/// provides the following features:
///   -# Given a vertex, it provides: the list of edges that are connected to 
/// it, the list of its neighbors, its link, its star, etc.
///   -# Given an edge, it provides: its vertices, its star, etc.
///   -# Given a triangle, its provides: its vertices, its edges, etc.
///   -# Given a tetrahedron, its provides: its vertices, its edges, its 
/// neighbor tetrahedra, etc.
///   -# Given a triangulation, it provides: its list of vertices, edges, 
/// triangles and tetrahedra.
///
/// Triangulation supports both explicit and implicit triangulations: 
///   -# Explicit triangulations: Given a list of points and a list of cells,
/// Triangulation provides time efficient accesses (requiring adequate 
/// pre-processing, see the documentation furtherdown).
///   -# Implicit triangulations: Given a regular grid (origin, spacings and
/// dimensions), Triangulation will perform an implicit triangulation of the 
/// grid, enabling both time and memory efficient traversals of triangulations 
/// of regular grids.
///
/// Apart from pre-processes, Triangulation requires no memory overhead in
/// addition to the input data. 
/// 
/// \note
/// Only pre-process the information you need! See the documentation further 
/// down.
/// \sa vtkTriangulation

#ifndef _TRIANGULATION_H
#define _TRIANGULATION_H

// base code includes
#include                  <AbstractTriangulation.h>
#include                  <ImplicitTriangulation.h>
#include                  <ExplicitTriangulation.h>

namespace wtfit{
  
  class Triangulation : public AbstractTriangulation{

    public:
        
      Triangulation();
      
      ~Triangulation();

      /// Get the \p localEdgeId-th edge of the \p cellId-th cell.
      /// 
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessCellEdges() needs to be called on this object prior to any 
      /// traversal, in a clearly distinct pre-processing step that involves no 
      /// traversal at all. An error will be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      ///
      /// \param cellId Input global cell identifier.
      /// \param localEdgeId Input local edge identifier, 
      /// in [0, getCellEdgeNumber()].
      /// \param edgeId Output global edge identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getCellEdgeNumber()
      inline int getCellEdge(const int &cellId,
        const int &localEdgeId, int &edgeId) const{
         
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedCellEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "CellEdge query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessCellEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getCellEdge(
            cellId, localEdgeId, edgeId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedCellEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "CellEdge query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessCellEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getCellEdge(
            cellId, localEdgeId, edgeId);
        }
        
        return 0;
      }

      /// Get the number of edges for the \p cellId-th cell.
      /// 
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// \pre For this function to behave correctly, preprocessCellEdges() 
      /// needs to be called on this object prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param cellId Input global cell identifier.
      /// \return Returns the number of cell edges.
      inline int getCellEdgeNumber(const int &cellId) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedCellEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "CellEdgeNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessCellEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getCellEdgeNumber(cellId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze          
          if(!implicitTriangulation_.hasPreprocessedCellEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "CellEdgeNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessCellEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getCellEdgeNumber(cellId);
        }
        return 0;
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of edges for all cells.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// The number of entries in this list is equal to the number of cells.
      /// Each entry is a vector of identifiers whose size is equal to the 
      /// number of edges for the corresponding cell.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessCellEdges() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the cell edge list.
      inline const vector<vector<int> > *getCellEdges(){
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedCellEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "CellEdges query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessCellEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return explicitTriangulation_.getCellEdges();
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedCellEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "CellEdges query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessCellEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return implicitTriangulation_.getCellEdges();        
        }
        return NULL;
      }
      
      /// Get the \p localNeighborId-th cell neighbor of the \p cellId-th cell.
      /// 
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessCellNeighbors() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      ///
      /// \param cellId Input global cell identifier.
      /// \param localNeighborId Input local neighbor identifier, 
      /// in [0, getCellNeighborNumber()].
      /// \param neighborId Output global neighbor cell identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getCellNeighborNumber()
      inline int getCellNeighbor(const int &cellId,
        const int &localNeighborId, int &neighborId) const{
         
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedCellNeighbors()){
            stringstream msg;
            msg << "[Triangulation] "
              << "CellNeighbor query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessCellNeighbors() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getCellNeighbor(
            cellId, localNeighborId, neighborId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedCellNeighbors()){
            stringstream msg;
                  msg << "[Triangulation] "
          << "CellNeighbor query without pre-process!"
          << endl;
                  msg << "[Triangulation] "
          << "Please call preprocessCellNeighbors() in a"
          << " pre-process." << endl;
                  dMsg(cerr, msg.str(), Debug::fatalMsg);
          }
#endif
          return implicitTriangulation_.getCellNeighbor(
            cellId, localNeighborId, neighborId);
        }
        return 0;
      }
      
      /// Get the number of cell neighbors for the \p cellId-th cell.
      /// 
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessCellNeighbors() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param cellId Input global cell identifier.
      /// \return Returns the number of cell neighbors.
      inline int getCellNeighborNumber(const int &cellId) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedCellNeighbors()){
            stringstream msg;
            msg << "[Triangulation] "
              << "CellNeighborNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessCellNeighbors() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getCellNeighborNumber(cellId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedCellNeighbors()){
            stringstream msg;
            msg << "[Triangulation] "
              << "CellNeighborNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessCellNeighbors() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getCellNeighborNumber(cellId);
        }
        return 0;
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of cell neighbors for all cells.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// The number of entries in this list is equal to the number of cells.
      /// Each entry is a vector of identifiers whose size is equal to the 
      /// number of neighbor cells for the corresponding cell.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessCellNeighbors() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the cell neighbor list.
      inline const vector<vector<int> > *getCellNeighbors(){
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedCellNeighbors()){
            stringstream msg;
            msg << "[Triangulation] "
              << "CellNeighbors query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessCellNeighbors() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return explicitTriangulation_.getCellNeighbors();
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedCellNeighbors()){
            stringstream msg;
            msg << "[Triangulation] "
              << "CellNeighbors query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessCellNeighbors() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return implicitTriangulation_.getCellNeighbors();        
        }
        return NULL;
      }
      
      /// Get the \p localTriangleId-th triangle id of the \p cellId-th cell.
      /// 
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// Also, the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessCellTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      ///
      /// \param cellId Input global cell identifier.
      /// \param localTriangleId Input local triangle identifier, 
      /// in [0, getCellTriangleNumber()].
      /// \param triangleId Output global triangle identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getCellTriangleNumber()
      inline int getCellTriangle(const int &cellId,
        const int &localTriangleId, int &triangleId) const{
         
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedCellTriangles()){
            stringstream msg;
            msg << "[Triangulation] "
              << "CellTriangle query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessCellTriangles() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getCellTriangle(
            cellId, localTriangleId, triangleId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedCellTriangles()){
            stringstream msg;
                  msg << "[Triangulation] "
          << "CellTriangle query without pre-process!"
          << endl;
                  msg << "[Triangulation] "
          << "Please call preprocessCellTriangles() in a"
          << " pre-process." << endl;
                  dMsg(cerr, msg.str(), Debug::fatalMsg);
          }
#endif
          return implicitTriangulation_.getCellTriangle(
            cellId, localTriangleId, triangleId);
        }
        return 0;
      }
      
      /// Get the number of triangles for the \p cellId-th cell.
      /// 
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// Also, the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessCellTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param cellId Input global cell identifier.
      /// \return Returns the number of cell triangles.
      inline int getCellTriangleNumber(const int &cellId) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedCellTriangles()){
            stringstream msg;
            msg << "[Triangulation] "
              << "CellTriangleNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessCellTriangles() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getCellTriangleNumber(cellId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedCellTriangles()){
            stringstream msg;
            msg << "[Triangulation] "
              << "CellTriangleNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessCellTriangles() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getCellTriangleNumber(cellId);
        }
        return 0;
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of triangles for all cells.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// Also, the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// The number of entries in this list is equal to the number of cells.
      /// Each entry is a vector of identifiers whose size is equal to the 
      /// number of triangles for the corresponding cell.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessCellTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the cell triangle list.
      inline const vector<vector<int> > *getCellTriangles(){
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedCellTriangles()){
            stringstream msg;
            msg << "[Triangulation] "
              << "CellTriangles query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessCellTriangles() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return explicitTriangulation_.getCellTriangles();
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedCellTriangles()){
            stringstream msg;
            msg << "[Triangulation] "
              << "CellTriangles query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessCellTriangles() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return implicitTriangulation_.getCellTriangles();        
        }
        return NULL;
      }
      
      /// Get the \p localVertexId-th vertex identifier of the \p cellId-th 
      /// cell.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      /// \param cellId Input global cell identifier.
      /// \param localVertexId Input local vertex identifier,
      /// in [0, getCellVertexNumber()].
      /// \param vertexId Ouput global vertex identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getCellVertexNumber()
      inline int getCellVertex(const int &cellId,
        const int &localVertexId, int &vertexId) const{
        if(!explicitTriangulation_.isEmpty()){
          return explicitTriangulation_.getCellVertex(
            cellId, localVertexId, vertexId);
        }
        else if(!implicitTriangulation_.isEmpty()){
          return implicitTriangulation_.getCellVertex(
            cellId, localVertexId, vertexId);
        }
        return 0;
      }
      
      /// Get the number of vertices in a cell.
      /// 
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      /// \param cellId Input global cell identifier.
      /// \returns Number of vertices in the cell.
      inline int getCellVertexNumber(const int &cellId) const{
        if(!explicitTriangulation_.isEmpty()){
          return explicitTriangulation_.getCellVertexNumber(cellId);
        } 
        else if(!implicitTriangulation_.isEmpty()){
          return implicitTriangulation_.getCellVertexNumber(cellId);
        }
        return 0;
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of edges of the triangulation.
      ///
      /// Here the notion of edge only makes sense if the triangulation has a 
      /// dimension greater than 1 (otherwise, use the cell information).
      ///
      /// The number of entries in this list is equal to the number of edges.
      /// Each entry is a pair of vertex identifiers.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      /// \pre For this function to behave correctly, 
      /// preprocessEdges() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the edge list.
      inline const vector<pair<int, int> > *getEdges(){
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "Edges query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return explicitTriangulation_.getEdges();
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "Edges query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return implicitTriangulation_.getEdges();
        }
        return NULL;
      }
      
      /// Get the \p localLinkId-th cell of the link of the \p edgeId-th 
      /// edge.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension IN THE LINK (3D: triangles, 2D: edges, 1D: vertices).
      ///
      /// The link cell is given as a vector of integers. The first one 
      /// refers to the number of vertices in the cell, while the following
      /// ones refer to the vertex identifiers of the cell.
      /// \pre For this function to behave correctly, 
      /// preprocessEdgeLinks() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param edgeId Input global edge identifier.
      /// \param localLinkId Input local link simplex identifier, 
      /// in [0, getEdgeLinkNumber()].
      /// \param link Output link cell.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getEdgeLinkNumber()
      /// \warning This function is not implemented in this version of the API 
      /// (it is a placeholder for a future version).
      inline int getEdgeLink(const int &edgeId, 
        const int &localLinkId, vector<long long int> &link) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedEdgeLinks()){
            stringstream msg;
            msg << "[Triangulation] "
              << "EdgeLink query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdgeLinks() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getEdgeLink(
            edgeId, localLinkId, link);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedEdgeLinks()){
            stringstream msg;
            msg << "[Triangulation] "
              << "EdgeLink query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdeLinks() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getEdgeLink(
            edgeId, localLinkId, link);          
        }
        return 0;
      }
      
      /// Get the number of cells in the link of the \p edgeId-th edge.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension IN THE LINK (3D: triangles, 2D: edges, 1D: vertices).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessEdgeLinks() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param edgeId Input global edge identifier.
      /// \return Returns the number of cells in the link of the edge. 
      /// \warning This function is not implemented in this version of the API 
      /// (it is a placeholder for a future version).
      inline int getEdgeLinkNumber(const int &edgeId) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedEdgeLinks()){
            stringstream msg;
            msg << "[Triangulation] "
              << "EdgeLinkNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdgeLinks() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getEdgeLinkNumber(edgeId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedEdgeLinks()){
            stringstream msg;
            msg << "[Triangulation] "
              << "EdgeLinkNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdgeLinks() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getEdgeLinkNumber(edgeId);
        }
        return 0;
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of link cells for all edges.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension IN THE LINK (3D: triangles, 2D: edges, 1D: vertices).
      ///
      /// The number of entries in this list is equal to the number of edges.
      /// Each entry is a vector of identifiers representing a cell.
      /// The first one refers to the number of vertices in the cell, while
      /// the following ones refer to the vertex identifiers of the cell.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessEdgeLinks() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the edge link list.
      /// \warning This function is not implemented in this version of the API 
      /// (it is a placeholder for a future version).
      inline const vector<vector<long long int> > *getEdgeLinks(){
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedEdgeLinks()){
            stringstream msg;
            msg << "[Triangulation] "
              << "EdgeLinks query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdgeLinks() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return explicitTriangulation_.getEdgeLinks();
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedEdgeLinks()){
            stringstream msg;
            msg << "[Triangulation] "
              << "EdgeLinks query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdgeLinks() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return implicitTriangulation_.getEdgeLinks();  
        }
        return 0;
      }
      
      /// Get the \p localStarId-th cell of the star of the \p edgeId-th edge.
      /// 
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// Also, the notion of edge only makes sense if the triangulation has a 
      /// dimension greater than 1 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessEdgeStars() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      ///
      /// \param edgeId Input global edge identifier.
      /// \param localStarId Input local star cell identifier,
      /// in [0, getEdgeStarNumber()].
      /// \param starId Output global star cell identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getEdgeStarNumber()
      inline int getEdgeStar(const int &edgeId,
        const int &localStarId, int &starId) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedEdgeStars()){
            stringstream msg;
            msg << "[Triangulation] "
              << "EdgeStar query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdgeStars() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getEdgeStar(
            edgeId, localStarId, starId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedEdgeStars()){
            stringstream msg;
            msg << "[Triangulation] "
              << "EdgeStar query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdgeStars() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getEdgeStar(
            edgeId, localStarId, starId);         
        }
        return 0;
      }
      
      /// Get the number of star cells for the \p edgeId-th edge.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// Also, the notion of edge only makes sense if the triangulation has a 
      /// dimension greater than 1 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessEdgeStars() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param edgeId Input global edge identifier
      /// \return Returns the number of star cells.
      inline int getEdgeStarNumber(const int &edgeId) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedEdgeStars()){
            stringstream msg;
            msg << "[Triangulation] "
              << "EdgeStarNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdgeStars() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getEdgeStarNumber(edgeId);
        } 
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedEdgeStars()){
            stringstream msg;
            msg << "[Triangulation] "
              << "EdgeStarNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdgeStars() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getEdgeStarNumber(edgeId);          
        }
        return 0;
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of star cell identifiers for all edges.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// Also, the notion of edge only makes sense if the triangulation has a 
      /// dimension greater than 1 (otherwise, use the cell information).
      ///
      /// The number of entries in this list is equal to the number of edges.
      /// Each entry is a vector of identifiers whose size is equal to the 
      /// number of star cells for the corresponding edge.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessEdgeStars() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the edge star list.
      inline const vector<vector<int> > *getEdgeStars(){
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedEdgeStars()){
            stringstream msg;
            msg << "[Triangulation] "
              << "EdgeStars query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdgeStars() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return explicitTriangulation_.getEdgeStars();
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedEdgeStars()){
            stringstream msg;
            msg << "[Triangulation] "
              << "EdgeStars query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdgeStars() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return implicitTriangulation_.getEdgeStars();
        }
        return NULL;
      }
      
      /// Get the \p localTriangleId-th triangle id of the \p edgeId-th edge.
      /// 
      /// Here the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessEdgeTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      ///
      /// \param edgeId Input global edge identifier.
      /// \param localTriangleId Input local triangle identifier, 
      /// in [0, getEdgeTriangleNumber()].
      /// \param triangleId Output global triangle identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getEdgeTriangleNumber()
      inline int getEdgeTriangle(const int &edgeId,
        const int &localTriangleId, int &triangleId) const{
         
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedEdgeTriangles()){
            stringstream msg;
            msg << "[Triangulation] "
              << "EdgeTriangle query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdgeTriangles() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getEdgeTriangle(
            edgeId, localTriangleId, triangleId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedEdgeTriangles()){
            stringstream msg;
                  msg << "[Triangulation] "
          << "EdgeTriangle query without pre-process!"
          << endl;
                  msg << "[Triangulation] "
          << "Please call preprocessEdgeTriangles() in a"
          << " pre-process." << endl;
                  dMsg(cerr, msg.str(), Debug::fatalMsg);
          }
#endif
          return implicitTriangulation_.getEdgeTriangle(
            edgeId, localTriangleId, triangleId);
        }
        return 0;
      }
      
      /// Get the number of triangles for the \p edgeId-th edge.
      /// 
      /// Here the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessEdgeTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param edgeId Input global edge identifier.
      /// \return Returns the number of edge triangles.
      inline int getEdgeTriangleNumber(const int &edgeId) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedEdgeTriangles()){
            stringstream msg;
            msg << "[Triangulation] "
              << "EdgeTriangleNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdgeTriangles() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getEdgeTriangleNumber(edgeId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedEdgeTriangles()){
            stringstream msg;
            msg << "[Triangulation] "
              << "EdgeTriangleNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdgeTriangles() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getEdgeTriangleNumber(edgeId);
        }
        return 0;
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of triangles for all edges.
      ///
      /// Here the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// The number of entries in this list is equal to the number of edges.
      /// Each entry is a vector of identifiers whose size is equal to the 
      /// number of triangles for the corresponding edge.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessEdgeTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the edge triangle list.
      inline const vector<vector<int> > *getEdgeTriangles(){
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedEdgeTriangles()){
            stringstream msg;
            msg << "[Triangulation] "
              << "EdgeTriangles query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdgeTriangles() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return explicitTriangulation_.getEdgeTriangles();
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedEdgeTriangles()){
            stringstream msg;
            msg << "[Triangulation] "
              << "EdgeTriangles query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdgeTriangles() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return implicitTriangulation_.getEdgeTriangles();        
        }
        return NULL;
      }
      
      /// Get the \p localVertexId-th vertex identifier of the \p edgeId-th 
      /// edge.
      ///
      /// Here the notion of edge only makes sense if the triangulation has a 
      /// dimension greater than 1 (otherwise, use the cell information).
      /// \pre For this function to behave correctly, 
      /// preprocessEdges() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param edgeId Input global edge identifier.
      /// \param localVertexId Input local vertex identifier (0 or 1).
      /// \param vertexId Output global vertex identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      inline int getEdgeVertex(const int &edgeId, 
        const int &localVertexId, int &vertexId) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "EdgeVertex query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getEdgeVertex(
            edgeId, localVertexId, vertexId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "EdgeVertex query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getEdgeVertex(
            edgeId, localVertexId, vertexId);          
        }
        return 0;
      }
      
      /// Get the number of cells in the triangulation.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      /// \return Returns the number of cells.
      inline int getNumberOfCells() const{
        if(!explicitTriangulation_.isEmpty()){
          return explicitTriangulation_.getNumberOfCells();
        }
        else if(!implicitTriangulation_.isEmpty()){
          return implicitTriangulation_.getNumberOfCells();
        }
        return 0;
      }
      
      /// Get the number of edges in the triangulation.
      ///
      /// Here the notion of edge only makes sense if the triangulation has a 
      /// dimension greater than 1 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessEdges() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns the number of edges.
      inline int getNumberOfEdges() const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "NumberOfEdges query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getNumberOfEdges();
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "NumberOfEdges query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getNumberOfEdges();
        }
        return 0;
      }
      
      /// Get the number of triangles in the triangulation.
      ///
      /// Here the notion of triangle only makes sense if the triangulation has
      /// a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns the number of triangles.
      inline int getNumberOfTriangles() const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedTriangles()){
            stringstream msg;
            msg << "[Triangulation] "
              << "NumberOfTriangles query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangles() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getNumberOfTriangles();
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedTriangles()){
            stringstream msg;
            msg << "[Triangulation] "
              << "NumberOfTriangles query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangles() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getNumberOfTriangles();          
        }
        return 0;
      }
      
      /// Get the number of vertices in the triangulation.
      /// \return Returns the number of vertices.
      inline int getNumberOfVertices() const{
        if(!explicitTriangulation_.isEmpty()){
          return explicitTriangulation_.getNumberOfVertices();
        }
        else if(!implicitTriangulation_.isEmpty()){
          return implicitTriangulation_.getNumberOfVertices();
        }
        return 0;
      }
     
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of triangles of the triangulation.
      ///
      /// Here the notion of triangle only makes sense if the triangulation has 
      /// a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// The number of entries in this list is equal to the number of 
      /// triangles.
      /// Each entry is a vector of vertex identifiers.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      /// \pre For this function to behave correctly, 
      /// preprocessTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the triangle list.
      inline const vector<vector<int> > *getTriangles(){
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedTriangles()){
            stringstream msg;
            msg << "[Triangulation] "
              << "Triangles query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangles() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return explicitTriangulation_.getTriangles();
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedTriangles()){
            stringstream msg;
            msg << "[Triangulation] "
              << "Triangles query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangles() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return implicitTriangulation_.getTriangles();
        }
        return NULL;
      }
      
      /// Get the \p localEdgeId-th edge of the \p triangleId-th triangle.
      ///
      /// Here the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessTriangleEdges() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param triangleId Input global triangle identifier.
      /// \param localEdgeId Input local edge identifier, 
      /// in [0, getTriangleEdgeNumber()].
      /// \param edgeId Output global edge identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getTriangleEdgeNumber()
      inline int getTriangleEdge(const int &triangleId, 
        const int &localEdgeId, int &edgeId) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedTriangleEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "TriangleEdge query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangleEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getTriangleEdge(
            triangleId, localEdgeId, edgeId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedTriangleEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "TriangleEdge query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangleEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getTriangleEdge(
            triangleId, localEdgeId, edgeId);
        }
        return 0;
      }
      
      /// Get the number of edges of the \p triangleId-th triangle.
      ///
      /// Here, the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessTriangleEdges() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param triangleId Input global triangle identifier.
      /// \return Returns the number of cells in the link of the triangle. 
      inline int getTriangleEdgeNumber(const int &triangleId) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedTriangleEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "TriangleEdgeNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangleEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getTriangleEdgeNumber(triangleId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedTriangleEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "TriangleEdgeNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangleEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getTriangleEdgeNumber(triangleId);
        }
        return 0;
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of edges for all triangles.
      ///
      /// Here the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// The number of entries in this list is equal to the number of 
      /// triangles. Each entry is a vector of identifiers representing the
      /// edges connected to the triangle (3).
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessTriangleEdges() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the triangle edge list.
      inline const vector<vector<int> > *getTriangleEdges(){
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedTriangleEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "TriangleEdges query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangleEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return explicitTriangulation_.getTriangleEdges();
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedTriangleEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "TriangleEdges query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangleEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return implicitTriangulation_.getTriangleEdges();  
        }
        return 0;
      }
      
      /// Get the \p localLinkId-th cell of the link of the \p triangleId-th 
      /// triangle.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension IN THE LINK (3D: triangles, 2D: edges, 1D: vertices).
      ///
      /// Also, the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// The link cell is given as a vector of integers. The first one 
      /// refers to the number of vertices in the cell, while the following
      /// ones refer to the vertex identifiers of the cell.
      /// \pre For this function to behave correctly, 
      /// preprocessTriangleLinks() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param triangleId Input global triangle identifier.
      /// \param localLinkId Input local link simplex identifier, 
      /// in [0, getTriangleLinkNumber()].
      /// \param link Output link cell.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getTriangleLinkNumber()
      /// \warning This function is not implemented in this version of the API 
      /// (it is a placeholder for a future version).
      inline int getTriangleLink(const int &triangleId, 
        const int &localLinkId, vector<long long int> &link) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedTriangleLinks()){
            stringstream msg;
            msg << "[Triangulation] "
              << "TriangleLink query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangleLinks() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getTriangleLink(
            triangleId, localLinkId, link);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedTriangleLinks()){
            stringstream msg;
            msg << "[Triangulation] "
              << "TriangleLink query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangleLinks() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getTriangleLink(
            triangleId, localLinkId, link);          
        }
        return 0;
      }
      
      /// Get the number of cells in the link of the \p triangleId-th triangle.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension IN THE LINK (3D: triangles, 2D: edges, 1D: vertices).
      ///
      /// Also, the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessTriangleLinks() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param triangleId Input global triangle identifier.
      /// \return Returns the number of cells in the link of the triangle. 
      /// \warning This function is not implemented in this version of the API 
      /// (it is a placeholder for a future version).
      inline int getTriangleLinkNumber(const int &triangleId) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedTriangleLinks()){
            stringstream msg;
            msg << "[Triangulation] "
              << "TriangleLinkNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangleLinks() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getTriangleLinkNumber(triangleId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedTriangleLinks()){
            stringstream msg;
            msg << "[Triangulation] "
              << "TriangleLinkNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangleLinks() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getTriangleLinkNumber(triangleId);
        }
        return 0;
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of link cells for all triangle.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension IN THE LINK (3D: triangles, 2D: edges, 1D: vertices).
      ///
      /// Also, the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// The number of entries in this list is equal to the number of 
      /// triangles.
      /// Each entry is a vector of identifiers representing a cell.
      /// The first one refers to the number of vertices in the cell, while
      /// the following ones refer to the vertex identifiers of the cell.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessTriangleLinks() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the triangle link list.
      /// \warning This function is not implemented in this version of the API 
      /// (it is a placeholder for a future version).
      inline const vector<vector<long long int> > *getTriangleLinks(){
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedTriangleLinks()){
            stringstream msg;
            msg << "[Triangulation] "
              << "TriangleLinks query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangleLinks() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return explicitTriangulation_.getTriangleLinks();
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedTriangleLinks()){
            stringstream msg;
            msg << "[Triangulation] "
              << "TriangleLinks query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangleLinks() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return implicitTriangulation_.getTriangleLinks();  
        }
        return 0;
      }
     
      /// Get the \p localStarId-th cell of the star of the \p triangleId-th 
      /// triangle.
      /// 
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// Also, the notion of triangle only makes sense if the triangulation has
      /// a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessTriangleStars() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      ///
      /// \param triangleId Input global triangle identifier.
      /// \param localStarId Input local star cell identifier,
      /// in [0, getTriangleStarNumber()].
      /// \param starId Output global star cell identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getTriangleStarNumber()
      inline int getTriangleStar(const int &triangleId,
        const int &localStarId, int &starId) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedTriangleStars()){
            stringstream msg;
            msg << "[Triangulation] "
              << "TriangleStar query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangleStars() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getTriangleStar(
            triangleId, localStarId, starId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedTriangleStars()){
            stringstream msg;
            msg << "[Triangulation] "
              << "TriangleStar query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangleStars() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getTriangleStar(
            triangleId, localStarId, starId);          
        }
        return 0;
      }
     
      /// Get the number of star cells for the \p triangleId-th triangle.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// Also,the notion of triangle only makes sense if the triangulation has 
      /// a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessTriangleStars() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param triangleId Input global triangle identifier.
      /// \return Returns the number of star cells.
      inline int getTriangleStarNumber(const int &triangleId) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedTriangleStars()){
            stringstream msg;
            msg << "[Triangulation] "
              << "TriangleStarNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangleStars() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getTriangleStarNumber(triangleId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedTriangleStars()){
            stringstream msg;
            msg << "[Triangulation] "
              << "TriangleStarNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangleStars() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getTriangleStarNumber(triangleId);          
        }
        return 0;
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of star cell identifiers for all triangles.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// Also, the notion of triangle only makes sense if the triangulation
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// The number of entries in this list is equal to the number of 
      /// triangles.
      /// Each entry is a vector of identifiers whose size is equal to the 
      /// number of star cells for the corresponding triangle.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessTriangleStars() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the triangle star list.
      inline const vector<vector<int> > *getTriangleStars(){
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedTriangleStars()){
            stringstream msg;
            msg << "[Triangulation] "
              << "TriangleStars query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangleStars() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return explicitTriangulation_.getTriangleStars();
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedTriangleStars()){
            stringstream msg;
            msg << "[Triangulation] "
              << "TriangleStars query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangleStars() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return implicitTriangulation_.getTriangleStars();         
        }
        return NULL;
      }
     
      /// Get the \p localVertexId-th vertex identifier of the \p triangleId-th 
      /// triangle.
      ///
      /// Here the notion of triangle only makes sense if the triangulation has 
      /// a dimension greater than 2 (otherwise, use the cell information).
      /// \pre For this function to behave correctly, 
      /// preprocessTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param triangleId Input global edge identifier.
      /// \param localVertexId Input local vertex identifier (in [0, 2]).
      /// \param vertexId Output global vertex identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      inline int getTriangleVertex(const int &triangleId,
        const int &localVertexId, int &vertexId) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedTriangles()){
            stringstream msg;
            msg << "[Triangulation] "
              << "TriangleVertex query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangles() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getTriangleVertex(
            triangleId, localVertexId, vertexId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedTriangles()){
            stringstream msg;
            msg << "[Triangulation] "
              << "TriangleVertex query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessTriangles() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getTriangleVertex(
            triangleId, localVertexId, vertexId);          
        }
        return 0;
      }
     
      /// Get the \p localEdgeId-th edge identifier connected to the 
      /// \p vertexId-th 
      /// vertex.
      ///
      /// Here the notion of edge only makes sense if the triangulation has a 
      /// dimension greater than 1 (otherwise, use the cell information).
      /// \pre For this function to behave correctly, 
      /// preprocessVertexEdges() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param vertexId Input global vertex identifier.
      /// \param localEdgeId Input local edge identifier,
      /// in [0, getVertexEdgeNumber()].
      /// \param edgeId Output global edge identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getVertexEdgeNumber()
      inline int getVertexEdge(const int &vertexId, 
        const int &localEdgeId, int &edgeId) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedVertexEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexEdge query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getVertexEdge(
            vertexId, localEdgeId, edgeId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedVertexEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexEdge query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getVertexEdge(
            vertexId, localEdgeId, edgeId);          
        }
        return 0;
      }
     
      /// Get the number of edges connected to the \p vertexId-th vertex.
      ///
      /// Here,the notion of edge only makes sense if the triangulation has 
      /// a dimension greater than 1 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexEdges() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param vertexId Input global vertex identifier.
      /// \return Returns the number of edges connected to the vertex.
      inline int getVertexEdgeNumber(const int &vertexId) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedVertexEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexEdgeNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getVertexEdgeNumber(vertexId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedVertexEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexEdgeNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getVertexEdgeNumber(vertexId);
        }
        return 0;
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of edge identifiers for all vertices.
      ///
      /// Here, the notion of edge only makes sense if the triangulation
      /// has a dimension greater than 1 (otherwise, use the cell information).
      ///
      /// The number of entries in this list is equal to the number of 
      /// vertices.
      /// Each entry is a vector of identifiers whose size is equal to the 
      /// number of edges connected to the corresponding vertex.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexEdges() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the vertex edge list.
      inline const vector<vector<int> > *getVertexEdges(){
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedVertexEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexEdges query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return explicitTriangulation_.getVertexEdges();
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedVertexEdges()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexEdges query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexEdges() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return implicitTriangulation_.getVertexEdges();          
        }
        return 0;
      }
      
      /// Get the \p localLinkId-th cell of the link of the \p vertexId-th 
      /// vertex.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension IN THE LINK (3D: triangles, 2D: edges, 1D: vertices).
      ///
      /// The link cell is given as a vector of integers. The first one 
      /// refers to the number of vertices in the cell, while the following
      /// ones refer to the vertex identifiers of the cell.
      /// \pre For this function to behave correctly, 
      /// preprocessVertexLinks() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param vertexId Input global vertex identifier.
      /// \param localLinkId Input local link simplex identifier, 
      /// in [0, getVertexLinkNumber()].
      /// \param link Output link cell.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getVertexLinkNumber()
      inline int getVertexLink(const int &vertexId, 
        const int &localLinkId, vector<long long int> &link) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedVertexLinks()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexLink query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexLinks() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getVertexLink(
            vertexId, localLinkId, link);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedVertexLinks()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexLink query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexLinks() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getVertexLink(
            vertexId, localLinkId, link);          
        }
        return 0;
      }
      
      /// Get the number of cells in the link of the \p vertexId-th vertex.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension IN THE LINK (3D: triangles, 2D: edges, 1D: vertices).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexLinks() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param vertexId Input global vertex identifier.
      /// \return Returns the number of cells in the link of the vertex. 
      inline int getVertexLinkNumber(const int &vertexId) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedVertexLinks()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexLinkNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexLinks() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getVertexLinkNumber(vertexId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedVertexLinks()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexLinkNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexLinks() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getVertexLinkNumber(vertexId);          
        }
        return 0;
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of link cells for all vertices.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension IN THE LINK (3D: triangles, 2D: edges, 1D: vertices).
      ///
      /// The number of entries in this list is equal to the number of 
      /// vertices.
      /// Each entry is a vector of identifiers representing a cell.
      /// The first one refers to the number of vertices in the cell, while
      /// the following ones refer to the vertex identifiers of the cell.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexLinks() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the vertex link list.
      inline const vector<vector<long long int> > *getVertexLinks(){
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedVertexLinks()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexLinks query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexLinks() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return explicitTriangulation_.getVertexLinks();
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedVertexLinks()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexLinks query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexLinks() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return implicitTriangulation_.getVertexLinks();  
        }
        return 0;
      }
      
      /// Get the \p localNeighborId-th vertex neighbor of the \p vertexId-th 
      /// vertex.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexNeighbors() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param vertexId Input global vertex identifier.
      /// \param localNeighborId Input local neighbor identifier,
      /// in [0, getVertexNeighborNumber()].
      /// \param neighborId Output global neighbor vertex identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getVertexNeighborNumber()
      inline int getVertexNeighbor(const int &vertexId, 
        const int &localNeighborId, int &neighborId) const{
          
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedVertexNeighbors()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexNeighbor query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexNeighbors() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getVertexNeighbor(
            vertexId, localNeighborId, neighborId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedVertexNeighbors()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexNeighbor query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexNeighbors() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getVertexNeighbor(
            vertexId, localNeighborId, neighborId);
        }
        return 0;
      }
      
      /// Get the number of vertex neighbors for the \p vertexId-th vertex.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexNeighbors() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param vertexId Input global vertex identifier.
      /// \return Returns the number vertex neighbors.
      inline int getVertexNeighborNumber(const int &vertexId) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedVertexNeighbors()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexNeighborNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexNeighbors() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif     
          return explicitTriangulation_.getVertexNeighborNumber(vertexId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedVertexNeighbors()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexNeighborNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexNeighbors() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif     
          return implicitTriangulation_.getVertexNeighborNumber(vertexId);          
        }
        return 0;
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of vertex neighbor identifiers for all vertices.
      ///
      /// The number of entries in this list is equal to the number of 
      /// vertices.
      /// Each entry is a vector of identifiers whose size is equal to the 
      /// number of vertex neighbors for the corresponding vertex.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexNeighbors() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the vertex neighbor list.
      inline const vector<vector<int> > *getVertexNeighbors(){
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedVertexNeighbors()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexNeighbors query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexNeighbors() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif 
          return explicitTriangulation_.getVertexNeighbors();
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedVertexNeighbors()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexNeighbors query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexNeighbors() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif 
          return implicitTriangulation_.getVertexNeighbors();          
        }
        return 0;
      }
      
      /// Get the point (3D coordinates) for the \p vertexId-th vertex.
      /// \param vertexId Input global vertex identifier.
      /// \param x Output x coordinate.
      /// \param y Output y coordinate.
      /// \param z Output z coordinate.
      /// \return Returns 0 upon success, negative values otherwise.
      inline int getVertexPoint(const int &vertexId, 
        float &x, float &y, float &z) const{
        if(!explicitTriangulation_.isEmpty()){
          return explicitTriangulation_.getVertexPoint(vertexId, x, y, z);
        }
        else if(!implicitTriangulation_.isEmpty()){
          return implicitTriangulation_.getVertexPoint(vertexId, x, y, z);
        }
        return 0;
      }
      
      /// Get the \p localStarId-th cell of the star of the \p vertexId-th 
      /// vertex.
      /// 
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexStars() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      ///
      /// \param vertexId Input global vertex identifier.
      /// \param localStarId Input local star cell identifier,
      /// in [0, getVertexStarNumber()].
      /// \param starId Output global star cell identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getVertexStarNumber()
      inline int getVertexStar(const int &vertexId, 
        const int &localStarId, int &starId) const{
        
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedVertexStars()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexStar query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexStar() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif          
          return explicitTriangulation_.getVertexStar(
            vertexId, localStarId, starId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedVertexStars()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexStar query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexStar() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif          
          return implicitTriangulation_.getVertexStar(
            vertexId, localStarId, starId);          
        }
        return 0;
      }
        
      /// Get the number of star cells for the \p vertexId-th vertex.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexStars() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param vertexId Input global vertex identifier
      /// \return Returns the number of star cells.
      inline int getVertexStarNumber(const int &vertexId) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedVertexStars()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexStarNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexStars() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getVertexStarNumber(vertexId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedVertexStars()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexStarNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexStars() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getVertexStarNumber(vertexId);         
        }
        return 0;
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of star cell identifiers for all vertices.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// The number of entries in this list is equal to the number of vertices.
      /// Each entry is a vector of identifiers whose size is equal to the 
      /// number of star cells for the corresponding vertex.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexStars() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the vertex star list.
      inline const vector<vector<int> > *getVertexStars(){
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedVertexStars()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexStars query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexStars() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif         
          return explicitTriangulation_.getVertexStars();
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedVertexStars()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexStars query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexStars() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif         
          return implicitTriangulation_.getVertexStars();          
        }
        return 0;
      }
      
      /// Get the \p localTriangleId-th triangle id of the 
      /// \p vertexId-th vertex.
      /// 
      /// Here the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      ///
      /// \param vertexId Input global vertex identifier.
      /// \param localTriangleId Input local triangle identifier, 
      /// in [0, getVertexTriangleNumber()].
      /// \param triangleId Output global triangle identifier.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getVertexTriangleNumber()
      /// \warning This function is not implemented in this version of the API 
      /// (it is a placeholder for a future version).
      inline int getVertexTriangle(const int &vertexId,
        const int &localTriangleId, int &triangleId) const{
         
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedVertexTriangles()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexTriangle query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexTriangles() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getVertexTriangle(
            vertexId, localTriangleId, triangleId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedVertexTriangles()){
            stringstream msg;
                  msg << "[Triangulation] "
          << "VertexTriangle query without pre-process!"
          << endl;
                  msg << "[Triangulation] "
          << "Please call preprocessVertexTriangles() in a"
          << " pre-process." << endl;
                  dMsg(cerr, msg.str(), Debug::fatalMsg);
          }
#endif
          return implicitTriangulation_.getVertexTriangle(
            vertexId, localTriangleId, triangleId);
        }
        return 0;
      }
      
      /// Get the number of triangles for the \p vertexId-th vertex.
      /// 
      /// Here the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \param vertexId Input global vertex identifier.
      /// \return Returns the number of vertex triangles.
      /// \warning This function is not implemented in this version of the API 
      /// (it is a placeholder for a future version).
      inline int getVertexTriangleNumber(const int &vertexId) const{
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedVertexTriangles()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexTriangleNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexTriangles() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return explicitTriangulation_.getVertexTriangleNumber(vertexId);
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedVertexTriangles()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexTriangleNumber query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexTriangles() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return -1;
          }
#endif
          return implicitTriangulation_.getVertexTriangleNumber(vertexId);
        }
        return 0;
      }
      
      /// \warning
      /// YOU SHOULD NOT CALL THIS FUNCTION UNLESS YOU REALLY KNOW WHAT YOU ARE
      /// DOING.
      ///
      /// Get the list of triangles for all vertices.
      ///
      /// Here the notion of triangle only makes sense if the triangulation 
      /// has a dimension greater than 2 (otherwise, use the cell information).
      ///
      /// The number of entries in this list is equal to the number of vertices.
      /// Each entry is a vector of identifiers whose size is equal to the 
      /// number of triangles for the corresponding vertex.
      ///
      /// In implicit mode, this function will force the creation of such a 
      /// list (which will be time and memory consuming). 
      /// THIS IS USUALLY A BAD IDEA.
      ///
      /// \pre For this function to behave correctly, 
      /// preprocessVertexTriangles() needs to be called
      /// on this object prior to any traversal, in a clearly distinct 
      /// pre-processing step that involves no traversal at all. An error will 
      /// be returned otherwise.
      /// \note It is recommended to exclude such a pre-processing step 
      /// from any time performance measurement.
      /// \return Returns a pointer to the vertex triangle list.
      /// \warning This function is not implemented in this version of the API 
      /// (it is a placeholder for a future version).
      inline const vector<vector<int> > *getVertexTriangles(){
        if(!explicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!explicitTriangulation_.hasPreprocessedVertexTriangles()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexTriangles query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexTriangles() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return explicitTriangulation_.getVertexTriangles();
        }
        else if(!implicitTriangulation_.isEmpty()){
#ifndef withKamikaze
          if(!implicitTriangulation_.hasPreprocessedVertexTriangles()){
            stringstream msg;
            msg << "[Triangulation] "
              << "VertexTriangles query without pre-process!"
              << endl;
            msg << "[Triangulation] "
              << "Please call preprocessVertexTriangles() in a"
              << " pre-process." << endl;
            dMsg(cerr, msg.str(), Debug::fatalMsg);
            return NULL;
          }
#endif
          return implicitTriangulation_.getVertexTriangles();        
        }
        return NULL;
      }
      
      /// Check if the data structure is empty or not.
      /// \return Returns true if empty, false otherwise.
      inline bool isEmpty() const {
         return (explicitTriangulation_.isEmpty() 
           && implicitTriangulation_.isEmpty());
      }
      
      /// Pre-process the cell edges.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getCellEdge()
      ///   - getCellEdgeNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getCellEdge()
      /// \sa getCellEdgeNumber()
      inline int preprocessCellEdges(){
        if(!explicitTriangulation_.isEmpty()){
          return explicitTriangulation_.preprocessCellEdges();
        }
        else if(!implicitTriangulation_.isEmpty()){
          return implicitTriangulation_.preprocessCellEdges();
        }
        return 0;
      }

      /// Pre-process the cell neighbors.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getCellNeighbor()
      ///   - getCellNeighbors()
      ///   - getCellNeighborNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getCellNeighbor()
      /// \sa getCellNeighbors()
      /// \sa getCellNeighborNumber()
      inline int preprocessCellNeighbors(){
        if(!explicitTriangulation_.isEmpty()){
          return explicitTriangulation_.preprocessCellNeighbors();
        }
        else if(!implicitTriangulation_.isEmpty()){
          return implicitTriangulation_.preprocessCellNeighbors();
        }
        return 0;
      }
      
      /// Pre-process the cell triangles.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getCellTriangle()
      ///   - getCellTriangles()
      ///   - getCellTriangleNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getCellTriangle()
      /// \sa getCellTriangles()
      /// \sa getCellTriangleNumber()
      inline int preprocessCellTriangles(){
        if(!explicitTriangulation_.isEmpty()){
          return explicitTriangulation_.preprocessCellTriangles();
        }
        else if(!implicitTriangulation_.isEmpty()){
          return implicitTriangulation_.preprocessCellTriangles();
        }
        return 0;
      }
      
      /// Pre-process the edges.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getEdges()
      ///   - getEdgeVertex()
      ///   - getNumberOfEdges()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getEdges()
      /// \sa getEdgeVertex()
      /// \sa getNumberOfEdges()
      inline int preprocessEdges(){
        if(!explicitTriangulation_.isEmpty()){
          return explicitTriangulation_.preprocessEdges();
        }
        else if(!implicitTriangulation_.isEmpty()){
          return implicitTriangulation_.preprocessEdges();
        }
        return 0;
      }
      
      /// Pre-process the edge links.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getEdgeLink()
      ///   - getEdgeLinks()
      ///   - getEdgeLinkNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getEdgeLink()
      /// \sa getEdgeLinks()
      /// \sa getEdgeLinkNumber()
      inline int preprocessEdgeLinks(){
        if(!explicitTriangulation_.isEmpty()){
          return explicitTriangulation_.preprocessEdgeLinks();
        }
        else if(!implicitTriangulation_.isEmpty()){
          return implicitTriangulation_.preprocessEdgeLinks();
        }
        return 0;
      }
      
      /// Pre-process the edge stars.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getEdgeStar()
      ///   - getEdgeStars()
      ///   - getEdgeStarNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getEdgeStar()
      /// \sa getEdgeStars()
      /// \sa getEdgeStarNumber()
      inline int preprocessEdgeStars(){
        if(!explicitTriangulation_.isEmpty()){
          return explicitTriangulation_.preprocessEdgeStars();
        }
        else if(!implicitTriangulation_.isEmpty()){
          return implicitTriangulation_.preprocessEdgeStars();
        }
        return 0;
      }
      
      /// Pre-process the edge triangles.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getEdgeTriangle()
      ///   - getEdgeTriangles()
      ///   - getEdgeTriangleNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getEdgeTriangle()
      /// \sa getEdgeTriangles()
      /// \sa getEdgeTriangleNumber()
      inline int preprocessEdgeTriangles(){
        if(!explicitTriangulation_.isEmpty()){
          return explicitTriangulation_.preprocessEdgeTriangles();
        }
        else if(!implicitTriangulation_.isEmpty()){
          return implicitTriangulation_.preprocessEdgeTriangles();
        }
        return 0;
      }
      
      /// Pre-process the triangles.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getNumberOfTriangles()
      ///   - getTriangles()
      ///   - getTriangleVertex()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getNumberOfTriangles()
      /// \sa getTriangles()
      /// \sa getTriangleVertex()
      inline int preprocessTriangles(){
        if(!explicitTriangulation_.isEmpty()){
          return explicitTriangulation_.preprocessTriangles();
        }
        else if(!implicitTriangulation_.isEmpty()){
          return implicitTriangulation_.preprocessTriangles();
        }
        return 0;
      }
      
      /// Pre-process the triangle edges.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getTriangleEdge()
      ///   - getTriangleEdges()
      ///   - getTriangleEdgeNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getTriangleEdge()
      /// \sa getTriangleEdges()
      /// \sa getTriangleEdgeNumber()
      inline int preprocessTriangleEdges(){
        if(!explicitTriangulation_.isEmpty()){
          return explicitTriangulation_.preprocessTriangleEdges();
        }
        else if(!implicitTriangulation_.isEmpty()){
          return implicitTriangulation_.preprocessTriangleEdges();
        }
        return 0;
      }
      
      /// Pre-process the triangle links.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getTriangleLink()
      ///   - getTriangleLinks()
      ///   - getTriangleLinkNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getTriangleLink()
      /// \sa getTriangleLinks()
      /// \sa getTriangleLinkNumber()
      inline int preprocessTriangleLinks(){
        if(!explicitTriangulation_.isEmpty()){
          return explicitTriangulation_.preprocessTriangleLinks();
        }
        else if(!implicitTriangulation_.isEmpty()){
          return implicitTriangulation_.preprocessTriangleLinks();
        }
        return 0;
      }
      
      /// Pre-process the triangle stars.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getTriangleStar()
      ///   - getTriangleStars()
      ///   - getTriangleStarNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getTriangleStar()
      /// \sa getTriangleStars()
      /// \sa getTriangleStarNumber()
      inline int preprocessTriangleStars(){
        if(!explicitTriangulation_.isEmpty()){
          return explicitTriangulation_.preprocessTriangleStars();
        }
        else if(!implicitTriangulation_.isEmpty()){
          return implicitTriangulation_.preprocessTriangleStars();
        }
        return 0;
      }
      
      /// Pre-process the vertex edges.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getVertexEdge()
      ///   - getVertexEdges()
      ///   - getVertexEdgeNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getVertexEdge()
      /// \sa getVertexEdges()
      /// \sa getVertexEdgeNumber()
      inline int preprocessVertexEdges(){
        if(!explicitTriangulation_.isEmpty()){
          return explicitTriangulation_.preprocessVertexEdges();
        }
        else if(!implicitTriangulation_.isEmpty()){
          return implicitTriangulation_.preprocessVertexEdges();
        }
        return 0;
      }
      
      /// Pre-process the vertex links.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getVertexLink()
      ///   - getVertexLinks()
      ///   - getVertexLinkNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getVertexLink()
      /// \sa getVertexLinks()
      /// \sa getVertexLinkNumber()
      inline int preprocessVertexLinks(){
        if(!explicitTriangulation_.isEmpty()){
          return explicitTriangulation_.preprocessVertexLinks();
        }
        else if(!implicitTriangulation_.isEmpty()){
          return implicitTriangulation_.preprocessVertexLinks();
        }
        return 0;
      }
      
      /// Pre-process the vertex neighbors.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getVertexNeighbor()
      ///   - getVertexNeighbors()
      ///   - getVertexNeighborNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getVertexNeighbor()
      /// \sa getVertexNeighbors()
      /// \sa getVertexNeighborNumber()
      inline int preprocessVertexNeighbors(){
        if(!explicitTriangulation_.isEmpty()){
          return explicitTriangulation_.preprocessVertexNeighbors();
        }
        else if(!implicitTriangulation_.isEmpty()){
          return implicitTriangulation_.preprocessVertexNeighbors();
        }
        return 0;
      }
      
      /// Pre-process the vertex stars.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getVertexStar()
      ///   - getVertexStars()
      ///   - getVertexStarNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getVertexStar()
      /// \sa getVertexStars()
      /// \sa getVertexStarNumber()
      inline int preprocessVertexStars(){
        if(!explicitTriangulation_.isEmpty()){
          return explicitTriangulation_.preprocessVertexStars();
        }
        else if(!implicitTriangulation_.isEmpty()){
          return implicitTriangulation_.preprocessVertexStars();
        }
        return 0;
      }
      
      /// Pre-process the vertex triangles.
      ///
      /// This function should ONLY be called as a pre-condition to the 
      /// following functions:
      ///   - getVertexTriangle()
      ///   - getVertexTriangles()
      ///   - getVertexTriangleNumber()
      ///
      /// \pre This function should be called prior to any traversal, in a 
      /// clearly distinct pre-processing step that involves no traversal at 
      /// all. An error will be returned otherwise.
      /// \note It is recommended to exclude this pre-processing function from
      /// any time performance measurement.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \sa getVertexTriangle()
      /// \sa getVertexTriangles()
      /// \sa getVertexTriangleNumber()
      inline int preprocessVertexTriangles(){
        if(!explicitTriangulation_.isEmpty()){
          return explicitTriangulation_.preprocessVertexTriangles();
        }
        else if(!implicitTriangulation_.isEmpty()){
          return implicitTriangulation_.preprocessVertexTriangles();
        }
        return 0;
      }
     
      /// Tune the debug level (default: 0)
      inline const int setDebugLevel(const int &debugLevel){
        explicitTriangulation_.setDebugLevel(debugLevel);
        implicitTriangulation_.setDebugLevel(debugLevel);
        debugLevel_ = debugLevel;
        return 0;
      }


      /// Set the input cells for the triangulation.
      ///
      /// Here the notion of cell refers to the simplicices of maximal 
      /// dimension (3D: tetrahedra, 2D: triangles, 1D: edges).
      ///
      /// \param cellNumber Number of input cells.
      /// \param cellArray Pointer to the input cells. This pointer should point
      /// to an array of long long int where cells are stored one after the 
      /// other. In particular, each cell starts by the number of vertices in 
      /// it, followed by the identifiers of its vertices. This corresponds to 
      /// the default cell array representation in VTK.
      /// \return Returns 0 upon success, negative values otherwise.
      ///
      /// \note This function does not need to be called if the current object 
      /// is a vtkTriangulation (this function is automatically called 
      /// if needed through vtkTriangulation::setInputData()).
      inline int setInputCells(const int &cellNumber,
        const long long int *cellArray){
        
        return explicitTriangulation_.setInputCells(cellNumber, cellArray);
      }
      
      /// Set the specifications of the input grid to implicitly represent as a
      /// triangulation.
      /// \param xOrigin Input x coordinate of the grid origin.
      /// \param yOrigin Input y coordinate of the grid origin.
      /// \param zOrigin Input z coordinate of the grid origin.
      /// \param xSpacing Input spacing along the x dimension.
      /// \param ySpacing Input spacing along the y dimension.
      /// \param zSpacing Input spacing along the z dimension.
      /// \param xDim Input number of vertices along the x dimension.
      /// \param yDim Input number of vertices along the y dimension.
      /// \param zDim Input number of vertices along the z dimension.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \note This function does not need to be called if the current object 
      /// is a vtkTriangulation (this function is automatically called 
      /// if needed through vtkTriangulation::setInputData()).
      inline int setInputGrid(
        const float &xOrigin, const float &yOrigin, const float &zOrigin,
        const float &xSpacing, const float &ySpacing, const float &zSpacing,
        const int &xDim, const int &yDim, const int &zDim){
      
        return implicitTriangulation_.setInputGrid(
          xOrigin, yOrigin, zOrigin,
          xSpacing, ySpacing, zSpacing,
          xDim, yDim, zDim);
        return 0;
      }
      
      /// Set the input 3D points of the triangulation.
      /// \param pointNumber Number of input vertices.
      /// \param pointSet Pointer to the 3D points. This pointer should point to
      /// an array of float where points are stored one after the other. 
      /// In particular, each point is represented by X-Y-Z coordinates (one 
      /// after the other). This corresponds to the default point set 
      /// representation in VTK.
      /// \return Returns 0 upon success, negative values otherwise.
      /// \note This function does not need to be called if the current object 
      /// is a vtkTriangulation (this function is automatically called 
      /// if needed through vtkTriangulation::setInputData()).
      inline int setInputPoints(const int &pointNumber, const float *pointSet){
        return explicitTriangulation_.setInputPoints(pointNumber, pointSet);
      }

      /// Tune the number of active threads (default: number of logical cores)
      inline int setThreadNumber(const int &threadNumber){
        explicitTriangulation_.setThreadNumber(threadNumber);
        implicitTriangulation_.setThreadNumber(threadNumber);
        threadNumber_ = threadNumber;
        return 0;
      }


      /// Internal usage. Pass the execution context (debug level, number of 
      /// threads, etc.) to the implementing classes.
      inline int setWrapper(const Wrapper *wrapper){
        explicitTriangulation_.setWrapper(wrapper);
        implicitTriangulation_.setWrapper(wrapper);
        return 0;
      }
      
    protected:
    
      ExplicitTriangulation 
                          explicitTriangulation_;
      ImplicitTriangulation
                          implicitTriangulation_;
  };
}

// if the package is not a template, comment the following line
// #include                  <Triangulation.cpp>

#endif // _TRIANGULATION_H
