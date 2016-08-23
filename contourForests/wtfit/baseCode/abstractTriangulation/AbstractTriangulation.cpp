#include                  <AbstractTriangulation.h>

AbstractTriangulation::AbstractTriangulation(){

  hasPreprocessedCellEdges_ = false;
  hasPreprocessedCellNeighbors_ = false;
  hasPreprocessedCellTriangles_ = false;
  hasPreprocessedEdges_ = false;
  hasPreprocessedEdgeLinks_ = false;
  hasPreprocessedEdgeStars_ = false;
  hasPreprocessedEdgeTriangles_ = false;
  hasPreprocessedTriangles_ = false;
  hasPreprocessedTriangleEdges_ = false;
  hasPreprocessedTriangleLinks_ = false;
  hasPreprocessedTriangleStars_ = false;
  hasPreprocessedVertexEdges_ = false;
  hasPreprocessedVertexLinks_ = false;
  hasPreprocessedVertexNeighbors_ = false;
  hasPreprocessedVertexStars_ = false;
  hasPreprocessedVertexTriangles_ = false;
}

AbstractTriangulation::~AbstractTriangulation(){

}
