#include <vtkContourForests.h>

vtkStandardNewMacro(vtkContourForests);

vtkContourForests::vtkContourForests()
    :  /// Base ///
      FieldId{},
      isChecked_{},
      isLoaded_{},
      calculSegmentation_{true},
      lessPartition_{false},
      tree_{},
      // Here the given number of core only serve for preprocess,
      // a clean tree append before the true process and re-set
      // the good number of threads
      contourTree_{new ParallelContourTree(OsCall::getNumberOfCores())},
      skeletonNodes_{vtkPolyData::New()},
      skeletonArcs_{vtkPolyData::New()},
      segmentation_{},
      CTPersistenceCurve_{vtkTable::New()},
      MTPersistenceCurve_{vtkTable::New()},
      STPersistenceCurve_{vtkTable::New()},
      simplificationThresholdCurve_{vtkTable::New()},
      CTPersistenceDiagram_{vtkUnstructuredGrid::New()},
      MTPersistenceDiagram_{vtkUnstructuredGrid::New()},
      STPersistenceDiagram_{vtkUnstructuredGrid::New()},

      /// Void ///
      voidUnstructuredGrid_{vtkUnstructuredGrid::New()},
      voidPolyData_{vtkPolyData::New()},
      voidTable_{vtkTable::New()},

      /// Configuration ///
      checkConformance_{},
      varyingMeshGeometry_{true},
      varyingMeshConnectivity_{true},
      varyingDataValues_{true},
      treeType_{TreeType::CONTOUR_TREE},
      showMin_{},
      showMax_{},
      showSaddle1_{},
      showSaddle2_{},
      showArc_{},
      arcResolution_{1},
      partitionNum_{-1},
      skeletonSmoothing_{},
      simplificationType_{},
      simplificationThreshold_{},
      simplificationThresholdBuffer_{},
      simplifyPersistenceDiagram_{},

      /// Computation handles ///
      toUpdateVertexSoSoffsets_{true},
      toComputeContourTree_{true},
      toComputeSimplification_{true},
      toUpdateTree_{true},
      toComputeSkeleton_{true},
      toComputeSegmentation_{true},
      toShowPersistenceCurve_{true},
      toShowPersistenceDiagram_{false},

      /// Convenient storage ///
      deltaScalar_{},
      originalNumberOfCriticalPoints_{},
      numberOfVertices_{},
      triangulation_{new vtkTriangulation},
      vertexPositions_{new vector<vector<double>>},
      vertexSoSoffsets_{new vector<int>},
      criticalPoints_{new vector<int>},
      vertexScalars_{nullptr},
      inputScalars_{new vector<vector<double>>},
      inputScalarsName_{new vector<string>},
      mergePairs_{new vector<pair<pair<int, int>, double>>},
      splitPairs_{new vector<pair<pair<int, int>, double>>},
      pairs_{new vector<pair<pair<int, int>, double>>},
      samples_{new vector<vector<vector<vector<int>>>>},
      barycenters_{new vector<vector<vector<vector<double>>>>}
{
   contourTree_->setWrapper(this);
   contourTree_->setDebugLevel(debugLevel_);
   UseAllCores = false;

   /// VTK Interface ///
   SetNumberOfInputPorts(1);
   SetNumberOfOutputPorts(3);

}

vtkContourForests::~vtkContourForests()
{
   /// Base ///
   delete contourTree_;
   skeletonNodes_->Delete();
   skeletonArcs_->Delete();
   // segmentation_->Delete();
   CTPersistenceCurve_->Delete();
   MTPersistenceCurve_->Delete();
   STPersistenceCurve_->Delete();
   simplificationThresholdCurve_->Delete();
   CTPersistenceDiagram_->Delete();
   MTPersistenceDiagram_->Delete();
   STPersistenceDiagram_->Delete();

   /// Void ///
   voidUnstructuredGrid_->Delete();
   voidPolyData_->Delete();
   voidTable_->Delete();

   /// Convenient storage ///
   delete triangulation_;
   delete vertexPositions_;
   delete vertexSoSoffsets_;
   delete criticalPoints_;
   delete inputScalars_;
   delete inputScalarsName_;
   delete mergePairs_;
   delete splitPairs_;
   delete pairs_;
   delete samples_;
   delete barycenters_;
}

void vtkContourForests::clearSkeleton()
{
   samples_->clear();
   barycenters_->clear();

   skeletonNodes_->Delete();
   skeletonNodes_ = vtkPolyData::New();
   skeletonArcs_->Delete();
   skeletonArcs_ = vtkPolyData::New();
}

void vtkContourForests::clearSegmentation()
{
   segmentation_->Delete();
}

void vtkContourForests::clearTree()
{
   tree_ = nullptr;

   int nbThread = threadNumber_;

   if(lessPartition_){
    if(nbThread%2){
        cout << "Less partition need a pair number of thread, increment " << endl;
        nbThread++;
    }

    nbThread /= 2;

   }

   delete contourTree_;
   contourTree_ = new ParallelContourTree(nbThread);
   contourTree_->setWrapper(this);
   contourTree_->setDebugLevel(debugLevel_);
   contourTree_->setLessPartition(lessPartition_);
   contourTree_->setThreadNumber(nbThread);

   mergePairs_->clear();
   splitPairs_->clear();
   pairs_->clear();
}

// transmit abort signals -- to copy paste in other wrappers
bool vtkContourForests::needsToAbort()
{
   return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int vtkContourForests::updateProgress(const float& progress)
{
   {
      stringstream msg;
      msg << "[vtkContourForests] " << progress * 100 << "% processed...." << endl;
      dMsg(cout, msg.str(), advancedInfoMsg);
   }

   UpdateProgress(progress);
   return 0;
}

int vtkContourForests::FillInputPortInformation(int port, vtkInformation* info)
{
   if (port == 0)
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
   return 1;
}

int vtkContourForests::FillOutputPortInformation(int port, vtkInformation* info)
{
   switch (port) {
      case 0:
         info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
         break;
      case 1:
         info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
         break;
      case 2:
         info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
         break;
      //case 3:
         //info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
         //break;
      //case 4:
         //info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
         //break;
      //case 5:
         //info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
         //break;
   }

   return 1;
}

void vtkContourForests::SetThreads()
{
   if (!UseAllCores)
      threadNumber_ = ThreadNumber;
   else
      threadNumber_ = OsCall::getNumberOfCores();

   Modified();
}

void vtkContourForests::SetThreadNumber(int threadNumber)
{
   ThreadNumber = threadNumber;
   SetThreads();
}

void vtkContourForests::SetDebugLevel(int d)
{
   Debug::setDebugLevel(d);
   debugLevel_ = d;
}

void vtkContourForests::SetUseAllCores(bool onOff)
{
   UseAllCores = onOff;
   SetThreads();
}

void vtkContourForests::SetUseInputOffsetScalarField(bool onOff)
{
   toUpdateVertexSoSoffsets_ = true;
   toComputeContourTree_     = true;
   toUpdateTree_             = true;
   toComputeSkeleton_        = true;
   toComputeSegmentation_    = true;

   useInputOffsetScalarField_ = onOff;
   Modified();
}

void vtkContourForests::CheckConformance(int checkConformance)
{
   checkConformance_ = checkConformance;
   Modified();
}

void vtkContourForests::SetScalarField(string scalarField)
{
   toComputeContourTree_  = true;
   toComputeSkeleton_     = true;
   toComputeSegmentation_ = true;

   scalarField_ = scalarField;
   Modified();
}

void vtkContourForests::EnableVaryingMeshGeometry(bool state)
{
   varyingMeshGeometry_ = state;
   Modified();
}

void vtkContourForests::EnableVaryingMeshConnectivity(bool state)
{
   varyingMeshConnectivity_ = state;
   Modified();
}

void vtkContourForests::EnableVaryingDataValues(bool state)
{
   varyingDataValues_ = state;
   Modified();
}

void vtkContourForests::SetTreeType(int treeType)
{
   if (treeType >= 0 && treeType <= 2) {
      toUpdateTree_ = true;
      if (treeType_ == TreeType::CONTOUR_TREE ||
          static_cast<TreeType>(treeType) == TreeType::CONTOUR_TREE) {
         toComputeContourTree_ = true;
      }

      toComputeSkeleton_     = true;
      toComputeSegmentation_ = true;

      treeType_ = static_cast<TreeType>(treeType);

      Modified();
   }
}

void vtkContourForests::ShowMin(bool state)
{
   toComputeSkeleton_ = true;

   showMin_ = state;
   Modified();
}

void vtkContourForests::ShowMax(bool state)
{
   toComputeSkeleton_ = true;

   showMax_ = state;
   Modified();
}

void vtkContourForests::ShowSaddle1(bool state)
{
   toComputeSkeleton_ = true;

   showSaddle1_ = state;
   Modified();
}

void vtkContourForests::ShowSaddle2(bool state)
{
   toComputeSkeleton_ = true;

   showSaddle2_ = state;
   Modified();
}

void vtkContourForests::ShowArc(bool state)
{
   toComputeSkeleton_ = true;

   showArc_ = state;
   Modified();
}

void vtkContourForests::ShowSegmentation(bool segm)
{
   calculSegmentation_ = segm;
   if (segm) {
      toComputeSkeleton_ = true;
      Modified();
   }
}

void vtkContourForests::SetArcResolution(int arcResolution)
{
   if (arcResolution >= 0) {
      toComputeSkeleton_ = true;

      arcResolution_ = arcResolution;
      Modified();
   }
}

void vtkContourForests::SetPartitionNumber(int partitionNum)
{
   partitionNum_ = partitionNum;
   Modified();
}

void vtkContourForests::SetLessPartition(bool l)
{
    lessPartition_ = l;
    Modified();
}

void vtkContourForests::SetSkeletonSmoothing(double skeletonSmoothing)
{
   if (skeletonSmoothing >= 0) {
      toComputeSkeleton_ = true;

      skeletonSmoothing_ = skeletonSmoothing;
      Modified();
   }
}

void vtkContourForests::ShowPersistenceCurve(bool state)
{
   toShowPersistenceCurve_ = state;
   Modified();
}

void vtkContourForests::ShowPersistenceDiagram(bool state)
{
   toShowPersistenceDiagram_ = state;
   Modified();
}

void vtkContourForests::SetSimplificationType(int type)
{
    simplificationType_ = type;
    Modified();
}

void vtkContourForests::SetSimplificationThreshold(double simplificationThreshold)
{
   if (simplificationThreshold >= 0 && simplificationThreshold <= 1) {
      simplificationThresholdBuffer_ = simplificationThreshold;
      toComputeSimplification_ = true;
      toUpdateTree_            = true;
      toComputeSkeleton_       = true;
      toComputeSegmentation_   = true;

      Modified();
   }
}

void vtkContourForests::SetSimplifyPersistenceDiagram(bool state)
{
   simplifyPersistenceDiagram_ = state;
   Modified();
}

int vtkContourForests::vtkDataSetToStdVector(vtkDataSet* input)
{
   // init
   if (varyingMeshGeometry_ || !numberOfVertices_)
      numberOfVertices_ = input->GetNumberOfPoints();

   // scalars
   if (varyingDataValues_ || !inputScalarsName_->size()) {
      if (input->GetPointData()) {
         int numberOfArrays = input->GetPointData()->GetNumberOfArrays();
         int numberOfScalarArrays{};

         for (int i = 0; i < numberOfArrays; ++i) {
            vtkDataArray* inputArray = input->GetPointData()->GetArray(i);
            if (inputArray) {
               if (inputArray->GetNumberOfTuples() == numberOfVertices_ &&
                   inputArray->GetNumberOfComponents() == 1) {
                  ++numberOfScalarArrays;
               }
            }
         }

         inputScalars_->resize(numberOfScalarArrays);
         inputScalarsName_->resize(numberOfScalarArrays);

         int k{};
         for (int i = 0; i < numberOfArrays; ++i) {
            vtkDataArray* inputArray = input->GetPointData()->GetArray(i);
            if (inputArray) {
               if (inputArray->GetNumberOfTuples() == numberOfVertices_ &&
                   inputArray->GetNumberOfComponents() == 1) {
                  (*inputScalars_)[k].resize(numberOfVertices_);
                  (*inputScalarsName_)[k] = inputArray->GetName();

                  for (unsigned int j = 0; j < numberOfVertices_; ++j) {
                     (*inputScalars_)[k][j] = inputArray->GetTuple1(j);
                  }

                  ++k;
               }
            }
         }
      }
   }

   if (scalarField_.size() == 0) {
      vertexScalars_ = &((*inputScalars_)[FieldId]);
      scalarField_ = (*inputScalarsName_)[FieldId];

   } else {

      for (unsigned int i = 0; i < inputScalarsName_->size(); ++i) {
         if ((*inputScalarsName_)[i] == scalarField_)
            vertexScalars_ = &((*inputScalars_)[i]);
      }
   }

   auto result      = std::minmax_element(vertexScalars_->begin(), vertexScalars_->end());
   double scalarMin = *result.first;
   double scalarMax = *result.second;
   deltaScalar_           = (scalarMax - scalarMin);

   // neighbors
   if (varyingMeshConnectivity_ || triangulation_->isEmpty()) {
      if (triangulation_)
         delete triangulation_;

      triangulation_ = new vtkTriangulation;
      triangulation_->setWrapper(this);
      //triangulation_->setDebugLevel(debugLevel_);
      triangulation_->setDebugLevel(0);
      triangulation_->setThreadNumber(threadNumber_);
      triangulation_->setInputData(input);
      //triangulation_->initData();
      triangulation_->preprocessEdges();
      triangulation_->preprocessVertexNeighbors();
   }

   // positions
   if (varyingMeshGeometry_ || !vertexPositions_->size()) {
      vertexPositions_->resize(numberOfVertices_);

      double point[3];
      for (unsigned int i = 0; i < vertexPositions_->size(); ++i) {
         (*vertexPositions_)[i].resize(3);
         input->GetPoint(i, point);
         for (unsigned int j = 0; j < 3; ++j)
            (*vertexPositions_)[i][j] = point[j];
      }
   }

   // offsets
   if (varyingMeshConnectivity_ || toUpdateVertexSoSoffsets_ || !vertexSoSoffsets_->size()) {
      vertexSoSoffsets_->resize(numberOfVertices_);

      if (useInputOffsetScalarField_ and inputOffsetScalarFieldName_.length()) {
         auto offsets = input->GetPointData()->GetArray(inputOffsetScalarFieldName_.data());
         for (unsigned int i = 0; i < vertexSoSoffsets_->size(); ++i)
            (*vertexSoSoffsets_)[i] = offsets->GetTuple1(i);
      } else {
         for (unsigned int i = 0; i < vertexSoSoffsets_->size(); ++i)
            (*vertexSoSoffsets_)[i] = i;
      }
      toUpdateVertexSoSoffsets_ = false;
   }

   if (varyingMeshConnectivity_ || varyingDataValues_ || !isLoaded_) {
      stringstream msg;
      msg << "[vtkContourForests] Convenient data storage has been loaded." << endl;
      msg << "[vtkContourForests]   Number of input scalars: " << inputScalars_->size() << endl;
      msg << "[vtkContourForests]   Input scalars name:" << endl;
      for (unsigned int i = 0; i < inputScalarsName_->size(); ++i)
         msg << "[vtkContourForests]     " << (*inputScalarsName_)[i] << endl;
      msg << "[vtkContourForests]   Active scalar name: " << scalarField_ << endl;
      msg << "[vtkContourForests]   Number of tuples: " << vertexScalars_->size() << endl;
      msg << "[vtkContourForests]   [min max]: [" << scalarMin << " " << scalarMax << "]" << endl;
      msg << "[vtkContourForests]   Number of vertices: " << numberOfVertices_ << endl;
      msg << "[vtkContourForests]   Vertex positions: " << boolalpha
          << (bool)vertexPositions_->size() << endl;
      msg << "[vtkContourForests]   Vertex offsets: " << boolalpha << (bool)vertexSoSoffsets_->size()
          << endl;
      dMsg(cout, msg.str(), detailedInfoMsg);
   }

  stringstream msg;
  msg << "[vtkContourForests] Launching computation for field '"
    << scalarField_ << "'..." << endl;
  dMsg(cout, msg.str(), timeMsg);

   isLoaded_ = true;
   return 0;
}

bool vtkContourForests::isCoincident(double p1[], double p2[])
{
   double sPrev[3];
   double sNext[3];
   for (unsigned int k = 0; k < 3; k++) {
      sPrev[k] = p2[k] - p1[k];
      sNext[k] = sPrev[k];
   }

   return (vtkMath::Normalize(sNext) == 0.0);
}

void vtkContourForests::getSkeletonArcs()
{
   vtkSmartPointer<vtkAppendPolyData> app = vtkSmartPointer<vtkAppendPolyData>::New();

   vtkDoubleArray* scalars{};
   vtkIntArray*    identifierScalars{};
   vtkIntArray*    typeScalars{};
   vtkIntArray*    sizeScalars{};
   vtkDoubleArray* spanScalars{};
   int type = static_cast<int>(TreeComponent::ARC);

   double point1[3];
   vector<double> point2(3);
   // get skeleton scalars
   vector<vector<vector<double>>> skeletonScalars(inputScalars_->size());
   for (unsigned int f = 0; f < inputScalars_->size(); ++f)
      getSkeletonScalars((*inputScalars_)[f], skeletonScalars[f]);

   double    inputScalar;
   SuperArc* a;
   int       regionSize;
   double    regionSpan;
   int       currentZone = 0;
   int       regionId;

   for (unsigned int i = 0; i < tree_->getNumberOfSuperArcs(); ++i) {
      a = tree_->getSuperArc(i);

      if (a->isVisible()) {
         int      upNodeId   = tree_->getSuperArc(i)->getUpNodeId();
         int      upVertex   = tree_->getNode(upNodeId)->getVertexId();
         float    coordUp[3];
         triangulation_->getVertexPoint(upVertex, coordUp[0], coordUp[1], coordUp[2]);

         int      downNodeId   = tree_->getSuperArc(i)->getDownNodeId();
         int      downVertex   = tree_->getNode(downNodeId)->getVertexId();
         float    coordDown[3];
         triangulation_->getVertexPoint(downVertex, coordDown[0], coordDown[1], coordDown[2]);

         regionSize = tree_->getNumberOfVisibleRegularNode(i);
         regionSpan = Geometry::distance(coordUp, coordDown, 3);
         regionId       = currentZone++;

         /// Line ///
         if ((*barycenters_)[static_cast<int>(treeType_)][i].size()) {
            // init: min
            int downNodeVId;
            if (treeType_ == TreeType::SPLIT_TREE)
               downNodeVId = tree_->getNode(a->getUpNodeId())->getVertexId();
            else
               downNodeVId = tree_->getNode(a->getDownNodeId())->getVertexId();

            for (unsigned int k = 0; k < 3; ++k) {
               point1[k] = (*vertexPositions_)[downNodeVId][k];
            }
            vtkSmartPointer<vtkLineSource> line = vtkSmartPointer<vtkLineSource>::New();
            line->SetPoint1(point1);

            const auto nbBarycenter = (*barycenters_)[static_cast<int>(treeType_)][i].size();
            for (unsigned int j = 0; j < nbBarycenter; ++j) {
               point2 = (*barycenters_)[static_cast<int>(treeType_)][i][j];
               line->SetPoint2(point2.data());

               if (!isCoincident(point1, point2.data())) {
                  line->Update();
                  vtkPolyData* lineData = line->GetOutput();

                  /// Point data ///
                  for (unsigned int f = 0; f < inputScalars_->size(); ++f) {
                     inputScalar = skeletonScalars[f][i][j];

                     scalars = vtkDoubleArray::New();
                     scalars->SetName((*inputScalarsName_)[f].data());
                     for (unsigned int k = 0; k < 2; ++k)
                        scalars->InsertTuple1(k, inputScalar);
                     lineData->GetPointData()->AddArray(scalars);
                     scalars->Delete();
                  }

                  /// Cell data ///
                  // Identifier
                  identifierScalars = vtkIntArray::New();
                  identifierScalars->SetName("SegmentationId");
                  for (unsigned int k = 0; k < 2; ++k)
                     identifierScalars->InsertTuple1(k, regionId);
                  lineData->GetCellData()->AddArray(identifierScalars);
                  identifierScalars->Delete();
                  // Type
                  typeScalars = vtkIntArray::New();
                  typeScalars->SetName("Type");
                  for (unsigned int k = 0; k < 2; ++k)
                     typeScalars->InsertTuple1(k, type);
                  lineData->GetCellData()->AddArray(typeScalars);
                  typeScalars->Delete();
                  // Size
                  sizeScalars = vtkIntArray::New();
                  sizeScalars->SetName("RegionSize");
                  for (unsigned int k = 0; k < 2; ++k)
                     sizeScalars->InsertTuple1(k, regionSize);
                  lineData->GetCellData()->AddArray(sizeScalars);
                  sizeScalars->Delete();
                  // Span
                  spanScalars = vtkDoubleArray::New();
                  spanScalars->SetName("RegionSpan");
                  for (unsigned int k = 0; k < 2; ++k)
                     spanScalars->InsertTuple1(k, regionSpan);
                  lineData->GetCellData()->AddArray(spanScalars);
                  spanScalars->Delete();

                  app->AddInputData(lineData);
               }

               line = vtkSmartPointer<vtkLineSource>::New();
               line->SetPoint1(point2.data());
            }

            // end: max
            int upNodeVId;
            if (treeType_ == TreeType::SPLIT_TREE)
               upNodeVId = tree_->getNode(a->getDownNodeId())->getVertexId();
            else
               upNodeVId = tree_->getNode(a->getUpNodeId())->getVertexId();
            for (unsigned int k = 0; k < 3; ++k)
               point2[k] = (*vertexPositions_)[upNodeVId][k];
            line->SetPoint2(point2.data());

            if (!isCoincident(point1, point2.data())) {
               line->Update();
               vtkPolyData* lineData = line->GetOutput();

               /// Point data ///
               for (unsigned int f = 0; f < inputScalars_->size(); ++f) {
                  inputScalar =
                      skeletonScalars[f][i][(*barycenters_)[static_cast<int>(treeType_)][i].size()];

                  scalars = vtkDoubleArray::New();
                  scalars->SetName((*inputScalarsName_)[f].data());
                  for (unsigned int k = 0; k < 2; ++k)
                     scalars->InsertTuple1(k, inputScalar);
                  lineData->GetPointData()->AddArray(scalars);
                  scalars->Delete();
               }

               /// Cell data ///
               // Identifier
               identifierScalars = vtkIntArray::New();
               identifierScalars->SetName("SegmentationId");
               for (unsigned int k = 0; k < 2; ++k)
                  identifierScalars->InsertTuple1(k, regionId);
               lineData->GetCellData()->AddArray(identifierScalars);
               identifierScalars->Delete();
               // Type
               typeScalars = vtkIntArray::New();
               typeScalars->SetName("Type");
               for (unsigned int k = 0; k < 2; ++k)
                  typeScalars->InsertTuple1(k, type);
               lineData->GetCellData()->AddArray(typeScalars);
               typeScalars->Delete();
               // Size
               sizeScalars = vtkIntArray::New();
               sizeScalars->SetName("RegionSize");
               for (unsigned int k = 0; k < 2; ++k)
                  sizeScalars->InsertTuple1(k, regionSize);
               lineData->GetCellData()->AddArray(sizeScalars);
               sizeScalars->Delete();
               // Span
               spanScalars = vtkDoubleArray::New();
               spanScalars->SetName("RegionSpan");
               for (unsigned int k = 0; k < 2; ++k)
                  spanScalars->InsertTuple1(k, regionSpan);
               lineData->GetCellData()->AddArray(spanScalars);
               spanScalars->Delete();

               app->AddInputData(lineData);
            }
         } else {
            vtkSmartPointer<vtkLineSource> line = vtkSmartPointer<vtkLineSource>::New();

            int downNodeVId = tree_->getNode(a->getDownNodeId())->getVertexId();
            for (unsigned int k = 0; k < 3; ++k)
               point1[k] = (*vertexPositions_)[downNodeVId][k];
            line->SetPoint1(point1);

            int upNodeVId = tree_->getNode(a->getUpNodeId())->getVertexId();
            for (unsigned int k = 0; k < 3; ++k)
               point2[k] = (*vertexPositions_)[upNodeVId][k];
            line->SetPoint2(point2.data());

            if (!isCoincident(point1, point2.data())) {
               line->Update();
               vtkPolyData* lineData = line->GetOutput();

               /// Point data ///
               for (unsigned int f = 0; f < inputScalars_->size(); ++f) {
                  inputScalar = skeletonScalars[f][i][0];

                  scalars = vtkDoubleArray::New();
                  scalars->SetName((*inputScalarsName_)[f].data());
                  for (unsigned int k = 0; k < 2; ++k)
                     scalars->InsertTuple1(k, inputScalar);
                  lineData->GetPointData()->AddArray(scalars);
                  scalars->Delete();
               }

               /// Cell data ///
               // Identifier
               identifierScalars = vtkIntArray::New();
               identifierScalars->SetName("SegmentationId");
               for (int k = 0; k < 2; ++k)
                  identifierScalars->InsertTuple1(k, regionId);
               lineData->GetCellData()->AddArray(identifierScalars);
               identifierScalars->Delete();
               // Type
               typeScalars = vtkIntArray::New();
               typeScalars->SetName("Type");
               for (unsigned int k = 0; k < 2; ++k)
                  typeScalars->InsertTuple1(k, type);
               lineData->GetCellData()->AddArray(typeScalars);
               typeScalars->Delete();
               // Size
               sizeScalars = vtkIntArray::New();
               sizeScalars->SetName("RegionSize");
               for (unsigned int k = 0; k < 2; ++k)
                  sizeScalars->InsertTuple1(k, regionSize);
               lineData->GetCellData()->AddArray(sizeScalars);
               sizeScalars->Delete();
               // Span
               spanScalars = vtkDoubleArray::New();
               spanScalars->SetName("RegionSpan");
               for (unsigned int k = 0; k < 2; ++k)
                  spanScalars->InsertTuple1(k, regionSpan);
               lineData->GetCellData()->AddArray(spanScalars);
               spanScalars->Delete();

               app->AddInputData(lineData);
            }
         }
      } else {
         //cout << " pruned _ :" << tree_->getNode(a->getDownNodeId())->getVertexId() << " -  "
              //<< tree_->getNode(a->getUpNodeId())->getVertexId() << endl;
      }
   }

   app->Update();
   skeletonArcs_->ShallowCopy(app->GetOutput());
}

int vtkContourForests::getSkeletonScalars(const vector<double>& scalars,
                                         vector<vector<double>>& skeletonScalars) const
{
   skeletonScalars.clear();
   skeletonScalars.resize(tree_->getNumberOfSuperArcs());

   int nodeId;
   int vertexId;

   double f;
   double f0;
   double f1;
   double fmin;
   double fmax;
   int nodeMinId;
   int nodeMaxId;
   int nodeMinVId;
   int nodeMaxVId;
   const SuperArc* a;
   for (unsigned int i = 0; i < tree_->getNumberOfSuperArcs(); ++i) {
      a = tree_->getSuperArc(i);

      if (!a->isPruned()) {
         if (treeType_ == TreeType::SPLIT_TREE) {
            nodeMinId = a->getUpNodeId();
            nodeMaxId = a->getDownNodeId();
         } else {
            nodeMaxId = a->getUpNodeId();
            nodeMinId = a->getDownNodeId();
         }

         nodeMaxVId = tree_->getNode(nodeMaxId)->getVertexId();
         nodeMinVId = tree_->getNode(nodeMinId)->getVertexId();

         fmax = scalars[nodeMaxVId];
         fmin = scalars[nodeMinVId];

         // init: min
         f0 = fmin;

         // iteration
         for (unsigned int j = 0; j < (*samples_)[static_cast<int>(treeType_)][i].size(); ++j) {
            const vector<int>& sample = (*samples_)[static_cast<int>(treeType_)][i][j];

            f = 0;
            for (unsigned int k = 0; k < sample.size(); ++k) {
               nodeId   = sample[k];
               vertexId = nodeId;
               f += scalars[vertexId];
            }
            if (sample.size()) {
               f /= sample.size();

               f1 = f;
               // update the arc
               skeletonScalars[i].push_back((f0 + f1) / 2);
               f0 = f1;
            }
         }

         // end: max
         f1 = fmax;

         // update the arc
         skeletonScalars[i].push_back((f0 + f1) / 2);
      }
   }

   return 0;
}

void vtkContourForests::getSkeletonNodes()
{
   vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
   double point[3];

   double scalar{};
   vector<vtkDoubleArray*> scalars(inputScalars_->size());
   for (unsigned int f = 0; f < inputScalars_->size(); ++f) {
      scalars[f] = vtkDoubleArray::New();
      scalars[f]->SetName((*inputScalarsName_)[f].data());
   }

   vtkIntArray* nodeIdentifierScalars = vtkIntArray::New();
   nodeIdentifierScalars->SetName("NodeIdentifier");

   vtkIntArray* vertexIdentifierScalars = vtkIntArray::New();
   vertexIdentifierScalars->SetName("VertexIdentifier");

   int type{};
   vtkIntArray* nodeTypeScalars = vtkIntArray::New();
   nodeTypeScalars->SetName("NodeType");

   vector<double> persistence(numberOfVertices_, -1);
   if (treeType_ == TreeType::MERGE_TREE) {
      for (unsigned int i = 0; i < mergePairs_->size(); ++i) {
         persistence[(*mergePairs_)[i].first.first] =
             std::max(persistence[(*mergePairs_)[i].first.first], (*mergePairs_)[i].second);
         persistence[(*mergePairs_)[i].first.second] =
             std::max(persistence[(*mergePairs_)[i].first.second], (*mergePairs_)[i].second);
      }
   } else if (treeType_ == TreeType::SPLIT_TREE) {
      for (unsigned int i = 0; i < splitPairs_->size(); ++i) {
         persistence[(*splitPairs_)[i].first.first] =
             std::max(persistence[(*splitPairs_)[i].first.first], (*splitPairs_)[i].second);
         persistence[(*splitPairs_)[i].first.second] =
             std::max(persistence[(*splitPairs_)[i].first.second], (*splitPairs_)[i].second);
      }
   } else {
      for (unsigned int i = 0; i < pairs_->size(); ++i) {
         persistence[(*pairs_)[i].first.first] =
             std::max(persistence[(*pairs_)[i].first.first], (*pairs_)[i].second);
         persistence[(*pairs_)[i].first.second] =
             std::max(persistence[(*pairs_)[i].first.second], (*pairs_)[i].second);
      }
   }
   vtkDoubleArray* nodePersistenceScalars = vtkDoubleArray::New();
   nodePersistenceScalars->SetName("NodePersistence");

   int identifier{};
   for (unsigned i = 0; i < criticalPoints_->size(); ++i) {
      int nodeId        = (*criticalPoints_)[i];
      if(tree_->getNode(nodeId)->isHidden()) continue;
      int vertexId      = tree_->getNode(nodeId)->getVertexId();
      NodeType nodeType = getNodeType(nodeId);

      if ((nodeType == NodeType::LOCAL_MINIMUM and showMin_) or
          (nodeType == NodeType::LOCAL_MAXIMUM and showMax_) or
          (nodeType == NodeType::SADDLE1 and showSaddle1_) or
          (nodeType == NodeType::SADDLE2 and showSaddle2_)) {
         // Positions
         for (unsigned k = 0; k < 3; ++k)
            point[k] = (*vertexPositions_)[vertexId][k];
         points->InsertPoint(identifier, point);

         // Scalars
         for (unsigned int f = 0; f < inputScalars_->size(); ++f) {
            scalar = (*inputScalars_)[f][vertexId];
            scalars[f]->InsertTuple1(identifier, scalar);
         }

         // NodeIdentifier
         nodeIdentifierScalars->InsertTuple1(identifier, nodeId);

         // VertexIdentifier
         vertexIdentifierScalars->InsertTuple1(identifier, vertexId);

         // Type
         type = static_cast<int>(nodeType);
         nodeTypeScalars->InsertTuple1(identifier, type);

         // Persistence
         nodePersistenceScalars->InsertTuple1(identifier, persistence[vertexId]);

         ++identifier;
      }
   }
   skeletonNodes_->SetPoints(points);
   for (unsigned int f = 0; f < inputScalars_->size(); ++f)
      skeletonNodes_->GetPointData()->AddArray(scalars[f]);
   skeletonNodes_->GetPointData()->AddArray(nodeIdentifierScalars);
   skeletonNodes_->GetPointData()->AddArray(vertexIdentifierScalars);
   skeletonNodes_->GetPointData()->AddArray(nodeTypeScalars);
   skeletonNodes_->GetPointData()->AddArray(nodePersistenceScalars);

   for (unsigned int f = 0; f < inputScalars_->size(); ++f)
      scalars[f]->Delete();
   nodeIdentifierScalars->Delete();
   vertexIdentifierScalars->Delete();
   nodeTypeScalars->Delete();
   nodePersistenceScalars->Delete();
}

NodeType vtkContourForests::getNodeType(int id)
{
   return getNodeType(id, treeType_, tree_);
}

NodeType vtkContourForests::getNodeType(int id, TreeType type, MergeTree* tree)
{
   int upDegree{};
   int downDegree{};
   if (type == TreeType::MERGE_TREE || type == TreeType::CONTOUR_TREE) {
      upDegree   = tree->getNode(id)->getUpValence();
      downDegree = tree->getNode(id)->getDownValence();
   } else {
      downDegree   = tree->getNode(id)->getUpValence();
      upDegree = tree->getNode(id)->getDownValence();
   }
   int degree = upDegree + downDegree;

   // saddle point
   if (degree > 1) {
      if (upDegree > 1)
         return NodeType::SADDLE2;
      else if(downDegree >1)
         return NodeType::SADDLE1;
      else
         return NodeType::REGULAR;
   }
   // local extremum
   else {
      if (upDegree)
         return NodeType::LOCAL_MINIMUM;
      else
         return NodeType::LOCAL_MAXIMUM;
   }
}

void vtkContourForests::getCriticalPoints()
{
   vector<bool> isCriticalPoint(numberOfVertices_);

   criticalPoints_->clear();
   for (unsigned int i = 0; i < numberOfVertices_; ++i)
      isCriticalPoint[i] = false;

   //const int nbVert = triangulation_->getNumberOfOriginalVertices();

   // looking for critical points
   for (unsigned int i = 0; i < tree_->getNumberOfSuperArcs(); ++i) {
      auto a = tree_->getSuperArc(i);

      if (!a->isPruned()) {
         int upId   = a->getUpNodeId();
         int up_vId = tree_->getNode(upId)->getVertexId();
         if (!isCriticalPoint[up_vId]) {
            isCriticalPoint[up_vId] = true;
            criticalPoints_->push_back(upId);
         }

         int downId   = a->getDownNodeId();
         int down_vId = tree_->getNode(downId)->getVertexId();
         if (!isCriticalPoint[down_vId]) {
            isCriticalPoint[down_vId] = true;
            criticalPoints_->push_back(downId);
         }
      }
   }
   //{
      //stringstream msg;
      //msg << "[vtkContourForests] List of critical points :" << endl;
      //for (unsigned int it = 0; it < criticalPoints_->size(); ++it)
         //msg << "[vtkContourForests]   NodeId:" << (*criticalPoints_)[it]
             //<< ", VertexId:" << tree_->getNode(it)->getVertexId() << endl;
      //dMsg(cout, msg.str(), advancedInfoMsg);
   //}
}

void vtkContourForests::getSimplificationThresholdCurve()
{
   vtkSmartPointer<vtkDoubleArray> thresholdXScalars = vtkSmartPointer<vtkDoubleArray>::New();
   thresholdXScalars->SetName("Threshold X");

   vtkSmartPointer<vtkDoubleArray> thresholdYScalars = vtkSmartPointer<vtkDoubleArray>::New();
   thresholdYScalars->SetName("Threshold Y");

   vtkSmartPointer<vtkTable> thresholdCurve = vtkSmartPointer<vtkTable>::New();

   thresholdCurve->SetNumberOfRows(2);

   if (originalNumberOfCriticalPoints_) {
      thresholdXScalars->SetNumberOfTuples(2);
      thresholdYScalars->SetNumberOfTuples(2);

      thresholdXScalars->SetTuple1(0, simplificationThreshold_);
      thresholdYScalars->SetTuple1(0, 1);
      thresholdXScalars->SetTuple1(1, simplificationThreshold_);
      thresholdYScalars->SetTuple1(1, originalNumberOfCriticalPoints_ / 2 + 1);

      thresholdCurve->AddColumn(thresholdXScalars);
      thresholdCurve->AddColumn(thresholdYScalars);

      simplificationThresholdCurve_->ShallowCopy(thresholdCurve);
   }
}

void vtkContourForests::getPersistenceCurve(TreeType type)
{
   vtkSmartPointer<vtkDoubleArray> persistenceScalars = vtkSmartPointer<vtkDoubleArray>::New();
   persistenceScalars->SetName("Persistence");
   vtkSmartPointer<vtkIntArray> numberOfPairsScalars = vtkSmartPointer<vtkIntArray>::New();
   numberOfPairsScalars->SetName("Number Of Pairs");

   vtkSmartPointer<vtkTable> persistenceCurve = vtkSmartPointer<vtkTable>::New();

   vector<pair<double, int>> plot;
   switch (type) {
      case TreeType::MERGE_TREE:
         contourTree_->getJoinTree()->computePersistencePlot<double>(plot, mergePairs_);
         break;
      case TreeType::SPLIT_TREE:
         contourTree_->getSplitTree()->computePersistencePlot<double>(plot, splitPairs_);
         break;
      case TreeType::CONTOUR_TREE:
      default:
         contourTree_->computePersistencePlot<double>(plot, mergePairs_, splitPairs_, pairs_);
         break;
   }

   unsigned int N = plot.size();
   persistenceCurve->SetNumberOfRows(N);

   if (N) {
      persistenceScalars->SetNumberOfTuples(N);
      numberOfPairsScalars->SetNumberOfTuples(N);

      for (unsigned int i = 0; i < N; ++i) {
         persistenceScalars->SetTuple1(i, plot[i].first);
         numberOfPairsScalars->SetTuple1(i, plot[i].second);
      }

      persistenceCurve->AddColumn(persistenceScalars);
      persistenceCurve->AddColumn(numberOfPairsScalars);

      switch (type) {
         case TreeType::MERGE_TREE:
            MTPersistenceCurve_->ShallowCopy(persistenceCurve);
            break;
         case TreeType::SPLIT_TREE:
            STPersistenceCurve_->ShallowCopy(persistenceCurve);
            break;
         case TreeType::CONTOUR_TREE:
            CTPersistenceCurve_->ShallowCopy(persistenceCurve);
            break;
      }
   }
}

void vtkContourForests::getCurves()
{
   if (treeType_ == TreeType::CONTOUR_TREE) {
      getPersistenceCurve(TreeType::CONTOUR_TREE);
   } else {
      getPersistenceCurve(TreeType::MERGE_TREE);
      getPersistenceCurve(TreeType::SPLIT_TREE);
   }

   {
      stringstream msg;
      msg << "[vtkContourForests] Persistence curves built." << endl;
      dMsg(cout, msg.str(), infoMsg);
   }
}

void vtkContourForests::getPersistenceDiagramInfo(TreeType type, int pair, bool extremum,
                                                 int& vertexId, int& nodeId, int& nodeType)
{
   _persistenceCmp2 cmp(vertexScalars_);
   auto mergePairs = *mergePairs_;
   std::sort(mergePairs.begin(), mergePairs.end(), cmp);
   auto splitPairs = *splitPairs_;
   std::sort(splitPairs.begin(), splitPairs.end(), cmp);
   auto pairs = *pairs_;
   std::sort(pairs.begin(), pairs.end(), cmp);

   if (extremum) {
      switch (type) {
         case TreeType::MERGE_TREE:
            vertexId = mergePairs[pair].first.second;
            nodeId   = contourTree_->getJoinTree()->getCorrespondingNode(vertexId);
            nodeType = static_cast<int>(getNodeType(nodeId, type, contourTree_->getJoinTree()));
            break;
         case TreeType::SPLIT_TREE:
            vertexId = splitPairs[pair].first.first;
            nodeId   = contourTree_->getSplitTree()->getCorrespondingNode(vertexId);
            nodeType = static_cast<int>(getNodeType(nodeId, type, contourTree_->getSplitTree()));
            break;
         case TreeType::CONTOUR_TREE:
         default:
            vertexId = pairs[pair].first.second;
            nodeId   = contourTree_->getCorrespondingNode(vertexId);
            nodeType = static_cast<int>(getNodeType(nodeId, type, contourTree_));
            break;
      }
   } else {
      switch (type) {
         case TreeType::MERGE_TREE:
            vertexId = mergePairs[pair].first.first;
            nodeId   = contourTree_->getJoinTree()->getCorrespondingNode(vertexId);
            nodeType = static_cast<int>(getNodeType(nodeId, type, contourTree_->getJoinTree()));
            break;
         case TreeType::SPLIT_TREE:
            vertexId = splitPairs[pair].first.second;
            nodeId   = contourTree_->getSplitTree()->getCorrespondingNode(vertexId);
            nodeType = static_cast<int>(getNodeType(nodeId, type, contourTree_->getSplitTree()));
            break;
         case TreeType::CONTOUR_TREE:
         default:
            vertexId = pairs[pair].first.first;
            nodeId   = contourTree_->getCorrespondingNode(vertexId);
            nodeType = static_cast<int>(getNodeType(nodeId, type, contourTree_));
            break;
      }
   }
}

void vtkContourForests::getPersistenceDiagram(TreeType type)
{
   vtkSmartPointer<vtkPoints> points                       = vtkSmartPointer<vtkPoints>::New();
   vtkSmartPointer<vtkUnstructuredGrid> persistenceDiagram = vtkUnstructuredGrid::New();

   vtkSmartPointer<vtkDoubleArray> birthScalars = vtkSmartPointer<vtkDoubleArray>::New();
   birthScalars->SetName("Birth");

   vtkSmartPointer<vtkDoubleArray> deathScalars = vtkSmartPointer<vtkDoubleArray>::New();
   deathScalars->SetName("Death");

   int vertexId{};
   vtkSmartPointer<vtkIntArray> vertexIdentifierScalars = vtkSmartPointer<vtkIntArray>::New();
   vertexIdentifierScalars->SetName("VertexIdentifier");

   int nodeId{};
   vtkSmartPointer<vtkIntArray> nodeIdentifierScalars = vtkSmartPointer<vtkIntArray>::New();
   nodeIdentifierScalars->SetName("NodeIdentifier");

   int nodeType{};
   vtkSmartPointer<vtkIntArray> nodeTypeScalars = vtkSmartPointer<vtkIntArray>::New();
   nodeTypeScalars->SetName("NodeType");

   int pairIdentifier{};
   vtkSmartPointer<vtkIntArray> pairIdentifierScalars = vtkSmartPointer<vtkIntArray>::New();
   pairIdentifierScalars->SetName("PairIdentifier");

   vtkSmartPointer<vtkDoubleArray> scalars = vtkSmartPointer<vtkDoubleArray>::New();
   scalars->SetName(scalarField_.data());

   vector<pair<double, double>> diagram;
   switch (type) {
      case TreeType::MERGE_TREE:
         contourTree_->getJoinTree()->computePersistenceDiagram<double>(diagram, mergePairs_);
         break;
      case TreeType::SPLIT_TREE:
         contourTree_->getSplitTree()->computePersistenceDiagram<double>(diagram, splitPairs_);
         break;
      case TreeType::CONTOUR_TREE:
      default:
         contourTree_->computePersistenceDiagram<double>(diagram, mergePairs_, splitPairs_, pairs_);
         break;
   }

   unsigned int N = diagram.size();
   vtkIdType ids[3];
   double p[3] = {0, 0, 0};
   int point_cpt{};
   int cell_cpt{};
   if (N) {
      birthScalars->SetNumberOfComponents(1);
      deathScalars->SetNumberOfComponents(1);
      vertexIdentifierScalars->SetNumberOfComponents(1);
      nodeIdentifierScalars->SetNumberOfComponents(1);
      nodeTypeScalars->SetNumberOfComponents(1);
      scalars->SetNumberOfComponents(1);
      pairIdentifierScalars->SetNumberOfComponents(1);

      birthScalars->InsertTuple1(point_cpt, diagram[0].first);
      deathScalars->InsertTuple1(point_cpt, diagram[0].first);
      getPersistenceDiagramInfo(type, 0, false, vertexId, nodeId, nodeType);
      vertexIdentifierScalars->InsertTuple1(point_cpt, vertexId);
      nodeIdentifierScalars->InsertTuple1(point_cpt, nodeId);
      nodeTypeScalars->InsertTuple1(point_cpt, nodeType);
      scalars->InsertTuple1(point_cpt, (*vertexScalars_)[vertexId]);
      p[0]   = diagram[0].first;
      p[1]   = diagram[0].first;
      ids[0] = points->InsertNextPoint(p);
      ++point_cpt;

      birthScalars->InsertTuple1(point_cpt, diagram[0].first);
      deathScalars->InsertTuple1(point_cpt, diagram[0].second);
      getPersistenceDiagramInfo(type, 0, true, vertexId, nodeId, nodeType);
      vertexIdentifierScalars->InsertTuple1(point_cpt, vertexId);
      nodeIdentifierScalars->InsertTuple1(point_cpt, nodeId);
      nodeTypeScalars->InsertTuple1(point_cpt, nodeType);
      scalars->InsertTuple1(point_cpt, (*vertexScalars_)[vertexId]);

      p[0]   = diagram[0].first;
      p[1]   = diagram[0].second;
      ids[1] = points->InsertNextPoint(p);
      ++point_cpt;

      birthScalars->InsertTuple1(point_cpt, diagram[0].first);
      deathScalars->InsertTuple1(point_cpt, diagram[0].first);
      getPersistenceDiagramInfo(type, 0, false, vertexId, nodeId, nodeType);
      vertexIdentifierScalars->InsertTuple1(point_cpt, vertexId);
      nodeIdentifierScalars->InsertTuple1(point_cpt, nodeId);
      nodeTypeScalars->InsertTuple1(point_cpt, nodeType);
      scalars->InsertTuple1(point_cpt, (*vertexScalars_)[vertexId]);

      p[0]   = diagram[0].first;
      p[1]   = diagram[0].first;
      ids[2] = points->InsertNextPoint(p);
      ++point_cpt;

      pairIdentifierScalars->InsertTuple1(cell_cpt, pairIdentifier);
      persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);
      ++cell_cpt;
      pairIdentifierScalars->InsertTuple1(cell_cpt, pairIdentifier);
      persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids + 1);
      ++cell_cpt;

      for (unsigned int i = 1; i < N; ++i) {
         ids[0] = ids[2];

         if (type == TreeType::CONTOUR_TREE or
             (simplifyPersistenceDiagram_ and
              diagram[i].second - diagram[i].first > simplificationThreshold_)) {
            birthScalars->InsertTuple1(point_cpt, diagram[i].first);
            deathScalars->InsertTuple1(point_cpt, diagram[i].first);
            getPersistenceDiagramInfo(type, i, false, vertexId, nodeId, nodeType);
            vertexIdentifierScalars->InsertTuple1(point_cpt, vertexId);
            nodeIdentifierScalars->InsertTuple1(point_cpt, nodeId);
            nodeTypeScalars->InsertTuple1(point_cpt, nodeType);
            scalars->InsertTuple1(point_cpt, (*vertexScalars_)[vertexId]);

            p[0]   = diagram[i].first;
            p[1]   = diagram[i].first;
            ids[1] = points->InsertNextPoint(p);
            ++point_cpt;

            pairIdentifierScalars->InsertTuple1(cell_cpt, -1);
            persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);
            ++cell_cpt;

            ids[0] = ids[1];

            birthScalars->InsertTuple1(point_cpt, diagram[i].first);
            deathScalars->InsertTuple1(point_cpt, diagram[i].second);
            getPersistenceDiagramInfo(type, i, true, vertexId, nodeId, nodeType);
            vertexIdentifierScalars->InsertTuple1(point_cpt, vertexId);
            nodeIdentifierScalars->InsertTuple1(point_cpt, nodeId);
            nodeTypeScalars->InsertTuple1(point_cpt, nodeType);
            scalars->InsertTuple1(point_cpt, (*vertexScalars_)[vertexId]);

            p[0]   = diagram[i].first;
            p[1]   = diagram[i].second;
            ids[1] = points->InsertNextPoint(p);
            ++point_cpt;

            pairIdentifierScalars->InsertTuple1(cell_cpt, i);
            persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);
            ++cell_cpt;

            birthScalars->InsertTuple1(point_cpt, diagram[i].first);
            deathScalars->InsertTuple1(point_cpt, diagram[i].first);
            getPersistenceDiagramInfo(type, i, false, vertexId, nodeId, nodeType);
            vertexIdentifierScalars->InsertTuple1(point_cpt, vertexId);
            nodeIdentifierScalars->InsertTuple1(point_cpt, nodeId);
            nodeTypeScalars->InsertTuple1(point_cpt, nodeType);
            scalars->InsertTuple1(point_cpt, (*vertexScalars_)[vertexId]);

            p[0]   = diagram[i].first;
            p[1]   = diagram[i].first;
            ids[2] = points->InsertNextPoint(p);
            ++point_cpt;

            pairIdentifierScalars->InsertTuple1(cell_cpt, i);
            persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids + 1);
            ++cell_cpt;
         }
      }

      persistenceDiagram->SetPoints(points);
      persistenceDiagram->GetPointData()->AddArray(birthScalars);
      persistenceDiagram->GetPointData()->AddArray(deathScalars);
      persistenceDiagram->GetPointData()->AddArray(vertexIdentifierScalars);
      persistenceDiagram->GetPointData()->AddArray(nodeIdentifierScalars);
      persistenceDiagram->GetPointData()->AddArray(nodeTypeScalars);
      persistenceDiagram->GetPointData()->AddArray(scalars);
      persistenceDiagram->GetCellData()->AddArray(pairIdentifierScalars);

      switch (type) {
         case TreeType::MERGE_TREE:
            MTPersistenceDiagram_->ShallowCopy(persistenceDiagram);
            break;
         case TreeType::SPLIT_TREE:
            STPersistenceDiagram_->ShallowCopy(persistenceDiagram);
            break;
         case TreeType::CONTOUR_TREE:
            CTPersistenceDiagram_->ShallowCopy(persistenceDiagram);
            break;
      }
   }
}

void vtkContourForests::getDiagrams()
{
   if (treeType_ == TreeType::CONTOUR_TREE) {
      getPersistenceDiagram(TreeType::CONTOUR_TREE);
   } else {
      getPersistenceDiagram(TreeType::MERGE_TREE);
      getPersistenceDiagram(TreeType::SPLIT_TREE);
   }

   {
      stringstream msg;
      msg << "[vtkContourForests] Persistence diagrams built." << endl;
      dMsg(cout, msg.str(), infoMsg);
   }
}

int vtkContourForests::sample(unsigned int samplingLevel)
{
   samples_->resize(3);
   (*samples_)[static_cast<int>(treeType_)].resize(tree_->getNumberOfSuperArcs());
   vector<vector<int>> sampleList(samplingLevel);

   SuperArc* a;
   for (unsigned int i = 0; i < tree_->getNumberOfSuperArcs(); ++i) {
      a = tree_->getSuperArc(i);

      if (!a->isPruned()) {
         for (unsigned int j = 0; j < samplingLevel; ++j) {
            sampleList[j].clear();
         }

         double fmax, fmin;
         int nodeMaxId, nodeMinId;
         int nodeMaxVId, nodeMinVId;
         double delta;
         if (a->getNumberOfRegularNodes()) {
            if (treeType_ == TreeType::SPLIT_TREE) {
               nodeMaxId = a->getDownNodeId();
               nodeMinId = a->getUpNodeId();
            } else {
               nodeMaxId = a->getUpNodeId();
               nodeMinId = a->getDownNodeId();
            }

            nodeMaxVId = tree_->getNode(nodeMaxId)->getVertexId();
            nodeMinVId = tree_->getNode(nodeMinId)->getVertexId();

            fmax = (*vertexScalars_)[nodeMaxVId];
            fmin = (*vertexScalars_)[nodeMinVId];

            delta = (fmax - fmin) / samplingLevel;

            double f;
            int nodeId;
            int vertexId;
            for (int j = 0; j < a->getNumberOfRegularNodes(); ++j) {
               nodeId   = a->getRegularNodeId(j);
               if(a->isMasqued(j)) continue;
               vertexId = nodeId;
               f        = (*vertexScalars_)[vertexId];

               for (unsigned int k = 0; k < samplingLevel; ++k) {
                  if (f <= (k + 1) * delta + fmin) {
                     sampleList[k].push_back(nodeId);
                     break;
                  }
               }
            }

            // update the arc
            for (unsigned int j = 0; j < sampleList.size(); ++j)
               (*samples_)[static_cast<int>(treeType_)][i].push_back(sampleList[j]);
         }
      }
   }

   return 0;
}

int vtkContourForests::computeBarycenters()
{
   barycenters_->resize(3);
   (*barycenters_)[static_cast<int>(treeType_)].resize(tree_->getNumberOfSuperArcs());
   vector<double> barycenter(3);
   int vertexId;

   const SuperArc* a;
   for (unsigned int i = 0; i < tree_->getNumberOfSuperArcs(); ++i) {
      a = tree_->getSuperArc(i);
      if (!a->isPruned()) {
         for (unsigned int j = 0; j < (*samples_)[static_cast<int>(treeType_)][i].size(); ++j) {
            vector<int>& sample = (*samples_)[static_cast<int>(treeType_)][i][j];

            for (unsigned int k = 0; k < 3; ++k)
               barycenter[k] = 0;

            for (unsigned int k = 0; k < sample.size(); ++k) {
               vertexId = sample[k];


               for (unsigned int l = 0; l < 3; ++l)
                  barycenter[l] += (*vertexPositions_)[vertexId][l];
            }
            if (sample.size()) {
               for (unsigned int k = 0; k < 3; ++k)
                  barycenter[k] /= sample.size();

               // update the arc
               unsigned int nbBar = (*barycenters_)[static_cast<int>(treeType_)][i].size();
               (*barycenters_)[static_cast<int>(treeType_)][i].resize(nbBar + 1);
               (*barycenters_)[static_cast<int>(treeType_)][i][nbBar].resize(3);

               for (unsigned int k = 0; k < 3; ++k)
                  (*barycenters_)[static_cast<int>(treeType_)][i][nbBar][k] = barycenter[k];
            }
         }
      }
   }

   return 0;
}

void vtkContourForests::computeSkeleton(unsigned int arcRes)
{
   sample(arcRes);
   computeBarycenters();
}

void vtkContourForests::smoothSkeleton(unsigned int skeletonSmoothing)
{
   for (unsigned int i = 0; i < skeletonSmoothing; i++) {
      for (unsigned int j = 0; j < tree_->getNumberOfSuperArcs(); j++) {
         if (!tree_->getSuperArc(j)->isPruned()) {
            smooth(j, !(treeType_ == TreeType::SPLIT_TREE));
         }
      }
   }
}

void vtkContourForests::smooth(const int idArc, bool order)
{
   int N = (*barycenters_)[static_cast<int>(treeType_)][idArc].size();
   if (N) {
      /// init ///
      vector<vector<double>> barycenterList(N);
      for (unsigned int i = 0; i < barycenterList.size(); ++i)
         barycenterList[i].resize(3);

      int up_vId;
      int down_vId;
      if (order) {
         up_vId   = tree_->getNode(tree_->getSuperArc(idArc)->getUpNodeId())->getVertexId();
         down_vId = tree_->getNode(tree_->getSuperArc(idArc)->getDownNodeId())->getVertexId();
      } else {
         down_vId = tree_->getNode(tree_->getSuperArc(idArc)->getUpNodeId())->getVertexId();
         up_vId   = tree_->getNode(tree_->getSuperArc(idArc)->getDownNodeId())->getVertexId();
      }

      double p0[3];
      double p1[3];
      for (unsigned int k = 0; k < 3; ++k) {
         p0[k] = (*vertexPositions_)[down_vId][k];
         p1[k] = (*vertexPositions_)[up_vId][k];
      }

      /// filtering ///
      if (N > 1) {
         // first
         for (unsigned int k = 0; k < 3; ++k)
            barycenterList[0][k] =
                (p0[k] + (*barycenters_)[static_cast<int>(treeType_)][idArc][1][k]) * 0.5;

         // main
         for (int i = 1; i < N - 1; ++i) {
            for (unsigned int k = 0; k < 3; ++k)
               barycenterList[i][k] =
                   ((*barycenters_)[static_cast<int>(treeType_)][idArc][i - 1][k] +
                    (*barycenters_)[static_cast<int>(treeType_)][idArc][i + 1][k]) *
                   0.5;
         }
         // last
         for (unsigned int k = 0; k < 3; ++k)
            barycenterList[N - 1][k] =
                (p1[k] + (*barycenters_)[static_cast<int>(treeType_)][idArc][N - 1][k]) * 0.5;
      } else {
         for (unsigned int k = 0; k < 3; ++k)
            barycenterList[0][k] = (p0[k] + p1[k]) * 0.5;
      }

      /// copy ///
      for (int i = 0; i < N; ++i) {
         for (unsigned int k = 0; k < 3; ++k)
            (*barycenters_)[static_cast<int>(treeType_)][idArc][i][k] = barycenterList[i][k];
      }
   }
}

void vtkContourForests::getSkeleton()
{
   Timer t;
   computeSkeleton(arcResolution_);
   smoothSkeleton(skeletonSmoothing_);

   // nodes
   if (showMin_ || showMax_ || showSaddle1_ || showSaddle2_)
      getSkeletonNodes();
   else
      skeletonNodes_->ShallowCopy(voidUnstructuredGrid_);

   // arcs
   if (showArc_)
      getSkeletonArcs();
   else
      skeletonArcs_->ShallowCopy(voidPolyData_);

   // ce qui est fait n'est plus à faire
   toComputeSkeleton_ = false;

   {
      stringstream msg;
      msg << "[vtkContourForests] Topological skeleton built in " << t.getElapsedTime()
          << "s :" << endl;
      msg << "[vtkContourForests]   Arc - Resolution: " << arcResolution_ << endl;
      msg << "[vtkContourForests]   Smoothing: " << boolalpha << skeletonSmoothing_ << endl;
      dMsg(cout, msg.str(), timeMsg);
   }
}

void vtkContourForests::getSegmentation(vtkDataSet* input)
{
   Timer t;

   // field
   int regionId{};
   vtkSmartPointer<vtkIntArray> scalarsRegionId = vtkSmartPointer<vtkIntArray>::New();
   scalarsRegionId->SetName("SegmentationId");
   scalarsRegionId->SetNumberOfTuples(vertexScalars_->size());

   int regionType{};
   vtkSmartPointer<vtkIntArray> scalarsRegionType = vtkSmartPointer<vtkIntArray>::New();
   scalarsRegionType->SetName("RegionType");
   scalarsRegionType->SetNumberOfTuples(vertexScalars_->size());

   int regionSize{};
   vtkSmartPointer<vtkIntArray> scalarsRegionSize = vtkSmartPointer<vtkIntArray>::New();
   scalarsRegionSize->SetName("RegionSize");
   scalarsRegionSize->SetNumberOfTuples(vertexScalars_->size());

   double regionSpan{};
   vtkSmartPointer<vtkDoubleArray> scalarsRegionSpan = vtkSmartPointer<vtkDoubleArray>::New();
   scalarsRegionSpan->SetName("RegionSpan");
   scalarsRegionSpan->SetNumberOfTuples(vertexScalars_->size());

   int currentZone{};

   if (!segmentation_) {
      segmentation_ = input->NewInstance();
      segmentation_->ShallowCopy(input);
   }

   for (unsigned int i = 0; i < numberOfVertices_; i++) {
      scalarsRegionId->SetTuple1(i, -1);
   }

   // nodes
   for (unsigned int it = 0; it < criticalPoints_->size(); ++it) {
      int nodeId   = (*criticalPoints_)[it];
      int vertexId = tree_->getNode(nodeId)->getVertexId();

      // RegionType
      regionType = -1;
      scalarsRegionType->SetTuple1(vertexId, regionType);
   }

   // arcs
   for (unsigned int i = 0; i < tree_->getNumberOfSuperArcs(); ++i) {
      auto a = tree_->getSuperArc(i);
      if (a->isVisible()) {
         int      upNodeId   = tree_->getSuperArc(i)->getUpNodeId();
         NodeType upNodeType = getNodeType(upNodeId);
         int      upVertex   = tree_->getNode(upNodeId)->getVertexId();
         float    coordUp[3];
         triangulation_->getVertexPoint(upVertex, coordUp[0], coordUp[1], coordUp[2]);

         int      downNodeId   = tree_->getSuperArc(i)->getDownNodeId();
         NodeType downNodeType = getNodeType(downNodeId);
         int      downVertex   = tree_->getNode(downNodeId)->getVertexId();
         float    coordDown[3];
         triangulation_->getVertexPoint(downVertex, coordDown[0], coordDown[1], coordDown[2]);

         regionSize = tree_->getNumberOfVisibleRegularNode(i);
         regionSpan = Geometry::distance(coordUp, coordDown);
         regionId   = currentZone++;
         //regionId = i;

         //cout << "arc : " << tree_->printArc(i);
         //cout << " span : " << regionSpan;
         //cout << " coords : ";
         //cout << coordDown[0] << ",";
         //cout << coordDown[1] << ",";
         //cout << coordDown[2] << " || ";
         //cout << coordUp[0] << ",";
         //cout << coordUp[1] << ",";
         //cout << coordUp[2] << endl;

         scalarsRegionId->SetTuple1(tree_->getNode(downNodeId)->getVertexId(), regionId);
         scalarsRegionId->SetTuple1(tree_->getNode(upNodeId)->getVertexId(), regionId);

         scalarsRegionSize->SetTuple1(tree_->getNode(downNodeId)->getVertexId(), regionSize);
         scalarsRegionSize->SetTuple1(tree_->getNode(upNodeId)->getVertexId(), regionSize);

         scalarsRegionSpan->SetTuple1(tree_->getNode(downNodeId)->getVertexId(), regionSpan);
         scalarsRegionSpan->SetTuple1(tree_->getNode(upNodeId)->getVertexId(), regionSpan);

         for (int j = 0; j < tree_->getSuperArc(i)->getNumberOfRegularNodes(); ++j) {
            int nodeId   = tree_->getSuperArc(i)->getRegularNodeId(j);
            int vertexId = nodeId;
            // cout << vertexId << ", ";
            if(tree_->getSuperArc(i)->isMasqued(j)) {
               // cout << vertexId << ", ";
               continue;
            }

            //cout << vertexId << ", ";
            scalarsRegionId->SetTuple1(vertexId, regionId);
            scalarsRegionSize->SetTuple1(vertexId, regionSize);
            scalarsRegionSpan->SetTuple1(vertexId, regionSpan);
         }
         //cout << endl;

         // RegionType
         if (upNodeType == NodeType::LOCAL_MINIMUM && downNodeType == NodeType::LOCAL_MAXIMUM)
            regionType = static_cast<int>(ArcType::MIN_ARC);
         else if (upNodeType == NodeType::LOCAL_MINIMUM || downNodeType == NodeType::LOCAL_MINIMUM)
            regionType = static_cast<int>(ArcType::MIN_ARC);
         else if (upNodeType == NodeType::LOCAL_MAXIMUM || downNodeType == NodeType::LOCAL_MAXIMUM)
            regionType = static_cast<int>(ArcType::MAX_ARC);
         else if (upNodeType == NodeType::SADDLE1 && downNodeType == NodeType::SADDLE1)
            regionType = static_cast<int>(ArcType::SADDLE1_ARC);
         else if (upNodeType == NodeType::SADDLE2 && downNodeType == NodeType::SADDLE2)
            regionType = static_cast<int>(ArcType::SADDLE2_ARC);
         else
            regionType = static_cast<int>(ArcType::SADDLE1_SADDLE2_ARC);

         for (int j = 0; j < tree_->getSuperArc(i)->getNumberOfRegularNodes(); ++j) {
            int nodeId   = tree_->getSuperArc(i)->getRegularNodeId(j);
            if(tree_->getSuperArc(i)->isMasqued(j)){
                // Ignore masqued ones
                continue;
            }
            int vertexId = nodeId;
            scalarsRegionType->SetTuple1(vertexId, regionType);
         }
      }
   }

   // output
   segmentation_->GetPointData()->AddArray(scalarsRegionId);
   segmentation_->GetPointData()->AddArray(scalarsRegionType);
   segmentation_->GetPointData()->AddArray(scalarsRegionSize);
   segmentation_->GetPointData()->AddArray(scalarsRegionSpan);

   // ce qui est fait n'est plus à faire
   toComputeSegmentation_ = false;

   {
      stringstream msg;
      msg << "[vtkContourForests] Topological segmentation built in " << t.getElapsedTime()
          << "s :" << endl;
      msg << "[vtkContourForests]   RegionType: " << boolalpha
          << (bool)scalarsRegionType->GetNumberOfTuples() << endl;
      msg << "[vtkContourForests]   SegmentationId: " << boolalpha
          << (bool)scalarsRegionId->GetNumberOfTuples() << endl;
      dMsg(cout, msg.str(), timeMsg);
   }
}

void vtkContourForests::getTree()
{
   // set ContourTree input and call
   setDebugLevel(debugLevel_);
   contourTree_->setTriangulation(triangulation_);
   contourTree_->setVertexScalars(vertexScalars_->data());
   contourTree_->setPartitionNum(partitionNum_);
   contourTree_->setSimplificationMethod(simplificationType_);
   contourTree_->build<double>(treeType_ == TreeType::CONTOUR_TREE, calculSegmentation_,
                               simplificationThreshold_);
   setDebugLevel(1);

   // ce qui est fait n'est plus à faire
   toComputeContourTree_    = false;
   toComputeSimplification_ = false;
}

void vtkContourForests::updateTree()
{
   // polymorphic tree
   switch (treeType_) {
      case TreeType::MERGE_TREE:
         tree_ = contourTree_->getJoinTree();
         break;
      case TreeType::SPLIT_TREE:
         tree_ = contourTree_->getSplitTree();
         break;
      case TreeType::CONTOUR_TREE:
         tree_ = contourTree_;
         break;
   }

   getCriticalPoints();
   if (!originalNumberOfCriticalPoints_)
      originalNumberOfCriticalPoints_ = criticalPoints_->size();
   //getSimplificationThresholdCurve();

   toUpdateTree_ = false;
}

int vtkContourForests::check(vtkDataSet* input)
{
   // Check point data
   bool hasPointData{};
   if (!input->GetPointData()) {
      vtkGenericWarningMacro(<< "Input has no point data.");
      return -1;
   } else
      hasPointData = true;

   // Check scalar data
   bool hasScalarData{};
   int N = input->GetPointData()->GetNumberOfArrays();
   for (int i = 0; i < N; ++i) {
      vtkDataArray* inputArray = input->GetPointData()->GetArray(i);
      if (inputArray->GetNumberOfComponents() == 1)
         hasScalarData = true;
   }
   if (!hasScalarData) {
      vtkGenericWarningMacro(<< "Input has no scalar data.");
      return -2;
   }

   /// Check triangulation
   bool isTriangle{true};
   bool isTetra{true};
   vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
   for (int i = 0; i < input->GetNumberOfCells() && (isTetra || isTriangle); ++i) {
      input->GetCell(i, cell);

      if (cell->GetCellType() != VTK_TRIANGLE)
         isTriangle = false;
      if (cell->GetCellType() != VTK_TETRA)
         isTetra = false;
   }
   if (!isTriangle && !isTetra) {
      vtkGenericWarningMacro(<< "Input mesh is not tetrahedralized.");
      return -3;
   }

   /// Check connectivity
   vtkSmartPointer<vtkConnectivityFilter> connectivityFilter =
       vtkSmartPointer<vtkConnectivityFilter>::New();
   connectivityFilter->SetInputData(input);
   connectivityFilter->Update();

   if (connectivityFilter->GetNumberOfExtractedRegions() > 1) {
      vtkGenericWarningMacro(<< "Input mesh has multiple connected components.");
      return -4;
   }

   isChecked_ = true;

   {
      stringstream msg;
      msg << "[vtkContourForests] Input has been checked:" << endl;
      msg << "[vtkContourForests]   Point data: " << boolalpha << hasPointData << endl;
      msg << "[vtkContourForests]   Scalar data: " << boolalpha << hasScalarData << endl;
      msg << "[vtkContourForests]   Tetrahedralized: " << boolalpha << (isTriangle || isTetra)
          << endl;
      msg << "[vtkContourForests]   Number of points: " << input->GetNumberOfPoints() << endl;
      msg << "[vtkContourForests]   Number of cells: " << input->GetNumberOfCells() << endl;
      dMsg(cout, msg.str(), detailedInfoMsg);
   }

   return 0;
}

int vtkContourForests::doIt(vtkDataSet* input, vtkPolyData* outputSkeletonNodes,
                           vtkPolyData* outputSkeletonArcs, vtkDataSet* outputSegmentation,
                           vtkTable* outputPersistenceCurve,
                           vtkTable* outputSimplificationThresholdCurve,
                           vtkUnstructuredGrid* outputPersistenceDiagram)
{
   if (checkConformance_ and !isChecked_)
      check(input);
   else {

      // conversion
      vtkDataSetToStdVector(input);

       if(simplificationType_ == 0){
          simplificationThreshold_ = simplificationThresholdBuffer_ * deltaScalar_;
       } else if (simplificationType_ == 1){
          double coord0[3], coord1[3], spanTotal;
          double* bounds = input->GetBounds();
          coord0[0] = bounds[0];
          coord1[0] = bounds[1];
          coord0[1] = bounds[2];
          coord1[1] = bounds[3];
          coord0[2] = bounds[4];
          coord1[2] = bounds[5];
          spanTotal = Geometry::distance(coord0,coord1);
          simplificationThreshold_ =
              simplificationThresholdBuffer_ * spanTotal;
       } else if (simplificationType_ == 2) {
          simplificationThreshold_ =
              simplificationThresholdBuffer_ * triangulation_->getNumberOfVertices();
       }

       //cout << "doit threshold : " << simplificationThreshold_ << endl;

      /// ContourTree ///
      if (varyingMeshGeometry_ || varyingMeshConnectivity_ || varyingDataValues_ ||
          toComputeContourTree_) {
         clearTree();
         getTree();
         //getCurves();
      }

      // persistence diagram
      //if (varyingMeshGeometry_ || varyingMeshConnectivity_ || varyingDataValues_ ||
          //toComputeContourTree_ || toComputeSimplification_) {
         //getDiagrams();
      //}

      // update the trees
      if (varyingMeshGeometry_ || varyingMeshConnectivity_ || varyingDataValues_ || toUpdateTree_)
         updateTree();

      /// Skeleton ///
      if (varyingMeshGeometry_ || varyingMeshConnectivity_ || varyingDataValues_ ||
          toComputeSkeleton_) {
         clearSkeleton();
         getSkeleton();
      }

      /// Segmentation ///
      if (varyingMeshGeometry_ || varyingMeshConnectivity_ || varyingDataValues_ ||
          toComputeSegmentation_) {
         getSegmentation(input);
      }

      /// Output ///
      // skeleton
      outputSkeletonNodes->ShallowCopy(skeletonNodes_);
      outputSkeletonArcs->ShallowCopy(skeletonArcs_);
      // segmentation
      outputSegmentation->ShallowCopy(segmentation_);
      // persistence curve
      if (false and toShowPersistenceCurve_) {
         switch (treeType_) {
            case TreeType::MERGE_TREE:
               outputPersistenceCurve->ShallowCopy(MTPersistenceCurve_);
               break;
            case TreeType::SPLIT_TREE:
               outputPersistenceCurve->ShallowCopy(STPersistenceCurve_);
               break;
            case TreeType::CONTOUR_TREE:
            default:
               outputPersistenceCurve->ShallowCopy(CTPersistenceCurve_);
               break;
         }

         // simplification threshold curve
         outputSimplificationThresholdCurve->ShallowCopy(simplificationThresholdCurve_);
      } else {
         //outputPersistenceCurve->ShallowCopy(voidTable_);
         //outputSimplificationThresholdCurve->ShallowCopy(voidTable_);
      }

      // persistence diagram
      if ( false and toShowPersistenceDiagram_) {
         switch (treeType_) {
            case TreeType::MERGE_TREE:
               outputPersistenceDiagram->ShallowCopy(MTPersistenceDiagram_);
               break;
            case TreeType::SPLIT_TREE:
               outputPersistenceDiagram->ShallowCopy(STPersistenceDiagram_);
               break;
            case TreeType::CONTOUR_TREE:
            default:
               outputPersistenceDiagram->ShallowCopy(CTPersistenceDiagram_);
               break;
         }
      } else {
          //outputPersistenceDiagram->ShallowCopy(voidUnstructuredGrid_);
      }
   }

   return 0;
}

int vtkContourForests::RequestData(vtkInformation* request, vtkInformationVector** inputVector,
                                  vtkInformationVector* outputVector)
{
   Memory m;
   vtkDataSet* input = vtkDataSet::GetData(inputVector[0]);

   /// Skeleton Nodes ///
   vtkInformation* outInfo = outputVector->GetInformationObject(0);
   vtkPolyData* outputSkeletonNodes =
       vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

   /// Skeleton Arcs ///
   outInfo = outputVector->GetInformationObject(1);
   vtkPolyData* outputSkeletonArcs =
       vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

   /// Segmentation ///
   outInfo = outputVector->GetInformationObject(2);
   vtkDataSet* outputSegmentation =
       vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

   /// Persistence Curve ///
   //outInfo = outputVector->GetInformationObject(3);
   vtkTable* outputPersistenceCurve = nullptr;
       //vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

   /// Simplification Threshold Curve ///
   //outInfo = outputVector->GetInformationObject(4);
   vtkTable* outputSimplificationThresholdCurve = nullptr;
       //vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

   /// Persistence Diagram ///
   //outInfo = outputVector->GetInformationObject(5);
   vtkUnstructuredGrid* outputPersistenceDiagram = nullptr;
       //vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

   // processing
   doIt(input, outputSkeletonNodes, outputSkeletonArcs, outputSegmentation, outputPersistenceCurve,
        outputSimplificationThresholdCurve, outputPersistenceDiagram);

   {
      stringstream msg;
      msg << "[vtkContourForests] Memory usage: " << m.getElapsedUsage() << " MB." << endl;
      dMsg(cout, msg.str(), memoryMsg);
   }

   return 1;
}
