#include "CompositeToGridTransformCLP.h"

// MRML includes
#include <vtkMRMLTransformNode.h>
#include <vtkMRMLScalarVolumeNode.h>
#include <vtkMRMLVolumeArchetypeStorageNode.h>
#include <vtkMRMLTransformStorageNode.h>

// VTK includes
#include <vtkGeneralTransform.h>
#include <vtkImageData.h>
#include <vtkMatrix4x4.h>
#include <vtkOrientedGridTransform.h>
#include <vtkOrientedGridTransform.h>

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{
void GetTransformedPointSamplesAsVectorImage(vtkImageData* vectorImage, vtkMRMLTransformNode* inputTransformNode, vtkMatrix4x4* ijkToRAS)
{
  vtkNew<vtkGeneralTransform> inputTransform;
  inputTransformNode->GetTransformFromWorld(inputTransform.GetPointer());

  // The orientation of the volume cannot be set in the image
  // therefore the volume will not appear in the correct position
  // if the direction matrix is not identity.
  vectorImage->AllocateScalars(VTK_FLOAT, 3);

  double point_RAS[4] = { 0, 0, 0, 1 };
  double transformedPoint_RAS[4] = { 0, 0, 0, 1 };
  double point_IJK[4] = { 0, 0, 0, 1 };
  float* voxelPtr = static_cast<float*>(vectorImage->GetScalarPointer());
  int* extent = vectorImage->GetExtent();
  int* dim = vectorImage->GetDimensions();
  float numberOfVoxels = dim[0] * dim[1] * dim[2];
  unsigned int voxelCount = 0;
  for (point_IJK[2] = extent[4]; point_IJK[2] <= extent[5]; point_IJK[2]++)
  {
    for (point_IJK[1] = extent[2]; point_IJK[1] <= extent[3]; point_IJK[1]++)
    {
      for (point_IJK[0] = extent[0]; point_IJK[0] <= extent[1]; point_IJK[0]++)
      {
        ijkToRAS->MultiplyPoint(point_IJK, point_RAS);

        inputTransform->TransformPoint(point_RAS, transformedPoint_RAS);

        // store the pointDislocationVector_RAS components in the image
        *(voxelPtr++) = static_cast<float>(transformedPoint_RAS[0] - point_RAS[0]);
        *(voxelPtr++) = static_cast<float>(transformedPoint_RAS[1] - point_RAS[1]);
        *(voxelPtr++) = static_cast<float>(transformedPoint_RAS[2] - point_RAS[2]);

        voxelCount++;
        if (voxelCount % 10000 == 0)
        {
          std::cout << "<filter-progress>" << (voxelCount / numberOfVoxels) << "</filter-progress>" << std::endl << std::flush;
        }
      }
    }
  }

}

} // end of anonymous namespace

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  // INPUT TRANSFORM 1

  bool transform1FromNode = inputTransform1Node.compare("None") != 0;
  bool transform1FromFile = !inputTransform1File.empty();

  if ((transform1FromNode && transform1FromFile) || (!transform1FromNode && !transform1FromFile))
  {
    std::cerr << "Specify either a MRML transform node or a transform file name for input transform 1." << std::endl;
    return EXIT_FAILURE;
  }

  vtkNew<vtkMRMLTransformNode> transform1Node;
  vtkNew<vtkMRMLTransformStorageNode> storageNode1;
  storageNode1->SetFileName(transform1FromNode ? inputTransform1Node.c_str() : inputTransform1File.c_str());
  if (!storageNode1->ReadData(transform1Node))
  {
    std::cerr << "Failed to read input transform 1 from file " << std::endl;
    return EXIT_FAILURE;
  }

  // INPUT TRANSFORM 2

  bool transform2FromNode = inputTransform2Node.compare("None") != 0;
  bool transform2FromFile = !inputTransform2File.empty();

  if ((transform2FromNode && transform2FromFile) || (!transform2FromNode && !transform2FromFile))
  {
    std::cerr << "Specify either a MRML transform node or a transform file name for input transform 2." << std::endl;
    return EXIT_FAILURE;
  }

  vtkNew<vtkMRMLTransformNode> transform2Node;
  vtkNew<vtkMRMLTransformStorageNode> storageNode2;
  storageNode2->SetFileName(transform2FromNode ? inputTransform2Node.c_str() : inputTransform2File.c_str());
  if (!storageNode2->ReadData(transform2Node))
  {
    std::cerr << "Failed to read input transform 2 from file " << std::endl;
    return EXIT_FAILURE;
  }

  // REFERENCE VOLUME

  bool referenceVolumeFromNode = inputReferenceVolumeNode.compare("None") != 0;
  bool referenceVolumeFromFile = !inputReferenceVolumeFile.empty();

  if ((referenceVolumeFromNode && referenceVolumeFromFile) || (!referenceVolumeFromNode && !referenceVolumeFromFile))
  {
    std::cerr << "Specify either a MRML volume node or a volume file name for reference volume." << std::endl;
    return EXIT_FAILURE;
  }

  vtkNew<vtkMRMLScalarVolumeNode> referenceVolumeNode;
  vtkNew<vtkMRMLVolumeArchetypeStorageNode> storageNode3;
  storageNode3->SetFileName(referenceVolumeFromNode ? inputReferenceVolumeNode.c_str() : inputReferenceVolumeFile.c_str());
  if (!storageNode3->ReadData(referenceVolumeNode))
  {
    std::cerr << "Failed to read reference volume from file " << std::endl;
    return EXIT_FAILURE;
  }

  // OUTPUT

  bool saveToNode = outputDisplacementField.compare("None") != 0;
  bool saveToFile = !outputFileName.empty();

  if ((saveToNode && saveToFile) || (!saveToNode && !saveToFile))
  {
    std::cerr << "Specify either a MRML output node or an output file name" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "<filter-comment>" << "Set up" << "</filter-comment>" << std::endl << std::flush;

  // Create a grid transform
  vtkSmartPointer<vtkMRMLTransformNode> outputGridTransformNode;
  outputGridTransformNode = vtkSmartPointer<vtkMRMLTransformNode>::New();

  vtkOrientedGridTransform* outputGridTransform = vtkOrientedGridTransform::SafeDownCast(
    outputGridTransformNode->GetTransformToParentAs("vtkOrientedGridTransform",
    false /* don't report conversion error */,
    true /* we would like to modify the transform */));
  if (outputGridTransform == nullptr)
  {
    // we cannot reuse the existing transform, create a new one
    vtkNew<vtkOrientedGridTransform> newOutputGridTransform;
    outputGridTransform = newOutputGridTransform.GetPointer();
    outputGridTransformNode->SetAndObserveTransformFromParent(outputGridTransform);
  }
  // Create/get displacement field image
  vtkImageData* outputVolume = outputGridTransform->GetDisplacementGrid();
  if (outputVolume == nullptr)
  {
    vtkNew<vtkImageData> newOutputVolume;
    outputVolume = newOutputVolume.GetPointer();
    outputGridTransform->SetDisplacementGridData(outputVolume);
  }
  // Update geometry based on reference image
  vtkNew<vtkMatrix4x4> ijkToRas; // RAS refers to world
  if (referenceVolumeNode != nullptr)
  {
    referenceVolumeNode->GetIJKToRASMatrix(ijkToRas.GetPointer());
    vtkNew<vtkMatrix4x4> rasToWorld;
    if (vtkMRMLTransformNode::GetMatrixTransformBetweenNodes(referenceVolumeNode->GetParentTransformNode(), nullptr /* world */, rasToWorld.GetPointer()))
    {
      vtkMatrix4x4::Multiply4x4(rasToWorld.GetPointer(), ijkToRas.GetPointer(), ijkToRas.GetPointer());
    }
    else
    {
      std::cerr << "vtkSlicerTransformLogic::ConvertToGridTransform: non-linearly transformed reference volume" \
       " is not supported. Harden or remove the transform from of the reference volume." << std::endl;
      return EXIT_FAILURE;
    }
    double origin[3] = { 0, 0, 0 };
    double spacing[3] = { 1, 1, 1 };
    vtkNew<vtkMatrix4x4> ijkToRasDirection; // normalized direction matrix
    for (int c = 0; c < 3; c++)
    {
      origin[c] = ijkToRas->GetElement(c, 3);
      spacing[c] = sqrt(ijkToRas->Element[0][c] * ijkToRas->Element[0][c]
        + ijkToRas->Element[1][c] * ijkToRas->Element[1][c]
        + ijkToRas->Element[2][c] * ijkToRas->Element[2][c]);
      if (spacing[c] == 0)
      {
        // Prevent division by zero in case there is a projection matrix is in the transform chain
        spacing[c] = 1.0;
      }
      for (int r = 0; r < 3; r++)
      {
        ijkToRasDirection->SetElement(r, c, ijkToRas->GetElement(r, c) / spacing[c]);
      }
    }
    outputVolume->SetExtent(referenceVolumeNode->GetImageData()->GetExtent());
    outputVolume->SetOrigin(origin);
    outputVolume->SetSpacing(spacing);
    // vtkImageData cannot store directions, therefore that has to be set in the grid transform
    outputGridTransform->SetGridDirectionMatrix(ijkToRasDirection.GetPointer());
  }

  // RUN

  vtkNew<vtkGeneralTransform> hardeningTransform;
  transform2Node->GetTransformToWorld(hardeningTransform.GetPointer());
  transform1Node->ApplyTransform(hardeningTransform.GetPointer());

  std::cout << "<filter-comment>" << "Computing" << "</filter-comment>" << std::endl << std::flush;
  GetTransformedPointSamplesAsVectorImage(outputVolume, transform1Node, ijkToRas.GetPointer());

  std::cout << "<filter-comment>" << "Writing" << "</filter-comment>" << std::endl << std::flush;
  vtkNew<vtkMRMLTransformStorageNode> storageNode;
  storageNode->SetFileName(saveToNode ? outputDisplacementField.c_str() : outputFileName.c_str());
  if (!storageNode->WriteData(outputGridTransformNode))
  {
    std::cerr << "Failed to write output transform" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
