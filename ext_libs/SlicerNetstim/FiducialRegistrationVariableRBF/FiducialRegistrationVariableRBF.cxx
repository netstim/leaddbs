// The logic of this computation is taken from plastimatch and modified to use variable RBF Radius.
// Instead of using plastimatch classes, here ITK and SimpleITK abstractions are used.
// https://gitlab.com/plastimatch/plastimatch

#include "FiducialRegistrationVariableRBFCLP.h"

// ITK includes
#include <itkImage.h>

#include "SimpleITK.h"

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

namespace sitk = itk::simple;

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{
  // Function to convert a point from std::vector to itk::Point
  itk::Point<double, 3> convertStdVectorToITKPoint(const std::vector<float> & vec)
  {
    itk::Point<double, 3> p;
    p[0] = vec[0];
    p[1] = vec[1];
    p[2] = vec[2];
    return p;
  }

  static float RBFValue (itk::Point<double, 3> * rbf_center, itk::Point<double, 3> * loc, float radius)
  {
    float r = rbf_center->EuclideanDistanceTo(*loc) / radius;
    float val = exp( -r*r );   
    return val;
  }

  void RBFGaussUpdateVectorField( 
    sitk::Image *vectorField, 
    float *coeff, 
    std::vector< itk::Point<double, 3> > * fixedLandmarks,
    float * adaptRadius)
  {
    unsigned int i, landmarkIndex, imageLinearIndex;
    float rbf;
    unsigned int numLandmarks = fixedLandmarks->size();

    itk::Point<double, 3> physicalPointITK;
    std::vector<double> physicalPoint;

    std::vector<unsigned int> size = vectorField->GetSize();
    float * buffer = vectorField->GetBufferAsFloat();
    std::vector<itk::int64_t> ijk {0,0,0};

    for(ijk[0]=0; ijk[0]<size[0]; ijk[0]=ijk[0]+1){
    for(ijk[1]=0; ijk[1]<size[1]; ijk[1]=ijk[1]+1){
    for(ijk[2]=0; ijk[2]<size[2]; ijk[2]=ijk[2]+1){

      physicalPoint = vectorField->TransformIndexToPhysicalPoint(ijk);
      for (i=0; i<3; i++){
        physicalPointITK[i] = physicalPoint[i];
      }

      imageLinearIndex = ijk[0] + (size[0] * (ijk[1] + size[1] * ijk[2]));

      for (landmarkIndex=0; landmarkIndex < numLandmarks; landmarkIndex++) {
        
        rbf = RBFValue(&fixedLandmarks->at(landmarkIndex), &physicalPointITK, adaptRadius[landmarkIndex]);

        for (i=0; i<3; i++){
          buffer[3*imageLinearIndex+i] += coeff[3*landmarkIndex+i] * rbf;
        }
      }

    }
    }
    }

  }

  float * BSplineRBFFindCoeffs(
    std::vector< itk::Point<double, 3> > * fixedLandmarks, 
    std::vector< itk::Point<double, 3> > * movingLandmarks,
    float * adaptRadius,
    float stiffness)
  {
    float rbfv1, rbfv2;
    int i, j, k, d;
    float rbf_prefactor, reg_term, r2, tmp;
    unsigned int numLandmarks = fixedLandmarks->size();

    float * coeff = (float*) malloc (3 * numLandmarks * sizeof(float));

    typedef vnl_matrix <double> Vnl_matrix;
    typedef vnl_svd <double> SVDSolverType;
    Vnl_matrix A, b;

    float RBFRadius = 0.0;
    for (i=0; i<numLandmarks; i++){
      RBFRadius += adaptRadius[i];
    }
    RBFRadius = RBFRadius / numLandmarks;

    A.set_size (3 * numLandmarks, 3 * numLandmarks);
    A.fill(0.);

    b.set_size (3 * numLandmarks, 1);
    b.fill (0.0);

    // right-hand side
    for (i=0; i<numLandmarks; i++) {
	  for (j=0; j<numLandmarks; j++) {
	    rbfv1 = RBFValue (&fixedLandmarks->at(i), &fixedLandmarks->at(j), adaptRadius[j]);
		
	    for (d=0;d<3;d++) {
		    b (3*i +d, 0) -= rbfv1 * (fixedLandmarks->at(j)[d] - movingLandmarks->at(j)[d]);
	    }
	  }
    }

    // matrix
    for (i = 0; i < numLandmarks; i++) {
	  for (j = 0; j < numLandmarks; j++) {
	    tmp = 0;
	    for (k = 0; k < numLandmarks; k++) {

		  rbfv1 = RBFValue (&fixedLandmarks->at(k), &fixedLandmarks->at(i), adaptRadius[k]);
		  rbfv2 = RBFValue (&fixedLandmarks->at(k), &fixedLandmarks->at(j), adaptRadius[k]);

		  tmp += rbfv1*rbfv2;
	    }
	    for (d=0;d<3;d++){
		    A(3*i+d, 3*j+d) = tmp;
	    }
	  }
    }

    //add regularization terms to the matrix
    rbf_prefactor = sqrt(M_PI/2.)*sqrt(M_PI/2.)*sqrt(M_PI/2.)/RBFRadius;
    for (d=0;d<3;d++) {
	  for (i=0;i<numLandmarks;i++) {
	    for (j=0;j<numLandmarks;j++) {
        tmp = A(3*i+d, 3*j+d);
        reg_term = 0.;			
        if (i==j) {
            reg_term = rbf_prefactor * 15.;
        }
        else
        {
          // r2 = sq distance between landmarks i,j in mm
          float d = fixedLandmarks->at(i).EuclideanDistanceTo(fixedLandmarks->at(j));
          r2 = (d * d) / (adaptRadius[i] * adaptRadius[j]);
          reg_term = rbf_prefactor * exp(-r2/2.) * (-10 + (r2-5.)*(r2-5.));
        }
        A (3*i+d,3*j+d) = tmp + reg_term * stiffness;
      } 
	  }
    }

    SVDSolverType svd (A, 1e-6);
    Vnl_matrix x = svd.solve (b);

    for (i=0; i<3*numLandmarks; i++) {
	    coeff[i] = x(i,0);
    }
    return coeff;
  }


} // end of anonymous namespace

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  // Checking conditions.

  if( fixedFiducials.size() <= 0 || movingFiducials.size() <= 0 ||
    fixedFiducials.size() != movingFiducials.size() )
  {
    std::cerr << "Fixed and moving landmark lists must be of the same size "
      << "and contain at least one point" << std::endl;
  }

  unsigned long numFiducials = fixedFiducials.size();

  if(outputDisplacementField.empty())
  {
    std::cerr << "An output transform must be specified" << std::endl;
  }

  unsigned long numRBFRadius = RBFRadius.size();

  if (numRBFRadius > 1 && numRBFRadius != numFiducials){
    std::cerr << "The number of RBF radius specified is more than one but does not match the number of fiducials.\n"
              << "Specify one RBF radius to use as a global parameter or multiple comma separated values corresponding to each fiducial.\n";
    return EXIT_FAILURE;
  }

  float * adaptRadius = (float *)malloc(numFiducials*sizeof(float));
  for (unsigned long i = 0; i<numFiducials; i++){
    if (numRBFRadius == 1){
      adaptRadius[i] = RBFRadius[0];
    }else{
      adaptRadius[i] = RBFRadius[i];
    }
  }

  typedef std::vector< itk::Point<double, 3> > PointList;

  PointList fixedPoints(fixedFiducials.size());
  PointList movingPoints(movingFiducials.size());

  // Convert both points lists to ITK points

  std::transform(fixedFiducials.begin(), fixedFiducials.end(),
    fixedPoints.begin(),
    convertStdVectorToITKPoint);

  std::transform(movingFiducials.begin(), movingFiducials.end(),
    movingPoints.begin(),
    convertStdVectorToITKPoint);

  // BSpline coeff

  float * coeff = BSplineRBFFindCoeffs(&fixedPoints, &movingPoints, adaptRadius, stiffness);

  // Create vector field
  sitk::ImageFileReader reader;
  reader.SetFileName(referenceVolume);
  sitk::Image &&referenceImage = reader.Execute();

  sitk::Image output(referenceImage.GetSize(), sitk::sitkVectorFloat32);
  output.SetOrigin(referenceImage.GetOrigin());
  output.SetSpacing(referenceImage.GetSpacing());
  output.SetDirection(referenceImage.GetDirection());

  RBFGaussUpdateVectorField(&output, coeff, &fixedPoints, adaptRadius);

  sitk::WriteImage(output, outputDisplacementField);

  return EXIT_SUCCESS;
}
