%% applyANTSTransformToPolyCoeffs - transform polynomial coefficient matrix in world coordinats of sourceRefImage to
% world coordinates in targetRefImage applying a ANTS affine transformation 
%
% Parameters: coeffs - Nx3 matrix, coeffs(end,:) is the bias term (x^0)
%             antsTransformFileStrings  -  '[antsTransformFile, useInverseFlag]' String with tuples of
%                  ANTS transformation file(s) names (binary
%                  .mat) strings and useInverseFlags, if multiple transformations should be concattenated
%                  pass a cell of transformation tuples. Note that the order is
%                  ANTs-like, i.e. stacked (last is applied first), and that
%                  ANTs specifies point transformations in the *inverse*
%                  direction than non-point transformations, thus useInverseFlag must be 1 for 
%                  a transform from A->B and 0 for a Transform A<-B.
%                  Refer to ANTs documentation for details.
%
% Returns: transCoeffs - Nx3 matrix
%
% Example: transCoeffs = applyANTSTransformToPolyCoeffs(coeffs,{'[trans_rigid_ct_post_to_t10GenericAffine.mat,1]','[trans_rigid_t1_to_mni0GenericAffine.mat,1]'})
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine
% 2017
% mail@andreashusch.de, husch.andreas@chl.lu

function transCoeffs = applyANTSTransformToPolyCoeffs(coeffs, antsTransformFileStrings)
LtoR = [-1 -1 1];
coeffs = coeffs * diag(LtoR); % LPI / RPI

trans = applyANTSTransformToPoints(coeffs,antsTransformFileStrings);
translation = applyANTSTransformToPoints([0 0 0],antsTransformFileStrings); % valid if center of rotation equals origin of space
transTranslation = applyANTSTransformToPoints(translation,antsTransformFileStrings); % valid if center of rotation equals origin of space

transBias = trans(end,:);% .* LtoR;
transCoef = trans(1:end-1,:);% .* LtoR;
% fix translation of coefs ( not to be translated)
transCoef = transCoef - repmat(translation,length(transCoef),1);

%FIXME document this
transCoeffs = [transCoef; transBias];
transCoeffs = transCoeffs * diag(LtoR);
end