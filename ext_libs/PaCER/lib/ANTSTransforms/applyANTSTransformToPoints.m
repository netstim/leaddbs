%% applyANTSTransformToPoints - transform points in world coordinats of a sourceRefImage to
% world coordinates in targetRefImage applying a ANTS affine transformation
%
% Parameters: transPoints - Nx3 matrix of points in world cooridnates
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
% Returns: Transformed Points
%
% Example: applyANTSTransformToPoints(rand(20,3), {'[trans_1_0GenericAffine.mat,1]','[trans_2_0GenericAffine.mat,1]'})
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine
% (c) 2017
% mail@andreashusch.de, husch.andreas@chl.lu
function transPoints = applyANTSTransformToPoints(points, antsTransformFileStrings)
try
[~, t] = system('echo $ANTSPATH');
antspath = t(1:end-1); % remove line ending
end
if(~isempty(t))
    applyTransformsToPointsCmd = [antspath 'antsApplyTransformsToPoints'];
else
    applyTransformsToPointsCmd = 'antsApplyTransformsToPoints'; % FIXME: Make multi-os compatabile and find bins automatically
end

tempFile = [tempname() '.csv'];
tempFileOut = [tempname() '.csv'];
fileHandle = fopen(tempFile, 'w');
fprintf(fileHandle, 'x,y,z,t,label,comment\n');
fclose(fileHandle);
dlmwrite(tempFile,  [points zeros(size(points,1),3)], '-append');

cmd = [applyTransformsToPointsCmd ' -d 3 -p 1 -i ' tempFile ];

if(iscell(antsTransformFileStrings))
    for i=1:length(antsTransformFileStrings)
        cmd = [cmd ' -t ' antsTransformFileStrings{i}]; %#ok<AGROW>
    end
else
    cmd = [cmd ' -t ' antsTransformFileStrings];
end

cmd = [cmd ' -o ' tempFileOut];

ret = ea_runcmd(cmd);

assert(ret == 0);
points = readmatrix(tempFileOut);
delete(tempFileOut);
delete(tempFile);

transPoints = points(:,1:3);
