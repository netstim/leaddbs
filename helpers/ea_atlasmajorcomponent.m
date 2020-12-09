function atlases=ea_atlasmajorcomponent(atlases,toclean)
% from specified component, only preserve largest component

for c=toclean
   for side=1:2
       % figure, patch(atlases.roi{c,side}.fv,'FaceColor','none')

       fvsplit=splitFV(atlases.roi{c,side}.fv.faces,atlases.roi{c,side}.fv.vertices);
       for comp=1:length(fvsplit)
           lengths(comp)=length(fvsplit(comp).faces);
       end
       [~,ix]=max(lengths);
       atlases.roi{c,side}.fv.faces=fvsplit(ix).faces;
       atlases.roi{c,side}.fv.vertices=fvsplit(ix).vertices;
   end
end



function fvOut = splitFV( f, v )
%SPLITFV Splits faces and vertices into connected pieces
%   FVOUT = SPLITFV(F,V) separates disconnected pieces inside a patch defined by faces (F) and
%   vertices (V). FVOUT is a structure array with fields "faces" and "vertices". Each element of
%   this array indicates a separately connected patch.
%
%   FVOUT = SPLITFV(FV) takes in FV as a structure with fields "faces" and "vertices"
%
%   For example:
%     fullpatch.vertices = [2 4; 2 8; 8 4; 8 0; 0 4; 2 6; 2 2; 4 2; 4 0; 5 2; 5 0];
%     fullpatch.faces = [1 2 3; 1 3 4; 5 6 1; 7 8 9; 11 10 4];
%     figure, subplot(2,1,1), patch(fullpatch,'facecolor','r'), title('Unsplit mesh');
%     splitpatch = splitFV(fullpatch);
%     colours = lines(length(splitpatch));
%     subplot(2,1,2), hold on, title('Split mesh');
%     for i=1:length(splitpatch)
%         patch(splitpatch(i),'facecolor',colours(i,:));
%     end
%
%   Note: faces and vertices should be defined such that faces sharing a coincident vertex reference
%   the same vertex number, rather than having a separate vertice defined for each face (yet at the
%   same vertex location). In other words, running the following command: size(unique(v,'rows') ==
%   size(v) should return TRUE. An explicit test for this has not been included in this function so
%   as to allow for the deliberate splitting of a mesh at a given location by simply duplicating
%   those vertices.
%
%   See also PATCH

%   Copyright Sven Holcombe
%   $Date: 2010/05/19 $

%% Extract f and v
if nargin==1 && isstruct(f) && all(isfield(f,{'faces','vertices'}))
    v = f.vertices;
    f = f.faces;
elseif nargin==2
    % f and v are already defined
else
    error('splitFV:badArgs','splitFV takes a faces/vertices structure, or these fields passed individually')
end

%% Organise faces into connected fSets that share nodes
fSets = zeros(size(f,1),1,'uint32');
currentSet = 0;

while any(fSets==0)
    currentSet = currentSet + 1;
    fprintf('Connecting set #%d vertices...',currentSet);
    nextAvailFace = find(fSets==0,1,'first');
    openVertices = f(nextAvailFace,:);
    while ~isempty(openVertices)
        availFaceInds = find(fSets==0);
        [availFaceSub, ~] = find(ismember(f(availFaceInds,:), openVertices));
        fSets(availFaceInds(availFaceSub)) = currentSet;
        openVertices = f(availFaceInds(availFaceSub),:);
    end
    fprintf(' done! Set #%d has %d faces.\n',currentSet,nnz(fSets==currentSet));
end
numSets = currentSet;

%% Create separate faces/vertices structures for each fSet
fvOut = repmat(struct('faces',[],'vertices',[]),numSets,1);

for currentSet = 1:numSets
    setF = f(fSets==currentSet,:);
    [unqVertIds, ~, newVertIndices] = unique(setF);
    fvOut(currentSet).faces = reshape(newVertIndices,size(setF));
    fvOut(currentSet).vertices = v(unqVertIds,:);
end
