%% NiftiSeg - Class representing a Nifti Segmentation File
%  No frills full volume variant of the NifitSegmentationSubvolume class by Andreas Husch / Florian
%  Bernard
%
%  Andreas Husch 2017
classdef NiftiSeg < NiftiMod & matlab.mixin.Copyable  & configurable & plotable3D
    properties (SetAccess = protected, GetAccess=public, SetObservable = true) %allow access via specific methods
        % Inherited:
        %         filepath = '';  % full qualified path + filename (normaly a nifti file)
        %         data = [];      % raw nifti data isLoaded in ram as matlab matrix TODO:
        %                           maybe make this dependent (can be calculated from
        %                           binaryCroppedVolume and offset!)
        %         header;            % the nifti header
        %         voxdim = NaN(3,1); % voxel dimensions, x,y,z
        %         isLoaded = false;    % true = in memory, false = on disk
        labelNames = {};                    % names of all existing labels
        presentLabelIntensities = []; % index Numbers of the labels present in the loaded segmentations
        presentLabelNames = []; % names of the labels present in the loaded segmentation
        presentAndRelevantLabelNames = []; % relevant label names that are present in this file
        LABEL_FILE_COLORS = [];
        perObjectfilters = {};
        segmentationBinaryIndices = {};
        surfaceData = [];   %variable to cache surface data for electrode intersection calculation
        logicalSurfaceData = [];
        plotHandle3D = {};
        labelFilePath = '';                %; 'c:/SVNRepos/CNM/resources/Region_Labels.lbl'
        plotClearLabel = []; % init in arg parser
        faceVertexStruct = {};     
    end
    
    properties (SetAccess = public, GetAccess=public, SetObservable) % set
        relevantLabelNames = {}; % label names configured to be relevant   %e.g.{'SNr+STN_R'}
        COLORS = []; %array of colors to plot the contained structures
        PLOT_TYPE_3D = 1; % type of 3D plot 1==isosurface, 2==voxel cubes
        FACE_ALPHA_3D = 0.6; % alpha of 3D plots
        SHOW_EDGES_3D = false; % show edges in 3D plots?
        SMOOTH_PLOT = true;
    end
    
    properties (Access=public, SetObservable = true)
        selectedStructureNameOrIdx = 1;  % Index of the currently selected structures. All methods of this class are applied to this structure if no labelNameOrIdx is given %TODO set via dependend method and to sanity check
        structuresToPlot = {};     % initializied with structuresToPlot = relevantLabelNAme => plot ALL structures      
    end
    
    methods
        function this = NiftiSeg(filepath, varargin)% filename, minSubvolumeOffset, maxSubvolumeOffset
            
            this = this@NiftiMod(filepath, varargin{:});
            
            this.argParser.addOptional('smoothPlot',  true);
            this.argParser.addOptional('showAllLabels',  true); % ignore relevant labels but show everything thats in the file
            this.argParser.addOptional('relevantLabelNames', {'SNr+STN_L', 'SNr+STN_R', 'NR_L', 'NR_R'...
			'Th_L', 'Th_R', 'Acc_L', 'Acc_R', 'Put+GP_L', 'Put+GP_R', 'STN_L', 'STN_R', 'SN_L', 'SN_R'});
            this.argParser.addOptional('labelFilePath', '');%Default_Region_Labels.label		
            this.argParser.addOptional('allowMissingLabels', true);
            this.argParser.addOptional('plotClearLabel', false);
            this.argParser.addOptional('preCompSurfacesFilePath','');
            this.argParser.addOptional('preComputeSurfaces', true);

            this.argParser.parse(varargin{:});
            Args = this.argParser.Results;
            
            this.relevantLabelNames = Args.relevantLabelNames;
            this.labelFilePath = Args.labelFilePath;
            this.plotClearLabel = Args.plotClearLabel;
            
            
            % select all contained and relevant structures to be plotted initially
            this.presentLabelIntensities = uint32(unique(this.img));
            
            if(~isempty(this.labelFilePath))
                [this.labelNames, allLblColors, labelFileIntensities] = readLabels(this.labelFilePath); %evtl im objekt "hardcoden" (in methode damit sch�n versoniert und sicher funktionierend)      
            else
                allLblColors = parula(length(this.presentLabelIntensities));
                labelFileIntensities = 1:length(this.presentLabelIntensities); 
                this.labelNames{1} = 'Clear Label';
                for i = 2:length(this.presentLabelIntensities)
                    this.labelNames{i}= ['No name given ' num2str(i-1) ];
                end
            end
            getLabelFileIdxForIntensity = @(intensity)(find(labelFileIntensities==intensity));
            getImageIdxForIntensity = @(intensity)(find(this.presentLabelIntensities==intensity));

            presImgAndLblFileIntensities  = intersect(this.presentLabelIntensities, labelFileIntensities);
            
            if(length(this.presentLabelIntensities) ~= length(labelFileIntensities)) 
                warning('NiftiSeg: No of. label intensities in image file does not match no. of labels specified in label file! Names might be wrong. Did you specifiy the correct label file for the given data?');
            end
            this.presentLabelNames = this.labelNames(this.presentLabelIntensities + 1);
            
            if(Args.showAllLabels)
                this.presentAndRelevantLabelNames  = this.presentLabelNames;
            else
                this.presentAndRelevantLabelNames = intersect(this.relevantLabelNames, ...
                    this.presentLabelNames, 'stable'); % stable is important here, as without this argument the order is changed and some code relies on the order of labels
                missingLabels = setdiff(this.relevantLabelNames, this.presentAndRelevantLabelNames);
                if (  numel(missingLabels) > 0 )
                    warning(['For the labels ' strjoin(missingLabels, ',') ' there is no segmentation!'])
                end
            end
            
            this.structuresToPlot = this.presentAndRelevantLabelNames;
            this.LABEL_FILE_COLORS = allLblColors;
            this.COLORS = parula(length(this.presentLabelIntensities));
            for i = 1: length(presImgAndLblFileIntensities)
                this.COLORS(getImageIdxForIntensity(presImgAndLblFileIntensities(i)),:) = ...
                    this.LABEL_FILE_COLORS(getLabelFileIdxForIntensity(presImgAndLblFileIntensities(i)),:);
            end
            this.SMOOTH_PLOT = Args.smoothPlot;
            
            for i=1:numel(this.presentAndRelevantLabelNames)
                labelIntensity = ...
                    getIntensityForLabel(this.labelNames, this.presentAndRelevantLabelNames{i});
                for j = 1:length(labelIntensity) % if(length(labelIntensity) > 2)  => label defined multiple times (e.g. unnamed structure)
                    % (binary) preprocessing of subvolume
                    
                    binVol = (this.img==labelIntensity(j));
                    if ( isempty(binVol) )
                        warning(['No segmentation for object ' ...
                            this.relevantLabelNames{i} ' available in file ' filepath]);
                    end
                    
                    this.segmentationBinaryIndices{i} = find(binVol);	                
                end  
            end
            
            %% Precompute FVs
            if(~exist(Args.preCompSurfacesFilePath,2) && Args.preComputeSurfaces)
                
                disp('Precomputing Surfaces...')
                
                for i=1:length(this.presentLabelNames)
                    this.faceVertexStruct{i} = this.computeIsoSurfaceFVWorld(this.presentLabelNames{i});
                end
            elseif(Args.preComputeSurfaces)
                disp('Loading precomputed surfaces from file...')
                this.faceVertexStruct = load(Args.preCompSurfacesFilePath);
            end
            
            %% Add listeners to model variables that should update the plot
            addlistener(this, 'COLORS', 'PostSet', @this.updatePlot3D);
            addlistener(this, 'PLOT_TYPE_3D', 'PostSet', @this.updatePlot3D);
            addlistener(this, 'SHOW_EDGES_3D', 'PostSet', @this.updatePlot3D);
            addlistener(this, 'FACE_ALPHA_3D', 'PostSet', @this.updatePlot3D); 
        end
        
        function set.COLORS(this, newColors)
            if(isequal([length(this.presentLabelIntensities) 3], size(newColors)))
                this.COLORS = newColors;
            else
                error('setCOLORS:invalidColorVector', 'Given Color Vector is invalid!');
            end
        end
                
        function pointCloud = pointCloud(this, labelNameOrIdx) % pointClound in intrinsic Coordinates 
            if ( ~isnumeric(labelNameOrIdx) )
                segBinInd=find(ismember(this.relevantLabelNames,labelNameOrIdx)==1);
            else
                segBinInd = labelNameOrIdx;
            end
            [i,j,k] = ind2sub(size(this.img),this.segmentationBinaryIndices{segBinInd});
            pointCloud = [i j k];
        end
           
        function pointCloudWorld = pointCloudWorld(this, labelNameOrIdx) % pointClounds in nifti world coordinate system
            if(nargin < 2)
                labelNameOrIdx = this.selectedStructureNameOrIdx;
            end
            voxelList = this.pointCloud(labelNameOrIdx);
            voxelList(:,[1 2 3]) = voxelList(:,[2 1 3]); % Swap i,j to X,Y 
            
            pointCloudWorld = this.getNiftiWorldCoordinatesFromMatlabIdx(voxelList')';
        end     
        
        function subvolumeCell=getBinaryVolumeCell(this, labelNamesCell)
            if ( ~exist('labelNamesCell', 'var') )
                labelNamesCell = this.presentAndRelevantLabelNames;
            end
            if ( ~iscell(labelNamesCell) ) % single argument is given as string -> convert to cell
                labelNamesCell = {labelNamesCell};
            end

            subvolumeCell = cell(1,numel(labelNamesCell));
            for i=1:numel(labelNamesCell)
                currVolume = false(this.voxdim);

                segBinInd=find(ismember(this.relevantLabelNames,labelNamesCell{i})==1);
                currVolume(this.segmentationBinaryIndices{segBinInd}) = true;

                subvolumeCell{i} = currVolume;
            end
         end
         
        function volume=getBinaryVolume(this, labelNameOrIdx)
            segBinInd = [];
            volume = false(this.voxdim');
            
            if(~exist('labelNameOrIdx', 'var')) % HA: return binary volume of selectedStructureNameOrIdx of the object
               labelNameOrIdx = this.selectedStructureNameOrIdx;
            end
            
            if ( ~iscell(labelNameOrIdx) )
                labelNameOrIdx = {labelNameOrIdx};
            end
%           labelNames = readLabels(Const.labelFile);
            
            for i=1:numel(labelNameOrIdx)
                if ( ~isnumeric(labelNameOrIdx{i}) ) % label name given
                    if(sum(ismember(this.presentAndRelevantLabelNames,labelNameOrIdx{i})) == 0)
                       disp(['Given labelname ' labelNameOrIdx{i} ' is not contained in relevant label names!']); 
                    end
                    segBinInd=union(segBinInd, ...
                        find(ismember(this.presentAndRelevantLabelNames,labelNameOrIdx{i})==1));
                else % index given
                    segBinInd = union(segBinInd, labelNameOrIdx{i});
                end
                volume(this.segmentationBinaryIndices{segBinInd(i)}) = true;
            end  
        end
               
        function faceVertexStruct = getIsoSurfaceFVWorld(this, labelNameOrIdx) %
            if(nargin < 2)
                labelNameOrIdx = this.selectedStructureNameOrIdx;
            end
            if(~isempty(this.faceVertexStruct)) %cached
                if(~isnumeric(labelNameOrIdx))
                    idx = this.getLabelIdxByName(labelNameOrIdx);
                else
                    idx = labelNameOrIdx;
                end
                faceVertexStruct = this.faceVertexStruct{idx};
            else
                faceVertexStruct = this.computeIsoSurfaceFVWorld(this, labelNameOrIdx);
            end
        end
        
        function faceVertexStruct = computeIsoSurfaceFVWorld(this, labelNameOrIdx) %
            if(nargin < 2)
                labelNameOrIdx = this.selectedStructureNameOrIdx;
            end
            if (this.SMOOTH_PLOT)
                faceVertexStruct = isosurface(smooth3(this.getBinaryVolume(labelNameOrIdx), 'gaussian', 3), 0.5);
            else
                faceVertexStruct = isosurface(this.getBinaryVolume(labelNameOrIdx), 0.5);
            end
            voxelList(:,[1 2 3]) = faceVertexStruct.vertices(:,[2 1 3]); % Swap i,j to X,Y
            
            faceVertexStruct.vertices = this.getNiftiWorldCoordinatesFromMatlabIdx(voxelList')';
        end
        
        function str = toString(this)
            [~, str] = fileparts(this.filepath);
        end
         
        function graphicsHandle = initPlot3D(this, parentAxes) %TODO: Config Panel 
            graphicsHandle = hggroup('Parent', parentAxes);
            
            for i=1:length(this.structuresToPlot)
                labelIntensity =  getIntensityForLabel(this.labelNames, this.structuresToPlot{i});
                labelIdx = this.getLabelIdxByName(this.structuresToPlot{i});
                
                if(strcmp(this.structuresToPlot{i}, 'Clear Label') || (labelIntensity == 0)  && ~this.plotClearLabel)
                    continue;
                end
                if(~isempty(labelIdx)) %make sure structure exists                
                    translatedIsosurface = this.getIsoSurfaceFVWorld(this.structuresToPlot{i});             
                
                    % todo enable changing the number of faces (configpanel)
                    %                     nFaces = 200;
                    %                     reducedIsosurface = reducepatch(translatedIsosurface, nFaces);
                    reducedIsosurface = translatedIsosurface;
                    this.plotHandle3D {end + 1} = patch(reducedIsosurface,...
                        'Parent', graphicsHandle); %offset is already mm!
                    if(this.SHOW_EDGES_3D)
                        edgeColor = 'k';
                    else
                        edgeColor = 'none';
                    end
                    set(this.plotHandle3D{end} , 'FaceColor', this.COLORS(labelIdx,:), ...
                        'EdgeColor', edgeColor, 'FaceAlpha', this.FACE_ALPHA_3D, 'Parent', graphicsHandle, ...
                        'Tag', this.structuresToPlot{i});
                end
            end
            graphicsHandle.Tag = 'NiftiSeg';
            daspect([1 1 1]);
            lighting gouraud;
            % add light if none is present
            if(isempty(findobj(parentAxes, 'Type', 'light')))
                camlight headlight;
            end
            material(graphicsHandle, 'dull');
            this.plotHandleHashMap3D(double(parentAxes)) = double(graphicsHandle);
        end
            

        function createVolumeFromVertFaces(this,~,~)
            binVol = zeros(this.voxdim');
                        for i=1:numel(this.structuresToPlot)
                            VertFaceStruct = findobj(this.plotHandle3D{i});
                            fv.faces = get(VertFaceStruct,'Faces');
                            fv.vertices = get(VertFaceStruct,'Vertices');
                            
                            binVol = binVol | isosurface2volumeFast(fv);
                        end
                        binVol = permute(binVol,[2 1 3]);
                        this.volumeData = binVol;                       
        end
           
        function idx = getLabelIdxByName(this, labelName)
            idx=find(ismember(this.presentLabelNames,labelName)==1);
        end
        
        function setLabelColorByName(this, labelName, rgbTriple)
            if(iscell(labelName)) % set of names given, set same color
                for i = 1:length(labelName)
                    idx = ismember(this.presentLabelNames,labelName{i})==1;
                    this.COLORS(idx,:) = rgbTriple;
                end
            else % single label name given
                idx = ismember(this.presentLabelNames,labelName)==1;
                this.COLORS(idx,:) = rgbTriple;
            end
        end
        
        function setLabelColorByIdx(this, labelIdx, rgbTriple)
            if(iscell(labelIdx)) % set of idxes given, set same color
                for i = 1:length(labelIdx)
                    idx = labelIdx{i};
                    this.COLORS(idx,:) = rgbTriple;
                end
            else % single label idx given
                this.COLORS(labelIdx,:) = rgbTriple;
            end
        end
    end

    
end