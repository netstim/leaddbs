classdef LeG_vox2atlaslabel < handle
    
    properties
        NMMDir; %default atlas (required)
        AtlasDir; AtlasNames; 
        AtlasLabels; AtlasVals; 
        AtlasData; AtlasMat; AtlasCnt; AtlasImg;
        DefFile; DefMat; Def; %ptvox2mnimm deformation (iy_MR.nii)
        VoxDef; %ptvox2atlasvox deformation (iy_MR.nii transformed using inv(atlasmat) -> mm2vox) 
        ProbRadius; %radius of sphere in mm for probability calculation
    end
    
    methods
        
        function obj = LeG_vox2atlaslabel(varargin)
            %Initialize as follows. Include NMMDir, AtlasDir (contains
            %other atlases of interest), and DefFile.
            %
            % AObj = LeG_vox2atlaslabel('NMMDir','./nmm','AtlasDir','./atlases','DefFile','./iy_MR.nii');
            
            obj.parseInput(varargin{:});
            if isempty(obj.NMMDir)
                error('NMMDir must be provided to initialize class!')          
            end                                              
            
            NMMName = dir(fullfile(obj.NMMDir,'*.nii'));
            AtlasNames = dir(fullfile(obj.AtlasDir,'*.nii'));
            obj.AtlasNames = [{NMMName.name};setdiff({AtlasNames.name}',NMMName.name)]; %make nmm 1st and remove any duplicates
            obj.AtlasCnt = length(obj.AtlasNames);
            obj.AtlasData = cell(obj.AtlasCnt,1);
            obj.AtlasImg = cell(obj.AtlasCnt,1); %atlas labeled image warped to patient space
            obj.AtlasMat = cell(obj.AtlasCnt,1);
            obj.AtlasVals = cell(obj.AtlasCnt,1);
            obj.AtlasLabels = cell(obj.AtlasCnt,1);
            obj.VoxDef = cell(obj.AtlasCnt,1); %vox2vox deformation (populated in findLabel)
            for k=1:length(obj.AtlasNames)
                if k==1 %NMM
                    AtlasObj = nifti(fullfile(obj.NMMDir,obj.AtlasNames{k}));
                    AtlasFID = fopen(fullfile(obj.NMMDir,regexprep(obj.AtlasNames{k},'\.nii$','.txt')));
                else %all other atlases
                    AtlasObj = nifti(fullfile(obj.AtlasDir,obj.AtlasNames{k}));
                    AtlasFID = fopen(fullfile(obj.AtlasDir,regexprep(obj.AtlasNames{k},'\.nii$','.txt')));
                end
                obj.AtlasData(k) = {double(AtlasObj.dat)};
                obj.AtlasMat(k) = {double(AtlasObj.mat)};
                AtlasTXT = fread(AtlasFID,[1,Inf],'*char');
                AtlasCell = regexp(AtlasTXT,'\r\n','split');
                AtlasCell(cellfun(@isempty,AtlasCell)) = [];
                AtlasCell = regexp(AtlasCell,'\t','split'); AtlasCell = reshape([AtlasCell{:}],length(AtlasCell{1}),[])';
                obj.AtlasVals(k) = {cellfun(@str2double,AtlasCell(:,1))};
                AtlasLabels = regexp(AtlasCell(:,2:end),'\.','split');
                AtlasLabelsExp = {};
                b = all(cellfun(@length,AtlasLabels)>1);
                for m=1:length(b)
                    if b(m)
                        atlaslabels = AtlasLabels(:,m);
                        AtlasLabelsExp = [AtlasLabelsExp,reshape([atlaslabels{:}],length(AtlasLabels{1,m}),[])'];
                    else
                        AtlasLabelsExp = [AtlasLabelsExp,cellfun(@cell2mat,AtlasLabels(:,m),'UniformOutput',false)];
                    end
                end
                obj.AtlasLabels(k) = {regexprep(regexprep(AtlasLabelsExp,'*','Unknown'),'''','')};
                fclose(AtlasFID);
            end
        end %end constructor
        
        function parseInput(obj, varargin)
            %Parses name-value input pairs and saves the value to the
            %matching class property.
            PropNames = properties(obj); 
            InputNames = varargin(1:2:end);
            InputVals = varargin(2:2:end);           
            for k=1:length(InputNames)
                if any(strcmp(InputNames{k},PropNames))
                    if isa(InputVals{k},'string')
                        obj.(InputNames{k}) = char(InputVals{k});
                    else
                        obj.(InputNames{k}) = InputVals{k};
                    end
                end
            end
        end
        
        function loadDefFile(obj, varargin)
            obj.parseInput(varargin{:});
            if isempty(obj.DefFile)
                error('Path to deformation file (DefFile) must be provided!')
            end
            
            Nii = nifti(obj.DefFile);
            obj.Def = single(Nii.dat(:,:,:,1,:));
            d = size(obj.Def);
            if d(4)~=1 || d(5)~=3, error('Deformation field is wrong!'); end
            obj.Def = reshape(obj.Def,[d(1:3) d(5)]); %imgvox2mnimm image deformation
            obj.DefMat = Nii.mat;
        end
        
        function [Label,XYZMNImm] = findLabel(obj, XYZIdx, AtlasName, varargin)
            %[Label,XYZMNImm] = findLabel([1,1,1],'NMM');
            %XYZIdx are voxel indices (not mm) in original image (not mni)
            %associated with deformation file (DefFile)
            if isempty(obj.Def)
                obj.loadDefFile(varargin{:});
            end
            
            Idx = strcmp(regexprep(obj.AtlasNames,'\.nii$',''),AtlasName);            
            if isempty(obj.VoxDef{Idx})
                mm2vox_atlas = inv(obj.AtlasMat{Idx}); %atlas transform
                obj.VoxDef(Idx) = {affine(obj.Def,mm2vox_atlas)}; %vox2vox transform (patient to mni)
            end
                        
            %%%%%%%%% Faster linear index version %%%%%%%%%%%%%%%%
            sz_vox = size(obj.VoxDef{Idx});
            n_vox = prod(sz_vox(1:3));
            
            lidx = calcLinIdx(sz_vox(1:3),XYZIdx);
            
            %MNI mm locations for each voxel in pt space (XYZIdx)
            MNImm = nan(length(lidx),3);
            MNImm(:,1) = obj.Def(lidx);
            MNImm(:,2) = obj.Def(lidx+n_vox);
            MNImm(:,3) = obj.Def(lidx+2*n_vox);
            XYZMNImm = nanmean(MNImm);
                        
            %Atlas voxel indices for each voxel in pt space (XYZIdx)
            atlasvox = nan(length(lidx),3);
            atlasvox(:,1) = round((obj.VoxDef{Idx}(lidx)));
            atlasvox(:,2) = round((obj.VoxDef{Idx}(lidx+n_vox)));
            atlasvox(:,3) = round((obj.VoxDef{Idx}(lidx+2*n_vox)));
            
            sz_atlas = size(obj.AtlasData{Idx});
            
            lidx_atlas = calcLinIdx(sz_atlas(1:3),atlasvox);
            nan_idx = isnan(lidx_atlas);
            
            %Atlas label indices (corresponds to atlas labels in text file)
            AtlasVal = zeros(length(lidx_atlas),1);
            AtlasVal(~nan_idx) = round(obj.AtlasData{Idx}(lidx_atlas(~nan_idx))); %label indices should be integers, but rounding just in case
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Label = getLabels(AtlasVal, obj.AtlasVals{Idx}, obj.AtlasLabels{Idx});
        end %end findLabel
        
        function LabelProbs = findLabelProbabilities(obj, COMIdx, AtlasName, varargin)
            %Uses voxel COM (COMIdx) to create a sphere of specified radius
            %for finding label probabilities for specified atlas
            %(AtlasName). COMIdx are voxel indices in patient space.
                        
            if isempty(obj.Def)
                obj.loadDefFile(varargin{:});
            else
                obj.parseInput(varargin{:});
            end

            RadiusMM = 10; %listening radius (mm) for finding labels
            if ~isempty(obj.ProbRadius)
                RadiusMM = obj.ProbRadius;
            end
            
            Idx = strcmp(regexprep(obj.AtlasNames,'\.nii$',''),AtlasName);            
            if isempty(obj.VoxDef{Idx})
                mm2vox_atlas = inv(obj.AtlasMat{Idx}); %atlas transform
                obj.VoxDef(Idx) = {affine(obj.Def,mm2vox_atlas)}; %vox2vox transform (patient to mni)
            end
            
            voxdef = obj.VoxDef{Idx};
            atlasdata = obj.AtlasData{Idx};
            atlasvals = obj.AtlasVals{Idx};
            atlaslabels = obj.AtlasLabels{Idx};
            
            XYZScale = sqrt(sum(obj.DefMat(1:3,1:3).^2)); %pt mm/vox scaling
            RadiusVox = ceil(RadiusMM./XYZScale);
            
            LabelProbs = cell(size(COMIdx,1),1);
            for m=1:size(COMIdx,1)
                
                % full index values centered on click point (XYZ)
                xidx = COMIdx(m,1)-RadiusVox(1):COMIdx(m,1)+RadiusVox(1);
                yidx = COMIdx(m,2)-RadiusVox(2):COMIdx(m,2)+RadiusVox(2);
                zidx = COMIdx(m,3)-RadiusVox(3):COMIdx(m,3)+RadiusVox(3);
                [YIdx,XIdx,ZIdx] = meshgrid(yidx,xidx,zidx);
                
                com = RadiusVox+1;
                
                [my,mx,mz] = meshgrid(1:length(yidx),1:length(xidx),1:length(zidx));
                mx = (mx - com(1)).*XYZScale(1);
                my = (my - com(2)).*XYZScale(2);
                mz = (mz - com(3)).*XYZScale(3);
                
                elecimg = (sqrt(mx.^2+my.^2+mz.^2)<=RadiusMM); %1s are sphere around center, 0s everywhere else
                
                fullidx = [XIdx(elecimg),YIdx(elecimg),ZIdx(elecimg)]; %full index of electrode sphere
                
                %%%%%%%%% Faster linear index version %%%%%%%%%%%%%%%%
                sz_vox = size(voxdef);
                n_vox = prod(sz_vox(1:3));
                
                lidx = calcLinIdx(sz_vox(1:3),fullidx);
                
                atlasvox = nan(length(lidx),3);
                atlasvox(:,1) = round((voxdef(lidx)));
                atlasvox(:,2) = round((voxdef(lidx+n_vox)));
                atlasvox(:,3) = round((voxdef(lidx+2*n_vox)));
                
                sz_atlas = size(atlasdata);
                
                lidx_atlas = calcLinIdx(sz_atlas(1:3),atlasvox);
                nan_idx = isnan(lidx_atlas);    
                
                AtlasVal = zeros(length(lidx_atlas),1);
                AtlasVal(~nan_idx) = round(atlasdata(lidx_atlas(~nan_idx)));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                [~,labelprobs] = getLabels(AtlasVal, atlasvals, atlaslabels);
                LabelProbs(m) = {labelprobs};
            end %end for
            
        end %end function

        function LabelProbs = findLabelProbabilities2(obj, COMIdx, AtlasName, varargin)
            %Uses voxel COM (COMIdx) to create a sphere of specified radius
            %for finding label probabilities for specified atlas
            %(AtlasName). COMIdx are voxel indices in patient space. This
            %version (#2) uses the atlas label image warped to pt space
            %(AtlasImg) to calc the probablilities. Slightly faster than
            %version #1.
            
            RadiusMM = 10; %search radius (mm) for finding labels
            
            Idx = strcmp(regexprep(obj.AtlasNames,'\.nii$',''),AtlasName); 
            if isempty(obj.AtlasImg{Idx})
                obj.warpAtlas2PatientSpace(AtlasName,varargin{:});
            end
            
            voxdef = obj.VoxDef{Idx};
            atlasvals = obj.AtlasVals{Idx};
            atlaslabels = obj.AtlasLabels{Idx};
            atlasimg = obj.AtlasImg{Idx};
            
            XYZScale = sqrt(sum(obj.DefMat(1:3,1:3).^2)); %pt mm/vox scaling
            RadiusVox = ceil(RadiusMM./XYZScale);
            
            LabelProbs = cell(size(COMIdx,1),1);
            for m=1:size(COMIdx,1)
                
                % full index values centered on click point (XYZ)
                xidx = COMIdx(m,1)-RadiusVox(1):COMIdx(m,1)+RadiusVox(1);
                yidx = COMIdx(m,2)-RadiusVox(2):COMIdx(m,2)+RadiusVox(2);
                zidx = COMIdx(m,3)-RadiusVox(3):COMIdx(m,3)+RadiusVox(3);
                [YIdx,XIdx,ZIdx] = meshgrid(yidx,xidx,zidx);
                
                com = RadiusVox+1;
                
                [my,mx,mz] = meshgrid(1:length(yidx),1:length(xidx),1:length(zidx));
                mx = (mx - com(1)).*XYZScale(1);
                my = (my - com(2)).*XYZScale(2);
                mz = (mz - com(3)).*XYZScale(3);
                
                elecimg = (sqrt(mx.^2+my.^2+mz.^2)<=RadiusMM); %1s are sphere around center, 0s everywhere else
                fullidx = [XIdx(elecimg),YIdx(elecimg),ZIdx(elecimg)]; %full index of electrode sphere
                
                sz_vox = size(voxdef);
                lidx = calcLinIdx(sz_vox(1:3),fullidx);
                AtlasVal = round(atlasimg(lidx));
                
                [~,labelprobs] = getLabels(AtlasVal, atlasvals, atlaslabels);
                LabelProbs(m) = {labelprobs};
            end %end for
            
        end %end function

        function varargout = warpAtlas2PatientSpace(obj, AtlasName, varargin)
            if isempty(obj.Def)
                obj.loadDefFile(varargin{:});
            end
            
            Idx = strcmp(regexprep(obj.AtlasNames,'\.nii$',''),AtlasName);
            if isempty(obj.VoxDef{Idx})
                mm2vox_atlas = inv(obj.AtlasMat{Idx}); %atlas transform
                obj.VoxDef(Idx) = {affine(obj.Def,mm2vox_atlas)}; %vox2vox transform (patient to mni)
            end
            
            atlasfile = fullfile(fileparts(obj.DefFile),['lw',AtlasName,'.nii']);
            if isempty(obj.AtlasImg{Idx})
                if isfile(atlasfile)
                    obj.AtlasImg(Idx) = {spm_read_vols(spm_vol(atlasfile))};
                    if nargout
                        varargout(1) = obj.AtlasImg(Idx);
                    end
                    return;
                end
            else
                if nargout
                    varargout(1) = obj.AtlasImg(Idx);
                end
                return;
            end
                                                
            sz_vox = size(obj.VoxDef{Idx},1:3);
            n_vox = prod(sz_vox);
            sz_atlas = size(obj.AtlasData{Idx});

            %vox2vox transform (rows = linear index into pt image, cols = vox indices into atlas)
            %nans exist in this transform and values could be out of range
            %wrt to atlas size
            atlasvox = reshape(round(double(obj.VoxDef{Idx})),[],3); %these are voxels so need to round 
            
            lidx_vox = (1:n_vox);            
            lidx_atlas = calcLinIdx(sz_atlas(1:3),atlasvox);
            
            nanidx = isnan(lidx_atlas);
            lidx_atlas(nanidx) = [];
            lidx_vox(nanidx) = [];
            
            atlasimg = zeros(sz_vox(1:3));
            atlasimg(lidx_vox) = round(obj.AtlasData{Idx}(lidx_atlas));
                      
            %Saving
            %Pinfo scale/offset factors should be 1 and 0 to avoid rounding
            %errors with atlas labeled data. If pinfo is not provided,
            %scale/offset will be calculated based on the data values. This
            %results in the data being scaled when saving and then
            %converted back to integers when loading, which adds a rounding
            %error (i.e. decimal) to each integer.
            info = spm_vol(obj.DefFile);
            info = rmfield(info,'private');
%             info = rmfield(info,'pinfo');
            info.pinfo = [1;0;352]; %pinfo(1:2) = scale/offset factors and pinfo(3) = 352 for NIfTI-1 or 544 for NIfTI-2 (see spm_create_vol.m)
            info.dt = [4,0]; %[int16 datatype, little-endian] (16 is float32, 1 is big-endian)
            info.descrip = 'atlas label deformation to patient space';
            info.fname = fullfile(fileparts(obj.DefFile),['lw',AtlasName,'.nii']);
            spm_write_vol(info,int16(atlasimg)); %write back to nii (with "lw" prefix)  
            
            obj.AtlasImg(Idx) = {atlasimg};
            if nargout
                varargout(1) = obj.AtlasImg(Idx);
            end
        end %end function

    end %end methods
    
end %end classdef

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Def = affine(y,M)
Def          = zeros(size(y),'single');
Def(:,:,:,1) = y(:,:,:,1)*M(1,1) + y(:,:,:,2)*M(1,2) + y(:,:,:,3)*M(1,3) + M(1,4);
Def(:,:,:,2) = y(:,:,:,1)*M(2,1) + y(:,:,:,2)*M(2,2) + y(:,:,:,3)*M(2,3) + M(2,4);
Def(:,:,:,3) = y(:,:,:,1)*M(3,1) + y(:,:,:,2)*M(3,2) + y(:,:,:,3)*M(3,3) + M(3,4);
end

function lidx = calcLinIdx(sz,idx)
%Calculates the linear index given the row/col/page indices (idx) and size
%of the matrix (sz). Only works for 3 dimensions. If values are out of
%range, they are clipped to fit sz.
k = cumprod(sz);

idx(idx<1) = 1;
idx(idx(:,1)>sz(1),1) = sz(1);
idx(idx(:,2)>sz(2),2) = sz(2);
idx(idx(:,3)>sz(3),3) = sz(3);

lidx = idx(:,1);
lidx = lidx + (idx(:,2)-1).*k(1);
lidx = lidx + (idx(:,3)-1).*k(2);
end

function [label,labelprobs] = getLabels(atlasval, atlasvals, atlaslabels, varargin)
label = {'Unknown'};
labelprobs = {};
if ~isempty(atlasval)
    [uAtlasVal,~,uAtlasIdx] = unique(atlasval);
    uAtlasCnt = sum((1:length(uAtlasVal))==uAtlasIdx)';
    uAtlasVal = sortrows([uAtlasVal,uAtlasCnt],2,'descend');
    [uVoxB,uVoxVal] = ismember(uAtlasVal(:,1),atlasvals);
    labels = atlaslabels(uVoxVal(uVoxB),:);
    if ~isempty(labels)
        idx = find(cellfun(@isempty,regexp(labels,'White|Ventricle|Vent|CSF|vessel|Chiasm|Unknown|Fluid')),1); %these are not gray matter so skip if gray matter exists nearby
        label = labels(1,:); %most common label
        labelprobs = [labels,num2cell(uAtlasVal(uVoxB,2)/length(atlasval))];
        if ~isempty(idx)
            if (idx>=1) && (idx<=size(labels,1))
                label = labels(idx,:); %most common non-white label (biases toward gray)
            end
        end
    end
end
end


