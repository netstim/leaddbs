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
        
%         function loadDefFile(obj, varargin)
%             obj.parseInput(varargin{:});
%             if isempty(obj.DefFile)
%                 error('Path to deformation file (DefFile) must be provided!')
%             end
%             
%             Nii = nifti(obj.DefFile);
%             obj.Def = single(Nii.dat(:,:,:,1,:));
%             d = size(obj.Def);
%             if d(4)~=1 || d(5)~=3, error('Deformation field is wrong!'); end
%             obj.Def = reshape(obj.Def,[d(1:3) d(5)]); %imgvox2mnimm image deformation
%             obj.DefMat = Nii.mat;
%         end
        function loadDefFile(obj, varargin)
            % Deformation loader using ea_load_nii (SPM-based).
            % Expects a field with 3 spatial components (x,y,z) spread across volumes
            % or a trailing dim of size 3.
            %
            % Sets:
            %   obj.Def    -> single [X Y Z 3]
            %   obj.DefMat -> 4x4 affine from NIfTI header
        
            obj.parseInput(varargin{:});
            if isempty(obj.DefFile)
                error('Path to deformation file (DefFile) must be provided!');
            end
            nii = ea_load_nii(obj.DefFile);
            data = nii.img;          
            obj.Def    = single(nii.img);
            obj.DefMat = nii.mat;
        end
        
%         function [Label,XYZMNImm] = findLabel(obj, XYZIdx, AtlasName, varargin)
%             %[Label,XYZMNImm] = findLabel([1,1,1],'NMM');
%             %XYZIdx are voxel indices (not mm) in original image (not mni)
%             %associated with deformation file (DefFile)
%             if isempty(obj.Def)
%                 obj.loadDefFile(varargin{:});
%             end
%             
%             Idx = strcmp(regexprep(obj.AtlasNames,'\.nii$',''),AtlasName);            
%             if isempty(obj.VoxDef{Idx})  
%                 mm2vox_atlas = inv(obj.AtlasMat{Idx}); %atlas transform
%                 obj.VoxDef(Idx) = {affine(obj.Def,mm2vox_atlas)}; %vox2vox transform (patient to mni)
%             end
%                         
%             %%%%%%%%% Faster linear index version %%%%%%%%%%%%%%%%
%             sz_vox = size(obj.VoxDef{Idx});
%             n_vox = prod(sz_vox(1:3));
%             
%             lidx = calcLinIdx(sz_vox(1:3),XYZIdx);
%             
%             %MNI mm locations for each voxel in pt space (XYZIdx)
%             MNImm = nan(length(lidx),3);
%             MNImm(:,1) = obj.Def(lidx);
%             MNImm(:,2) = obj.Def(lidx+n_vox);
%             MNImm(:,3) = obj.Def(lidx+2*n_vox);
%             XYZMNImm = nanmean(MNImm);
%                         
%             %Atlas voxel indices for each voxel in pt space (XYZIdx)
%             atlasvox = nan(length(lidx),3);
%             atlasvox(:,1) = round((obj.VoxDef{Idx}(lidx)));
%             atlasvox(:,2) = round((obj.VoxDef{Idx}(lidx+n_vox)));
%             atlasvox(:,3) = round((obj.VoxDef{Idx}(lidx+2*n_vox)));
%             
%             sz_atlas = size(obj.AtlasData{Idx});
%             
%             lidx_atlas = calcLinIdx(sz_atlas(1:3),atlasvox);
%             nan_idx = isnan(lidx_atlas);
%             
%             %Atlas label indices (corresponds to atlas labels in text file)
%             AtlasVal = zeros(length(lidx_atlas),1);
%             AtlasVal(~nan_idx) = round(obj.AtlasData{Idx}(lidx_atlas(~nan_idx))); %label indices should be integers, but rounding just in case
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             
%             Label = getLabels(AtlasVal, obj.AtlasVals{Idx}, obj.AtlasLabels{Idx});
%         end %end findLabel

%         function [Label, XYZMNImm, AtlasVox, AtlasVal] = findLabel(obj, XYZIdx, AtlasName, varargin)
%         % Map patient-space voxel indices -> MNI mm via ea_map_coords, then read atlas labels.
%         % Returns a label PER VOXEL plus the mean MNI mm of the set.
%         %
%         % Required name-values:
%         %   'SrcImg'      : source image path whose grid XYZIdx refer to (e.g., MR or CT)
%         %   'InvXformDir' : subject's inverseTransform directory
%         %   'Template'    : template (MNI) image path used for normalization
%         %   'NormMethod'  : 'ants' or 'spm'
%         
%             % ---- parse inputs
%             p = inputParser;
%             addParameter(p,'SrcImg','',@(s)ischar(s)||isstring(s));
%             addParameter(p,'InvXformDir','',@(s)ischar(s)||isstring(s));
%             addParameter(p,'Template','',@(s)ischar(s)||isstring(s));
%             addParameter(p,'NormMethod','',@(s)ischar(s)||isstring(s));
%             parse(p,varargin{:});
%             srcImg   = char(p.Results.SrcImg);
%             invDir   = char(p.Results.InvXformDir);
%             template = char(p.Results.Template);
%             normmtd  = char(p.Results.NormMethod);
%         
%             if any(cellfun(@isempty,{srcImg,invDir,template,normmtd}))
%                 error('Provide SrcImg, InvXformDir, Template, and NormMethod.');
%             end
%             if size(XYZIdx,2)~=3, error('XYZIdx must be N×3 voxel indices.'); end
%         
%             % ---- pick atlas
%             atlasClean = regexprep(obj.AtlasNames,'\.nii(\.gz)?$','');
%             targetName = regexprep(AtlasName,'\.nii(\.gz)?$','');
%             aIdx = find(strcmp(atlasClean, targetName),1);
%             if isempty(aIdx), error('Atlas "%s" not found.', AtlasName); end
%         
%             Avox2mm      = obj.AtlasMat{aIdx};        % atlas vox→mm
%             mm2vox_atlas = inv(Avox2mm);
%             atlasVol     = obj.AtlasData{aIdx};
%             szA          = size(atlasVol);
%         
%             % ---- vox -> native mm (vectorized)
%             Vsrc      = spm_vol(srcImg);
%             hom_vox   = [double(XYZIdx), ones(size(XYZIdx,1),1)];   % N×4
%             mm_native = (Vsrc.mat * hom_vox.').';                   % N×4
%             mm_native = mm_native(:,1:3);                           % N×3 (mm)
%         
%             % ---- native mm -> MNI mm (vectorized if supported)
%             try
%                 mm_mni = ea_map_coords(mm_native', srcImg, invDir, template, normmtd);  % N×3
%             catch
%                 % fallback: loop if your ea_map_coords only accepts 1×3
%                 mm_mni = nan(size(mm_native));
%                 for i = 1:size(mm_native,1)
%                     mm_mni(i,:) = ea_map_coords(mm_native(i,:), srcImg, invDir, template, normmtd);
%                 end
%             end
%             mm_mni = mm_mni';
%             % report mean MNI mm for the set
%             XYZMNImm = mean(mm_mni, 1, 'omitnan');
%         
%             % ---- MNI mm -> atlas voxels, sample labels
%             hom_mni   = [mm_mni, ones(size(mm_mni,1),1)];
%             AtlasVox  = (mm2vox_atlas * hom_mni.').';               % N×4
%             AtlasVox  = round(AtlasVox(:,1:3));                     % nearest voxel
%         
%             inb = AtlasVox(:,1)>=1 & AtlasVox(:,1)<=szA(1) & ...
%                   AtlasVox(:,2)>=1 & AtlasVox(:,2)<=szA(2) & ...
%                   AtlasVox(:,3)>=1 & AtlasVox(:,3)<=szA(3);
%         
%             AtlasVal = zeros(size(AtlasVox,1),1);
%             if any(inb)
%                 lin = sub2ind(szA(1:3), AtlasVox(inb,1), AtlasVox(inb,2), AtlasVox(inb,3));
%                 AtlasVal(inb) = round(atlasVol(lin));
%             end
%         
%             Label = getLabels(AtlasVal, obj.AtlasVals{aIdx}, obj.AtlasLabels{aIdx});  % 1 label per voxel
%         end

        function [Label, XYZMNImm, AtlasVox, AtlasVal] = findLabel(obj, XYZMNImm, AtlasName, varargin)
        % FINDLABEL  Get atlas label(s) from MNI mm coordinates.
        %
        %   [Label, XYZMNImm, AtlasVox, AtlasVal] = obj.findLabel(XYZMNImm, AtlasName)
        %
        % Inputs
        %   XYZMNImm : N×3 MNI coordinates (mm) in template space
        %   AtlasName: name of the atlas in obj.AtlasNames (with or without .nii/.nii.gz)
        %
        % Outputs
        %   Label    : cell array of label strings, one per row of XYZMNImm
        %   XYZMNImm : echo of the input MNI coords (N×3)
        %   AtlasVox : N×3 voxel indices in atlas space
        %   AtlasVal : N×1 numeric atlas values (indices into AtlasVals/AtlasLabels)
        
            %#ok<*INUSD>  % varargin unused now but kept for backwards compatibility
        
            % ---- basic checks ----
            if isempty(XYZMNImm)
                Label    = {};
                AtlasVox = [];
                AtlasVal = [];
                return;
            end
        
            XYZMNImm = double(XYZMNImm);
            % Allow 3×1 vector
            if isvector(XYZMNImm) && numel(XYZMNImm)==3
                XYZMNImm = XYZMNImm(:).';   % 1×3
            end
            if size(XYZMNImm,2) ~= 3
                error('XYZMNImm must be N×3 MNI mm coordinates.');
            end
        
            % ---- locate atlas ----
            atlasClean = regexprep(obj.AtlasNames,'\.nii(\.gz)?$','');
            targetName = regexprep(AtlasName,'\.nii(\.gz)?$','');
            aIdx = find(strcmp(atlasClean, targetName),1);
            if isempty(aIdx)
                error('Atlas "%s" not found in obj.AtlasNames.', AtlasName);
            end
        
            Avox2mm  = obj.AtlasMat{aIdx};    % atlas vox -> mm
            if isequal(size(Avox2mm),[3 4])
                Avox2mm = [Avox2mm; 0 0 0 1]; % make it 4×4 if stored as 3×4
            end
            mm2vox_atlas = inv(Avox2mm);
        
            atlasVol = obj.AtlasData{aIdx};
            szA      = size(atlasVol);
        
            % ---- MNI mm -> atlas vox ----
            N       = size(XYZMNImm,1);
            hom_mni = [XYZMNImm, ones(N,1)];            % N×4
            vox4    = (mm2vox_atlas * hom_mni.').';     % N×4
            AtlasVox = round(vox4(:,1:3));              % N×3
        
            % ---- bounds + sampling ----
            inb = AtlasVox(:,1)>=1 & AtlasVox(:,1)<=szA(1) & ...
                  AtlasVox(:,2)>=1 & AtlasVox(:,2)<=szA(2) & ...
                  AtlasVox(:,3)>=1 & AtlasVox(:,3)<=szA(3);
        
            AtlasVal = zeros(N,1);
            if any(inb)
                lin = sub2ind(szA(1:3), ...
                              AtlasVox(inb,1), ...
                              AtlasVox(inb,2), ...
                              AtlasVox(inb,3));
                AtlasVal(inb) = round(atlasVol(lin));
            end
        
            % ---- convert numeric labels → text labels ----
            Label = getLabels(AtlasVal, obj.AtlasVals{aIdx}, obj.AtlasLabels{aIdx});
        end
        
%         function LabelProbs = findLabelProbabilities(obj, COMIdx, AtlasName, varargin)
%             %Uses voxel COM (COMIdx) to create a sphere of specified radius
%             %for finding label probabilities for specified atlas
%             %(AtlasName). COMIdx are voxel indices in patient space.
%                         
%             if isempty(obj.Def)
%                 obj.loadDefFile(varargin{:});
%             else
%                 obj.parseInput(varargin{:});
%             end
% 
%             RadiusMM = 10; %listening radius (mm) for finding labels
%             if ~isempty(obj.ProbRadius)
%                 RadiusMM = obj.ProbRadius;
%             end
%             
%             Idx = strcmp(regexprep(obj.AtlasNames,'\.nii$',''),AtlasName);            
%             if isempty(obj.VoxDef{Idx})
%                 mm2vox_atlas = inv(obj.AtlasMat{Idx}); %atlas transform
%                 obj.VoxDef(Idx) = {affine(obj.Def,mm2vox_atlas)}; %vox2vox transform (patient to mni)
%             end
%             
%             voxdef = obj.VoxDef{Idx};
%             atlasdata = obj.AtlasData{Idx};
%             atlasvals = obj.AtlasVals{Idx};
%             atlaslabels = obj.AtlasLabels{Idx};
%             
%             XYZScale = sqrt(sum(obj.DefMat(1:3,1:3).^2)); %pt mm/vox scaling
%             RadiusVox = ceil(RadiusMM./XYZScale);
%             
%             LabelProbs = cell(size(COMIdx,1),1);
%             for m=1:size(COMIdx,1)
%                 
%                 % full index values centered on click point (XYZ)
%                 xidx = COMIdx(m,1)-RadiusVox(1):COMIdx(m,1)+RadiusVox(1);
%                 yidx = COMIdx(m,2)-RadiusVox(2):COMIdx(m,2)+RadiusVox(2);
%                 zidx = COMIdx(m,3)-RadiusVox(3):COMIdx(m,3)+RadiusVox(3);
%                 [YIdx,XIdx,ZIdx] = meshgrid(yidx,xidx,zidx);
%                 
%                 com = RadiusVox+1;
%                 
%                 [my,mx,mz] = meshgrid(1:length(yidx),1:length(xidx),1:length(zidx));
%                 mx = (mx - com(1)).*XYZScale(1);
%                 my = (my - com(2)).*XYZScale(2);
%                 mz = (mz - com(3)).*XYZScale(3);
%                 
%                 elecimg = (sqrt(mx.^2+my.^2+mz.^2)<=RadiusMM); %1s are sphere around center, 0s everywhere else
%                 
%                 fullidx = [XIdx(elecimg),YIdx(elecimg),ZIdx(elecimg)]; %full index of electrode sphere
%                 
%                 %%%%%%%%% Faster linear index version %%%%%%%%%%%%%%%%
%                 sz_vox = size(voxdef);
%                 n_vox = prod(sz_vox(1:3));
%                 
%                 lidx = calcLinIdx(sz_vox(1:3),fullidx);
%                 
%                 atlasvox = nan(length(lidx),3);
%                 atlasvox(:,1) = round((voxdef(lidx)));
%                 atlasvox(:,2) = round((voxdef(lidx+n_vox)));
%                 atlasvox(:,3) = round((voxdef(lidx+2*n_vox)));
%                 
%                 sz_atlas = size(atlasdata);
%                 
%                 lidx_atlas = calcLinIdx(sz_atlas(1:3),atlasvox);
%                 nan_idx = isnan(lidx_atlas);    
%                 
%                 AtlasVal = zeros(length(lidx_atlas),1);
%                 AtlasVal(~nan_idx) = round(atlasdata(lidx_atlas(~nan_idx)));
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 
%                 [~,labelprobs] = getLabels(AtlasVal, atlasvals, atlaslabels);
%                 LabelProbs(m) = {labelprobs};
%             end %end for
%             
%         end %end function
        function LabelProbs = findLabelProbabilities(obj, COMIdx, AtlasName, varargin)
        % Uses voxel COM(s) (COMIdx, N×3 in source/patient grid) to create a sphere
        % (RadiusMM) and compute label probabilities from the chosen atlas in MNI space.
        %
        % Name-Value (required):
        %   'SrcImg'      : path to the source image whose grid COMIdx refer to (e.g., MR/CT)
        %   'Template'    : path to the atlas/template image in MNI space
        %   'InvXformDir' : subject directory containing 'inverseTransform' (Lead-DBS)
        %   'NormMethod'  : normalization method string, e.g. 'ants','spm','fsl'
        %
        % Name-Value (optional):
        %   'RadiusMM'    : scalar spherical radius in mm (default = 10 or obj.ProbRadius if set)
        %
        % Output:
        %   LabelProbs : cell(N,1). Each cell is the label-probability map returned by:
        %                [~, labelprobs] = getLabels(AtlasVal, atlasvals, atlaslabels);
        
            % ---------- parse inputs ----------
            p = inputParser;
            addParameter(p,'SrcImg','',@(s)ischar(s)||isstring(s));
            addParameter(p,'Template','',@(s)ischar(s)||isstring(s));
            addParameter(p,'InvXformDir','',@(s)ischar(s)||isstring(s));
            addParameter(p,'NormMethod','',@(s)ischar(s)||isstring(s));
            addParameter(p,'RadiusMM',[],@(x)isnumeric(x)&&isscalar(x)&&x>0);
            parse(p,varargin{:});
        
            srcImg   = char(p.Results.SrcImg);
            template = char(p.Results.Template);
            invDir   = char(p.Results.InvXformDir);
            normmtd  = upper(char(p.Results.NormMethod));
        
            if any(cellfun(@isempty,{srcImg,template,invDir,normmtd}))
                error('findLabelProbabilities:MissingInputs', ...
                      'Provide SrcImg, Template, InvXformDir, and NormMethod.');
            end
            if size(COMIdx,2)~=3
                error('findLabelProbabilities:BadCOMIdx','COMIdx must be N×3 voxel indices.');
            end
        
            % radius selection: NV pair > obj.ProbRadius > default 10
            if ~isempty(p.Results.RadiusMM)
                RadiusMM = p.Results.RadiusMM;
            elseif ~isempty(obj.ProbRadius)
                RadiusMM = obj.ProbRadius;
            else
                RadiusMM = 10;
            end
        
            % ---------- resolve atlas (affine + volume + lut) ----------
            atlasClean = regexprep(obj.AtlasNames,'\.nii(\.gz)?$','');
            targetName = regexprep(AtlasName,'\.nii(\.gz)?$','');
            aIdx = find(strcmp(atlasClean, targetName), 1);
            if isempty(aIdx), error('Atlas "%s" not found.', AtlasName); end
        
            Avox2mm = obj.AtlasMat{aIdx};                % atlas vox -> mm
            if isequal(size(Avox2mm),[3 4]), Avox2mm=[Avox2mm; 0 0 0 1]; end
            if ~isequal(size(Avox2mm),[4 4]), error('Unexpected atlas affine size.'); end
            mm2vox_atlas = inv(Avox2mm);
        
            atlasVol    = obj.AtlasData{aIdx};
            atlasvals   = obj.AtlasVals{aIdx};
            atlaslabels = obj.AtlasLabels{aIdx};
            szA         = size(atlasVol);
        
            % ---------- source image geometry ----------
            Vsrc  = spm_vol(srcImg);
            vxmm  = sqrt(sum(Vsrc.mat(1:3,1:3).^2));     % [dx dy dz] mm/vox
            isz   = Vsrc.dim(1:3);
            Rvox  = ceil(RadiusMM ./ vxmm);              % radius in voxels
        
            % inverseTransform stub for ea_map_coords
            transStub = invDir;
            if ~endsWith(transStub,'inverseTransform')
                transStub = fullfile(transStub,'inverseTransform');
            end
        
            % ---------- per-COM neighborhood ----------
            LabelProbs = cell(size(COMIdx,1),1);
        
            for m = 1:size(COMIdx,1)
                c = double(COMIdx(m,:));  % [i j k] in src grid
        
                % neighborhood ranges (clipped to bounds)
                xidx = (c(1)-Rvox(1)):(c(1)+Rvox(1)); xidx = xidx(xidx>=1 & xidx<=isz(1));
                yidx = (c(2)-Rvox(2)):(c(2)+Rvox(2)); yidx = yidx(yidx>=1 & yidx<=isz(2));
                zidx = (c(3)-Rvox(3)):(c(3)+Rvox(3)); zidx = zidx(zidx>=1 & zidx<=isz(3));
                if isempty(xidx) || isempty(yidx) || isempty(zidx)
                    LabelProbs{m} = containers.Map(); %#ok<AGROW>
                    continue
                end
        
                % grid the subvolume and form a spherical mask in mm
                [YIdx,XIdx,ZIdx] = meshgrid(yidx,xidx,zidx);              % note (y,x,z) order for meshgrid
                % find the COM index within the subvolume
                cx = find(xidx==c(1),1); cy = find(yidx==c(2),1); cz = find(zidx==c(3),1);
                [my,mx,mz] = meshgrid(1:numel(yidx),1:numel(xidx),1:numel(zidx));
                mx = (mx - cx) .* vxmm(1);
                my = (my - cy) .* vxmm(2);
                mz = (mz - cz) .* vxmm(3);
                sphereMask = (mx.^2 + my.^2 + mz.^2) <= (RadiusMM^2);
        
                % voxels inside sphere (N×3 in src grid)
                fullidx = [XIdx(sphereMask), YIdx(sphereMask), ZIdx(sphereMask)];
        
                % ---- map patient vox -> MNI mm (batch; ea_map_coords prefers 3×N) ----
                XYZ_vx_3xN = fullidx.';                             % 3×N
                mm_mni_3xN = ea_map_coords(XYZ_vx_3xN, srcImg, transStub, template, normmtd); % 3×N
        
                % ---- MNI mm -> atlas vox; sample labels ----
                hom_mni     = [mm_mni_3xN; ones(1,size(mm_mni_3xN,2))];   % 4×N
                vox_atlas_4 = mm2vox_atlas * hom_mni;                      % 4×N
                AtlasVox    = round(vox_atlas_4(1:3,:)).';                 % N×3
        
                inb = AtlasVox(:,1)>=1 & AtlasVox(:,1)<=szA(1) & ...
                      AtlasVox(:,2)>=1 & AtlasVox(:,2)<=szA(2) & ...
                      AtlasVox(:,3)>=1 & AtlasVox(:,3)<=szA(3);
        
                AtlasVal = zeros(size(AtlasVox,1),1,'like',atlasVol);
                if any(inb)
                    linA = sub2ind(szA(1:3), AtlasVox(inb,1), AtlasVox(inb,2), AtlasVox(inb,3));
                    AtlasVal(inb) = round(atlasVol(linA));
                end
        
                % aggregate to probabilities
                [~, labelprobs] = getLabels(AtlasVal, atlasvals, atlaslabels);
                LabelProbs{m} = labelprobs;                              %#ok<AGROW>
            end
        end

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


