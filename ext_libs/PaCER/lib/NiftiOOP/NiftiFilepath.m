classdef NiftiFilepath
    properties (Access = public) 
        patientId = [];
        modality = [];
        referenceSystem = [];
        warningAttribute = [];
        dataRootFolder = [];
    end
    
    properties (Access = protected)
        directFilepath = [];
        openDialogIfNotExists = [];
        
    end
    
    properties (Constant) 
        supportedReferenceSystems = {'rigidIntraCT', 'affineMni', 'affineMniDeprecated', ...
                                     'nonlinearRefMni', 'raw', 'rigidIntra', 'rigidMni' };
        supportedModalities = {'ct', 'swan', 'swan-brain', 'swan-sift', 't1', 't1-brain', 't1-sift', 't2', 't2-brain',...
            't2-sift', 'ct-preop', 'ct-postop', 'xml-plan', 'extracted-electrodes', 't1-vess', 't1-vess-dilate', 'swan-vess', 'swan-vess-dilate' }
        supportedFileEndings = {'.nii', '.nii.gz'};
%         , 'segmentation_swan', 'segmentation', 'segmentation-manual' ...
%             'segmentation-auto', 'segmentation-morel', 'DTI'};%... %             'segmentation_asmFromLm', ...
%             'segmentation_asm', 'segmentation_by_registration'};%, 'segmentation_asm'};
        % TODO: segmentation-asm, segmentation-asm-lm, dti
        noWarning = 1;
        showWarning = 2;
        error = 3;
    end
    
    properties (Dependent)
        filepath;
        filename;
        folder;
        patientFolder;
    end
    
    methods
        function this = NiftiFilepath(patientIdOrFilepath, modality, referenceSystem, ...
                openDialogIfNotExists, warningAttribute, dataRootFolder)
            if (nargin == 0)
                
            elseif ( nargin == 1 )
                if ( exist(patientIdOrFilepath, 'file') == 2 ) % first argument is full filepath
%                    warning('NiftiFilepath:staticPath', 'NiftiFilepath object has been created with static path, therefore this instance does NOT support dynamic filepath generation!');
                    this.directFilepath = patientIdOrFilepath;
                else % try to parse filepath and check if file exists in new dataRootFolder (useful if Const.dataRootFolder  has been changed)
                    try
                        [folder, name, ext] = fileparts(patientIdOrFilepath);
                        [folderWithoutRefSys, refSys] = fileparts(folder);
                        [~,patId] = fileparts(folderWithoutRefSys);
                    
                        this.directFilepath = [Const.dataRootFolder  filesep patId filesep refSys filesep name ext];
                        
                        if ( ~exist(this.directFilepath, 'file') )
                            throw(MException('',''));
                        end
                        
                    catch ME
                        warning off backtrace;
                        throw(MException('NiftiFilepath:NiftiFilepath:staticPathNotFound', ['File ''' patientIdOrFilepath ''' does not exist!']));
                    end
                end
            else
                if ( ~isempty(regexpi(modality, '*segmentation*') ) )
                    warning('Using NiftiFilepath with segmentations is deprecated. Use NiftiSegmentationFilepath instead!');
                end
                if ( ~exist('openDialogIfNotExists', 'var') )
                    openDialogIfNotExists = false;
                end
                if ( ~exist('warningAttribute', 'var') )
                    warningAttribute = NiftiFilepath.error;
                end
                this.patientId = patientIdOrFilepath;
                this.modality = modality;  
                this.referenceSystem = referenceSystem;  
                
                this.openDialogIfNotExists = openDialogIfNotExists;
                this.warningAttribute = warningAttribute;
                if(nargin == 6)
                    this.dataRootFolder = dataRootFolder;
                else
                    this.dataRootFolder = Const.dataRootFolder;
                end
                
            end
        end
        
        function fp=get.filepath(this)
            fp=this.getFilepath();
        end
        
        function fp = getFilepath(this)
            if ( ~isempty(this.directFilepath) )
                fp = this.directFilepath;
            else
                fp = NiftiFilepath.getFilepathStatic(this.patientId, ...
                                                    this.modality, ...
                                                    this.referenceSystem, ...
                                                    this.openDialogIfNotExists, ...
                                                    this.warningAttribute,...
                                                this.dataRootFolder);
            end
        end
        
        function fp=get.folder(this)
            fp = this.filepath;
            fp = fileparts(fp);
        end
        
        function f=get.patientFolder(this)
            f=fileparts(this.folder);
        end
        
        function fname=get.filename(this)
            fp = this.filepath;
            
            
            [~,fname,~]=fileparts(this.filepath);
            [~,fname,~]=fileparts(fname);
%             feIdx = regexpi(fname,[filesep '.']); % '\' is not valid on windows and leads to error!
%             if (~isempty(feIdx))
%                 fname(feIdx:end) = [];
%             end    
        end
    end
    
    methods (Static)
        function niftiFilepath = getFilepathStatic(patientId, modality, referenceSystem, openDialogIfNotExists, warningAttribute, dataRootFolder)
            if ( ~exist('openDialogIfNotExists', 'var') )
                openDialogIfNotExists = false;
            end
            if ( ~exist('warningAttribute', 'var') )
                warningAttribute = NiftiFilepath.error;
            end
            if ( ~exist('openDialogIfNotExists', 'var') )
                dataRootFolder = Const.dataRootFolder;
            end
            if ( isnumeric(patientId) )
                patientId = num2str(patientId);
            end
            patientDir = getFilepathFromWildcard(dataRootFolder , ['*' patientId '*']);
            
            % patient does not exist
            if ( isempty(patientDir) )
                if ( openDialogIfNotExists )
                    niftiFilepath = NiftiFilepath.openFilepathDialog();      
                    return;
                else
                    warnErrorStr = 'Patient ID invalid. Set argument openDialogIfNotExists to true in order to show a file dialog for manual filepath selection!';
                    if ( warningAttribute == NiftiFilepath.noWarning )
                        niftiFilepath = [];
                        return;
                    elseif ( warningAttribute == NiftiFilepath.showWarning )
                        warning(warnErrorStr);
                        niftiFilepath = [];
                        return;
                    else
                        throw(MException('NiftiFilepath:getFilepathStatic:InvalidPatientId', warnErrorStr));
                    end
                end
            end
            
            % unsupported referenceSystem
            if ( ~ismember(lower(referenceSystem), lower(NiftiSegmentationFilepath.supportedReferenceSystems)) )
                if ( openDialogIfNotExists )
                    niftiFilepath = NiftiFilepath.openFilepathDialog();      
                    return;
                else
                    suppRefSysStr = cell2str(NiftiFilepath.supportedReferenceSystems, ', ');
                    warnErrorStr = ['Reference System must be one of {' suppRefSysStr '}. Set argument openDialogIfNotExists to true in order to show a file dialog for manual filepath selection!'];
                    if ( warningAttribute == NiftiFilepath.noWarning )
                        niftiFilepath = [];
                        return;
                    elseif ( warningAttribute == NiftiFilepath.showWarning )
                        warning(warnErrorStr);
                        niftiFilepath = [];
                        return;
                    else
                        throw(MException('NiftiFilepath:getFilepathStatic:UnkownReferenceSystem', warnErrorStr));
                    end
                end
            end
            
            refSystemDir = [patientDir filesep referenceSystem];
            
            % supported reference system but no data of patient exists in
            % this reference system
            if ( exist(refSystemDir, 'file') ~= 7 ) % not a folder
                if ( openDialogIfNotExists )
                    niftiFilepath = NiftiFilepath.openFilepathDialog();      
                    return;
                else
                    warnErrorStr = ['No data of patient ' patientId ' exists in reference system {' referenceSystem '}. Set argument openDialogIfNotExists to true in order to show a file dialog for manual filepath selection!'];
  
                    if ( warningAttribute == NiftiFilepath.noWarning )
                        niftiFilepath = [];
                        return;
                    elseif ( warningAttribute == NiftiFilepath.showWarning )
                        warning(warnErrorStr);
                        niftiFilepath = [];
                        return;
                    else
                        throw(MException('NiftiFilepath:getFilepathStatic:NoDataForPatient', warnErrorStr));
                    end
                end
            end
            
            switch modality
                case {'xml-plan', 'XML-PLAN'}
                    modalityWildcard = '*\.xml$';
                case {'ct-preop', 'CT-PREOP', 'ct_preop', 'CT_PREOP'}
                    modalityWildcard = ['*CT_PREOP*' '*nii*'];
                case {'ct-postop', 'CT-POSTOP', 'ct_postop', 'CT_POSTOP'}
                    modalityWildcard = ['*CT_POSTOP*' '*nii*'];
                case {'DTI', 'dti'}
                    modalityWildcard = ['*DTI*' '*nii*']; %TODO adapt as soon as DTI data is available
                case {'SWAN', 'Swan', 'swan'}
                    modalityWildcard = ['*_SWAN*' '*nii*'];
                case {'SWAN-BRAIN', 'Swan-brain', 'swan-brain'}
                    modalityWildcard = ['*_SWAN*' Const.brainFileEndingRegexp '*nii*'];
                case {'T1', 't1'}
                    modalityWildcard = ['*3D_T1*' '*nii*'];
                case {'T1-brain', 't1-brain'}
                    modalityWildcard = ['*3D_T1*' Const.brainFileEndingRegexp '*nii*'];
                case {'T2', 't2'}
                    modalityWildcard =  ['*AX_*T2*' '*nii*'];
                case {'T2-brain', 't2-brain'}
                    modalityWildcard =  ['*AX_*T2*' Const.brainFileEndingRegexp '*nii*'];
                case {'t1-sift', 'T1-SIFT'}
                    modalityWildcard = ['*3D_T1*' Const.siftFeaturesFileEndingRegexp];
                case {'swan-sift', 'SWAN_SIFT'}
                    modalityWildcard = ['*MR_*_3D_SWAN*' Const.siftFeaturesFileEndingRegexp];
                case {'swan-vess', 'SWAN-VESS'}
                    modalityWildcard = ['*_SWAN_VESS*' '*nii*'];
                case {'swan-vess-dilate', 'SWAN-VESS-DILATE'}
                    modalityWildcard = ['*_SWAN_VESS_DILATE*' '*nii*'];
                case {'T1-vess', 'T1-VESS'}
                    modalityWildcard = ['*_T1_VESS*' '*nii*'];
                case {'T1-vess-dilate', 'T1-VESS-DILATE'}
                    modalityWildcard = ['*_T1_VESS_DILATE*' '*nii*'];
                case {'t2-sift', 'T2-SIFT'}
                    modalityWildcard = ['*AX_T2*' Const.siftFeaturesFileEndingRegexp];
                case {'extracted-electrodes', 'EXTRACTED-ELECTRODE'}
                    modalityWildcard = '*\.ele$';
                    
                    %                 case {'segmentation_swan', 'SEGMENTATION_SWAN', 'segmentation', 'SEGMENTATION', 'segmentation-manual'}
                    %                     modalityWildcard =  ['*' Const.segmentationFilePrefix '*' '*nii*'];
                    %                 case 'segmentation_asmFromLm'
                    %                     modalityWildcard =  ['*' Const.segmentationAsmFromLmFilePrefix '*' '*nii*'];
% %                 case 'segmentation_asm'
% %                     modalityWildcard =  ['*' Const.segmentationAsmFilePrefix '*.nii'];     
% %                 case 'segmentation_by_registration'
% %                     modalityWildcard =  ['*' Const.segmentationByRegistrationPrefix '*.nii'];
%                 case 'segmentation-auto'
%                     modalityWildcard =  ['*' Const.segmentationAutoFilePrefix '*' '*nii*'];
%                 case 'segmentation-morel'
%                     modalityWildcard =  ['*' Const.segmentationMorelFilePrefix '*' '*nii*'];
                otherwise
                    suppModStr = cell2str(NiftiFilepath.supportedModalities, ', ');
                    warnErrorStr = ['Modality must be one of {' suppModStr '}. Set argument openDialogIfNotExists to true in order to show a file dialog for manual filepath selection!'];
                    if ( warningAttribute == NiftiFilepath.noWarning )
                        niftiFilepath = [];
                        return;
                    elseif ( warningAttribute == NiftiFilepath.showWarning )
                        warning(warnErrorStr);
                        niftiFilepath = [];
                        return;
                    else
                        disp(warnErrorStr)
                        throw(MException('NiftiFilepath:getFilepathStatic:UnkonwnModality', warnErrorStr));
                    end
            end
                    
            niftiFilepath = getFilepathFromWildcard(refSystemDir, modalityWildcard);
        end
        
        function niftiFilepath = openFilepathDialog(this)
            [matfile, matpath, ~]=uigetfile({'*.nii.gz;*.nii'}, 'Please select .nii/.nii.gz file', ...
                                        this.dataRootFolder )
            niftiFilepath = [matpath matfile]; 
        end
        
        function supportedModalities = getSupportedModalities()
            supportedModalities = NiftiFilepath.supportedModalities;
        end
        
        function supportedReferenceSystems = getSupportedReferenceSystems()
            supportedReferenceSystems = NiftiFilepath.supportedReferenceSystems;
        end
        
        function patientsCell = getAvailablePatients(rootDir)
            if ( ~exist('rootDir', 'var') )
                rootDir = Const.dataRootFolder ;
            end
            patDir = dir(rootDir);
            
            patientsCell = {};
            for i=1:numel(patDir)
                folderName = patDir(i).name;
                if ( ~(strcmp(folderName, '.') || strcmp(folderName, '..') || strcmp(folderName(1), '.')) && exist([rootDir filesep folderName],'dir') == 7)
                    patientsCell{end+1} = folderName;
                end
            end
        end
        
    end
end