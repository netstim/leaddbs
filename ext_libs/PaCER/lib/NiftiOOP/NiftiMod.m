%% NiftiyMod -  Base class for all Nifti Objects. Wrapping the nifti Toolbox,
% rewritten no frills version of the more powerful NiftiModality Class
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicne
% 2016 - 2018
% mail@andreashusch.de
%
% Original NiftiModiality in collaboration Florian Bernard.

classdef NiftiMod < id & configurable
    properties (Access = protected) %references to data objects    
        header = [];         % the nifti header (real data) 
        data = [];      % raw nifti data isLoaded in ram as matlab matrix�
        filepathSuffix = [];
    end   
    properties (SetAccess = protected, GetAccess = public) %references to data objects
        voxdim = NaN(3,1); % voxel dimensions, P,L,I
        voxsize = NaN(3,1); % voxel size, P,L,I in mm
        transformationMatrix = NaN(4,4); % matrix that transforms voxel indices to world coordinates in mm, as a default this is set to the qform
        sform = NaN(4,4); % sform matrix that transforms voxel indices to world coordinates in mm
        qform = NaN(4,4); % qform (dervied) matrix that transforms voxel indices to world coordinates in mm

        isLoaded;  % true = in memory, false = on disk    
        argParser = []; 
        
        niftiFilepathObject = [];
    end
    properties (Access = public, SetObservable = true) % public access via specific set/get methods
        isToBeCached=true;  % true => data is kept in memory as long as possible
    end
    properties (Dependent = true, SetObservable = true, Hidden = true) % create "nifti structure like" interface with "virtual" properties
         hdr;
         img;
    end
    
    properties (Access = public, SetObservable = true, Transient = true, Dependent) % public access via specific set/get methods
       %HA: FIXME: Dependent!?
        filepath = '';  % full qualified path + filename (normaly a nifti file) 
    end
    
    methods
        function this = NiftiMod(niftiFilepathObject, varargin)
            if( ~isa(niftiFilepathObject, 'NiftiFilepath'))
                niftiFilepathObject = NiftiFilepath(niftiFilepathObject);
            end
            
            if ( ~isa(niftiFilepathObject, 'NiftiFilepath') )
                error('First argument must be a NiftiFilepath object!');
            end
                            
            this = this@id();                     % call superclass constructor to get an ID 
            
            this.isLoaded = false;
            this.data = [];
            this.niftiFilepathObject = niftiFilepathObject;
            
            %% Load hdr from disk
            this.header = load_untouch_header_only(niftiFilepathObject.filepath);
            
            %% Assign Header Fields (tranformationMatrix, voxdim, voxsize etc.)
            this.assignHeaderFields();
            
            %% Check
            checkNiftiHdr(this);
            
            %% argParser
            this.argParser = inputParser();
            this.argParser.KeepUnmatched = true;
            this.argParser.addOptional('isToBeCached', true);
            this.argParser.parse(varargin{:});
            args = this.argParser.Results;
            this.isToBeCached = args.isToBeCached;
        end
        
        function str = toString(this)
           [~, str] = fileparts(this.filepath);
        end
        
        function value = get.filepath(this)
            value = this.niftiFilepathObject.filepath;
        end
        
        function [image, header] = load(this)
            if( ~this.isLoaded && ~strcmp(this.filepath,'')) % it is not already isLoaded and we were not called with empty string (by Matlab callingn our get functions without command after initiating!)
                disp(['Loading ' this.filepath ' from disk...']);
                % nifti = load_nii(this.filepath, [], [], [], [], [], 0.1, 'q'); % load with 0.1 tolerance, prefer the q-form over the sform!
                nifti = load_untouch_nii(this.filepath); % load with 0.1 tolerance, prefer the q-form over the sform!
                
                this.header  = nifti.hdr; % header should maybe kept all the time ..
                
                %% Assigning Header Fields (tranformationMatrix, voxdim, voxsize etc.)
                this.assignHeaderFields();

                %% caching
                if(this.isToBeCached) % keep data in attributes if caching is enabled
                    disp('Caching is enabeld.');
                    
                    this.data = nifti.img;
                    this.isLoaded = true;
                else
                    disp('Caching is disabled. NOT keeping original volume in memory! Set this.isToBeCached = true to change this behaviour.');
                    this.data = [];
                    this.isLoaded = false;
                end
                header = nifti.hdr; % return data all the time
                
                %% precision
                image = single(nifti.img);
            end
        end
        
        function worldCoordinates=getNiftiWorldCoordinatesFromMatlabIdx(...
                this, voxelIdxList,... 
                orientation)
            trans = this.transformationMatrix;
            
            if ( size(voxelIdxList,1) ~= 3)
                error('voxelIdxList needs to be a 3 x N matrix');
            end

            worldCoordinates = trans*[voxelIdxList-1; ones(1, size(voxelIdxList,2))]; % not the -1 for the matlab indexing!
            worldCoordinates = worldCoordinates(1:3,:);
        end
        
        function [voxelIdxList,voxelToNiftiWorldMatrix] = ...
                getMatlabIdxFromNiftiWorldCoordinates(this, worldCoordList,...
                orientation)
            trans = this.transformationMatrix;
            
            if ( ~isempty(worldCoordList) && size(worldCoordList,1) ~= 3)
                error('worldCoordList needs to be a 3 x N matrix');
            end

			Imatlab = eye(4);
			Imatlab(1:3,end) = -1;
			voxelToNiftiWorldMatrix = trans*Imatlab; % matlab index starts with 1, not 0
			
			voxelIdxList = [];
			if ( ~isempty(worldCoordList) )
				voxelIdxList = voxelToNiftiWorldMatrix\...
					[worldCoordList; ones(1, size(worldCoordList,2))];
				voxelIdxList = voxelIdxList(1:3,:);
			end
        end
        
        function passivate(this) % remove data from memory. will be automatically reisLoaded if needed
            disp(['Passivating ' this.filepath ' ...']);
           
            this.data = [];
            this.isLoaded = false;
        end
        
        function image = get.img(this)
            if(~this.isLoaded)  % load nifti from disk when needed
                [image, ~] = this.load();
            else
                image = this.data;
            end          
        end
        
        function save(this, newFilepath, newImg, isSegmentation, forceFloat)
            % this.hdr has been changed during load_nii call (in contrast to
            % load_untouch_nii function). however, here we need exactly the
            % original header, therefore it is loaded again
            nii.hdr = load_untouch_header_only(this.niftiFilepathObject.filepath);
            
            if ( exist('isSegmentation','var') && isSegmentation )
                nii.hdr.dime.datatype = 2; % uint8 is sufficient for segmentation
            end
            if ( exist('forceFloat','var') && forceFloat )
                nii.hdr.dime.datatype = 16;
            end
            
            % nii.hdr = this.hdr;
            if ( this.hdr.dime.pixdim(1) ~= nii.hdr.dime.pixdim(1) )
                warning('niftiModality:swapLR', 'Swapping left and right because original and touched header do not agree, make sure to check if the nifti file is as expected!');
                nii.img = newImg(end:-1:1,:,:);
            else
                nii.img = newImg;
            end
            nii.untouch = 1;
            
            disp(['Saving new image data in ' newFilepath ' to disk...']);
            
            [pathstr,filename,ext] = fileparts(newFilepath);
            if(isequal(ext,''))
                newFilepath = [newFilepath, '.nii'];
            elseif(isequal(ext,'.gz'))
                newFilepath = [pathstr filesep filename];
            elseif(~isequal(ext,'.nii'))
                error('NiftiModality:save:WrongFileending','Current Fileending is unkown');
            end
            
            save_untouch_nii(nii, newFilepath);
            
            
            % save zipped file
            gzip(newFilepath); %TODO make this configurable
            delete(newFilepath);
            %             this.filepath = newFilepath;
            %             newFilepath = [newFilepath, '.gz'];
            %             this.niftiFilepathObject = NiftiFilepath(newFilepath);
            
        end
        
        function header = get.hdr(this)
            if(isempty(this.header))  % load nifti from disk when needed
                %  [~, header] = this.load();
                
                %                 as the header is transformed during load it is not possible
                %                 to load the header only as done below (leads to equal header
                %                 as returned by load_untouch_nii()
                this.header = load_untouch_header_only(this.filepath);
                
                %% Assign Header Fields (tranformationMatrix, voxdim, voxsize etc.)
                this.assignHeaderFields();
                
                header = this.header;
            else
                header = this.header;
            end
        end
    end
    methods (Access = private)
        function assignHeaderFields(this)
            this.voxdim  = this.hdr.dime.dim(2:4)';
            this.voxsize = this.hdr.dime.pixdim(2:4)';
            
            %% getting sform
            this.sform = [this.hdr.hist.srow_x;...
                this.hdr.hist.srow_y;...
                this.hdr.hist.srow_z;...
                [0 0 0 1]];
            
            if(trace(this.sform(1:3,1:3)) == 0)
                disp('Notice: NiftiMod: nifti sform is empty');
            end
            
            %% getting qform by transforming the quaternion data to homogenous matrix
            h = this.hdr.hist;
            h.pixdim = this.hdr.dime.pixdim;
            
            this.qform = cbiQuaternionToHomogeneous(h);
            
            %% setting transformationMatrix to the qform
            this.transformationMatrix = this.qform;
        end
    end
end