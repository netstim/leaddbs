%% NiftiyModSPM - Base class for all Nifti Objects. Wrapping the SPM Toolbox.
%
% Andreas Husch, Florian Bernard
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicne
% 2016 - 2017
% mail@andreashusch.de, fbernardpi@gmail.com
%
% Andreas Horn
% Charite University Medicine Berlin - Movement Disorders Unit
% 2018
%
% Andreas Husch
% University of Luxembourg - Luxembourg Centre for Systems Biomedicne
% 2018

classdef NiftiModSPM <  id & configurable
    properties (Access = protected) %references to data objects    
        header = [];         % the nifti header (real data) 
        data = [];      % raw nifti data isLoaded in ram as matlab matrix???
        filepathSuffix = [];
    end   
    properties (SetAccess = protected, GetAccess = public) %references to data objects
        voxdim = NaN(3,1); % voxel dimensions, P,L,I
        voxsize = NaN(3,1); % voxel size, P,L,I in mm
        transformationMatrix = NaN(4,4); % matrix that transforms voxel indices to world coordinates in mm
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
        function this = NiftiModSPM(niftiFilepathObject, varargin)
            if( ~isa(niftiFilepathObject, 'NiftiFilepath'))
                 niftiFilepathObject = NiftiFilepath(niftiFilepathObject);
            end
            
            if ( ~isa(niftiFilepathObject, 'NiftiFilepath') )
                error('First argument must be a NiftiFilepath object!');
            end
                            
            this = this@id();                     % call superclass constructor to get an ID 
            
            this.isLoaded = false;
            lhdr = ea_open_vol(niftiFilepathObject.filepath); % make sure to load header on object conostruction
            this.voxdim = lhdr.dim';
            this.voxsize= lhdr.voxsize';
            this.transformationMatrix = lhdr.mat;
            A = lhdr.mat(1:3,1:3);
            if(A(~eye(3)))
                warning('Transformation contains off-diagnonal elements. Check carefully and consider reordering of data');
            end

            this.data = [];
            this.niftiFilepathObject = niftiFilepathObject;
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
                
                disp(['Loading ' this.filepath ' from disk using NiftiModSPM...']);
                %nifti = ea_load_untouch_nii(this.filepath);
                
                eanifti = ea_load_nii(this.filepath);
                
                this.voxdim = eanifti.dim';
                this.voxsize= eanifti.voxsize';
                this.header = eanifti; % header should maybe kept all the time ..
                this.transformationMatrix = eanifti.mat;

                if(this.isToBeCached) % keep data in attributes if caching is enabled
                    disp('Caching is enabeled.');
                    
                    this.data = eanifti.img;
                    this.isLoaded = true;
                else
                    disp('Caching is disabled. NOT keeping original volume in memory! Set this.isToBeCached = true to change this behaviour.');
                    this.data = [];
                    this.isLoaded = false;
                end
                header = eanifti; % return data all the time
                image = single(eanifti.img);
            end
        end
        
        function worldCoordinates=getNiftiWorldCoordinatesFromMatlabIdx(this, ...
                voxelIdxList,varargin)
                                                            % matlab index starts with 1, not 0
            worldCoordinates = this.transformationMatrix*[voxelIdxList; ones(1, size(voxelIdxList,2))];
            worldCoordinates = worldCoordinates(1:3,:);
        end
        
        function voxelIdxList = getMatlabIdxFromNiftiWorldCoordinates(this, worldCoordList,...
                varargin)
            if ( ~isempty(worldCoordList) && size(worldCoordList,1) ~= 3)
                error('worldCoordList needs to be a 3 x N matrix');
            end            
            
            voxelIdxList = this.transformationMatrix\[worldCoordList; ones(1, size(worldCoordList,2))];
            voxelIdxList = voxelIdxList(1:3,:); % matlab index starts with 1, not 0 however SPM already compensated for that by altering the transformationMatrix :-)
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
            warning('Saving files is currently not supported by NiftiModSPM');
            return;
             % this.hdr has been changed during load_nii call (in contrast to
            % load_untouch_nii function). however, here we need exactly the
            % original header, therefore it is loaded again
            nii = ea_load_nii(this.niftiFilepathObject.filepath);
            
            if ( exist('isSegmentation','var') && isSegmentation )
                nii.dt(1) = 2; % uint8 is sufficient for segmentation
            end
            if ( exist('forceFloat','var') && forceFloat )
                nii.dt(1) = 16;
            end
            
%             %% we won't need below part due to using SPM for loading
%             % nii.hdr = this.hdr;
%             if ( this.dim(1) ~= nii.dim(1) )
%                 warning('niftiModality:swapLR', 'Swapping left and right because original and touched header do not agree, make sure to check if the nifti file is as expected!');
%                 nii.img = newImg(end:-1:1,:,:);
%             else
%                 nii.img = newImg;
%             end
%             nii.untouch = 1;
%             
%             disp(['Saving new image data in ' newFilepath ' to disk...']);
%             
%             [pathstr,filename,ext] = fileparts(newFilepath);
%             if(isequal(ext,''))
%                 newFilepath = [newFilepath, '.nii'];
%             elseif(isequal(ext,'.gz'))
%                 newFilepath = [pathstr filesep filename];
%             elseif(~isequal(ext,'.nii'))
%                  error('NiftiModality:save:WrongFileending','Current Fileending is unkown');
%             end
%             
%             save_untouch_nii(nii, newFilepath); 
%             
% 
%             % save zipped file
%             gzip(newFilepath); %TODO make this configurable
%             delete(newFilepath);
% %             this.filepath = newFilepath;
% %             newFilepath = [newFilepath, '.gz'];
% %             this.niftiFilepathObject = NiftiFilepath(newFilepath);
% %%
        end
        
        function header = get.hdr(this)
            if(isempty(this.header))  % load nifti from disk when needed
                [~, header] = this.load(); 
               
%                 as the header is transformed during load it is not possible
%                 to load the header only as done below (leads to equal header
%                 as returned by load_untouch_nii()
%                 [this.header] = load_nii_hdr(this.filepath);
%                 
%                 this.voxdim = this.header.dime.dim(2:4);
%                 this.voxsize= this.header.dime.pixdim(2:4);   
%                 header = this.header;
            else
                header = this.header;
            end
        end              
    end
end

function test() %#ok<DEFNU>
%%
niiFile = '/Users/fb/cnm/svn/resources/directions_mni.nii.gz';
nii = NiftiModSPM(NiftiFilepath(niiFile));
nii.load();

idx = [1 1 1; 10 10 10; 12 17 99]';
nii.getNiftiWorldCoordinatesFromMatlabIdx(idx)
%%
end