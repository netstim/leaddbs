%% NiftiyMod -  Base class for all Nifti Objects. Wrapping the nifti Toolbox,
% no frills version of the more powerful NiftiModality Class
%
% Andreas Husch, Florian Bernard
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicne
% 2016 - 2017
% mail@andreashusch.de, fbernardpi@gmail.com

classdef NiftiMod <  id & configurable
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
        function this = NiftiMod(niftiFilepathObject, varargin)
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
%             if(trace(this.transformationMatrix(1:3,1:3)) == 0) % workaround defective nifti headers
%                 this.transformationMatrix(1,1) = lhdr.dime.pixdim(2);
%                 this.transformationMatrix(2,2) = lhdr.dime.pixdim(3);
%                 this.transformationMatrix(3,3) = lhdr.dime.pixdim(4);
%                 this.transformationMatrix(1,4) = lhdr.hist.qoffset_x;
%                 this.transformationMatrix(2,4) = lhdr.hist.qoffset_y;
%                 this.transformationMatrix(3,4) = lhdr.hist.qoffset_z;
%             end
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

                disp(['Loading ' this.filepath ' from disk...']);
                %nifti = ea_load_untouch_nii(this.filepath);

                eanifti = ea_load_nii(this.filepath);

                this.voxdim = eanifti.dim';
                this.voxsize= eanifti.voxsize';
                this.header = eanifti; % header should maybe kept all the time ..
                this.transformationMatrix = eanifti.mat;
%                 if(trace(this.transformationMatrix(1:3,1:3)) == 0) % workaround defective nifti headers
%                     this.transformationMatrix(1,1) = nifti.hdr.dime.pixdim(2);
%                     this.transformationMatrix(2,2) = nifti.hdr.dime.pixdim(3);
%                     this.transformationMatrix(3,3) = nifti.hdr.dime.pixdim(4);
%                     this.transformationMatrix(1,4) = nifti.hdr.hist.qoffset_x; %
%                     this.transformationMatrix(2,4) = nifti.hdr.hist.qoffset_y;
%                     this.transformationMatrix(3,4) = nifti.hdr.hist.qoffset_z;
%                 end
                if(this.isToBeCached) % keep data in attributes if caching is enabled
                    disp('Caching is enabled.');

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

        function worldCoordinates=getNiftiWorldCoordinatesFromMatlabIdx(...
                this, voxelIdxList,...
                orientation)
            trans = this.transformationMatrix;

            if ( ~exist('orientation', 'var') )
                if(det(trans) < 0)
                    orientation = 'right-handed';
                else
                    orientation = 'left-handed';
                end
            end

            if ( size(voxelIdxList,1) ~= 3)
                error('voxelIdxList needs to be a 3 x N matrix');
            end

            switch orientation
                case {'left-handed'} % left-handed coordinate system
                    % nothing to do
                case {'right-handed'} % right-handed coordinate system
                    disp('NiftiMod: Altering transformation hopefully yielding LPI- coordinates. Validate slice orientation.');
                    trans(1,1) = -trans(1,1);
                    trans(1,4) = -((this.voxdim(1)-1) * this.voxsize(1)- trans(1,4));
                    % transformation matrix transfroms from voxel to RPI- world,
                    % however, we want LPI- world
                   % coordinates
            end
            worldCoordinates = trans*[voxelIdxList-1; ones(1, size(voxelIdxList,2))];
            worldCoordinates = worldCoordinates(1:3,:);
        end

        function [voxelIdxList,voxelToNiftiWorldMatrix] = ...
                getMatlabIdxFromNiftiWorldCoordinates(this, worldCoordList,...
                orientation)
            trans = this.transformationMatrix;
            if ( ~exist('orientation', 'var') )
                if(det(trans) < 0)
                    orientation = 'right-handed';
                else
                    orientation = 'left-handed';
                end
            end

            if ( ~isempty(worldCoordList) && size(worldCoordList,1) ~= 3)
                error('worldCoordList needs to be a 3 x N matrix');
            end
            switch orientation
                case {'left-handed'} % left-handed coordinate system
                    % nothing to do
                case {'right-handed'} % right-handed coordinate system
                    % transformation matrix transfroms from voxel to RPI-
                    % world, however, but we want LPI- world
                    % coordinates
                    trans(1,1) = -trans(1,1);
                    trans(1,4) = -((this.voxdim(1)-1) * this.voxsize(1)- trans(1,4));
            end

% 			voxelIdxList = trans\[worldCoordList; ones(1, size(worldCoordList,2))];
%             voxelIdxList = voxelIdxList(1:3,:)+1;

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
nii = NiftiModality(NiftiFilepath(niiFile));
nii.load();

idx = [1 1 1; 10 10 10; 12 17 99]';
nii.getItkWorldCoordinatesFromMatlabIdx(idx)
%%
end
