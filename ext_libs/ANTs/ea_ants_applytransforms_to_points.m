function output=ea_ants_applytransforms_to_points(varargin)
% Wrapper for antsApplyTransformsToPoints

directory=fullfile(varargin{1},filesep);
input=varargin{2};
useinverse=varargin{3};

if useinverse
    istr='Inverse';
else
    istr='';
end

options.prefs=ea_prefs(fileparts(directory));

[~,glprebase]=fileparts(options.prefs.gprenii);
[~,lprebase]=fileparts(options.prefs.prenii);
% use 'gl' affix for tranforms
try
    if exist([directory,lprebase,'Composite.h5'],'file')
        movefile([directory,lprebase,'Composite.h5'],[directory,glprebase,'Composite.h5']);
        movefile([directory,lprebase,'InverseComposite.h5'],[directory,glprebase,'InverseComposite.h5']);
    end
end
try
    if exist([directory,lprebase,'0GenericAffine.mat'],'file')
        movefile([directory,lprebase,'0GenericAffine.mat'],[directory,glprebase,'0GenericAffine.mat']);
        try movefile([directory,lprebase,'1Warp.nii.gz'],[directory,glprebase,'1Warp.nii.gz']); end
        try movefile([directory,lprebase,'1InverseWarp.nii.gz'],[directory,glprebase,'1InverseWarp.nii.gz']); end
    end
end

if nargin>3
    transform=varargin{4};
    tstring=[' --transform [',transform, ',',num2str(useinverse),']']; % [transformFileName,useInverse]
else
    if useinverse
        if exist([directory,glprebase,'Composite.h5'],'file')
            tstring=[' -t [',ea_path_helper([directory,glprebase,istr,'Composite.h5']),',0]'];
        else
            tstring=[' -t [',ea_path_helper([directory,glprebase,'0GenericAffine.mat']),',',num2str(useinverse),']',...
                ' -t [',ea_path_helper([directory,glprebase,'1',istr,'Warp.nii.gz']),',0]',...
                ];
        end

    else
        if exist([directory,glprebase,'Composite.h5'],'file')
            tstring=[' -t [',ea_path_helper([directory,glprebase,istr,'Composite.h5']),',0]'];
        else
            tstring=[' -t [',ea_path_helper([directory,glprebase,'1',istr,'Warp.nii.gz']),',0]',...
                ' -t [',ea_path_helper([directory,glprebase,'0GenericAffine.mat']),',',num2str(useinverse),']'...
                ];
        end
    end
end

ea_libs_helper;

basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    applyTransformsToPoints = ea_path_helper([basedir, 'antsApplyTransformsToPoints.exe']);
else
    applyTransformsToPoints = [basedir, 'antsApplyTransformsToPoints.', computer('arch')];
end

    uuid=ea_generate_uuid;

cmd = [applyTransformsToPoints, ...
    ' --dimensionality 3' ...   % dimensionality
    ' --precision 1' ...    % double precision
    ' --input ', ea_path_helper([directory,'tmpin_',uuid,'.csv']) ...  % input csv file with x,y,z,t (at least) as the column header
    ' --output ', ea_path_helper([directory,'tmpout_',uuid,'.csv']) ...    % warped output csv file
    tstring];

ea_writecsv([directory,'tmpin_',uuid,'.csv'],input);

if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end

output=ea_readcsv([directory,'tmpout_',uuid,'.csv']);
delete([directory,'tmpout_',uuid,'.csv']);
delete([directory,'tmpin_',uuid,'.csv']);


function coord=ea_readcsv(pth)
fid=fopen(pth);
C=textscan(fid,'%f %f %f %f','commentStyle', '#','delimiter', ',','Headerlines',1);
fclose(fid);
coord=cell2mat(C(1:3));


function ea_writecsv(pth,input)
fid=fopen(pth,'w');
try
    fprintf(fid,'x,y,z,t \n');
catch
    ea_error(['Cannot open file for writing at ',pth,'.']);
end
fprintf(fid,'%f,%f,%f,0\n',input'); %transpose needed for 'fprintf': matrix column to file row
fclose(fid);
