function output=ea_ants_applytransforms_to_points(varargin)
% Wrapper for antsApplyTransformsToPoints

subdir=varargin{1};
input=varargin{2};
useinverse=varargin{3};
if useinverse
    istr='Inverse';
else
    istr='';
end

[~,ptname]=fileparts(subdir);
options.prefs=ea_prefs(ptname);
[~,lprebase]=fileparts(options.prefs.prenii);
if nargin>3
    transform=varargin{4};
    tstring=[' --transform [',transform, ',',num2str(useinverse),']']; % [transformFileName,useInverse]
else
    if useinverse
        if exist([subdir,lprebase,'Composite.h5'],'file')
            tstring=[' -t [',ea_path_helper([subdir,lprebase]),istr,'Composite.h5,0]'];
        else
            tstring=    [  ' -t [',ea_path_helper([subdir,lprebase]),'0GenericAffine.mat,',num2str(useinverse),']',...
                ' -t [',ea_path_helper([subdir,lprebase]),'1',istr,'Warp.nii.gz,',num2str(0),']',...
                ];
        end
        
    else
        if exist([subdir,lprebase,'Composite.h5'],'file')
            tstring=[' -t [',ea_path_helper([subdir,lprebase]),istr,'Composite.h5,0]'];
        else
            tstring=[' -t [',ea_path_helper([subdir,lprebase]),'1',istr,'Warp.nii.gz,',num2str(0),']',...
                ' -t [',ea_path_helper([subdir,lprebase]),'0GenericAffine.mat,',num2str(useinverse),']'...
                ];
        end
    end
end

ea_libs_helper;

basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    applyTransformsToPoints = [basedir, 'antsApplyTransformsToPoints.exe'];
elseif isunix
    applyTransformsToPoints = [basedir, 'antsApplyTransformsToPoints.', computer];
end


cmd = [applyTransformsToPoints, ...
    ' --dimensionality 3' ...   % dimensionality
    ' --precision 1' ...    % double precision
    ' --input ', ea_path_helper([subdir,'tmpin.csv']) ...  % input csv file with x,y,z,t (at least) as the column header
    ' --output ', ea_path_helper([subdir,'tmpout.csv']) ...    % warped output csv file
tstring];
    
ea_writecsv([subdir,'tmpin.csv'],input);

if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end

output=ea_readcsv([subdir,'tmpout.csv']);
delete([subdir,'tmpout.csv']);
delete([subdir,'tmpin.csv']);




function c=ea_readcsv(pth)
fid=fopen(pth);
C=textscan(fid,'%f %f %f %f','commentStyle', '#','delimiter', ',','Headerlines',1);
fclose(fid);
for coord=1:length(C{1})
   c(:,coord)=[C{1}(coord);C{2}(coord);C{3}(coord);1];
end


function c=ea_writecsv(pth,input)

fid=fopen(pth,'w');
fprintf(fid,'x,y,z,t \n');
for c=1:size(input,2)
   fprintf(fid,[num2str(input(1,c)),',',num2str(input(2,c)),',',num2str(input(3,c)),',0 \n']); 
end
fclose(fid);
