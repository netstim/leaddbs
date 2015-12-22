% configure file for reading DWI data (Dicom, Bruker or other formats)
% DTI Configuration file
% BATCH system (Volkmar Glauche)
%
% File created by Susanne Schnell

function read_data = dti_cfg_read

% ---------------------------------------------------------------------
% dicom Input dicom file names
% ---------------------------------------------------------------------
dicom         = cfg_files;
dicom.tag     = 'dicom';
dicom.name    = 'Dicom Images';
dicom.help    = {'These are the images that will be needed for tensor calculation.'
    'Choose more than one directory if you want to average the data (this is only possible if the datasets are acquired with identical sequence parameters)!'};
dicom.filter  = 'dir';
dicom.num     = [1 Inf];

% ---------------------------------------------------------------------
% exdicom Dicom conversion
% ---------------------------------------------------------------------
exdicom       = cfg_exbranch;
exdicom.tag   = 'exdicom';
exdicom.name  = 'Import Dicom Data';
exdicom.help  = {
        'The tools are tested for Siemens Dicom VB12 and VB13. All other Dicom files can generall be read in, but diffusion specific information might be provided manually.'
        };
exdicom.val   = {dicom};
exdicom.prog  = @(job)dti_readdata_ui('dicom',job);
exdicom.vout  = @vout;

% ---------------------------------------------------------------------
% bruker Input bruker file names
% ---------------------------------------------------------------------
bruker         = cfg_files;
bruker.tag     = 'bruker';
bruker.name    = 'Bruker Images';
bruker.help    = {'This is the directory with the images that will be needed for tensor calculation.'};
bruker.filter  = 'dir';
bruker.num     = [1 1];

% ---------------------------------------------------------------------
% exbruker Bruker conversion
% ---------------------------------------------------------------------
exbruker       = cfg_exbranch;
exbruker.tag   = 'exbruker';
exbruker.name  = 'Import Bruker Data';
exbruker.val   = {bruker};
exbruker.prog  = @(job)dti_readdata_ui('bruker',job);
exbruker.vout  = @vout;

% ---------------------------------------------------------------------
% binary Input binary file name
% ---------------------------------------------------------------------
binary         = cfg_files;
binary.tag     = 'binary';
binary.name    = 'Binary Data File';
binary.help    = {'These are the images that will be needed for tensor calculation.'};
binary.filter  = 'any';
binary.ufilter = '\.bin$';
binary.num     = [1 1];

% ---------------------------------------------------------------------
% exbinary Binary conversion
% ---------------------------------------------------------------------
exbinary       = cfg_exbranch;
exbinary.tag   = 'exbinary';
exbinary.name  = 'Import Binary Data';
exbinary.help  = {'Binary is an inhouse format, which consists of two files: "_raw.bin" and "_info.mat" (see examples and manual).'
    'mrstruct is also an inhouse format, but can be created by the user (see manual).'};
exbinary.val   = {binary};
exbinary.prog  = @(job)dti_readdata_ui('binary',job);
exbinary.vout  = @vout;

% ---------------------------------------------------------------------
% descheme Input File containg diffusion encoding directions
% ---------------------------------------------------------------------
descheme         = cfg_files;
descheme.tag     = 'descheme';
descheme.name    = 'DEscheme';
descheme.help    = {'Please chose a text file (m-file or mat-file also possible) with diffusion encoding directions.'
    'The scans without diffusion weighting (b = 0) should be in the beginning.'};
descheme.filter  = 'any';
descheme.ufilter = '\.(txt)|(mat)|m$';
descheme.num     = [1 1];

% ---------------------------------------------------------------------
% nob0s The amount of b=0 scans
% ---------------------------------------------------------------------
nob0s         = cfg_entry;
nob0s.tag     = 'nob0s';
nob0s.name    = 'nob0s';
nob0s.val     = {1};
nob0s.help    = {'Enter the amount of scans without diffusion weighting (b = 0).'};
nob0s.strtype = 'e';
nob0s.num     = [1 1];



% ---------------------------------------------------------------------
% threshold
% ---------------------------------------------------------------------
threshold         = cfg_entry;
threshold.tag     = 'threshold';
threshold.name    = 'threshold';
threshold.val     = {40};
threshold.help    = {'Enter the threshold where the fit is performed.'};
threshold.strtype = 'e';
threshold.num     = [1 1];


% ---------------------------------------------------------------------
% bvalue the b-weighting 
% ---------------------------------------------------------------------
bvalue         = cfg_entry;
bvalue.tag     = 'bvalue';
bvalue.name    = 'bvalue';
bvalue.val     = {1000};
bvalue.help    = {'Enter the bvalue other than b = 0.'};
bvalue.strtype = 'e';
bvalue.num     = [1 Inf];

% ---------------------------------------------------------------------
% volname: name of appended volume structure
% ---------------------------------------------------------------------
volname        = cfg_entry;
volname.tag     = 'volname';
volname.name    = 'volume structure name';
volname.val     = {'mrStruct'};
volname.help    = {'Enter the name of the mrstruct to be added'};
volname.strtype = 's';
volname.num     = [1 Inf];

% ---------------------------------------------------------------------
% dtdname for append volume
% ---------------------------------------------------------------------
dtdname         = cfg_files;
dtdname.tag     = 'dtdname';
dtdname.name    = 'DTD filename';
dtdname.help    = {'Select the dtd-mat file for appending mrstruct.'};
dtdname.filter  = '_DTD.mat';
dtdname.ufilter = '.*';
dtdname.num     = [1 1];

% ---------------------------------------------------------------------
% filename of mrstruct for append volume
% ---------------------------------------------------------------------
filename2         = cfg_files;
filename2.tag     = 'filename2';
filename2.name    = 'Mrstruct Data File';
filename2.help    = {'Select the mrstruct for appending to DTD-mat file.'};
filename2.filter  = 'mat';
filename2.ufilter = '.*';
filename2.num     = [1 1];

% ---------------------------------------------------------------------
% mrstruct Input mrstruct file name
% ---------------------------------------------------------------------
filename         = cfg_files;
filename.tag     = 'filename';
filename.name    = 'Mrstruct Data File';
filename.help    = {'These are the images saved as mrstruct that will be needed for tensor calculation.'};
filename.filter  = 'mat';
filename.ufilter = '.*';
filename.num     = [1 1];

% ---------------------------------------------------------------------
% mrstruct Input mrstruct file name
% ---------------------------------------------------------------------
filenameHARDI         = cfg_files;
filenameHARDI.tag     = 'filename';
filenameHARDI.name    = 'mrStruct Data File containg HARDI information';
filenameHARDI.help    = {'These are the images that will be needed for tensor calculation.'};
filenameHARDI.filter  = 'mat';
filenameHARDI.ufilter = '.*';
filenameHARDI.num     = [1 1];


% ---------------------------------------------------------------------
% append volume to DTD (mrstruct)
% ---------------------------------------------------------------------
append         = cfg_exbranch;
append.tag     = 'mrstruct';
append.name    = 'append volume to DTD';
append.help    = {'Select an mrstruct of the same volume size for adding and saving into the dtd-mat file.'};
append.val  = {dtdname filename2 volname};
append.prog = @(job)dti_readdata_ui('append',job);
append.vout = @vout;

% ---------------------------------------------------------------------
% mrstruct Input mrstruct file name
% ---------------------------------------------------------------------
mrstruct         = cfg_exbranch;
mrstruct.tag     = 'mrstruct';
mrstruct.name    = 'Mrstruct Depending Information';
mrstruct.help    = {'Select mrstruct and additional diffusion specific information.', ...
        'For the matlab structure (mrstruct) or binary file (two files needed: "_raw.bin" and "_info.mat").'};
mrstruct.val  = {filename descheme bvalue nob0s};
mrstruct.prog = @(job)dti_readdata_ui('mrstruct',job);
mrstruct.vout = @vout;
% ---------------------------------------------------------------------
% compute DTD from HARDI
% ---------------------------------------------------------------------


fnameDTD         = cfg_entry;
fnameDTD.tag     = 'fnameDTD';
fnameDTD.name    = 'File Name of the resulting DTD';
fnameDTD.help    = {'Type in the name of the new file containing the tensor.'};
fnameDTD.strtype = 's';
fnameDTD.num     = [1 Inf];

dir      = cfg_files;
dir.tag  = 'dir';
dir.name = 'Output directory';
dir.help = {'Select the output directory.'};
dir.filter = 'dir';
dir.num  = [1 1];



kyes      = cfg_const;
kyes.tag  = 'kyes';
kyes.name = 'Yes';
kyes.val  = {1};
kno     = cfg_const;
kno.tag  = 'kno';
kno.name = 'No';
kno.val  = {1};

kurto        = cfg_choice;
kurto.tag    = 'kurto';
kurto.name   = 'Compute Kurtosis Parameters';
kurto.values = {kyes kno};
kurto.help   = {'If data contains multiple b-shells kurtosis parameters can be estimated'};



out      = cfg_branch;
out.tag  = 'out';
out.name = 'User-specified output location';
out.val  = {dir fnameDTD};
out.help = {'Specify output directory and filename.'};

auto      = cfg_const;
auto.tag  = 'auto';
auto.name = 'Automatically determine output filename';
auto.val  = {1};

newfileDTD        = cfg_choice;
newfileDTD.tag    = 'newfile';
newfileDTD.name   = 'Select output location';
newfileDTD.values = {auto out};
newfileDTD.help   = {'Name and location can be specified explicitly, or they can be derived from the input filename.'};


computeDTD         = cfg_exbranch;
computeDTD.tag     = 'computeDTD';
computeDTD.name    = 'compute DTD from HARDI';
computeDTD.help    = {'Select an mrstruct (HARDI) for tensor estimation.'};
computeDTD.val  = { filenameHARDI kurto threshold newfileDTD};
computeDTD.prog = @(job)dti_readdata_ui('computeDTDfromHARDI',job);
computeDTD.vout = @vout;

% ---------------------------------------------------------------------
% read_data 
% ---------------------------------------------------------------------
read_data         = cfg_choice;
read_data.tag     = 'ReadData';
read_data.name    = 'Data Import';
read_data.values  = {exdicom exbruker exbinary mrstruct append computeDTD};
read_data.help    = {'Read Dicom, Bruker, Matlab structure (mrstruct) or binary file DWI data.'};




% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
function dep = vout(job)
dep            = cfg_dep;
dep.sname      = 'DWI Data';
dep.src_output = substruct('.','files');
dep.tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});
