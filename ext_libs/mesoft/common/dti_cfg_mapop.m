% configure file for calculation of tensors
% DTI Configuration file
% BATCH system (Volkmar Glauche)
%
% File created by Susanne Schnell


function mapop = dti_cfg_mapop

% ---------------------------------------------------------------------
% filename1 name of probabilistic map 1
% ---------------------------------------------------------------------
filename1         = cfg_files;
filename1.tag     = 'filename1';
filename1.name    = 'Select first Map';
filename1.help    = {'help..'};
filename1.filter  = 'mat';
filename1.ufilter = '.*';
filename1.num     = [1 1];

% ---------------------------------------------------------------------
% filename2 name of probabilistic map 2
% ---------------------------------------------------------------------
filename2         = cfg_files;
filename2.tag     = 'filename2';
filename2.name    = 'Select second Map';
filename2.help    = {'help..'};
filename2.filter  = 'mat';
filename2.ufilter = '.*';
filename2.num     = [1 1];

% ---------------------------------------------------------------------
% filenames names of probabilistic maps 1
% ---------------------------------------------------------------------
filenames         = cfg_files;
filenames.tag     = 'filenames';
filenames.name    = 'Select List of Maps';
filenames.help    = {'help..'};
filenames.filter  = 'mat';
filenames.ufilter = '.*';
filenames.num     = [2 inf];

% ---------------------------------------------------------------------
% newfilename name of resulting file
% ---------------------------------------------------------------------
newfilename         = cfg_entry;
newfilename.tag     = 'newfilename';
newfilename.name    = 'New File Name';
newfilename.val     = {'.mat'};
newfilename.help    = {'Type in the name of the new file, if you leave this empty a default file name is generated with the endings "_normMAP.mat", "_addMAP.mat" or "_multMAP.mat"'};
newfilename.strtype = 's';
newfilename.num     = [1 Inf];

% ---------------------------------------------------------------------
% connmtx connection matrix for multiple multiplications
% ---------------------------------------------------------------------
connmtx         = cfg_entry;
connmtx.tag     = 'connmtx';
connmtx.name    = 'Connection matrix';
connmtx.help    = {['Enter the indices of the maps you want to multiply in ' ...
                    'a #connections-by-2 array. The output filenames will ' ...
                    'be determined automatically from the inputs.']};
connmtx.strtype = 'n';
connmtx.num     = [Inf 2];

% ---------------------------------------------------------------------
% addition input of two maps
% ---------------------------------------------------------------------
addition         = cfg_exbranch;
addition.tag     = 'addition';
addition.name    = 'Addition';
addition.help    = {'Select the two maps you want to add.'};
addition.val     = {filename1 filename2 newfilename};
addition.prog    = @(job)dti_mapop_ui('addition',job);
addition.vout    = @vout;

% ---------------------------------------------------------------------
% multiplication input of two maps
% ---------------------------------------------------------------------
multiplication         = cfg_exbranch;
multiplication.tag     = 'multiplication';
multiplication.name    = 'Multiplication (2 Maps)';
multiplication.help    = {'Select the two maps you want to multiply. Multiplaction results depend on the probabilistic tracking method chosen. If the extended version was used you will have three maps: the merging fibre map, the connecting fibre map and all-in-one map.'};
multiplication.val     = {filename1 filename2 newfilename};
multiplication.prog    = @(job)dti_mapop_ui('multiplication',job);
multiplication.vout    = @vout;

% ---------------------------------------------------------------------
% connmulti input of map list and connection matrix
% ---------------------------------------------------------------------
connmulti         = cfg_exbranch;
connmulti.tag     = 'connmulti';
connmulti.name    = 'Multiplication (List of Maps)';
connmulti.help    = {['Enter a list of maps and a #maps-by-2 matrix ' ...
                    'indicating which map pairs you want to multiply. Multiplaction results depend on the probabilistic tracking method chosen. If the extended version was used you will have three maps: the merging fibre map, the connecting fibre map and all-in-one map.']};
connmulti.val     = {filenames connmtx};
connmulti.prog    = @(job)dti_mapop_ui('connmulti',job);
connmulti.check   = @(job)dti_mapop_ui('connmulti_check',job);
connmulti.vout    = @vout;

% ---------------------------------------------------------------------
% normalize input of one map and the new filename
% ---------------------------------------------------------------------
normalize         = cfg_exbranch;
normalize.tag     = 'normalize';
normalize.name    = 'Normalization';
normalize.help    = {'Select the map to normalize and chose a filename for the resulting map'};
normalize.val     = {filename1 newfilename};
normalize.prog    = @(job)dti_mapop_ui('normalize',job);
normalize.vout    = @vout;

% ---------------------------------------------------------------------
% mapop 
% ---------------------------------------------------------------------
mapop         = cfg_choice;
mapop.tag     = 'mapop';
mapop.name    = 'Operations with Probability Maps';
mapop.help    = {'Mathematical operations with probability maps.'};
mapop.values  = {normalize addition multiplication connmulti};

% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
function dep = vout(job)
dep            = cfg_dep;
dep.sname      = 'ProbMap After Operation';
dep.src_output = substruct('.','files');
dep.tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});