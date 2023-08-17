function resultfig=ea_showgroup(M,options,ptidx,planes,structures)
clc;
% set options
options.leadprod = 'group';
% set pt specific options
options.root=M.root;
options.patientname='';

options.expstatvat.do=M.ui.statvat;
options.native=0;
try
    options.numcontacts=size(M.elstruct(1).coords_mm{1},1);
catch
    warning('Localizations seem not properly defined.');
end
options.elmodel=M.elstruct(1).elmodel;
 options=ea_resolve_elspec(options);
options.prefs=ea_prefs(options.patientname);
options.d3.verbose='on';

options.d3.elrendering=M.ui.elrendering;

try options.d3.isomatrix=M.isomatrix; end
try options.d3.isomatrix_name=M.isomatrix_name; end

options.d2.write=0;

options.d2.atlasopacity=0.15;

options.d3.isovscloud=M.ui.isovscloudpopup;
options.d3.showisovolume=M.ui.showisovolumecheck;

options.d3.colorpointcloud=M.ui.colorpointcloudcheck;
options.d3.exportBB=0;	% don't export brainbrowser struct by default

options.normregressor=M.ui.normregpopup;

% Prepare isomatrix (includes a normalization step if M.ui.normregpopup
% says so:

for reg=1:length(options.d3.isomatrix)
    try options.d3.isomatrix{reg}=ea_reformat_isomatrix(options.d3.isomatrix{reg},M,options); end
end



options.groupmode=1;

% overwrite active contacts information with new one from S (if present).
try
    for pt=1:length(M.elstruct)
        M.elstruct(pt).activecontacts=M.S(pt).activecontacts;
    end
end

try
    for pt=1:length(M.elstruct)
        M.elstruct(pt).groupcolors=M.groups.color;
    end
end
options.groupmode=1;
options.modality=3; % use template image
options.patient_list=M.patient.list;

% mer development
if ~exist('ptidx','var')
    ptidx=1:length(M.patient.list);
end
vizstruct.elstruct=M.elstruct(ptidx);
uipatdirs=M.patient.list(ptidx);
npts=length(uipatdirs);


% amend .pt to identify which patient is selected (needed for isomatrix).
for pt=1:length(ptidx)
    M.elstruct(ptidx(pt)).pt=ptidx(pt);
end

resultfig=ea_elvis(options,M.elstruct(ptidx));
ea_keepatlaslabels(structures{:})
ea_setplanes(planes.x,planes.y,planes.z);






