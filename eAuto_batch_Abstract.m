clear
clc
%% autocoord batches

patdirs={'FleischmannN','WitkowskaJus','KieslichH','LachmannG','Kerber','SchoenfeldG','DreissigD','ReiskiMarian','LangeWalter','SchabackerLu','SchatzMartin'};
%patdirs={'NeubauerOswa'}%,'TheelkeW', 'BouchenafaIn','ZwickC'}; % 3mm ones

%% set options
options.root='/Volumes/EspionageMounts/andreashorn/1065433271/bg/out/';

options.endtolerance=3; % how many slices to use with zero signal until end of electrode estimate.
options.sprungwert=4; % how far electrode centroid may be (in xy axis) from last to current slice.
options.tra_stdfactor=0.9; % Default: 2.0 - the lower this factor, the lower the threshold (more included pixels in tra process).
options.cor_stdfactor=1.0; % Default: 1.0 - the higher this factor, the lower the threshold (more included pixels in cor process).
options.verbose=2; % 4: Show figures but close them 3: Show all but close all figs except resultfig 2: Show all and leave figs open, 1: Show displays only, 0: Show no feedback. 
options.sides=[1:2]; %side=1 -> left electrode, side=2 -> right electrode. both: [1:2]
options.maskwindow=10;
options.slow=0; % if true, there will be some pauses at critical points so that the process can be better visualized. Mainly for demonstration or debugging problems.
options.zheights=6; % if 1: use cor only, if 2: use mean of cor and tra, if 3: use multiplied cor * tra 4: use ^10 version of cor, 5: use ^10 version of cor.*tra 6: smoothed and then like 5. -> to determine heights of electrode contacts.
options.targetknown=0; % don't use this for now.
options.target='stn';
options.overlays={'edge',   'gp_h',     'gpi',  'morel'};
options.ocolors={'white',   'grey',    'red',   'multi'};       % use off or multi, red, blue, green, yellow, magenta, cyan, white, grey and black only.
options.oopacities=[1,    .2,         1,  	1];                 % use range 0-1
options.manualheightcorrection=1;


eldist=3; % distance of electrode contacts. 3 mm for large (VIM) electrodes, 2 mm for small (STN/GPi) electrodes.


%% run analysis

for pat=1:length(patdirs)
try
    results=autocoord_v2(patdirs{pat},eldist,options);
    patfit(pat)=results.fit;
    
    write_fiducials(results.coords_mm,fullfile(options.root,patdirs{pat},'L_auto.fcsv'))

    
    save([options.root,patdirs{pat},filesep,'automanual_stats'],'results');
catch
    disp([patdirs{pat},' failed.']);
end
end
 
%disp(['Mean fit was ',num2str(mean(patfit)),'.']);