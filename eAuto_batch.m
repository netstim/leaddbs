clear
clc
%% autocoord batches

patdirs={'MildeJoachim','KieslichH','LangeWalter','DreissigD','RegenbergIng','LehmannB','LachmannG','ReiskiMarian'};
patdirs={'LachmannG','DreissigD'};
%patdirs={'MildeJoachim'};




%% set options
options.root='/Volumes/EspionageMounts/andreashorn/1065433271/bg/out/';
options.normalize=0;

options.endtolerance=3; % how many slices to use with zero signal until end of electrode estimate.
options.sprungwert=4; % how far electrode centroid may be (in xy axis) from last to current slice.
options.tra_stdfactor=0.9; % Default: 2.0 - the lower this factor, the lower the threshold (more included pixels in tra process).
options.cor_stdfactor=1.0; % Default: 1.0 - the higher this factor, the lower the threshold (more included pixels in cor process).
options.verbose=3; % 4: Show figures but close them 3: Show all but close all figs except resultfig 2: Show all and leave figs open, 1: Show displays only, 0: Show no feedback. 
options.sides=[1:2]; %side=1 -> left electrode, side=2 -> right electrode. both: [1:2]
options.maskwindow=10;
options.slow=0; % if true, there will be some pauses at critical points so that the process can be better visualized. Mainly for demonstration or debugging problems.
options.zheights=2; % if 1: use cor only, 2: use smoothed version of cor only, if 3: use mean of cor and tra, if 4: use multiplied cor * tra 5: use ^10 version of cor, 6: use ^10 version of cor.*tra 7: smoothed and then like 5. -> to determine heights of electrode contacts.
options.zresolution=10; % voxels are being parcellated into this amount of portions.
options.targetknown=1; % if 0, prior of stn will be used. Set priorstrength to 0 to not use any prior at all.
options.target='stn';
options.priorstrength=0.1; % how strong the prior information will be weighted. Set to 0 to not use prior at all.
options.writeoutimages=0;
options.overlays={'edge',   'gp_h',     'gpi',  'morel'};
options.ocolors={'white',   'grey',    'red',   'multi'};       % use off or multi, red, blue, green, yellow, magenta, cyan, white, grey and black only.
options.oopacities=[1,    .2,         1,  	1];                 % use range 0-1
options.manualheightcorrection=0;


eldist=2; % distance of electrode contacts. 3 mm for large (VIM) electrodes, 2 mm for small (STN/GPi) electrodes.


%% run analysis
close all
dispbn;
for pat=1:length(patdirs)

    results=autocoord_v2(patdirs{pat},eldist,options);
try
    patfit(pat)=results.fit;
end 
    write_fiducials(results.coords_mm,fullfile(options.root,patdirs{pat},'L_auto.fcsv'))

    
    save([options.root,patdirs{pat},filesep,'automanual_stats'],'results');
  % disp([patdirs{pat},' failed.']);
 
end
 

%disp(['Mean fit was ',num2str(mean(patfit)),'.']);