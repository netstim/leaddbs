clear
clc
%% autocoord batches

patdirs={'BrandtH','DettmerW','DreissigD','FleischmannN','FrommerToni','GadM','GenschmerGer','HartmannDiet','HolwasJ','HornH','KeilMarco','Kerber','KieslichH','KnausOlga','KohrtW','KubischUrsul','LachmannG','LangeWalter','LehmannB','MildeJoachim','NarogluBaki','NickelH','OttmarManfre','PrahstMirjam','PuchertBrigi'};
patdirs={'RegenbergIng','ReiskiMarian','Reusche','RoloffBirgit','SchabackerLu','SchatzMartin','SchiekeI','SchlenzkaChr','SchneiderW','SchoenfeldG','TurkerHasan','WitkowskaJus','ZahnJ','ZwickC'};
%patdirs={'GenschmerGer','HartmannDiet','HolwasJ','HornH','KeilMarco','Kerber','KieslichH','KnausOlga','KohrtW','KubischUrsul','LachmannG','LangeWalter','LehmannB','MildeJoachim','NarogluBaki','NickelH','OttmarManfre','PrahstMirjam','PuchertBrigi'};

patdirs={'RadkeEdwin'};

%% set options
options.root=detroot_andy; %'/Volumes/EspionageMounts/andreashorn/1065433271/bg/out/';
options.normalize=0;

options.endtolerance=10; % how many slices to use with zero signal until end of electrode estimate.
options.sprungwert=4; % how far electrode centroid may be (in xy axis) from last to current slice.
options.refinesteps=0; % how often to re-iterate to reconstruct trajectory. More than 2 should usually not be beneficial. Use 0 to use the direct measurement.
options.tra_stdfactor=0.9;%0.9; % Default: 2.0 - the lower this factor, the lower the threshold (more included pixels in tra process).
options.cor_stdfactor=1.0; % Default: 1.0 - the higher this factor, the lower the threshold (more included pixels in cor process).
options.verbose=3; % 4: Show figures but close them 3: Show all but close all figs except resultfig 2: Show all and leave figs open, 1: Show displays only, 0: Show no feedback. 
options.sides=[1:2]; %side=1 -> left electrode, side=2 -> right electrode. both: [1:2]
options.maskwindow=15; % size of the window that follows the trajectory
options.slow=0; % if true, there will be some pauses at critical points so that the process can be better visualized. Mainly for demonstration or debugging problems.
options.axiscontrast=10; % if 8: use tra only but smooth it before. % if 9: use mean of cor and tra but smooth it. % if 10: use raw tra only.
options.zheights=4; % if 1: use cor only, 2: use smoothed version of cor only, if 3: use mean of cor and tra, if 4: use multiplied cor * tra 5: use ^10 version of cor, 6: use ^10 version of cor.*tra 7: smoothed and then like 5. -> to determine heights of electrode contacts.
options.zresolution=10; % voxels are being parcellated into this amount of portions.
options.targetknown=0; % if 0, prior of stn will be used. Set priorstrength to 0 to not use any prior at all.
options.target='gpi';
options.priorstrength=0.1; % how strong the prior information will be weighted. Set to 0 to not use prior at all.
options.showatlases=0;
options.showfibres=0;
options.writeoutimages=0;
options.overlays={'edge',   'gp_h',     'gpi',  'morel'};
options.ocolors={'white',   'grey',    'red',   'multi'};       % use off or multi, red, blue, green, yellow, magenta, cyan, white, grey and black only.
options.oopacities=[1,    .2,         1,  	1];                 % use range 0-1
options.manualheightcorrection=1;


options.eldist=3; % distance of electrode contacts. 3 mm for large (VIM) electrodes, 2 mm for small (STN/GPi) electrodes.



options=defaultoptions(options);

%% run analysis
close all
dispbn;
for pat=1:length(patdirs)
    try
        results=autocoord_v2(patdirs{pat},options);
        
        try
            patfit(pat)=results.fit;
        end
        write_fiducials(results.coords_mm,fullfile(options.root,patdirs{pat},'L_auto.fcsv'))
        
        try
            save([options.root,patdirs{pat},filesep,'automanual_stats'],'results');
        end
    catch
        disp([patdirs{pat},' failed.']);
    end
    % disp([patdirs{pat},' failed.']);
    
end


%disp(['Mean fit was ',num2str(mean(patfit)),'.']);