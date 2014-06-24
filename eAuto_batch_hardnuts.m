clear
clc
%% autocoord batches

patdirs={'NickelH','HolwasJ','ZwickC','TheelkeW'}; 
%,'SchoenfeldG','KohrtW','ReiskiMarian','LangeWalter','SchiekeI','DreissigD','Fintel','GadM','KeilMarco','RoloffBirgit','NeubauerOswa','TeschendorfH','SchneiderW','MindingB','KrummowH','KrauseH','GrefrathG','BrandtH','MascheskiAgn','ZahnJ','NickelH','HolwasJ','ZwickC','TheelkeW'}; %{'SchoenfeldG','KohrtW','ReiskiMarian','LangeWalter','SchiekeI','DreissigD','Fintel','GadM','KeilMarco','RoloffBirgit','NeubauerOswa','TeschendorfH','SchneiderW','MindingB','KrummowH','KrauseH','GrefrathG','BrandtH','MascheskiAgn','ZahnJ','NickelH'}; %,'SchoenfeldG','KohrtW','ReiskiMarian','LangeWalter','SchiekeI','DreissigD','Fintel','GadM','KeilMarco','RoloffBirgit','NeubauerOswa','TeschendorfH','SchneiderW','MindingB','KrummowH','KrauseH','GrefrathG','BrandtH','MascheskiAgn','ZahnJ','NickelH','HolwasJ','ZwickC','TheelkeW'};
eldist=[2,2,2,2,2,2,2,2,2]; % distance of electrode contacts. 3 mm for large (VIM) electrodes, 2 mm for small (STN/GPi) electrodes.

%% set options
options.root='/Volumes/EspionageMounts/andreashorn/1065433271/bg/out/';

options.endtolerance=3; % how many slices to use with zero signal until end of electrode estimate.
options.sprungwert=4; % how far electrode centroid may be (in xy axis) from last to current slice.
options.tra_stdfactor=0.9; % Default: 2.0 - the lower this factor, the lower the threshold (more included pixels in tra process).
options.cor_stdfactor=1.0; % Default: 1.0 - the higher this factor, the lower the threshold (more included pixels in cor process).
options.verbose=3; % 4: Show figures but close them 3: Show all but close all figs except resultfig 2: Show all and leave figs open, 1: Show displays only, 0: Show no feedback. 
options.sides=[1:2]; %side=1 -> right electrode, side=2 -> left electrode. both: [1:2]
options.maskwindow=12;
options.slow=0; % if true, there will be some pauses at critical points so that the process can be better visualized. Mainly for demonstration or debugging problems.




%% run analysis

for pat=1:length(patdirs)

    results=autocoord_v2(patdirs{pat},eldist(pat),options);
    patfit(pat)=results.fit;
    
    save([options.root,patdirs{pat},filesep,'automanual_stats'],'results');
end
 
showdis(['Mean fit was ',num2str(mean(patfit)),'.'],options.verbose);