clear
clc
%% autocoord batches

patdirs={'SchoenfeldG'}; %,'KohrtW','ReiskiMarian','LangeWalter','SchiekeI','DreissigD','Fintel','GadM','KeilMarco','RoloffBirgit','NeubauerOswa','TeschendorfH','SchneiderW','MindingB','KrummowH','KrauseH','GrefrathG','BrandtH','MascheskiAgn','ZahnJ','NickelH','HolwasJ','ZwickC','TheelkeW'};


%% set options
options.root='/Volumes/EspionageMounts/andreashorn/1065433271/bg/out/';

options.endtolerance=1; % how many slices to use with zero signal until end of electrode estimate.
options.sprungwert=4; % how far electrode centroid may be (in xy axis) from last to current slice.
options.tra_stdfactor=2.2; % Default: 2.0 - the lower this factor, the lower the threshold (more included pixels in tra process).
options.cor_stdfactor=1.5; % Default: 1.0 - the higher this factor, the lower the threshold (more included pixels in cor process).
options.verbose=0; % 2: Show figures and displays, 1: Show displays only, 0: Show no feedback. 

eldist=1.7;


%% run analysis

for pat=1:length(patdirs)
cnt=1;
    for cor_stdfactor=1:0.2:3
    for tra_stdfactor=1:0.2:3
    
        options.cor_stdfactor=cor_stdfactor;
        options.tra_stdfactor=tra_stdfactor;
        try
        results=autocoord_v2(patdirs{pat},eldist,options);
        
        resvec(cnt)=[cor_stdfactor,tra_stdfactor,sum(results.distances),sum(results.optimdistances)];
        display(num2str(resvec(cnt)));
        cnt=cnt+1;
        catch
            display(['Electrode estimation did not work out with cor_stdfactor=',num2str(cor_stdfactor),' and tra_stdfactor=',num2str(tra_stdfactor),'.']);
        end
    end
    end
    save([patdirs{pat},'_resvec'],'resvec');
end


