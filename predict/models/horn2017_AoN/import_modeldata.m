function import_modeldata
% This script imported the raw patient data for the 2017 AoN model. It was
% only needed to build the original model and is included for
% reproducibility purposes. It does not need to be run in general.
lead path
[pth,fn]=fileparts(which(mfilename));
cd(pth)
pts=1:95;
improvements=[-3.44827586206897;43.9024390243902;33.3333333333333;33.3333333333333;48.6486486486487;61.5384615384615;44.2307692307692;15.3846153846154;50;40.9090909090909;38.8888888888889;34.7826086956522;39.2857142857143;-7.40740740740741;53.5714285714286;34.4827586206897;62.5000000000000;75.3846153846154;25.9259259259259;-4.16666666666667;81.8181818181818;78.7234042553192;46.4285714285714;64.1509433962264;39.1304347826087;53.5714285714286;56.6037735849057;50;67.9245283018868;76.2711864406780;16.6666666666667;67.6923076923077;20.8333333333333;69.0476190476191;36.8421052631579;63.8297872340426;61.1111111111111;37.9310344827586;57.1428571428571;75;72;66.6666666666667;43.6619718309859;41.7910447761194;32.5000000000000;23.0769230769231;68.7500000000000;25.7142857142857;36.3636363636364;71.1111111111111;-11.1111111111111;53.8461538461538;47.1698113207547;65.1162790697674;38.7755102040816;8.33333333333333;41.8604651162791;46.5116279069767;27.2727272727273;32.1428571428571;28.5714285714286;69.8630136986301;43.6619718309859;70.2127659574468;-16.1290322580645;48;76.2711864406780;76.9230769230769;69.7674418604651;71.6666666666667;59.7014925373134;-10.3448275862069;57.4468085106383;22.5000000000000;54.9019607843137;80.7692307692308;63.0434782608696;48.8372093023256;70;63.7681159420290;74.5454545454545;49.1228070175439;25;-2.27272727272727;64.0625000000000;81.5789473684211;50;40;45.4545454545455;40.5405405405405;41.7910447761194;63.4615384615385;84.7826086956522;-21;58.0645161290323]; % clinical improvements, in this case percent updrs-iii
rootfolder='/Users/andreashorn/Desktop/horn2017/recos_only/'; % specify root folder for data
gs='gs_20170903150355'; % specify lead group identity (generic stimulation name)
connectomes={'HCP_MGH_30fold_groupconnectome (Horn 2017)','GSP 1000 (Yeo 2011)_Full Set (Yeo 2011)','PPMI_90 (Ewert 2017)','PPMI 74_15 (Horn 2017)_Patients'}; % specify which connectomes were run / which ones to import
types={'dMRI','fMRI','dMRI','fMRI'}; % specify types for each connectome (dMRI/fMRI)
load modeldata
options.native = 0;
for c=1:length(connectomes)
    connectome=connectomes{c};
    type=types{c};
    for pt=pts
        [~, subPrefix] = fileparts([rootfolder,num2str(pt), '_']);
        subDir = fullfile(rootfolder,num2str(pt),'stimulations',ea_nt(options),gs,connectome);
        connName = ea_getConnLabel(connectome);
        switch type
            case 'dMRI'
                fis{pt} = fullfile(subDir, [subPrefix, 'sim-binary_model-simbio_seed-dMRI_conn-', connName, '_strucmap.nii']);
            case 'fMRI'
                fis{pt} = fullfile(subDir, [subPrefix, 'sim-binary_model-simbio_seed-fMRI_conn-', connName, '_desc-AvgRFz_funcmap.nii']);
        end
    end
    mkdir([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'combined_maps',filesep,type,filesep,connectome]);
    ea_Cmap(fis,improvements,...
        [ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'combined_maps',filesep,type,filesep,connectome,filesep,type,'_optimal.nii'],...
        modeldata.mask,getsk(type));

    if strcmp(type,'dMRI')
        delete([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'combined_maps',filesep,type,filesep,connectome,filesep,'s',type,'_optimal.nii']);
    end
    modelmap=ea_load_nii([ea_getearoot,'predict',filesep,'models',filesep,'horn2017_AoN',filesep,'combined_maps',filesep,type,filesep,connectome,filesep,type,'_optimal.nii']);
    modelvals=modelmap.img(modeldata.mask);

    for pt=pts
        if strcmp(type,'dMRI')
            ptmap=ea_load_nii(ea_dosk(fis{pt},modeldata.mask));
        else
            ptmap=ea_load_nii(fis{pt});
        end
        ptvals=ptmap.img(modeldata.mask);
        infs=isinf(modelvals);
        infs=logical(infs+isinf(ptvals));
        modelvals(infs)=[];
        ptvals(infs)=[];
        sim(pt)=corr(modelvals,ptvals,'type',getcorrtype(type),'rows','pairwise');
    end

    clear fis

    % add values to modeldata:
    modeldata.connectomes.(rmbracketspace(connectome)).([type,'sims'])=sim;
    modeldata.updrs3percimprov=improvements;
    save('modeldata','modeldata');
end

function str=rmbracketspace(str)
str=strrep(str,' ','_');
str=strrep(str,'(','_');
str=strrep(str,')','_');
str=strrep(str,'>','_');

function ct=getcorrtype(type)
switch type
    case 'dMRI'
        ct='spearman'; % non-Gaussian data even if more or less normalized by default
    case 'fMRI'
        ct='pearson'; % Fisher-z-transformed Gaussian data
end

function sk=getsk(type)
switch type
    case 'dMRI'
        sk='sk';
    case 'fMRI'
        sk='';
end
