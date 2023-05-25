function X=ea_genfeaturesfromvta(uipatdirs, stimname, atlname)

%ea_run_mapper(uipatdirs, stimname);
options.prefs=ea_prefs('');
if exist('atlname','var')

    copyfile([ea_space([],'labeling'),atlname,'.nii'],[tempdir,'label.nii']);
    ea_conformspaceto(fullfile(ea_getearoot,'predict','spaces','222.nii'),[tempdir,'label.nii'],0);
    nii=ea_load_nii([tempdir,'label.nii']);
    ea_delete([tempdir,'label.nii']);
    cnt=1;
    for parc=unique(nii.img(nii.img>0))'
        allfeatsix{cnt}=find(nii.img==parc);
        cnt=cnt+1;
    end
else
    fts=load(fullfile(ea_getearoot,'predict','models','horn_fox','feature_idx.mat'));
    allfeatsix=cell2mat(fts.idx');
end

vtaType = options.prefs.lcm.vatseed

options.native=0;
for pt=1:length(uipatdirs)
    thispt=uipatdirs{pt};
    options=ea_getrootptname(thispt,options);
    coords_mm=ea_load_reconstruction(options);
    load([thispt,filesep,'stimulations',filesep,ea_nt(options),stimname,filesep,'stimparameters_right.mat']);
    for side=1:2
        acnt{side}=mean(coords_mm{side}(logical(S.activecontacts{side}),:),1);
    end

    [~, subPrefix] = fileparts([thispt, '_']);
    modelLabel = ea_simModel2Label(S.model);
    fConnName = ea_getConnLabel(options.lcm.func.connectome);
    dConnName = ea_getConnLabel(options.lcm.struc.connectome);
    funcnii=ea_load_nii(fullfile(thispt,'stimulations',ea_nt(options),stimname,[subPrefix, 'sim-', vtaType, '_model-', modelLabel, '_seed-fMRI_conn-', fConnName, '_desc-AvgRFz_funcmap.nii']));
    strucnii=ea_load_nii(fullfile(thispt,'stimulations',ea_nt(options),stimname,[subPrefix, 'sim-', vtaType, '_model-', modelLabel, '_seed-dMRI_conn-', dConnName, '_strucmap.nii']));

    % assign feature vector X:
    if iscell(allfeatsix)
       for c=1:length(allfeatsix)
          Fu(pt,c)=mean(funcnii.img(allfeatsix{c}));
          St(pt,c)=mean(strucnii.img(allfeatsix{c}));
       end
       X(pt,:)=[St(pt,:),Fu(pt,:)];

    else
        %   X(pt,:)=[cell2mat(acnt),strucnii.img(allfeatsix)',funcnii.img(allfeatsix)'];
        X(pt,:)=[strucnii.img(allfeatsix)',funcnii.img(allfeatsix)'];
    end
end


function options=ea_getrootptname(folder,options)
[options.root,options.patientname]=fileparts(folder);
options.root=[options.root,filesep];


function ea_run_mapper(uipatdirs, stimname)
% - Lead-DBS Job created on 13-Feb-2017:
% --------------------------------------

% Execute job:
% ---------------------------------------

options=getoptslocal;
options.uipatdirs = uipatdirs';
options.lcm.struc.connectome = 'HCP842_1mm (Yeh 2011)';
options.lcm.func.connectome = 'GSP 1000 Groupmatrix (Yeo 2011) > Full Set (Yeo 2011)';
options.lcm.seeds = stimname;

options.leadprod = 'mapper';

allpatdirs=options.uipatdirs;
for pat=1:length(allpatdirs)
    % set subject specific options:
    options.root=[fileparts(allpatdirs{pat}),filesep];
    [~,thispatdir]=fileparts(allpatdirs{pat});
    options.patientname=thispatdir;
    options.uipatdirs=allpatdirs{pat};
    ea_run('run',options);
end


function options=getoptslocal
options.endtolerance = 10;
options.sprungwert = 4;
options.refinesteps = 0;
options.tra_stdfactor = 0.9;
options.cor_stdfactor = 1;
options.earoot = ea_getearoot;
options.normalize.do = 0;
options.normalize.method = [];
options.checkreg = false;
options.coregmr.method = '';
options.coregct.do = 0;
options.coregctcheck = 0;
options.modality = 1;
options.verbose = 3;
options.sides = [1 2];
options.doreconstruction = 0;
options.autoimprove = 0;
options.axiscontrast = 8;
options.zresolution = 10;
options.atl.genpt = 0;
options.atl.can = 1;
options.atl.pt = 0;
options.atl.ptnative = 0;
options.native = 0;
options.d2.col_overlay = 1;
options.d2.con_overlay = 1;
options.d2.con_color = [1 1 1];
options.d2.lab_overlay = 1;
options.d2.bbsize = 25;
options.d2.backdrop = 'MNI152NLin2009bAsym T2 (Fonov 2011)';
options.d2.write = 0;
options.d2.atlasopacity = 0.15;
options.refinelocalization = 0;
options.d3.write = 0;
options.d3.prolong_electrode = 2;
options.d3.verbose = 'on';
options.d3.elrendering = 1;
options.d3.hlactivecontacts = 0;
options.d3.showactivecontacts = 1;
options.d3.showpassivecontacts = 1;
options.d3.showisovolume = 0;
options.d3.isovscloud = 0;
options.d3.mirrorsides = 0;
options.d3.autoserver = 0;
options.d3.expdf = 0;
options.numcontacts = 4;
options.writeoutpm = 0;
options.expstatvat.do = 0;
options.fiberthresh = 10;
options.writeoutstats = 1;
%%
options.colormap = [0.2081 0.1663 0.5292
                    0.21162380952381 0.189780952380952 0.57767619047619
                    0.212252380952381 0.213771428571429 0.626971428571429
                    0.2081 0.2386 0.677085714285714
                    0.195904761904762 0.264457142857143 0.7279
                    0.170728571428571 0.291938095238095 0.779247619047619
                    0.125271428571429 0.324242857142857 0.830271428571429
                    0.0591333333333334 0.359833333333333 0.868333333333333
                    0.0116952380952381 0.387509523809524 0.881957142857143
                    0.00595714285714286 0.408614285714286 0.882842857142857
                    0.0165142857142857 0.4266 0.878633333333333
                    0.032852380952381 0.443042857142857 0.871957142857143
                    0.0498142857142857 0.458571428571429 0.864057142857143
                    0.0629333333333333 0.473690476190476 0.855438095238095
                    0.0722666666666667 0.488666666666667 0.8467
                    0.0779428571428571 0.503985714285714 0.838371428571429
                    0.079347619047619 0.52002380952381 0.831180952380952
                    0.0749428571428571 0.537542857142857 0.826271428571429
                    0.0640571428571428 0.556985714285714 0.823957142857143
                    0.0487714285714286 0.577223809523809 0.822828571428571
                    0.0343428571428572 0.596580952380952 0.819852380952381
                    0.0265 0.6137 0.8135
                    0.0238904761904762 0.628661904761905 0.803761904761905
                    0.0230904761904762 0.641785714285714 0.791266666666667
                    0.0227714285714286 0.653485714285714 0.776757142857143
                    0.0266619047619048 0.664195238095238 0.760719047619048
                    0.0383714285714286 0.674271428571429 0.743552380952381
                    0.0589714285714286 0.683757142857143 0.725385714285714
                    0.0843 0.692833333333333 0.706166666666667
                    0.113295238095238 0.7015 0.685857142857143
                    0.145271428571429 0.709757142857143 0.664628571428571
                    0.180133333333333 0.717657142857143 0.642433333333333
                    0.217828571428571 0.725042857142857 0.619261904761905
                    0.258642857142857 0.731714285714286 0.595428571428571
                    0.302171428571429 0.737604761904762 0.571185714285714
                    0.348166666666667 0.742433333333333 0.547266666666667
                    0.395257142857143 0.7459 0.524442857142857
                    0.442009523809524 0.748080952380952 0.503314285714286
                    0.487123809523809 0.749061904761905 0.483976190476191
                    0.530028571428571 0.749114285714286 0.466114285714286
                    0.570857142857143 0.748519047619048 0.449390476190476
                    0.609852380952381 0.747314285714286 0.433685714285714
                    0.6473 0.7456 0.4188
                    0.683419047619047 0.743476190476191 0.404433333333333
                    0.718409523809524 0.741133333333333 0.39047619047619
                    0.752485714285714 0.7384 0.376814285714286
                    0.785842857142857 0.735566666666667 0.363271428571429
                    0.818504761904762 0.732733333333333 0.349790476190476
                    0.850657142857143 0.7299 0.336028571428571
                    0.882433333333333 0.727433333333333 0.3217
                    0.913933333333333 0.725785714285714 0.30627619047619
                    0.944957142857143 0.726114285714286 0.288642857142857
                    0.973895238095238 0.731395238095238 0.266647619047619
                    0.993771428571429 0.745457142857143 0.240347619047619
                    0.999042857142857 0.765314285714286 0.216414285714286
                    0.995533333333333 0.786057142857143 0.196652380952381
                    0.988 0.8066 0.179366666666667
                    0.978857142857143 0.827142857142857 0.163314285714286
                    0.9697 0.848138095238095 0.147452380952381
                    0.962585714285714 0.870514285714286 0.1309
                    0.958871428571429 0.8949 0.113242857142857
                    0.95982380952381 0.921833333333333 0.0948380952380953
                    0.9661 0.951442857142857 0.0755333333333333
                    0.9763 0.9831 0.0538];
%%
options.dolc = 0;
options.lcm.seeddef = 'vats';
options.lcm.odir = '';
options.lcm.omask = [];
options.lcm.struc.do = 1;
options.lcm.struc.espace = 1;
options.lcm.func.do = 1;
options.lcm.cmd = 1;
options.lc = [];
