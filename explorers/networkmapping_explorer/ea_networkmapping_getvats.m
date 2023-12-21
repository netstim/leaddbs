function vatlist = ea_networkmapping_getvats(obj)
% Return list of VATs

numPatient = length(obj.allpatients);
vatlist = cell(numPatient*2,1);

[modlist,modType] = ea_genmodlist;
[~,ix] = ismember(obj.connectome,modlist);
if ~ix
    ea_error('Something went wrong, connectome disappeared?');
end
modType = modType(ix);
modStr = {'dMRI', 'fMRI'};
modStr = modStr{modType};

disp('Construct VAT list...')
for sub=1:numPatient
    [~, subPrefix] = fileparts([obj.allpatients{sub}, '_']);
    % Original VAT E-field
    stimFolder = [obj.allpatients{sub}, filesep, 'stimulations', filesep, ea_nt(0), 'gs_', obj.M.guid];
    stimParams = ea_regexpdir(stimFolder, 'stimparameters\.mat$', 0);
    load(stimParams{1}, 'S');
    modelLabel = ea_simModel2Label(S.model);
    vatlist{sub,1} = fullfile(stimFolder, [subPrefix, 'sim-efield_model-', modelLabel, '_seed-', modStr, '.nii']);

    if 1 % for now always recreate compound VTA seed
        rungenlocalmapper(obj,sub)
    end

    % Mirrored VAT E-field
    vatlist{numPatient+sub,1} = fullfile(stimFolder, [subPrefix, 'sim-efield_model-', modelLabel, '_seed-', modStr, '_hemidesc-flipped.nii']);
    if ~isfile(vatlist{numPatient+sub,1})
        ea_flip_lr_nonlinear(vatlist{sub,1}, vatlist{numPatient+sub,1}, 1);
    end
end


function rungenlocalmapper(obj,sub)

ptdir=obj.allpatients{sub};
guid=['gs_',obj.M.guid];

options = getoptslocal;
options.prefs.lcm.vatseed = 'efield';
options.lcm.seeds = guid;
options.lcm.seeddef = 'vats';
options.lcm.odir = [];
options.lcm.omask = [];
options.lcm.struc.do = 0;
options.lcm.struc.connectome = obj.connectome;
options.lcm.struc.espace = 1;
options.lcm.func.do = 0;
options.lcm.func.exportgmtc = 0;
options.lcm.func.connectome = strrep(obj.connectome,' > ','>');
options.lcm.cmd = 1;
options.lcm.onlygenvats=1; % force to only create VTAs (do not create network maps)
options.ecog.extractsurface.do = 0;
options.uivatdirs = {ptdir};
options.uipatdirs = {''};

ea_run('run', options);


function options = getoptslocal
options.endtolerance = 10;
options.sprungwert = 4;
options.refinesteps = 0;
options.tra_stdfactor = 0.9;
options.cor_stdfactor = 1;
options.earoot = ea_getearoot;
options.normalize.do = 0;
options.normalize.method = [];
options.checkreg = false;
options.normalize.refine = 0;
options.coregmr.check = 0;
options.coregmr.do = 0;
options.coregmr.method = '';
options.coregct.do = 0;
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
options.d2.lab_overlay = 0;
options.d2.bbsize = 50;
options.d2.backdrop = 'MNI152NLin2009bAsym T1 (Fonov 2011)';
options.d2.fid_overlay = 1;
options.d2.write = 0;
options.d2.atlasopacity = 0.15;
options.refinelocalization = 0;
options.scrf.do = 0;
options.scrf.mask = 'Coarse mask (Schönecker 2008)';
options.d3.write = 0;
options.d3.prolong_electrode = 2;
options.d3.verbose = 'on';
options.d3.elrendering = 1;
options.d3.exportBB = 0;
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
options.elmodel = 'Medtronic 3389';
options.expstatvat.do = 0;
options.fiberthresh = 10;
options.writeoutstats = 1;
%%
options.colormap = [0.2422 0.1504 0.6603
    0.250390476190476 0.164995238095238 0.707614285714286
    0.257771428571429 0.181780952380952 0.751138095238095
    0.264728571428571 0.197757142857143 0.795214285714286
    0.270647619047619 0.21467619047619 0.836371428571429
    0.275114285714286 0.234238095238095 0.870985714285714
    0.2783 0.255871428571429 0.899071428571429
    0.280333333333333 0.278233333333333 0.9221
    0.281338095238095 0.300595238095238 0.941376190476191
    0.281014285714286 0.322757142857143 0.957885714285714
    0.279466666666667 0.344671428571429 0.971676190476191
    0.275971428571429 0.366680952380952 0.982904761904762
    0.269914285714286 0.3892 0.9906
    0.260242857142857 0.412328571428571 0.995157142857143
    0.244033333333333 0.435833333333333 0.998833333333333
    0.220642857142857 0.460257142857143 0.997285714285714
    0.196333333333333 0.484719047619048 0.989152380952381
    0.183404761904762 0.507371428571429 0.979795238095238
    0.178642857142857 0.528857142857143 0.968157142857143
    0.176438095238095 0.549904761904762 0.952019047619048
    0.168742857142857 0.570261904761905 0.935871428571428
    0.154 0.5902 0.9218
    0.146028571428571 0.609119047619048 0.907857142857143
    0.13802380952381 0.627628571428572 0.897290476190476
    0.124814285714286 0.645928571428571 0.888342857142857
    0.111252380952381 0.6635 0.876314285714286
    0.0952095238095238 0.679828571428571 0.859780952380952
    0.0688714285714285 0.694771428571429 0.839357142857143
    0.0296666666666667 0.708166666666667 0.816333333333333
    0.00357142857142857 0.720266666666667 0.7917
    0.00665714285714287 0.731214285714286 0.766014285714286
    0.0433285714285715 0.741095238095238 0.739409523809524
    0.096395238095238 0.75 0.712038095238095
    0.140771428571429 0.7584 0.684157142857143
    0.1717 0.766961904761905 0.655442857142857
    0.193766666666667 0.775766666666667 0.6251
    0.216085714285714 0.7843 0.5923
    0.246957142857143 0.791795238095238 0.556742857142857
    0.290614285714286 0.797290476190476 0.518828571428572
    0.340642857142857 0.8008 0.478857142857143
    0.3909 0.802871428571428 0.435447619047619
    0.445628571428572 0.802419047619048 0.390919047619048
    0.5044 0.7993 0.348
    0.561561904761905 0.794233333333333 0.304480952380953
    0.617395238095238 0.787619047619048 0.261238095238095
    0.671985714285714 0.779271428571429 0.2227
    0.7242 0.769842857142857 0.191028571428571
    0.773833333333333 0.759804761904762 0.164609523809524
    0.820314285714286 0.749814285714286 0.153528571428571
    0.863433333333333 0.7406 0.159633333333333
    0.903542857142857 0.733028571428571 0.177414285714286
    0.939257142857143 0.728785714285714 0.209957142857143
    0.972757142857143 0.729771428571429 0.239442857142857
    0.995647619047619 0.743371428571429 0.237147619047619
    0.996985714285714 0.765857142857143 0.219942857142857
    0.995204761904762 0.789252380952381 0.202761904761905
    0.9892 0.813566666666667 0.188533333333333
    0.978628571428571 0.838628571428572 0.176557142857143
    0.967647619047619 0.8639 0.164290476190476
    0.961009523809524 0.889019047619048 0.153676190476191
    0.959671428571429 0.913457142857143 0.142257142857143
    0.962795238095238 0.937338095238095 0.126509523809524
    0.969114285714286 0.960628571428571 0.106361904761905
    0.9769 0.9839 0.0805];
%%
options.dolc = 0;
options.leadprod = 'mapper';
options.prefs = ea_prefs;
options.lc.general.parcellation = 'ABI_atlas_reduced_V2';
options.lc.general.parcellationn = 2;
options.lc.graph.struc_func_sim = 0;
options.lc.graph.nodal_efficiency = 0;
options.lc.graph.eigenvector_centrality = 0;
options.lc.graph.degree_centrality = 0;
options.lc.graph.fthresh = NaN;
options.lc.graph.sthresh = NaN;
options.lc.func.compute_CM = 0;
options.lc.func.compute_GM = 0;
options.lc.func.prefs.TR = 2.69;
options.lc.struc.compute_CM = 0;
options.lc.struc.compute_GM = 0;
options.lc.struc.ft.method = 'ea_ft_gqi_yeh';
options.lc.struc.ft.methodn = 2;
options.lc.struc.ft.do = 1;
options.lc.struc.ft.normalize = 1;
options.lc.struc.ft.dsistudio.fiber_count = 5000;
options.lc.struc.ft.upsample.factor = 1;
options.lc.struc.ft.upsample.how = 0;
options.exportedJob = 1;
