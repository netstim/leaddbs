function resultfig=ea_mnifigure(atlasname)
options=getoptslocal;
if exist('atlasname','var')
    options.atlasset=atlasname;
else
    options.atlasset='';
    options.d3.writeatlases=0;
end
options.leadprod='dbs';
options.d3.elrendering=1;
options.d3.exportBB=0;
resultfig=ea_elvis(options);
colormap(gray)
hold on
ea_zoomcenter(resultfig.CurrentAxes, [0,0,0], 3);


function options=getoptslocal
options.patientname='No Patient Selected';
options.endtolerance = 10;
options.sprungwert = 4;
options.refinesteps = 0;
options.tra_stdfactor = 0.9;
options.cor_stdfactor = 1;
options.earoot = ea_getearoot;
options.normalize.do = false;
options.normalize.settings = [];
options.normalize.method = 'ANTs (Avants 2008)';
options.checkreg = false;
options.coregmr.method = 'SPM';
options.coregct.do = false;
options.coregct.method = 'ANTs (Avants 2008)';
options.coregct.coregthreshs = NaN;
options.coregctcheck = 0;
options.modality = 1;
options.verbose = 3;
options.sides = [1 2];
options.doreconstruction = false;
options.maskwindow = 10;
options.automask = 1;
options.autoimprove = 0;
options.axiscontrast = 8;
options.zresolution = 10;
options.atl.genpt = false;
options.atl.can = true;
options.atl.pt = 0;
options.atl.ptnative = false;
options.native = 0;
options.d2.col_overlay = 0;
options.d2.con_overlay = 1;
options.d2.con_color = [0 0 0];
options.d2.lab_overlay = 0;
options.d2.bbsize = 10;
options.d2.backdrop = 'Patient Pre-OP';
options.d2.write = false;
options.d2.atlasopacity = 0.15;
options.d2.writeatlases = 1;
options.refinelocalization = false;
options.scrf = 0;
options.d3.write = true;
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
options.d3.writeatlases = 1;
options.numcontacts = 4;
options.entrypoint = 'STN, GPi or ViM';
options.entrypointn = 1;
options.writeoutpm = 0;
options.elmodeln = 1;
options.elmodel = 'Medtronic 3389';
options.atlasset = '';
options.atlassetn = 34;
options.expstatvat.do = 0;
options.fiberthresh = 10;
options.writeoutstats = 1;
%%
options.colormap = [0 0 0.5625
                    0 0 0.625
                    0 0 0.6875
                    0 0 0.75
                    0 0 0.8125
                    0 0 0.875
                    0 0 0.9375
                    0 0 1
                    0 0.0625 1
                    0 0.125 1
                    0 0.1875 1
                    0 0.25 1
                    0 0.3125 1
                    0 0.375 1
                    0 0.4375 1
                    0 0.5 1
                    0 0.5625 1
                    0 0.625 1
                    0 0.6875 1
                    0 0.75 1
                    0 0.8125 1
                    0 0.875 1
                    0 0.9375 1
                    0 1 1
                    0.0625 1 1
                    0.125 1 0.9375
                    0.1875 1 0.875
                    0.25 1 0.8125
                    0.3125 1 0.75
                    0.375 1 0.6875
                    0.4375 1 0.625
                    0.5 1 0.5625
                    0.5625 1 0.5
                    0.625 1 0.4375
                    0.6875 1 0.375
                    0.75 1 0.3125
                    0.8125 1 0.25
                    0.875 1 0.1875
                    0.9375 1 0.125
                    1 1 0.0625
                    1 1 0
                    1 0.9375 0
                    1 0.875 0
                    1 0.8125 0
                    1 0.75 0
                    1 0.6875 0
                    1 0.625 0
                    1 0.5625 0
                    1 0.5 0
                    1 0.4375 0
                    1 0.375 0
                    1 0.3125 0
                    1 0.25 0
                    1 0.1875 0
                    1 0.125 0
                    1 0.0625 0
                    1 0 0
                    0.9375 0 0
                    0.875 0 0
                    0.8125 0 0
                    0.75 0 0
                    0.6875 0 0
                    0.625 0 0
                    0.5625 0 0];
%%
options.dolc = 0;
options.uipatdirs = {'No Patient Selected'};
options.prefs = ea_prefs;
options.lc.general.parcellation = 'AICHA reordered (Joliot 2015)';
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
options.lc.struc.ft.do = 0;
options.lc.struc.ft.normalize = 1;
