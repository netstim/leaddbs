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
options.dicomimp.do = 0;
options.assignnii = 0;
options.normalize.do = false;
options.normalize.settings = [];
options.normalize.method = 'ANTs (Avants 2008)';
options.normalize.check = false;
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
options.atl.normalize = 0;
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
options.manualheightcorrection = false;
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
options.writeoutpm = 1;
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
options.prefs.dev.profile = 'ah';
options.prefs.pp.do = 0;
options.prefs.pp.csize = 4;
options.prefs.pp.profile = 'local';
options.prefs.prenii_searchstring = 'anat_*.nii';
options.prefs.prenii_order = {
                              't1'
                              't2'
                              'pd'
                              }';
options.prefs.prenii_unnormalized = 'anat_t2.nii';
options.prefs.prenii_unnormalized_t1 = 'anat_t1.nii';
options.prefs.prenii_unnormalized_pd = 'anat_pd.nii';
options.prefs.tranii_unnormalized = 'postop_tra.nii';
options.prefs.sagnii_unnormalized = 'postop_sag.nii';
options.prefs.cornii_unnormalized = 'postop_cor.nii';
options.prefs.rawctnii_unnormalized = 'postop_ct.nii';
options.prefs.ctnii_coregistered = 'rpostop_ct.nii';
options.prefs.tp_ctnii_coregistered = 'tp_rpostop_ct.nii';
options.prefs.patientdir = '';
options.prefs.prenii = 'lanat.nii';
options.prefs.tranii = 'lpostop_tra.nii';
options.prefs.cornii = 'lpostop_cor.nii';
options.prefs.sagnii = 'lpostop_sag.nii';
options.prefs.ctnii = 'lpostop_ct.nii';
options.prefs.gprenii = 'glanat.nii';
options.prefs.gtranii = 'glpostop_tra.nii';
options.prefs.gcornii = 'glpostop_cor.nii';
options.prefs.gsagnii = 'glpostop_sag.nii';
options.prefs.gctnii = 'glpostop_ct.nii';
options.prefs.tp_gctnii = 'tp_glpostop_ct.nii';
options.prefs.rest_searchstring = 'rest*.nii';
options.prefs.rest = 'rest.nii';
options.prefs.lc.struc.maxdist = 2;
options.prefs.lc.struc.minlen = 3;
options.prefs.lc.graphsurfc = [0.2081 0.1663 0.5292];
options.prefs.lc.matsurfc = [0.8 0.7 0.4];
options.prefs.lc.seedsurfc = [0.8 0.1 0.1];
options.prefs.lc.func.regress_global = 1;
options.prefs.lc.func.regress_wmcsf = 1;
options.prefs.lc.func.bphighcutoff = 0.08;
options.prefs.lc.func.bplowcutoff = 0.009;
options.prefs.lc.datadir = '';
options.prefs.lcm.vatseed = 'efield';
options.prefs.b0 = 'b0.nii';
options.prefs.fa = 'fa.nii';
options.prefs.fa2anat = 'fa2anat.nii';
options.prefs.FTR_unnormalized = 'FTR.mat';
options.prefs.FTR_normalized = 'wFTR.mat';
options.prefs.DTD = 'DTD.mat';
options.prefs.HARDI = 'HARDI.mat';
options.prefs.dti = 'dti.nii';
options.prefs.bval = 'dti.bval';
options.prefs.bvec = 'dti.bvec';
options.prefs.sampledtidicom = 'sample_dti_dicom.dcm';
options.prefs.normmatrix = 'lmat.txt';
options.prefs.normalize.default = 'ea_normalize_ants';
options.prefs.normalize.inverse.warp = 'inverse';
options.prefs.normalize.inverse.customtpm = 0;
options.prefs.normalize.createwarpgrids = 0;
options.prefs.normalize.fsl.warpres = 8;
options.prefs.normalize.spm.resolution = 1;
options.prefs.normalize.coreg = 'auto';
options.prefs.ctcoreg.default = 'ea_coregpostopct_ants';
options.prefs.mrcoreg.default = 'spm';
options.prefs.mrcoreg.writeoutcoreg = 0;
options.prefs.scrf.auto = 'mask1';
options.prefs.hullmethod = 2;
options.prefs.hullsmooth = 5;
options.prefs.hullsimplify = 0.6;
options.prefs.lhullmethod = 2;
options.prefs.lhullsmooth = 7;
options.prefs.lhullsimplify = 'auto';
options.prefs.d2.useprepost = 'post';
options.prefs.d2.groupcolors = 'lead';
options.prefs.d2.isovolsmoothed = 's';
options.prefs.d2.isovolcolormap = 'jet';
options.prefs.d2.isovolsepcomb = 'combined';
options.prefs.d3.fiberstyle = 'tube';
options.prefs.d3.fiberdiameter = 0.1;
options.prefs.d3.maxfibers = 200;
options.prefs.d3.colorjitter = 0;
options.prefs.d3.cortexcolor = [0.65 0.65 0.65];
options.prefs.d3.cortexalpha = 0.5;
options.prefs.d3.corticalatlas = 'DKT';
options.prefs.video.path = [-90 10
                            -110 10
                            -180 80
                            -250 10
                            -360 10
                            -450 10];
options.prefs.video.opts.FrameRate = 24;
options.prefs.video.opts.Duration = 30;
options.prefs.video.opts.Periodic = true;
options.prefs.mer.offset = 2;
options.prefs.mer.length = 24;
options.prefs.mer.defaulttract = 1;
options.prefs.mer.tag.visible = 'off';
options.prefs.dicom.dicomfiles = 0;
options.prefs.dicom.tool='dcm2niix';
options.prefs.addfibers = {};
options.prefs.native.warp = 'inverse';
options.prefs.ls.autosave = 0;
options.prefs.ls.dir = '';
options.prefs.env.dev = 1;
options.prefs.ixi.meanage = 60;
options.prefs.ixi.dir = '';
options.prefs.ltx.pdfconverter = '';
options.prefs.prenii_t1 = 'lanat_t1.nii';
options.prefs.prenii_pd = 'lanat_pd.nii';
options.prefs.gprenii_t1 = 'glanat_t1.nii';
options.prefs.gprenii_pd = 'glanat_pd.nii';
options.prefs.vat.gm = 'atlas';
options.prefs.firstrun = 'off';
options.prefs.lc.defaultParcellation = 'AICHA reordered (Joliot 2015)';
options.prefs.machine.defaultatlas = 'STN-Subdivisions (Accolla 2014)';
options.prefs.machine.d2.col_overlay = 0;
options.prefs.machine.d2.con_overlay = 1;
options.prefs.machine.d2.con_color = [0 0 0];
options.prefs.machine.d2.lab_overlay = 0;
options.prefs.machine.d2.bbsize = 10;
options.prefs.machine.d2.backdrop = 'Patient Pre-OP';
options.prefs.machine.d2.fid_overlay = 1;
options.prefs.machine.space = 'MNI_ICBM_2009b_NLIN_ASYM';
options.prefs.machine.lc.graph.struc_func_sim = 0;
options.prefs.machine.lc.graph.nodal_efficiency = 0;
options.prefs.machine.lc.graph.eigenvector_centrality = 0;
options.prefs.machine.lc.graph.degree_centrality = 0;
options.prefs.machine.lc.graph.fthresh = NaN;
options.prefs.machine.lc.graph.sthresh = NaN;
options.prefs.machine.lc.func.compute_CM = 0;
options.prefs.machine.lc.func.compute_GM = 0;
options.prefs.machine.lc.func.prefs.TR = 2.69;
options.prefs.machine.lc.struc.compute_CM = 0;
options.prefs.machine.lc.struc.compute_GM = 0;
options.prefs.machine.lc.struc.ft.method = 'ea_ft_gqi_yeh';
options.prefs.machine.lc.struc.ft.do = 0;
options.prefs.machine.lc.struc.ft.normalize = 1;
options.prefs.machine.methods_show = 0;
options.prefs.machine.methods.show = 1;
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
