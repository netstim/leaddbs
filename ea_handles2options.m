function options=ea_handles2options(handles)
% main function converting GUI handles to options struct used in
% ea_autocoord & ea_write (i.e. the main lead batch functions).

%% some manual options that can be set:
options.endtolerance=10; % how many slices to use with zero signal until end of electrode estimate.
options.sprungwert=4; % how far electrode centroid may be (in xy axis) from last to current slice.
options.refinesteps=0; % how often to re-iterate to reconstruct trajectory. More than 2 should usually not be beneficial. Use 0 to use the direct measurement.
options.tra_stdfactor=0.9; % Default: 0.9 - the lower this factor, the lower the threshold (more included pixels in tra process).
options.cor_stdfactor=1.0; % Default: 1.0 - the higher this factor, the lower the threshold (more included pixels in cor process).

%% set options
options.earoot = ea_getearoot;

try
    options.importdcm.do = handles.dicom2bidscheckbox.Value;
catch
    options.importdcm.do = 0;
end

try
    options.importdcm.tool = handles.dcm2niitool.Value;
end

try
    options.importnii.do = handles.nifti2bidscheckbox.Value;
catch
    options.importnii.do = 0;
end

try % not working when calling from lead_anatomy
    options.normalize.do=(get(handles.normalize_checkbox,'Value') == get(handles.normalize_checkbox,'Max'));
    options.normalize.settings=getappdata(handles.normsettings,'settings');
catch
    options.normalize.do=0;
end

try
    options.normalize.method = handles.normmethod.String{handles.normmethod.Value};
end

try % not working when calling from lead_anatomy
    options.checkreg=get(handles.checkreg,'Value');
catch
    options.checkreg=0;
end

try % also open up checkreg in case of dMRI check registrations
    if get(handles.checkregdmri,'Value')
        options.checkreg=1;
    end
end

try % also open up checkreg in case of dMRI check registrations
    if get(handles.checkregfmri,'Value')
        options.checkreg=1;
    end
end

options.normalize.refine=0;
try
    options.normalize.refine=get(handles.refinefit,'Value');
end

try % not working when calling from lead_anatomy
    options.coregmr.check=(get(handles.coregmrcheck,'Value') == get(handles.coregmrcheck,'Max'));
catch
    options.coregmr.check=0;
end

try
    options.overwriteapproved = handles.overwriteapproved.Value;
end

try
    options.coregmr.method=get(handles.coregmrmethod,'String');
    options.coregmr.method=options.coregmr.method{get(handles.coregmrmethod,'Value')};
    options.coregmr.do=get(handles.coreg_checkbox,'Value');
catch
    options.coregmr.do=0;
    options.coregmr.method='';
end

try
    options.coregb0.addSyN = get(handles.addSyN, 'Value');
end

try % not working when calling from lead_connectome
    % coreg CT
    options.coregct.do = get(handles.coreg_checkbox,'Value') == get(handles.coreg_checkbox,'Max');
    options.coregct.method = handles.coregctmethod.String{handles.coregctmethod.Value};
catch
    options.coregct.do=0;
end

try
    % set modality (MR/CT) in options
    options.modality = get(handles.MRCT,'Value');
catch
    options.modality=1;
end

options.verbose=3; % 4: Show figures but close them 3: Show all but close all figs except resultfig 2: Show all and leave figs open, 1: Show displays only, 0: Show no feedback.

try
    options.sides=ea_assignsides(handles);
catch
    options.sides=1:2;
end
try
    options.doreconstruction=(get(handles.doreconstruction,'Value') == get(handles.doreconstruction,'Max'));
    if strcmp(get(handles.maskwindow_txt,'String'),'auto')
        options.maskwindow=10; % initialize at 10
        options.automask=1; % set automask flag
    else
        options.maskwindow=str2num(get(handles.maskwindow_txt,'String')); % size of the window that follows the trajectory
        options.automask=0; % unset automask flag
    end
catch
    options.doreconstruction=0;
end

options.autoimprove=0; % if true, templates will be modified.
options.axiscontrast=8; % if 8: use tra only but smooth it before. % if 9: use mean of cor and tra but smooth it. % if 10: use raw tra only.
options.zresolution=10; % voxels are being parcellated into this amount of portions.

try
    options.atl.genpt=get(handles.vizspacepopup,'Value')==2; % generate patient specific atlases
catch
    options.atl.genpt=0;
end

if isfield(handles,'vizspacepopup')
    if get(handles.vizspacepopup,'Value')==2 && strcmp(handles.atlassetpopup.String(handles.atlassetpopup.Value),'Use none')
        options.atl.genpt=0;
    end
end

try
    options.atl.can=get(handles.vizspacepopup,'Value')==1; % display canonical atlases
catch
    options.atl.can=1;
end

options.atl.pt=0; % display patient specific atlases. This is not done anymore for now.
try
    options.atl.ptnative=get(handles.vizspacepopup,'Value')==2; % show results in native space.
catch
    options.atl.ptnative=0;
end

try
    if options.atl.ptnative
        options.native=1;
    else
        options.native=0;
    end
catch
    options.native=0;
end

try
    options.d2.write=(get(handles.writeout2d_checkbox,'Value') == get(handles.writeout2d_checkbox,'Max'));
catch
    options.d2.write=0;
end
options.d2.atlasopacity=0.15;

try
    options.refinelocalization=(get(handles.refinelocalization,'Value') == get(handles.refinelocalization,'Max'));
catch
    options.refinelocalization=0;
end

try
    options.scrf.do=get(handles.scrf,'Value');
catch
    options.scrf.do=0;
end

try
    options.scrf.mask = handles.scrfmask.String{handles.scrfmask.Value};
catch
    options.scrf.mask = 'Coarse mask (Schönecker 2008)';
end

try
    options.d3.write=(get(handles.render_checkbox,'Value') == get(handles.render_checkbox,'Max'));
catch
    options.d3.write=0;
end

options.d3.prolong_electrode=2;
options.d3.verbose='on';
options.d3.elrendering=1;
options.d3.exportBB=0; % don't export brainbrowser struct by default
options.d3.hlactivecontacts=0;
options.d3.showactivecontacts=1;
options.d3.showpassivecontacts=1;
options.d3.showisovolume=0;
options.d3.isovscloud=0;
options.d3.mirrorsides=0;

try
    options.d3.autoserver=get(handles.exportservercheck,'Value');
catch
    options.d3.autoserver=0;
end

options.d3.expdf=0;
options.numcontacts=4;

try
    options.entrypointn = handles.targetpopup.Value;
    options.entrypoint = handles.targetpopup.String{options.entrypointn};
end

options.writeoutpm = 0;

try
    options.elmodeln = handles.electrode_model_popup.Value;
    options.elmodel = handles.electrode_model_popup.String{options.elmodeln};
catch
    elms = ea_resolve_elspec;
    options.elmodel = elms{1};
end

try
    options.atlassetn = handles.atlassetpopup.Value;
    options.atlasset = handles.atlassetpopup.String{options.atlassetn};
end

try
    if strcmp(options.atlasset,'Use none')
        options.d3.writeatlases=0;
        options.d2.writeatlases=1;
    else
        options.d3.writeatlases=1;
        options.d2.writeatlases=1;
    end
end

try
    options.reconmethod=get(handles.reconmethod,'String');
    options.reconmethod=options.reconmethod{get(handles.reconmethod,'Value')};
end

options.expstatvat.do=0;

options.fiberthresh=10;
options.writeoutstats=1;

options.colormap=parula(64); % default colormap, use this explicitly rather than 'colormap' to avoid popup window.
try % not working when calling from lead_anatomy
    options.dolc=get(handles.include_lead_connectome_subroutine,'Value');
catch
    options.dolc=0;
end

% lead connectome mapper options:
try
    seed = handles.seeddefpopup.String;
    if iscell(seed)
        seed = seed{handles.seeddefpopup.Value};
    end

    switch seed
        case 'Manually choose seeds'
            options.lcm.seeds = getappdata(handles.seedbutton,'seeds');
            options.lcm.seeddef = 'manual';
        case 'Manually choose parcellation'
            options.lcm.seeds = getappdata(handles.seedbutton,'seeds');
            options.lcm.seeddef = 'parcellation';
        otherwise
            stimname = erase(seed, 'Use VAT: ');
            options.lcm.seeds = stimname;
            options.lcm.seeddef = 'vats';
    end

    try
        options.lcm.odir=getappdata(handles.odirbutton,'odir');
    catch % called from predict module.
        options=rmfield(options,'lcm');
    end

    % if isempty(options.lcm.odir)
    %     if ~strcmp(options.lcm.seeddef,'vats')
    %         try
    %         options.lcm.odir=[fileparts(options.lcm.seeds{1}),filesep];
    %         end
    %     else
    %         options.lcm.odir='';
    %     end
    % end

    options.lcm.omask=getappdata(handles.omaskbutton,'omask');
    options.lcm.struc.do=get(handles.dostructural,'Value');
    options.lcm.func.do=get(handles.dofunctional,'Value');
    options.lcm.func.exportgmtc=get(handles.exportgmtc,'Value');
    options.lcm.cmd=get(handles.command,'Value');
    options.lcm.struc.connectome=get(handles.fiberspopup,'String');
    if iscell(options.lcm.struc.connectome)
        options.lcm.struc.connectome=options.lcm.struc.connectome{get(handles.fiberspopup,'Value')};
    end
    options.lcm.func.connectome=strrep(get(handles.fmripopup,'String'),' > ','>');

    if iscell(options.lcm.func.connectome)
        options.lcm.func.connectome=options.lcm.func.connectome{get(handles.fmripopup,'Value')};
    end
    options.lcm.struc.espace=get(handles.strucexportspace,'Value');
end

% lead predict options:
try
    includes={'Coords','VTA','dMRI','fMRI'};
    todel=[];
    if ~strcmp(handles.inccoordinate.Enable,'on') || ~handles.inccoordinate.Value
        todel=[todel,1];
    end
    if ~strcmp(handles.incvta.Enable,'on') || ~handles.incvta.Value
        todel=[todel,2];
    end
    if ~strcmp(handles.incstructural.Enable,'on') || ~handles.incstructural.Value
        todel=[todel,3];
    end
    if ~strcmp(handles.incfunctional.Enable,'on') || ~handles.incfunctional.Value
        todel=[todel,4];
    end
    includes(todel)=[];
    options.predict.includes=includes;

    % dMRI connectome
    if ~iscell(handles.fiberspopup.String)
        options.predict.dMRIcon{1}=handles.fiberspopup.String;
    else
        options.predict.dMRIcon=handles.fiberspopup.String;
    end
    options.predict.dMRIcon=options.predict.dMRIcon{handles.fiberspopup.Value};

    % fMRI connectome
    if ~iscell(handles.fmripopup.String)
        options.predict.fMRIcon{1}=handles.fmripopup.String;
    else
    	options.predict.fMRIcon=handles.fmripopup.String;
    end
    options.predict.fMRIcon=options.predict.fMRIcon{handles.fmripopup.Value};

    % Chosen prediction model
    mfiles=getappdata(handles.predictionmodel,'mfiles');
    if ~iscell(handles.predictionmodel.String)
    	options.predict.model{1}=handles.predictionmodel.String;
    else
    	options.predict.model=handles.predictionmodel.String;
    end
    options.predict.model=options.predict.model{handles.predictionmodel.Value};
    options.predict.model_mfile=mfiles{handles.predictionmodel.Value};

    % Chosen stimulation name
    if ~iscell(handles.seeddefpopup.String)
    	options.predict.stimulation{1}=handles.seeddefpopup.String;
    else
    	options.predict.stimulation=handles.seeddefpopup.String;
    end
    options.predict.stimulation=options.predict.stimulation{handles.seeddefpopup.Value};
end

try
    options.ecog.extractsurface.do=get(handles.extractsurface,'Value');
    options.ecog.extractsurface.method=get(handles.surfacemethod,'Value');
    %options.ecog.localize=get(handles.localizeecog,'Value');
catch
    options.ecog.extractsurface.do=0;
end


function sides=ea_assignsides(handles)
cnt=1;
elnum = sum(cellfun(@(f) ~isempty(f), regexp(fieldnames(handles),'^side\d+$','match')));
for el=1:elnum
    if get(handles.(['side',num2str(el)]),'Value')
        sides(cnt)=el;
        cnt=cnt+1;
    end
end
