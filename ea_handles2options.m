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

%uipatdir=get(handles.patdir_choosebox,'String');

options.earoot=[ea_getearoot];
try % not working when calling from lead_anatomy
    options.dicomimp=get(handles.dicomcheck,'Value');
    options.normalize.do=(get(handles.normalize_checkbox,'Value') == get(handles.normalize_checkbox,'Max'));
    options.normalize.settings=getappdata(handles.normsettings,'settings');
catch
    options.dicomimp=0;
    options.normalize.do=0;
end
try
options.normalize.method=getappdata(handles.leadfigure,'normmethod');
options.normalize.method=options.normalize.method{get(handles.normmethod,'Value')};
options.normalize.methodn=get(handles.normmethod,'Value');
end
try % not working when calling from lead_anatomy
    options.normalize.check=(get(handles.normcheck,'Value') == get(handles.normcheck,'Max'));
catch
    options.normalize.check=0;
end


try
    options.coregmr.method=get(handles.coregmrpopup,'String');
    options.coregmr.method=options.coregmr.method{get(handles.coregmrpopup,'Value')};
catch
    options.coregmr.method='';
end
try % not working when calling from lead_connectome
    % coreg CT
    options.coregct.do=(get(handles.coregct_checkbox,'Value') == get(handles.coregct_checkbox,'Max'));
    options.coregct.method=getappdata(handles.leadfigure,'coregctmethod');
    options.coregct.method=options.coregct.method{get(handles.coregctmethod,'Value')};
    options.coregct.methodn=get(handles.coregctmethod,'Value');
    options.coregct.coregthreshs= eval( [ '[', get(handles.coregthreshs,'String'), ']' ] );
    
    options.coregctcheck=get(handles.coregctcheck,'Value');
catch
    options.coregct.do=0;
    options.coregctcheck=0;
end

try
        % set modality (MR/CT) in options
    options.modality = get(handles.MRCT,'Value');
catch
        options.modality=1;
end


options.verbose=3; % 4: Show figures but close them 3: Show all but close all figs except resultfig 2: Show all and leave figs open, 1: Show displays only, 0: Show no feedback.

%sidelog=[get(handles.right_checkbox,'Value') == get(handles.right_checkbox,'Max'),get(handles.left_checkbox,'Value') == get(handles.left_checkbox,'Max')];
%sidepos=[1,2];

%options.sides=sidepos(logical(sidelog)); %side=1 -> left electrode, side=2 -> right electrode. both: [1:2]
try
    switch get(handles.sidespopup,'Value')
        case 1
            options.sides=1:2;
        case 2
            options.sides=1;
        case 3
            options.sides=1;
    end
catch
    options.sides=1:2;
end

try
    options.doreconstruction=(get(handles.doreconstruction_checkbox,'Value') == get(handles.doreconstruction_checkbox,'Max'));
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
options.atl.normalize=0; % normalize patient specific atlasset. This is not done anymore for now.
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
    options.manualheightcorrection=(get(handles.manualheight_checkbox,'Value') == get(handles.manualheight_checkbox,'Max'));
catch
    options.manualheightcorrection=0;
end

try
    options.scrf=get(handles.scrf,'Value');
catch
    options.scrf=0;
end

try
    options.d3.write=(get(handles.render_checkbox,'Value') == get(handles.render_checkbox,'Max'));
catch
    options.d3.write=0;
end
options.d3.prolong_electrode=2;
options.d3.verbose='on';
options.d3.elrendering=1;
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
    options.entrypoint=get(handles.targetpopup,'String');
    options.entrypoint=options.entrypoint{get(handles.targetpopup,'Value')};
    options.entrypointn=get(handles.targetpopup,'Value');
end
options.writeoutpm=1;
try
    options.elmodeln = get(handles.electrode_model_popup,'Value');
    string_list = get(handles.electrode_model_popup,'String');
    options.elmodel=string_list{options.elmodeln};
end
try
    options.atlasset=get(handles.atlassetpopup,'String'); %{get(handles.atlassetpopup,'Value')}
    options.atlasset=options.atlasset{get(handles.atlassetpopup,'Value')};
    options.atlassetn=get(handles.atlassetpopup,'Value');
end
try
    if strcmp(options.atlasset,'Use none');
        options.d3.writeatlases=0;
        options.d2.writeatlases=1;
    else
        options.d3.writeatlases=1;
        options.d2.writeatlases=1;
    end
end


options.expstatvat.do=0;

options.fiberthresh=10;
options.writeoutstats=1;


options.colormap=colormap;
try % not working when calling from lead_anatomy
    options.dolc=get(handles.include_lead_connectome_subroutine,'Value');
catch
    options.dolc=0;
end


% lead connectome mapper options:

try

    sdp=get(handles.seeddefpopup,'String');
    if iscell(sdp)
    sdp=sdp{get(handles.seeddefpopup,'Value')};
    end
    switch sdp
        case 'Manually choose seeds'
            options.lcm.seeds=getappdata(handles.seedbutton,'seeds');
            options.lcm.seeddef='manual';
        otherwise
            stimname=sdp(11:end);
            options.lcm.seeds=stimname;
            options.lcm.seeddef='vats';
    end
    options.lcm.odir=getappdata(handles.odirbutton,'odir');
    if isempty(options.lcm.odir)
        if ~strcmp(options.lcm.seeddef,'vats')
            try
            options.lcm.odir=[fileparts(options.lcm.seeds{1}),filesep];
            end
        else
            options.lcm.odir='';
        end
    end
    options.lcm.omask=getappdata(handles.omaskbutton,'omask');
    options.lcm.struc.do=get(handles.dostructural,'Value');
    options.lcm.func.do=get(handles.dofunctional,'Value');
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