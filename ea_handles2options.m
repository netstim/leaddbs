function options=ea_handles2options(handles)

%% some manual options that can be set:


options.endtolerance=10; % how many slices to use with zero signal until end of electrode estimate.
options.sprungwert=4; % how far electrode centroid may be (in xy axis) from last to current slice.
options.refinesteps=0; % how often to re-iterate to reconstruct trajectory. More than 2 should usually not be beneficial. Use 0 to use the direct measurement.
options.tra_stdfactor=0.9; % Default: 0.9 - the lower this factor, the lower the threshold (more included pixels in tra process).
options.cor_stdfactor=1.0; % Default: 1.0 - the higher this factor, the lower the threshold (more included pixels in cor process).




%% set options

%uipatdir=get(handles.patdir_choosebox,'String');

options.earoot=[ea_getearoot];
options.dicomimp=get(handles.dicomcheck,'Value');

options.normalize.do=(get(handles.normalize_checkbox,'Value') == get(handles.normalize_checkbox,'Max'));
options.normalize.method=getappdata(handles.leadfigure,'normmethod');

options.normalize.method=options.normalize.method{get(handles.normmethod,'Value')};
options.normalize.methodn=get(handles.normmethod,'Value');

options.normalize.check=(get(handles.normcheck,'Value') == get(handles.normcheck,'Max'));

% coreg CT
options.coregct.do=(get(handles.coregct_checkbox,'Value') == get(handles.coregct_checkbox,'Max'));
options.coregct.method=getappdata(handles.leadfigure,'coregctmethod');
options.coregct.method=options.coregct.method{get(handles.coregctmethod,'Value')};
options.coregct.methodn=get(handles.coregctmethod,'Value');
options.coregct.coregthreshs= eval( [ '[', get(handles.coregthreshs,'String'), ']' ] );

options.coregctcheck=get(handles.coregctcheck,'Value');


options.coregmr.method=get(handles.coregmrpopup,'Value');

% set modality (MR/CT) in options
options.modality = get(handles.MRCT,'Value');




options.verbose=3; % 4: Show figures but close them 3: Show all but close all figs except resultfig 2: Show all and leave figs open, 1: Show displays only, 0: Show no feedback.

%sidelog=[get(handles.right_checkbox,'Value') == get(handles.right_checkbox,'Max'),get(handles.left_checkbox,'Value') == get(handles.left_checkbox,'Max')];
%sidepos=[1,2];

%options.sides=sidepos(logical(sidelog)); %side=1 -> left electrode, side=2 -> right electrode. both: [1:2]
options.sides=1:2;

options.doreconstruction=(get(handles.doreconstruction_checkbox,'Value') == get(handles.doreconstruction_checkbox,'Max'));
if strcmp(get(handles.maskwindow_txt,'String'),'auto')
options.maskwindow=10; % initialize at 10
options.automask=1; % set automask flag
else
options.maskwindow=str2num(get(handles.maskwindow_txt,'String')); % size of the window that follows the trajectory
options.automask=0; % unset automask flag
end
options.autoimprove=0; % if true, templates will be modified.
options.axiscontrast=8; % if 8: use tra only but smooth it before. % if 9: use mean of cor and tra but smooth it. % if 10: use raw tra only.
options.zresolution=10; % voxels are being parcellated into this amount of portions.

options.atl.genpt=get(handles.vizspacepopup,'Value')==2; % generate patient specific atlases
options.atl.normalize=0; % normalize patient specific atlasset. This is not done anymore for now.
options.atl.can=get(handles.vizspacepopup,'Value')==1; % display canonical atlases
options.atl.pt=0; % display patient specific atlases. This is not done anymore for now.
options.atl.ptnative=get(handles.vizspacepopup,'Value')==2; % show results in native space.
if options.atl.ptnative
    options.native=1;
else
    options.native=0;
end

options.d2.write=(get(handles.writeout2d_checkbox,'Value') == get(handles.writeout2d_checkbox,'Max'));

options.d2.atlasopacity=0.15;


options.manualheightcorrection=(get(handles.manualheight_checkbox,'Value') == get(handles.manualheight_checkbox,'Max'));
options.d3.write=(get(handles.render_checkbox,'Value') == get(handles.render_checkbox,'Max'));
options.d3.prolong_electrode=2;
options.d3.verbose='on';
options.d3.elrendering=1;
options.d3.hlactivecontacts=0;
options.d3.showactivecontacts=1;
options.d3.showpassivecontacts=1;
options.d3.showisovolume=0;
options.d3.isovscloud=0;
options.d3.autoserver=get(handles.exportservercheck,'Value');
options.d3.expdf=0;
options.numcontacts=4;
options.entrypoint=get(handles.targetpopup,'String');
options.entrypoint=options.entrypoint{get(handles.targetpopup,'Value')};
options.entrypointn=get(handles.targetpopup,'Value');

options.writeoutpm=1;

options.elmodeln = get(handles.electrode_model_popup,'Value');
string_list = get(handles.electrode_model_popup,'String');
options.elmodel=string_list{options.elmodeln};
options.atlasset=get(handles.atlassetpopup,'String'); %{get(handles.atlassetpopup,'Value')}
options.atlasset=options.atlasset{get(handles.atlassetpopup,'Value')};
options.atlassetn=get(handles.atlassetpopup,'Value');

if strcmp(options.atlasset,'Use none');
    options.d3.writeatlases=0;
    options.d2.writeatlases=1;
else
    options.d3.writeatlases=1;
    options.d2.writeatlases=1;
end


options.expstatvat.do=0;

options.fiberthresh=10;
options.writeoutstats=1;


options.colormap=colormap;

options.dolc=get(handles.include_lead_connectome_subroutine,'Value');

