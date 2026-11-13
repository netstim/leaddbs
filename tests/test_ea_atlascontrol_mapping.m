function tests = test_ea_atlascontrol_mapping
% Regression for the FOCUS volume/tract basename collisions.
tests = functiontests(localfunctions);
end

function testDuplicateNamesAndReorderedSubgroups(testCase)
assumeTrue(testCase,usejava('jvm'),'Atlas control requires Java.');
scene=figure('Visible','off');
cleanupScene=onCleanup(@() delete(scene)); %#ok<NASGU>
ax=axes('Parent',scene);
toolbar=uitoolbar(scene);
names={'palhy.nii.gz','palhy.mat','CNIII.nii.gz','CNIII.mat', ...
    'LHonly.nii.gz','midline.nii.gz'};
atlases.names=names;
atlases.types=[3 3 3 3 2 5];
atlases.pixdim={1 1;'fibers' 'fibers';1 1;'fibers' 'fibers';[] 1;1 []};
atlases.colormap=lines(numel(names));
atlases.colors=1:numel(names);
atlases.labelnames={'Display labels'};
atlases.labels={{'Pallidofugal','Pallidofugal','Oculomotor','Oculomotor', ...
    'Left only','Midline'}};
atlases.subgroups(1)=struct('label','Tracts','entries',[2 4 5]);
atlases.subgroups(2)=struct('label','Volumes','entries',[1 3 6]);
atlasmap=[1 1;1 2;2 1;2 2;3 1;3 2;4 1;4 2;5 2;6 1];
surfs=cell(size(atlasmap,1),1);
buttons=gobjects(size(atlasmap,1),1);
labels=gobjects(size(atlasmap,1),1);
for k=1:size(atlasmap,1)
    [~,stem]=fileparts(names{atlasmap(k,1)});
    [~,stem]=fileparts(stem);
    if atlasmap(k,1)==6
        suffix='_midline';
    elseif atlasmap(k,2)==1
        suffix='_right';
    else
        suffix='_left';
    end
    tag=[stem,suffix];
    surfs{k}=patch(ax,'Faces',[1 2 3],'Vertices',[0 0 0;1 0 0;0 1 0], ...
        'Tag',tag,'Visible','on');
    buttons(k)=uitoggletool(toolbar,'Tag',tag,'State','on');
    labels(k)=text(ax,0,0,0,tag,'Visible','off');
end
labelbutton=uitoggletool(toolbar,'State','off');
setappdata(scene,'atlasmap',atlasmap);
options.prefs=ea_prefs;
options.native=0;
options.atlasset='FOCUS mapping regression';
options.d3.verbose='off';
control=ea_atlasselect(buttons,surfs,atlases,options,scene,labelbutton,labels);
cleanupControl=onCleanup(@() delete(control)); %#ok<NASGU>
handles=getappdata(control,'handles');
h=getappdata(control,'h');
verifyEqual(testCase,h.sgsubsurf{1}{1},[3 4]);
verifyEqual(testCase,h.sgsubsurf{1}{2},[7 8]);
verifyEqual(testCase,h.sgsubsurf{1}{3},9);
verifyEqual(testCase,h.sgsubsurf{2}{3},10);
verifyEqual(testCase,surfs{1}.Tag,surfs{3}.Tag);
verifyEqual(testCase,surfs{5}.Tag,surfs{7}.Tag);

% The tract preset must select the fibers while hiding same-named volumes.
preset=struct('label','Tracts only','show',[2 4 5],'hide',[1 3 6], ...
    'default','absolute');
applyPreset(handles,preset);
expected=repmat({'off'},size(buttons));
expected([3 4 7 8 9])={'on'};
verifyEqual(testCase,arrayfun(@(b) char(b.State),buttons,'UniformOutput',false),expected);
verifyEqual(testCase,cellfun(@(p) char(p.Visible),surfs,'UniformOutput',false),expected);

% Show/hide all must reach all four formerly orphaned tract patches.
preset.show=[];
preset.hide=1:numel(names);
applyPreset(handles,preset);
verifyTrue(testCase,all(cellfun(@(p) strcmp(p.Visible,'off'),surfs)));
preset.show=1:numel(names);
preset.hide=[];
applyPreset(handles,preset);
verifyTrue(testCase,all(cellfun(@(p) strcmp(p.Visible,'on'),surfs)));

% A checkbox click must hide only the chosen oculomotor tract hemisphere.
import com.mathworks.mwswing.checkboxtree.*
jtree=getappdata(control,'jtree');
set(h.sgsubside{1}{2}{2},'SelectionState',SelectionState.NOT_SELECTED);
event=java.awt.event.MouseEvent(jtree,502,0,0,0,0,1,false);
ea_atlasselect('mouseReleasedCallback',jtree,event,h);
verifyEqual(testCase,char(surfs{8}.Visible),'off');
verifyEqual(testCase,char(buttons(8).State),'off');
verifyEqual(testCase,char(surfs{6}.Visible),'on');
verifyEqual(testCase,char(surfs{7}.Visible),'on');

% Synchronizing an individual toolbar change must target that exact tract.
buttons(3).State='off';
ea_synctree(handles);
verifyEqual(testCase,char(surfs{3}.Visible),'off');
verifyEqual(testCase,char(surfs{1}.Visible),'on');
verifyEqual(testCase,char(h.sgsubside{1}{1}{1}.getSelectionState), ...
    'not selected');
verifyEqual(testCase,char(h.sgsubside{2}{1}{1}.getSelectionState), ...
    'selected');
end

function applyPreset(handles,preset)
setappdata(handles.presets,'presetactions',{preset});
handles.presets.String={preset.label};
handles.presets.Value=1;
ea_atlasselect('presets_Callback',handles.presets,[],handles);
end
