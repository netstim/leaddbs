function [cortexH,cortex] = ea_showcortex(varargin)
% This function shows a cortical reconstruction in the 3D-Scene viewer. It
% reads in cortex.mat found in the eAuto_root/templates/cortex folder, or
% in the patientdirectory on button press (Cortical Reconstruction
% Visualization) in the Lead-Anatomy scene. To view patient specific
% cortical reconstructions, load them to the patient directory using
% ea_importfs via the "Import FS" button in the main Lead DBS GUI.
% __________________________________________________________________________________
% Copyright (C) 2017 University of Pittsburgh (UPMC), Brain Modulation Lab
% Ari Kappel

if nargin==2
    options=varargin{2};
end
if nargin>2
    elstruct=varargin{2};
    options=varargin{3};
end

resultfig=varargin{1};
set(0, 'currentfigure', resultfig);  % for figures
cortexH=getappdata(resultfig,'cortex');
try delete(cortexH{1}); end
try delete(cortexH{2}); end

% Initial Opening
% if ~isfield(getappdata(resultfig),'cortex')
color = options.prefs.d3.cortexcolor; % default color is gray
alpha = options.prefs.d3.cortexalpha; % default alph is 0.333
% else
%     color = options.prefs.d3.cortexcolor;
%     awin = getappdata(resultfig,'awin');
%     appdata = getappdata(awin,'UsedByGUIData_m');
%     alpha = str2double(appdata.cortexalpha.String);
%     clear appdata
% end


nm=[0:2]; % native and mni
try
    nmind=[options.atl.pt,options.atl.can,options.atl.ptnative]; % which shall be performed?
catch
    nmind=[0 1 0];
end
nm=nm(logical(nmind)); % select which shall be performed.

% switch between patient and template cortex
slicecontroldata = getappdata(gcf);
if nm==1 && contains(slicecontroldata.templateused,'Patient')
    nm = 0;
end

% macaque brain currently not supported % mcr=ea_checkmacaque(options);
% mcr = '';

for nativemni=nm % switch between native and mni space.

    switch nativemni
        case 0 % patient cortex in mni space in future release
            % root=[options.root,options.patientname,filesep];
            % adir=[root,''];
            if ~exist([ea_space(options,'cortex'),'CortexHiRes.mat'],'file')
                ea_error('Missing Template Cortex')
                return
            end
            adir=ea_space(options,'cortex');
            casestr = 'Template';
            reslice='yes'; %future option to import patient brain and reslice to mni space
        case 1 % template cortex in mni space
            %root=[options.earoot];
            if ~exist([ea_space(options,'cortex'),'CortexHiRes.mat'],'file')
                ea_error('Missing Template Cortex')
                return
            end
            adir=ea_space(options,'cortex');
            casestr = 'Template';
            reslice='no';
        case 2 % patient cortex in native space
            root=[options.root,options.patientname,filesep];
            adir=[root,'cortex/'];
            casestr = 'Patient';
            reslice='no';
    end
end

% Check if CortexHiRes.mat and CortexLowRes_*.mat already exists
files = dir(adir);
files = files(cellfun(@(x) isempty(regexp(x, '^\.', 'once')), {files.name}));
files = files(~[files.isdir]);
files = {files(contains({files.name},'Cortex')).name};

    % % Use this to choose which file to load:
    % if size(files,2)>=2
    %     str = [{['More than one Cortex found in the ' casestr ' folder.']},...
    %         {'Please, select which Cortex you would like to view.'}];
    %     file = char(files(listdlg('PromptString',str,'Name',options.patientname,...
    %         'SelectionMode','single','ListString',files,...
    %         'ListSize',[250 150])));
    %     cortex = load([adir,file],'Vertices_rh','Faces_rh','Vertices_lh','Faces_lh');
    % elseif ~isempty(files)
    %     file = files{1};
    %     cortex = load([adir,file],'Vertices_rh','Faces_rh','Vertices_lh','Faces_lh');
    % else
    %     fprintf('No Cortex found in patient directory. \nRunning Import FS...\n')
    %     ea_importfs(options)
    %     try
    %         cortex = load([adir,'CortexHiRes.mat'],'Vertices_rh','Faces_rh','Vertices_lh','Faces_lh');
    %     catch
    %     end
    % end
% Only option is CortexHiRes.mat - can resample annot file for low res in
% future release
str = 'Select Cortex';
file = listdlg('PromptString',str,'SelectionMode','single',...
    'ListString','Cortex','ListSize',[250,150]);
if ~isempty(file)
    try
        cortex = load([adir,'CortexHiRes.mat'],'Vertices_rh','Faces_rh','Vertices_lh','Faces_lh');
    catch
    end
else
    return
end

% % reslice patient cortex to mni space in future release
% % for now always use template cortex in mni space
% if strcmp(reslice,'yes'); end

% try
%     load([adir,'/annot_',options.prefs.d3.corticalatlas,'.mat'])
% end

% Show cortex
tagstr = {'RH','LH'};
set(0, 'currentfigure', resultfig); hold on
cortexH{1} = patch('vertices',cortex.Vertices_rh,'faces',cortex.Faces_rh(:,[1 3 2]),...
    'FaceVertexCData',repmat([0.65 0.65 0.65],[size(cortex.Vertices_rh,1),1]),...
    'edgecolor','none','FaceAlpha',alpha,'FaceColor','interp',...
    'facelighting', 'gouraud', 'specularstrength', .25,'Tag',tagstr{1});

cortexH{2} = patch('vertices',cortex.Vertices_lh,'faces',cortex.Faces_lh(:,[1 3 2]),...
    'FaceVertexCData',repmat([0.65 0.65 0.65],[size(cortex.Vertices_lh,1),1]),...
    'edgecolor','none','FaceAlpha',alpha,'FaceColor','interp',...
    'facelighting', 'gouraud', 'specularstrength', .25,'Tag',tagstr{2});

% Load annotation file
files = dir(adir);
files = {files(~cellfun(@(x) isempty(regexp(x, 'annot', 'once')), {files.name})).name};
if ~isempty(files)
    atlas = files{~cellfun(@(x) isempty(regexp(x, ['annot_',options.prefs.d3.cortex_defaultatlas,'.mat'], 'match')), files)};
    if ~isempty(atlas)
        load(fullfile(adir,atlas));

        for side = 1:2
            % Choose gyri
            structures{side} = annot(side).colortable.struct_names;
            labelidx{side} = mat2cell((1:length(structures{side}))', ones(1,length(structures{side})));
            annot(side).cdat = get(cortexH{side},'FaceVertexCData');
            structidx = arrayfun(@(x) find(annot(side).label==annot(side).colortable.table(x,5)),[labelidx{side}{:}],'uni',0);
            labels = cell2mat(arrayfun(@(x) find(x==annot(side).label),annot(side).colortable.table(:,5),'uni',0));

            for i = 1:size(annot(side).colortable.table,1)
                colorindex = structidx{i};
                visindex=setdiff(labels,colorindex);
                %         invisindex=[invisindex;setdiff(labels,colorindex)];
                annot(side).cdat(colorindex,:) = repmat(annot(side).colortable.table(i,1:3)/256,[length(colorindex),1]);
                annot(side).adat(colorindex,:) = ones(length(colorindex),1) * alpha;
            end
            set(cortexH{side},'FaceVertexCData',annot(side).cdat)
            set(cortexH{side},'FaceVertexAlphaData',annot(side).adat,'FaceAlpha','interp')
        end
        setappdata(resultfig,'cortex',cortexH);
        setappdata(resultfig,'annot',annot);

        selstate.structures = structures;
        selstate.labelidx = labelidx;
        setappdata(resultfig, 'cortsurfs', selstate)

        atlas = strrep(atlas,'.mat','');
        cswin = ea_cortexselect(cortex, annot, atlas, colorindex, structures{1}, options, resultfig);
        set(cswin, 'visible', options.d3.verbose);
        setappdata(resultfig, 'aswin', cswin);
    end
end
