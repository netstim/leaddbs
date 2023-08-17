function [PL,thresh]=ea_cvshowfiberconnectivities(resultfig,fibersfile,seedfile,targetsfile,thresh,sides,options,stimparams,changedstates,mode,showregs,showlabs)
% This function shows fiber-connectivity from a volume defined by a nx3
% point-list (volume). If stimparams.showconnectivities is set, connected
% areas to the volume are also visualized. To do so, the function uses
% inhull.m which is covered by the BSD-license (see below).
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

set(0,'CurrentFigure',resultfig)
colornames=['rbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywk' ...
    'rbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywk' ...
    'rbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywk' ...
    'rbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywk' ...
    'rbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywk' ...
    'rbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywk' ...
    'rbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywkrbgcmywk'];

fprintf('\nCalculating Fibers/Connectivity...\n\n\n');

hold on

%% loadstate from file
% set defaults
if ~options.prefs.env.dev || ~isfield(options,'savefibers')
    options.savefibers.save = 0;
    options.savefibers.load = 0;
end

if ~options.savefibers.load
    % run as normal...
    %% load fibers (either from file or from figure and store in figure for next time).
    % get app data
    disp('Loading fiberset...');
    if ~changedstates(1) && ~isempty(getappdata(resultfig,'fibers')) % fibers file already loaded
        fibers=getappdata(resultfig,'fibers');
        idxv=getappdata(resultfig,'idxv');
        fibers=fibers(:,1:3);
        fibersidx=getappdata(resultfig,'fibersidx');
    else % load data
        if ischar(fibersfile)
            [fibers,fibersidx]=ea_loadfibertracts(fibersfile);
        else
            fibers=fibersfile.fibers;
            fibersidx=fibersfile.fibersidx;
            clear fibersfile
        end
        idxv=fibers(:,4);
        fibers=fibers(:,1:3);
        setappdata(resultfig,'fibersidx',fibersidx);
        setappdata(resultfig,'fibers',fibers);
        setappdata(resultfig,'idxv',idxv);
        %setappdata(resultfig,'fibersfile',fibersfile);
    end
    numtotalfibs=length(fibersidx);

    if strcmp(thresh,'auto')
        switch mode
            case 'vat'
                thresh=round(numtotalfibs/20000);
            case 'mat'
                thresh=round(numtotalfibs/5000);
        end
    else
        thresh=round(str2double(thresh));
    end
    fprintf('Done.\n\n\n');

    %% load seed definition
    if ~changedstates(2)
        seed=getappdata(resultfig,'seed');
    else % load data
        if iscell(seedfile)
            seed = cell(1, length(seedfile));
            for s=1:length(seedfile)
                seed{s}=ea_load_nii(seedfile{s});
            end
        else
            seed{1}=seedfile;
            clear seedfile
        end
        setappdata(resultfig,'seed',seed);
    end

    %% load targets definition
    if ~changedstates(3)
        targets=getappdata(resultfig,'targets');
    else % load data
        if ischar(targetsfile)
            targets=ea_load_nii(targetsfile);
        else
            targets=targetsfile;
            clear targetfile
        end
        setappdata(resultfig,'targets',targets);
    end
    origtargets=targets; % original targets map.

    %% prepare fibers
    % ea_dispercent(0,'Preparing fibers');
    % [idx,~]=cellfun(@size,fibers);
    %
    % fibers=cell2mat(fibers');
    % idxv=zeros(size(fibers,1),1);
    % lid=1; cnt=1;
    % for id=idx
    %     ea_dispercent(cnt/length(idx));
    %     idxv(lid:lid+id-1)=cnt;
    %     lid=lid+id;
    %     cnt=cnt+1;
    % end
    % ea_dispercent(1,'end');

    %% select fibers that traverse through seed voxels
    [seed_fv,volume]=ea_fvseeds(seed,options);

    ea_dispercent(0,'Selecting connecting fibers...');
    cnt=1;
    selectedfibs = cell(1, length(seed));
    for side=1:length(seed)
        % adhusch: pre-filtering fibers outside "euclidian hull" around seed for further speed-up
        [~,D] = knnsearch(mean(seed_fv{side}.vertices),fibers);
        maxD = max(pdist2(mean(seed_fv{side}.vertices),seed_fv{side}.vertices));
        prefiltIdx = (D <= maxD); % has length all fibers

        % adhusch: filter fiberes with respect to the real seed shape
        filtPrefiltIdx = inpolyhedron(seed_fv{side}, fibers(prefiltIdx,:), 'flipnormals', true); % massive speed up compared to ea_intriangulation (approx. factor 10?)

        % filtPrefiltIdx has only length of prefiltered fibre subset!
        filtIdx = prefiltIdx;
        filtIdx(prefiltIdx) = filtPrefiltIdx; % filter of all fibers

        selectedfibs{side}=unique(idxv(filtIdx)); % mapping individual points back to fibers
        ea_dispercent(cnt/length(sides));
        cnt=cnt+1;
    end

    ea_dispercent(1,'end');
    connectingfibs=cell(2,1);

    %% reformat fibers
    disp('Reformating fibers...');
    fibers=mat2cell(fibers,fibersidx,3)';
    fprintf('Done.\n\n\n');
    sideselectedfibs = cell(1, length(seed));
    for side=1:length(seed)
        try
            sideselectedfibs{side}=unique(cell2mat(selectedfibs(:,side)));
        end

        try
            connectingfibs{side}=fibers(sideselectedfibs{side});
        end
    end

    %% check which areas are connected to VAT by fibers:
    doubleconnectingfibs=cell(2,1);

    la = 1;
    todelete = cell(1, length(seed));
    howmanyfibs = cell(1, length(seed));
    tareas = cell(1, length(seed));
    contargets = cell(1, length(seed));
    for side=1:length(seed)

        todelete{side}=[];
        cnt=1; % reset cnt.

        %% extract areas connected by fibers.
        if options.writeoutpm && ~exist('pm','var')
            pm=targets;
            pm.img(:)=0;
        end

        [~,labelname]=fileparts(targets.fname);
        aID = fopen(fullfile(ea_space([],'labeling'),[labelname,'.txt']));
        atlas_lgnd=textscan(aID,'%d %s');
        allcareas=[];

        fibmax=length(connectingfibs{side});
        ea_dispercent(0,'Gathering region information');
        for fib=1:fibmax
            ea_dispercent(fib/fibmax);

            thisfibendpoints=[connectingfibs{side}{fib}(1,1:3);connectingfibs{side}{fib}(end,1:3)];
            thisfibendpoints=targets.mat\[thisfibendpoints,ones(2,1)]'; % mm 2 vox
            thisfibendpoints=double(thisfibendpoints(1:3,:));

            conareas=spm_sample_vol(targets,thisfibendpoints(1,:),thisfibendpoints(2,:),thisfibendpoints(3,:),0);
            if any(conareas)
                doubleconnectingfibs{side}{la,cnt}=connectingfibs{side}{fib};
                todelete{side}=[todelete{side},fib];
                cnt=cnt+1;
            end

            if options.writeoutpm
                try
                    pm.img(round(thisfibendpoints(1,1)),round(thisfibendpoints(2,1)),round(thisfibendpoints(3,1)))=...
                        pm.img(round(thisfibendpoints(1,1)),round(thisfibendpoints(2,1)),round(thisfibendpoints(3,1)))+1;
                    pm.img(round(thisfibendpoints(1,2)),round(thisfibendpoints(2,2)),round(thisfibendpoints(3,2)))=...
                        pm.img(round(thisfibendpoints(1,2)),round(thisfibendpoints(2,2)),round(thisfibendpoints(3,2)))+1;
                end
            end
            allcareas=[allcareas, conareas];
        end
        allcareas=round(allcareas);
        ea_dispercent(1, 'end');

        atlength=length(atlas_lgnd{1});
        howmanyfibs{side}=zeros(atlength,1);
        tareas{side}=[];

        tcnt=1;
        for reg=1:atlength
            howmanyfibs{side}(reg)=sum(allcareas==reg); % how many fibers connect VAT and anat. region.
            if howmanyfibs{side}(reg)>(thresh(1))
                tareas{side}(tcnt)=reg;
                tcnt=tcnt+1;
            end
        end
        tareas{side}=unique(tareas{side});


        % Write out connectivity stats
        if options.writeoutstats
            load(options.subj.stats, 'ea_stats');
            % assign the place where to write stim stats into struct
            if isfield(options,'groupmode')
                if options.groupmode
                    stimparams.label=['gs_',options.groupid];
                end
                
            end
            [ea_stats,thisstim]=ea_assignstimcnt(ea_stats,stimparams);
            ea_stats.stimulation(thisstim).ft(side).fibercounts{la}=howmanyfibs{side}/numtotalfibs;

            ea_stats.stimulation(thisstim).ft(side).nfibercounts{la}=ea_stats.stimulation(thisstim).ft(side).fibercounts{la}/volume{side};
            ea_stats.stimulation(thisstim).ft(side).labels{la}=atlas_lgnd{2};
            save(options.subj.stats, 'ea_stats');
        end

        contargets{side}=round(targets.img);
        otargets=contargets{side};
        contargets{side}(:)=0;
        for target=1:atlength
            contargets{side}(otargets==target)=howmanyfibs{side}(target);
        end
        %targets.img(targets.img<thresh)=0;
    end

    if options.savefibers.save && exist(options.savefibers.dir,'dir')
        wsp = whos;
        argin = {'resultfig','fibersfile','seedfile','targetsfile','sides','options','stimparams','changedstates','mode','showregs','showlabs'};
        vars = strcat('''',setdiff({wsp(:).name},argin),''',');
        vars = [vars{:}];
        str2eval = sprintf('save([''%s'',''%s''],%s)',options.savefibers.dir,'workspace.mat',vars(1:end-1));
        eval(str2eval)
    end

elseif options.savefibers.load
    load([options.savefibers.dir,'workspace.mat'])
end

addht = getappdata(resultfig,'addht');
if isempty(addht)
    addht = uitoolbar(resultfig);
    setappdata(resultfig, 'addht', addht);
else
    togglebtns = findobj(get(addht, 'Children'), 'Type', 'uitoggletool');
    toogletags = arrayfun(@(obj) get(obj, 'Tag'), togglebtns, 'Uni', 0);
    toggleset = startsWith(toogletags, {'seedbtn', 'regionbtn', 'labelbtn', 'fibbtn'});
    delete(togglebtns(toggleset))
end

for side=1:length(seed)
    % always show seed patch (i.e. VAT)
    PL.matseedsurf{side}=ea_showseedpatch(resultfig,seed{side},seed{side}.img,options);

    [PL.matsurf{side},PL.conlabels{side}]=ea_showconnectivitypatch(resultfig,targets,contargets{side},thresh,atlas_lgnd{2},tareas{side},[],showregs,showlabs);

    switch side
        case 1
            seedtooltip = 'Seed Regions - Right Side';
            regtooltip = 'Connected Regions - Right Side';
            labeltooltip = 'Region Labels - Right Side';
            seedtag = 'seedbtn_right';
            regtag = 'regionbtn_right';
            labeltag = 'labelbtn_right';
        case 2
            seedtooltip = 'Seed Regions - Left Side';
            regtooltip = 'Connected Regions - Left Side';
            labeltooltip = 'Region Labels - Left Side';
            seedtag = 'seedbtn_left';
            regtag = 'regionbtn_left';
            labeltag = 'labelbtn_left';
        otherwise
            seedtooltip = ['Seed Regions - ', num2str(side)];
            regtooltip = ['Connected Regions - ', num2str(side)];
            labeltooltip = ['Region Labels - ', num2str(side)];
            seedtag = ['seedbtn_', num2str(side)];
            regtag = ['regionbtn_', num2str(side)];
            labeltag = ['labelbtn_', num2str(side)];
    end

    seedbtn=uitoggletool(addht,'CData',ea_get_icn('vat'),...
                           'TooltipString',seedtooltip,...
                           'OnCallback',{@objvisible,PL.matseedsurf{side}},...
                           'OffCallback',{@objinvisible,PL.matseedsurf{side}},...
                           'State','on',...
                           'Tag',seedtag);

    regionbtn=uitoggletool(addht,'CData',ea_get_icn('connectivities'),...
                           'TooltipString',regtooltip,...
                           'OnCallback',{@objvisible,PL.matsurf{side}},...
                           'OffCallback',{@objinvisible,PL.matsurf{side}},...
                           'State','on',...
                           'Tag',regtag);
    labelbtn=uitoggletool(addht,'CData',ea_get_icn('labels'),...
                           'TooltipString',labeltooltip,...
                           'OnCallback',{@objvisible,PL.conlabels{side}},...
                           'OffCallback',{@objinvisible,PL.conlabels{side}},...
                           'State','on',...
                           'Tag',labeltag);

    clear allcareas conareas
end
clear tareas

% Write out probability map of fiber terminals
if options.writeoutpm
    pm.dt(1) = 16;
    pm.fname=[options.root,options.patientname,filesep,'ea_pm','.nii'];
    spm_write_vol(pm,pm.img);
end

% plot fibers that do connect to seed:
for iside=1:length(options.sides)
    side=options.sides(iside);
    if ~isempty(connectingfibs{side})
        % Remove single point
        single = cellfun(@(x) all(size(x)==[1,3]),connectingfibs{side});
        connectingfibs{side}(single)=[];

        fibmax=length(connectingfibs{side});

        if fibmax>options.prefs.d3.maxfibers % if too many fibers are selected, reduce amount of them.
            connectingfibs{side}=connectingfibs{side}(round(linspace(1,length(connectingfibs{side}),options.prefs.d3.maxfibers)));
            fibmax=options.prefs.d3.maxfibers;
        end

        ea_dispercent(0,'Plotting fibers that connect to seed');

        for fib=1:fibmax
            ea_dispercent(fib/fibmax);
            %for segment=1:length(connectingfibs{fib})-1
            connectingfibs{side}{la,fib}=connectingfibs{side}{la,fib}';

            if ~isfield(stimparams,'group')
                connectingfibs{side}{la,fib}(4:6,:)=detcolor(connectingfibs{side}{la,fib}); % add coloring information to the 4th-6th column.
            else % if more than one group is analyzed, coloring info will be off the group color.
                RGB=zeros(1,3);
                RGB(:,1)=stimparams(1).groupcolors(stimparams(1).group,1);
                RGB(:,2)=stimparams(1).groupcolors(stimparams(1).group,2);
                RGB(:,3)=stimparams(1).groupcolors(stimparams(1).group,3);
                connectingfibs{side}{la,fib}(4:6,:) = repmat(RGB, size(connectingfibs{side}{fib},2), 1)';
            end

            for dim=1:size(connectingfibs{side}{fib},1)
                thisfib(dim,:)=double(interp1q((1:size(connectingfibs{side}{fib},2))',connectingfibs{side}{fib}(dim,:)',(1:0.1:size(connectingfibs{side}{fib},2))')');
            end
            % plot fibers
            [PL.fib_plots.fibs(side,fib),fv(fib)]=ea_plot_fiber(thisfib,6,0,options);

            % store for webexport
            jetlist=jet;
            try
                PL.bbfibfv(fib).vertices=thisfib(1:3,:)';
                PL.bbfibfv(fib).faces=[1:size(thisfib,2)-1;2:size(thisfib,2)]';
                PL.bbfibfv(fib).normals=zeros(size(PL.bbfibfv(fib).vertices,1),3);
                PL.bbfibfv(fib).colors=[squeeze(ind2rgb(round(thisfib(4,:)),jetlist)),repmat(0.7,size(thisfib,2),1)];
            end

            clear thisfib
        end

        if strcmp(options.prefs.d3.fiberstyle,'tube')
            fv=ea_concatfv(fv);
            set(0,'CurrentFigure',resultfig);
            PL.fib_plots.fibs(side,1)=patch(fv,'Facecolor', 'interp', 'EdgeColor', 'none','FaceAlpha',0.2);
            set(PL.fib_plots.fibs(side,1),'FaceVertexCData', get(PL.fib_plots.fibs(side,1),'FaceVertexCData'));
            PL.fib_plots.fibs(:,2:end)=[];
        else
            ea_dispercent(1,'end');
        end

        if strcmp(options.prefs.d3.fiberstyle,'line')
                    fibInd = ishandle(PL.fib_plots.fibs(side,:));
            if verLessThan('matlab','8.4') % ML <2014b support
                set(PL.fib_plots.fibs(side,logical(PL.fib_plots.fibs(side,fibInd))),'EdgeAlpha',0.05);
            else
                try
                    set(PL.fib_plots.fibs(side,fibInd),'EdgeAlpha',0.2);
                    set(PL.fib_plots.fibs(side,fibInd),'FaceLighting','phong');
                    set(PL.fib_plots.fibs(side,fibInd),'MarkerSize',0.01);
                    set(PL.fib_plots.fibs(side,fibInd),'LineWidth',0.2);
                    set(PL.fib_plots.fibs(side,fibInd), 'SpecularColorReflectance', 0);
                    set(PL.fib_plots.fibs(side,fibInd), 'SpecularExponent', 5);
                    set(PL.fib_plots.fibs(side,fibInd), 'SpecularStrength', 0.5)
                    set(PL.fib_plots.fibs(side,fibInd),'FaceAlpha',0);
                    set(PL.fib_plots.fibs(side,fibInd),'Tag',sprintf('Fiber%d',side));
                end
            end
        end

        switch side
            case 1
                fibtooltip = 'Connected Fibers - Right Side';
                fibtag = 'fibbtn_right';
            case 2
                fibtooltip = 'Connected Fibers - Left Side';
                fibtag = 'fibbtn_left';
            otherwise
                fibtooltip = ['Connected Fibers - ', num2str(side)];
                fibtag = ['fibbtn_', num2str(side)];

        end

        fibbtn=uitoggletool(addht,'CData',ea_get_icn('fibers_vat'),...
                            'TooltipString',fibtooltip,...
                            'OnCallback',{@objvisible,PL.fib_plots.fibs(side,:)},...
                            'OffCallback',{@objinvisible,PL.fib_plots.fibs(side,:)},...
                            'State','on',...
                            'Tag',fibtag);
    end
end

% plot seed surface:
setappdata(resultfig, [mode,'PL'], PL);


function objvisible(hobj, evt, obj)
set(obj, 'Visible', 'on');


function objinvisible(hobj, evt, obj)
set(obj, 'Visible', 'off');


function [fv,volume]=ea_fvseeds(seed,options)
volume = cell(1, length(seed));
fv = cell(1, length(seed));
for s=1:length(seed)
    volume{s}=sum(seed{s}.img(:))*seed{s}.mat(1)*seed{s}.mat(6)*seed{s}.mat(11);
    fv{s}=isosurface(permute(seed{s}.img,[2,1,3]),0.3);
    fv{s}.vertices=[fv{s}.vertices,ones(size(fv{s}.vertices,1),1)]';
    fv{s}.vertices=seed{s}.mat*fv{s}.vertices;
    fv{s}.vertices=fv{s}.vertices(1:3,:)';
end


function indcol=detcolor(mat) % determine color based on traversing direction.

xyz=abs(diff(mat,1,2));
rgb=xyz/max(xyz(:));

rgb=[rgb,rgb(:,end)];
rgbim=zeros(1,size(rgb,2),3);
rgbim(1,:,:)=rgb';
try
    indcol=rgb2ind(rgbim,jet);
    indcol=squeeze(ind2rgb(indcol,jet))';
catch
    keyboard
end


function in = inhull(testpts,xyz,tess,tol)

% Copyright (c) 2009, John D'Errico
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% inhull: tests if a set of points are inside a convex hull
% usage: in = inhull(testpts,xyz)
% usage: in = inhull(testpts,xyz,tess)
% usage: in = inhull(testpts,xyz,tess,tol)
%
% arguments: (input)
%  testpts - nxp array to test, n data points, in p dimensions
%       If you have many points to test, it is most efficient to
%       call this function once with the entire set.
%
%  xyz - mxp array of vertices of the convex hull, as used by
%       convhulln.
%
%  tess - tessellation (or triangulation) generated by convhulln
%       If tess is left empty or not supplied, then it will be
%       generated.
%
%  tol - (OPTIONAL) tolerance on the tests for inclusion in the
%       convex hull. You can think of tol as the distance a point
%       may possibly lie outside the hull, and still be perceived
%       as on the surface of the hull. Because of numerical slop
%       nothing can ever be done exactly here. I might guess a
%       semi-intelligent value of tol to be
%
%         tol = 1.e-13*mean(abs(xyz(:)))
%
%       In higher dimensions, the numerical issues of floating
%       point arithmetic will probably suggest a larger value
%       of tol.
%
%       DEFAULT: tol = 0
%
% arguments: (output)
%  in  - nx1 logical vector
%        in(i) == 1 --> the i'th point was inside the convex hull.
%
% Example usage: The first point should be inside, the second out
%
%  xy = randn(20,2);
%  tess = convhulln(xy);
%  testpoints = [ 0 0; 10 10];
%  in = inhull(testpoints,xy,tess)
%
% in =
%      1
%      0
%
% A non-zero count of the number of degenerate simplexes in the hull
% will generate a warning (in 4 or more dimensions.) This warning
% may be disabled off with the command:
%
%   warning('off','inhull:degeneracy')
%
% See also: convhull, convhulln, delaunay, delaunayn, tsearch, tsearchn
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 3.0
% Release date: 10/26/06

% get array sizes
% m points, p dimensions
p = size(xyz,2);
[n,c] = size(testpts);
if p ~= c
    error 'testpts and xyz must have the same number of columns'
end
if p < 2
    error 'Points must lie in at least a 2-d space.'
end

% was the convex hull supplied?
if (nargin<3) || isempty(tess)
    tess = convhulln(xyz);
end
[nt,c] = size(tess);
if c ~= p
    error 'tess array is incompatible with a dimension p space'
end

% was tol supplied?
if (nargin<4) || isempty(tol)
    tol = 0;
end

% build normal vectors
switch p
    case 2
        % really simple for 2-d
        nrmls = (xyz(tess(:,1),:) - xyz(tess(:,2),:)) * [0 1;-1 0];

        % Any degenerate edges?
        del = sqrt(sum(nrmls.^2,2));
        degenflag = (del<(max(del)*10*eps));
        if sum(degenflag)>0
            warning('inhull:degeneracy',[num2str(sum(degenflag)), ...
                ' degenerate edges identified in the convex hull'])

            % we need to delete those degenerate normal vectors
            nrmls(degenflag,:) = [];
            nt = size(nrmls,1);
        end
    case 3
        % use vectorized cross product for 3-d
        ab = xyz(tess(:,1),:) - xyz(tess(:,2),:);
        ac = xyz(tess(:,1),:) - xyz(tess(:,3),:);
        nrmls = cross(ab,ac,2);
        degenflag = false(nt,1);
    otherwise
        % slightly more work in higher dimensions,
        nrmls = zeros(nt,p);
        degenflag = false(nt,1);
        for i = 1:nt
            % just in case of a degeneracy
            % Note that bsxfun COULD be used in this line, but I have chosen to
            % not do so to maintain compatibility. This code is still used by
            % users of older releases.
            %  nullsp = null(bsxfun(@minus,xyz(tess(i,2:end),:),xyz(tess(i,1),:)))';
            nullsp = null(xyz(tess(i,2:end),:) - repmat(xyz(tess(i,1),:),p-1,1))';
            if size(nullsp,1)>1
                degenflag(i) = true;
                nrmls(i,:) = NaN;
            else
                nrmls(i,:) = nullsp;
            end
        end
        if sum(degenflag)>0
            warning('inhull:degeneracy',[num2str(sum(degenflag)), ...
                ' degenerate simplexes identified in the convex hull'])

            % we need to delete those degenerate normal vectors
            nrmls(degenflag,:) = [];
            nt = size(nrmls,1);
        end
end

% scale normal vectors to unit length
nrmllen = sqrt(sum(nrmls.^2,2));
% again, bsxfun COULD be employed here...
%  nrmls = bsxfun(@times,nrmls,1./nrmllen);
nrmls = nrmls.*repmat(1./nrmllen,1,p);

% center point in the hull
center = mean(xyz,1);

% any point in the plane of each simplex in the convex hull
a = xyz(tess(~degenflag,1),:);

% ensure the normals are pointing inwards
% this line too could employ bsxfun...
%  dp = sum(bsxfun(@minus,center,a).*nrmls,2);
dp = sum((repmat(center,nt,1) - a).*nrmls,2);
k = dp<0;
nrmls(k,:) = -nrmls(k,:);

% We want to test if:  dot((x - a),N) >= 0
% If so for all faces of the hull, then x is inside
% the hull. Change this to dot(x,N) >= dot(a,N)
aN = sum(nrmls.*a,2);

% test, be careful in case there are many points
in = false(n,1);

% if n is too large, we need to worry about the
% dot product grabbing huge chunks of memory.
memblock = 1e6;
blocks = max(1,floor(n/(memblock/nt)));
aNr = repmat(aN,1,length(1:blocks:n));
for i = 1:blocks
    j = i:blocks:n;
    if size(aNr,2) ~= length(j)
        aNr = repmat(aN,1,length(j));
    end
    in(j) = all((nrmls*testpts(j,:)' - aNr) >= -tol,1)';
end


function str=getstate(val)
switch val
    case 1
        str='on';
    case 0
        str='off';
end


function hn=ea_arrow3(p1,p2,s,w,h,ip,alpha,beta)
% ARROW3 (R13)
%   ARROW3(P1,P2) draws lines from P1 to P2 with directional arrowheads.
%   P1 and P2 are either nx2 or nx3 matrices.  Each row of P1 is an
%   initial point, and each row of P2 is a terminal point.
%
%   ARROW3(P1,P2,S,W,H,IP,ALPHA,BETA) can be used to specify properties
%   of the line, initial point marker, and arrowhead.  S is a character
%   string made with one element from any or all of the following 3
%   columns:
%
%     Color Switches      LineStyle            LineWidth
%     ------------------  -------------------  --------------------
%     k  blacK (default)  -  solid (default)   0.5 points (default)
%     y  Yellow           :  dotted            0   no lines
%     m  Magenta          -. dashdot           /   LineWidthOrder
%     c  Cyan             -- dashed
%     r  Red              *  LineStyleOrder            _______ __
%     g  Green                                       ^        |
%     b  Blue                                       / \       |
%     w  White                        Arrowhead    /   \   Height
%     a  Asparagus                                /     \     |
%     d  Dark gray                               /       \    |
%     e  Evergreen                              /___   ___\ __|__
%     f  Firebrick                             |    | |    |
%     h  Hot pink                              |-- Width --|
%     i  Indigo                                |    | |    |
%     j  Jade                                       | |
%     l  Light gray                                 | |
%     n  Nutbrown                                   | |
%     p  Pear                                       | |
%     q  kumQuat                      Line       -->| |<--LineWidth
%     s  Sky blue                                   | |
%     t  Tawny                                      | |
%     u  bUrgundy                                   | |
%     v  Violet                                     | |
%     z  aZure                                      | |
%     x  random                       Initial      /   \
%     o  colorOrder                   Point    -->(     )<--IP
%     |  magnitude                    Marker       \_ _/
%
%     -------------------------------------------------------------
%                          Color Equivalencies
%     -------------------------------------------------------------
%     ColorOrder     Arrow3         |     Simulink       Arrow3
%     ----------     ----------     |     ----------     ----------
%     Color1         Blue           |     LightBlue      aZure
%     Color2         Evergreen      |     DarkGreen      Asparagus
%     Color3         Red            |     Orange         kumQuat
%     Color4         Sky blue       |     Gray           Light gray
%     Color5         Violet         |
%     Color6         Pear           |
%     Color7         Dark gray      |
%     -------------------------------------------------------------
%
%   The components of S may be specified in any order.  Invalid
%   characters in S will be ignored and replaced by default settings.
%
%   Prefixing the color code with '_' produces a darker shade, e.g.
%   '_t' is dark tawny; prefixing the color code with '^' produces a
%   lighter shade, e.g. '^q' is light kumquat.  The relative brightness
%   of light and dark color shades is controlled by the scalar parameter
%   BETA.  Color code prefixes do not affect black (k), white (w), or
%   the special color switches (xo|).
%
%   ColorOrder may be achieved in two fashions:  The user may either
%   set the ColorOrder property (using RGB triples) or define the
%   global variable ColorOrder (using a string of valid color codes).
%   If the color switch is specified with 'o', and the global variable
%   ColorOrder is a string of color codes (color switches less 'xo|',
%   optionally prefixed with '_' or '^'), then the ColorOrder property
%   will be set to the sequence of colors indicated by the ColorOrder
%   variable.  The color sequence 'bersvpd' matches the default
%   ColorOrder property.  If the color switch is specified with 'o', and
%   the global variable ColorOrder is empty or invalid, then the current
%   ColorOrder property will be used.  Note that the ColorOrder variable
%   takes precedence over the ColorOrder property.
%
%   The magnitude color switch is used to visualize vector magnitudes
%   in conjunction with a colorbar.  If the color switch is specified
%   with '|', colors are linearly interpolated from the current ColorMap
%   according to the length of the associated line.  This option sets
%   CLim to [MinM,MaxM], where MinM and MaxM are the minimum and maximum
%   magnitudes, respectively.
%
%   The current LineStyleOrder property will be used if LineStyle is
%   specified with '*'.  MATLAB cycles through the line styles defined
%   by the LineStyleOrder property only after using all colors defined
%   by the ColorOrder property.  If however, the global variable
%   LineWidthOrder is defined, and LineWidth is specified with '/',
%   then each line will be drawn with sequential color, linestyle, and
%   linewidth.
%
%   W (default = 1) is a vector of arrowhead widths; use W = 0 for no
%   arrowheads.  H (default = 3W) is a vector of arrowhead heights.  If
%   vector IP is neither empty nor negative, initial point markers will
%   be plotted with diameter IP; for default diameter W, use IP = 0.
%   The units of W, H and IP are 1/72 of the PlotBox diagonal.
%
%   ALPHA (default = 1) is a vector of FaceAlpha values ranging between
%   0 (clear) and 1 (opaque).  FaceAlpha is a surface (arrowhead and
%   initial point marker) property and does not affect lines.  FaceAlpha
%   is not supported for 2D rendering.
%
%   BETA (default = 0.4) is a scalar that controls the relative
%   brightness of light and dark color shades, ranging between 0 (no
%   contrast) and 1 (maximum contrast).
%
%   Plotting lines with a single color, linestyle, and linewidth is
%   faster than plotting lines with multiple colors and/or linestyles.
%   Plotting lines with multiple linewidths is slower still.  ARROW3
%   chooses renderers that produce the best screen images; exported
%   or printed plots may benefit from different choices.
%
%   ARROW3(P1,P2,S,W,H,'cone',...) will plot cones with bases centered
%   on P1 in the direction given by P2.  In this instance, P2 is
%   interpreted as a direction vector instead of a terminal point.
%   Neither initial point markers nor lines are plotted with the 'cone'
%   option.
%
%   HN = ARROW3(P1,P2,...) returns a vector of handles to line and
%   surface objects created by ARROW3.
%
%   ARROW3 COLORS will plot a table of named colors with default
%   brightness.  ARROW3('colors',BETA) will plot a table of named
%   colors with brightness BETA.
%
%   ARROW3 attempts to preserve the appearance of existing axes.  In
%   particular, ARROW3 will not change XYZLim, View, or CameraViewAngle.
%   ARROW3 does not, however, support stretch-to-fill scaling.  AXIS
%   NORMAL will restore the current axis box to full size and remove any
%   restrictions on the scaling of units, but will likely result in
%   distorted arrowheads and initial point markers.  See
%   (arrow3_messes_up_my_plots.html).
%
%   If a particular aspect ratio or variable limit is required, use
%   DASPECT, PBASPECT, AXIS, or XYZLIM commands before calling ARROW3.
%   Changing limits or aspect ratios after calling ARROW3 may alter the
%   appearance of arrowheads and initial point markers.  ARROW3 sets
%   XYZCLimMode to manual for all plots, sets DataAspectRatioMode to
%   manual for linear plots, and sets PlotBoxAspectRatioMode to manual
%   for log plots and 3D plots.  CameraViewAngleMode is also set to
%   manual for 3D plots.
%
%   ARROW3 UPDATE will restore the appearance of arrowheads and
%   initial point markers that have become corrupted by changes to
%   limits or aspect ratios.  ARROW3('update',SF) will redraw initial
%   point markers and arrowheads with scale factor SF.  If SF has one
%   element, SF scales W, H and IP.  If SF has two elements, SF(1)
%   scales W and IP, and SF(2) scales H.  If SF has three elements,
%   SF(1) scales W, SF(2) scales H, and SF(3) scales IP.  All sizes are
%   relative to the current PlotBox diagonal.
%
%   ARROW3 UPDATE COLORS will update the magnitude coloring of
%   arrowheads, initial point markers, and lines to conform to the
%   current ColorMap.
%
%   HN = ARROW3('update',...) returns a vector of handles to updated
%   objects.
%
%   EXAMPLES:
%
%     % 2D vectors
%     arrow3([0 0],[1 3])
%     arrow3([0 0],[1 2],'-.e')
%     arrow3([0 0],[10 10],'--x2',1)
%     arrow3(zeros(10,2),50*rand(10,2),'x',1,3)
%     arrow3(zeros(10,2),[10*rand(10,1),500*rand(10,1)],'u')
%     arrow3(10*rand(10,2),50*rand(10,2),'x',1,[],1)
%
%     % 3D vectors
%     arrow3([0 0 0],[1 1 1])
%     arrow3(zeros(20,3),50*rand(20,3),'--x1.5',2)
%     arrow3(zeros(100,3),50*rand(100,3),'x',1,3)
%     arrow3(zeros(10,3),[10*rand(10,1),500*rand(10,1),50*rand(10,1)],'a')
%     arrow3(10*rand(10,3),50*rand(10,3),'x',[],[],0)
%
%     % 3D animation
%     t=(0:pi/40:8*pi)'; u=cos(t); v=sin(t);
%     plot3(20*t,u,v)
%     axis([0,600,-1.5,1.5,-1.5,1.5])
%     grid on, view(35,25)
%     hold on
%     pbaspect([1.8,1.4,1])
%     arrow3(zeros(3),diag([500,1.5,1.5]),'l',0.7,[],0)
%     p=[20*t,u,v]; inc=4:1:length(t);
%     p2=p(inc,:); p1=p(inc-1,:);
%     hn=arrow3(p1(1,:),p2(1,:),'0_b',0.7);
%     for i=2:1:length(p1)
%       delete(hn)
%       hn=arrow3(p1(i,:),p2(i,:),'0_b',0.7);
%       pause(0.01)
%     end
%     hold off
%
%     % Cone plot
%     t=(pi/8:pi/8:2*pi)'; p1=[cos(t) sin(t) t]; p2=repmat([0 0 1],16,1);
%     arrow3(p1,p2,'x',2,4,'cone'), hold on
%     plot3(p1(:,1),p1(:,2),p1(:,3)), hold off
%     pause % change cone size
%     arrow3('update',[1,2])
%
%     % Just for fun
%     arrow3(zeros(100,3),50*rand(100,3),'x',8,4,[],0.95)
%     light('position',[-10 -10 -10],'style','local')
%     light('position',[60,60,60]), lighting gouraud
%
%     % ColorOrder variable, color code prefixes, and Beta
%     global ColorOrder, ColorOrder='^ui^e_hq^v';
%     theta=[0:pi/22:pi/2]';
%     arrow3(zeros(12,2),[cos(theta),sin(theta)],'1.5o',1.5,[],[],[],0.5)
%
%     % ColorOrder property, LineStyleOrder, and LineWidthOrder
%     global ColorOrder, ColorOrder=[];
%     set(gca,'ColorOrder',[1,0,0;0,0,1;0.25,0.75,0.25;0,0,0])
%     set(gca,'LineStyleOrder',{'-','--','-.',':'})
%     global LineWidthOrder, LineWidthOrder=[1,2,4,8];
%     w=[1,2,3,4]; h=[4,6,4,2];
%     arrow3(zeros(4,2),[10*rand(4,1),500*rand(4,1)],'o*/',w,h,0)
%
%     % Magnitude coloring
%     colormap spring
%     arrow3(20*randn(20,3),50*randn(20,3),'|',[],[],0)
%     set(gca,'color',0.7*[1,1,1])
%     set(gcf,'color',0.5*[1,1,1]), grid on, colorbar
%     pause % change the ColorMap and update colors
%     colormap hot
%     arrow3('update','colors')
%
%     % LogLog plot
%     set(gca,'xscale','log','yscale','log');
%     axis([1e2,1e8,1e-2,1e-1]); hold on
%     p1=repmat([1e3,2e-2],15,1);
%     q1=[1e7,1e6,1e5,1e4,1e3,1e7,1e7,1e7,1e7,1e7,1e7,1e6,1e5,1e4,1e3];
%     q2=1e-2*[9,9,9,9,9,7,5,4,3,2,1,1,1,1,1]; p2=[q1',q2'];
%     global ColorOrder, ColorOrder=[];
%     set(gca,'ColorOrder',rand(15,3))
%     arrow3(p1,p2,'o'), grid on, hold off
%
%     % SemiLogX plot
%     set(gca,'xscale','log','yscale','linear');
%     axis([1e2,1e8,1e-2,1e-1]); hold on
%     p1=repmat([1e3,0.05],15,1);
%     q1=[1e7,1e6,1e5,1e4,1e3,1e7,1e7,1e7,1e7,1e7,1e7,1e6,1e5,1e4,1e3];
%     q2=1e-2*[9,9,9,9,9,7,5,4,3,2,1,1,1,1,1]; p2=[q1',q2'];
%     arrow3(p1,p2,'x'), grid on, hold off
%
%     % SemiLogY plot
%     set(gca,'xscale','linear','yscale','log');
%     axis([2,8,1e-2,1e-1]); hold on
%     p1=repmat([3,2e-2],17,1);
%     q1=[7,6,5,4,3,7,7,7,7,7,7,7,7,6,5,4,3];
%     q2=1e-2*[9,9,9,9,9,8,7,6,5,4,3,2,1,1,1,1,1]; p2=[q1',q2'];
%     set(gca,'LineStyleOrder',{'-','--','-.',':'})
%     arrow3(p1,p2,'*',1,[],0), grid on, hold off
%
%     % Color tables
%     arrow3('colors')           % default color table
%     arrow3('colors',0.3)       % low contrast color table
%     arrow3('colors',0.5)       % high contrast color table
%
%     % Update initial point markers and arrowheads
%     % relative to the current PlotBox diagonal
%     arrow3('update')           % redraw same size
%     arrow3('update',2)         % redraw double size
%     arrow3('update',0.5)       % redraw half size
%     arrow3('update',[0.5,2,1]) % redraw W half size,
%                                %        H double size, and
%                                %        IP same size
%
%     See also (arrow3_examples.html), (arrow3_messes_up_my_plots.html).

%   Copyright(c)2002-2013 Version 5.15
%     Tom Davis (tdavis@metzgerwillard.com)
%     Jeff Chang

%   Revision History:
%
%     01/15/13 - Use AppData instead of UserData. (TD)
%     07/27/11 - Added animation example. (TD)
%     05/13/09 - Corrected spelling errors (TD)
%     03/16/08 - Updated contact information (TD)
%     10/23/07 - Corrected zero magnitude exclusion (TD)
%     09/08/07 - Added cone plot option; removed adaptive grid
%                spacing; corrected scale factor; removed "nearly"
%                tight limits (TD)
%     07/24/07 - Ignore zero-magnitude input (TD)
%     07/08/07 - Modified named colors to match named Simulink
%                colors; added light and dark shades for basic
%                colors (ymcrgb) (TD)
%     07/01/07 - Modified named colors to match default ColorOrder
%                colors (TD)
%     06/24/07 - Error checking for empty P1, P2 (TD)
%     06/17/07 - Trim colors,W,H,IP,ALPHA to LENGTH(P1) (TD)
%     05/27/07 - Magnitude coloring and documentation revision (TD)
%     03/10/07 - Improved code metrics (TD)
%     02/21/07 - Preserve existing axis appearance;
%                use relative sizes for W, H, and IP;
%                removed version checking; minor bug fixes (TD)
%     01/09/04 - Replaced calls to LINSPACE, INTERP1, and
%                COLORMAP (TD)
%     12/17/03 - Semilog examples, CAXIS support, magnitude
%                coloring, and color updating; use CData instead
%                of FaceColor; minor bug fixes (TD)
%     07/17/03 - Changed 2D rendering from OpenGL to ZBuffer;
%                defined HN for COLORS and UPDATE options (TD)
%     02/27/03 - Replaced calls to RANDPERM, VIEW, REPMAT, SPHERE,
%                and CYLINDER; added ZBuffer for log plots, RESET
%                for CLA and CLF, and ABS for W and H (TD)
%     02/01/03 - Added UPDATE scale factor and MATLAB version
%                checking, replaced call to CROSS (TD)
%     12/26/02 - Added UserData and UPDATE option (TD)
%     11/16/02 - Added more named colors, color code prefix,
%                global ColorOrder, ALPHA , and BETA (TD)
%     10/12/02 - Added global LineWidthOrder,
%                vectorized W, H and IP (TD)
%     10/05/02 - Changed CLF to CLA for subplot support,
%                added ColorOrder and LineStyleOrder support (TD)
%     04/27/02 - Minor log plot revisions (TD)
%     03/26/02 - Added log plot support (TD)
%     03/24/02 - Adaptive grid spacing control to trade off
%                appearance vs. speed based on size of matrix (JC)
%     03/16/02 - Added "axis tight" for improved appearance (JC)
%     03/12/02 - Added initial point marker (TD)
%     03/03/02 - Added aspect ratio support (TD)
%     03/02/02 - Enhanced program's user friendliness (JC)
%                (lump Color, LineStyle, and LineWidth together)
%     03/01/02 - Replaced call to ROTATE (TD)
%     02/28/02 - Modified line plotting,
%                added linewidth and linestyle (TD)
%     02/27/02 - Minor enhancements on 3D appearance (JC)
%     02/26/02 - Minor enhancements for speed (TD&JC)
%     02/26/02 - Optimize PLOT3 and SURF for speed (TD)
%     02/25/02 - Return handler, error handling, color effect,
%                generalize for 2D/3D vectors (JC)
%     02/24/02 - Optimize PLOT3 and SURF for speed (TD)
%     02/23/02 - First release (JC&TD)

%-------------------------------------------------------------------------
% Error Checking
global LineWidthOrder ColorOrder
if nargin<8 || isempty(beta)
    beta=0.4;
end
beta=abs(beta(1));
if nargout
    hn=[];
end
if strcmpi(p1,'colors')                            % plot color table
    if nargin>1
        beta=abs(p2(1));
    end
    ea_LocalColorTable(1,beta);
    return
end
fig=gcf; ax=gca;
if strcmpi(p1,'update'), ad=getappdata(ax,'Arrow3');    % update
    ea_LocalLogCheck(ax);
    if size(ad,2)<13
        error('Invalid AppData');
    end
    setappdata(ax,'Arrow3',[]); sf=[1,1,1]; flag=0;
    if nargin>1
        if strcmpi(p2,'colors'), flag=1;               % update colors
        elseif ~isempty(p2)                            % update surfaces
            sf=p2(1)*sf; n=length(p2(:));
            if n>1
                sf(2)=p2(2);
                if n>2
                    sf(3)=p2(3);
                end
            end
        end
    end
    H=ea_LocalUpdate(fig,ax,ad,sf,flag);
    if nargout
        hn=H;
    end
    return
end
InputError=['Invalid input, type HELP ',upper(mfilename),...
    ' for usage examples'];
if nargin<2
    error(InputError);
end
[r1,c1]=size(p1); [r2,c2]=size(p2);
if c1<2 || c1>3 || r1*r2==0
    error(InputError);
end
if r1~=r2
    error('P1 and P2 must have same number of rows');
end
if c1~=c2
    error('P1 and P2 must have same number of columns');
end
p=sum(abs(p2-p1),2)~=0; cone=0;
if nargin>5 && ~isempty(ip) && strcmpi(ip,'cone')  % cone plot
    cone=1; p=sum(p2,2)~=0;
    if ~any(p)
        error('P2 cannot equal 0');
    end
    set(ax,'tag','Arrow3ConePlot');
elseif ~any(p)
    error('P1 cannot equal P2');
end
if ~all(p)
    warning('Arrow3:ZeroMagnitude','Zero magnitude ignored')
    p1=p1(p,:); p2=p2(p,:); [r1,c1]=size(p1);
end
n=r1; Zeros=zeros(n,1);
if c1==2
    p1=[p1,Zeros];
    p2=[p2,Zeros];
elseif ~any([p1(:,3);p2(:,3)])
    c1=2;
end
L=get(ax,'LineStyleOrder'); C=get(ax,'ColorOrder');
ST=get(ax,'DefaultSurfaceTag'); LT=get(ax,'DefaultLineTag');
EC=get(ax,'DefaultSurfaceEdgeColor');
if strcmp(get(ax,'nextplot'),'add') && strcmp(get(fig,'nextplot'),'add')
    Xr=get(ax,'xlim'); Yr=get(ax,'ylim'); Zr=get(ax,'zlim');
    [xs,ys,xys]=ea_LocalLogCheck(ax); restore=1;
    if xys
        mode='auto';
        if any([p1(:,3);p2(:,3)])
            error('3D log plot not supported');
        end
        if (xs && ~all([p1(:,1);p2(:,1)]>0)) || ...
           (ys && ~all([p1(:,2);p2(:,2)]>0))
            error('Nonpositive log data not supported')
        end
    else
        mode='manual';
        if strcmp(get(ax,'WarpToFill'),'on')
            warning('Arrow3:WarpToFill',['Stretch-to-fill scaling not ',...
                'supported;\nuse DASPECT or PBASPECT before calling ARROW3.']);
        end
    end
    set(ax,'XLimMode',mode,'YLimMode',mode,'ZLimMode',mode,...
        'CLimMode','manual');
else
    restore=0;
    cla reset;
    xys=0;
    set(fig,'nextplot','add');
    if c1==2
        azel=[0,90];
    else
        azel=[-37.5,30];
    end
    setappdata(ax,'Arrow3',[]);
    set(ax,'nextplot','add','View',azel);
end

%-------------------------------------------------------------------------
% Style Control
[vc,cn]=ea_LocalColorTable(0); prefix=''; OneColor=0;
if nargin<3
    [c,ls,lw]=LocalValidateCLSW;% default color, linestyle/width
else
    [c,ls,lw]=LocalValidateCLSW(s);
    if length(c)>1
        if sum('_^'==c(1))
            prefix=c(1);
        end
        c=c(2);
    end
    if c=='x'                              % random named color (less white)
        [ignore,i]=sort(rand(1,23)); c=cn(i,:);        %#ok
    elseif c=='o'                                    % ColorOrder
        if ~isempty(ColorOrder)
            [c,failed]=ea_LocalColorMap(lower(ColorOrder),vc,cn,beta);
            if failed, ColorOrderWarning=['Invalid ColorOrder ',...
                    'variable, current ColorOrder property will be used'];
                warning('Arrow3:ColorOrder',ColorOrderWarning)
            else
                C=c;
            end
        end, c=C;
    elseif c=='|'
        map=get(fig,'colormap');          % magnitude coloring
        M=(p1-p2);
        M=sqrt(sum(M.*M,2));
        minM=min(M);
        maxM=max(M);
        if maxM-minM<1
            minM=0;
        end
        set(ax,'clim',[minM,maxM]); c=ea_LocalInterp(minM,maxM,map,M);
    elseif ~sum(vc==c)
        c='k';
        ColorWarning=['Invalid color switch, ',...
            'default color (black) will be used'];
        warning('Arrow3:Color',ColorWarning)
    end
end
if length(c)==1                                    % single color
    c=ea_LocalColorMap([prefix,c],vc,cn,beta); OneColor=1;
end
set(ax,'ColorOrder',c); c=ea_LocalRepmat(c,[ceil(n/size(c,1)),1]);
if ls~='*'
    set(ax,'LineStyleOrder',ls);
end    % LineStyleOrder
if lw=='/'                                         % LineWidthOrder
    if ~isempty(LineWidthOrder)
        lw=ea_LocalRepmat(LineWidthOrder(:),[ceil(n/length(LineWidthOrder)),1]);
    else
        lw=0.5;
        LineWidthOrderWarning=['Undefined LineWidthOrder, ',...
            'default width (0.5) will be used'];
        warning('Arrow3:LineWidthOrder',LineWidthOrderWarning)
    end
end
if nargin<4 || isempty(w)
    w=1;
end                % width
w=ea_LocalRepmat(abs(w(:)),[ceil(n/length(w)),1]);
if nargin<5 || isempty(h)
    h=3*w;
end              % height
h=ea_LocalRepmat(abs(h(:)),[ceil(n/length(h)),1]);
if nargin>5 && ~isempty(ip) && ~cone               % ip
    ip=ea_LocalRepmat(ip(:),[ceil(n/length(ip)),1]);
    i=find(ip==0); ip(i)=w(i);
else
    ip=-ones(n,1);
end
if nargin<7 || isempty(alpha)
    alpha=1;
end
a=ea_LocalRepmat(alpha(:),[ceil(n/length(alpha)),1]); % FaceAlpha

%-------------------------------------------------------------------------
% Log Plot
if xys
    units=get(ax,'units'); set(ax,'units','points');
    pos=get(ax,'position'); set(ax,'units',units);
    if strcmp(get(ax,'PlotBoxAspectRatioMode'),'auto')
        set(ax,'PlotBoxAspectRatio',[pos(3),pos(4),1]);
    end
    par=get(ax,'PlotBoxAspectRatio');
    set(ax,'DataAspectRatio',[par(2),par(1),par(3)]);
    % map coordinates onto unit square
    q=[p1;p2]; xr=Xr; yr=Yr;
    if xs
        xr=log10(xr);
        q(:,1)=log10(q(:,1));
    end
    if ys
        yr=log10(yr);
        q(:,2)=log10(q(:,2));
    end
    q=q-ea_LocalRepmat([xr(1),yr(1),0],[2*n,1]);
    dx=xr(2)-xr(1); dy=yr(2)-yr(1);
    q=q*diag([1/dx,1/dy,1]);
    q1=q(1:n,:); q2=q(n+1:end,:);
else
    xs=0; ys=0; dx=0; dy=0; xr=0; yr=0;
end

%-------------------------------------------------------------------------
% Line
if ~cone
    set(ax,'DefaultLineTag','arrow3');
    if length(lw)==1
        if lw>0
            if OneColor && ls(end)~='*' && n>1 % single color, linestyle/width
                P=zeros(3*n,3); i=1:n;
                P(3*i-2,:)=p1(i,:); P(3*i-1,:)=p2(i,:); P(3*i,1)=NaN;
                H1=plot3(P(:,1),P(:,2),P(:,3),'LineWidth',lw);
            else                               % single linewidth
                H1=plot3([p1(:,1),p2(:,1)]',[p1(:,2),p2(:,2)]',...
                    [p1(:,3),p2(:,3)]','LineWidth',lw);
            end
        else
            H1=[];
        end
    else                                   % use LineWidthOrder
        ls=ea_LocalRepmat(cellstr(L),[ceil(n/size(L,1)),1]);
        H1=Zeros;
        for i=1:n
            H1(i)=plot3([p1(i,1),p2(i,1)],[p1(i,2),p2(i,2)],...
                [p1(i,3),p2(i,3)],ls{i},'Color',c(i,:),'LineWidth',lw(i));
        end
    end
else                                     % cone plot
    P=zeros(3*n,3); i=1:n;
    P(3*i-2,:)=p1(i,:); P(3*i-1,:)=p1(i,:); P(3*i,1)=NaN;
    H1=plot3(P(:,1),P(:,2),P(:,3));
end

%-------------------------------------------------------------------------
% Scale
if ~restore
    axis tight;
end
ar=get(ax,'DataAspectRatio'); ar=sqrt(3)*ar/norm(ar);
set(ax,'DataAspectRatioMode','manual');
if xys
    sf=1;
else
    xr=get(ax,'xlim');
    yr=get(ax,'ylim');
    zr=get(ax,'zlim');
    sf=norm(diff([xr;yr;zr],1,2)./ar')/72;
end

%-------------------------------------------------------------------------
% AppData
c=c(1:n,:); w=w(1:n); h=h(1:n); ip=ip(1:n); a=a(1:n);
setappdata(ax,'Arrow3',[getappdata(ax,'Arrow3');p1,p2,c,w,h,ip,a]);

%-------------------------------------------------------------------------
% Arrowhead
whip=sf*[w,h,ip];
if xys
    whip=whip*sqrt(2)/72;
    p1=q1;
    p2=q2;
end
w=whip(:,1); h=whip(:,2); ip=whip(:,3);
if cone                                            % cone plot
    delete(H1), H1=[];
    p2=p2./ea_LocalRepmat(sqrt(sum(p2.*p2,2)),[1,3]);
    p2=p1+p2.*ea_LocalRepmat(ar,[n,1]).*ea_LocalRepmat(h,[1,3]);
end
W=(p1-p2)./ea_LocalRepmat(ar,[n,1]);
W=W./ea_LocalRepmat(sqrt(sum(W.*W,2)),[1,3]);         % new z direction
U=[-W(:,2),W(:,1),Zeros];
N=sqrt(sum(U.*U,2)); i=find(N<eps); j=length(i);
U(i,:)=ea_LocalRepmat([1,0,0],[j,1]); N(i)=ones(j,1);
U=U./ea_LocalRepmat(N,[1,3]);                         % new x direction
V=[W(:,2).*U(:,3)-W(:,3).*U(:,2),...               % new y direction
    W(:,3).*U(:,1)-W(:,1).*U(:,3),...
    W(:,1).*U(:,2)-W(:,2).*U(:,1)];

m=20;                               % surface grid spacing
set(ax,'DefaultSurfaceTag','arrow3','DefaultSurfaceEdgeColor','none');
r=[0;1]; theta=(0:m)/m*2*pi; Ones=ones(1,m+1);
x=r*cos(theta); y=r*sin(theta); z=r*Ones;
G=surface(x/2,y/2,z); dar=diag(ar);
X=get(G,'XData'); Y=get(G,'YData'); Z=get(G,'ZData');
H2=Zeros; [j,k]=size(X);
for i=1:n   % translate, rotate, and scale
    H2(i)=copyobj(G,ax);
    xyz=[w(i)*X(:),w(i)*Y(:),h(i)*Z(:)]*[U(i,:);V(i,:);W(i,:)]*dar;
    x=reshape(xyz(:,1),j,k)+p2(i,1);
    y=reshape(xyz(:,2),j,k)+p2(i,2);
    z=reshape(xyz(:,3),j,k)+p2(i,3);
    ea_LocalSetSurface(xys,xs,ys,dx,dy,xr,yr,...
        x,y,z,a(i),c(i,:),H2(i),2,m+1);
end
delete(G);

%-------------------------------------------------------------------------
% Initial Point Marker
if any(ip>0)
    theta=(-m:2:m)/m*pi; phi=(-m:2:m)'/m*pi/2; cosphi=cos(phi);
    x=cosphi*cos(theta); y=cosphi*sin(theta); z=sin(phi)*Ones;
    G=surface(x*ar(1)/2,y*ar(2)/2,z*ar(3)/2);
    X=get(G,'XData'); Y=get(G,'YData'); Z=get(G,'ZData');
    H3=zeros(n,1);
    for i=1:n   % translate
        if ip(i)>0
            H3(i)=copyobj(G,ax);
            x=p1(i,1)+X*ip(i); y=p1(i,2)+Y*ip(i); z=p1(i,3)+Z*ip(i);
            ea_LocalSetSurface(xys,xs,ys,dx,dy,xr,yr,...
                x,y,z,a(i),c(i,:),H3(i),m+1,m+1);
        end
    end, delete(G);
else
    H3=[];
end

%-------------------------------------------------------------------------
% Finish
if restore
    xr=Xr;
    yr=Yr;
    zr=Zr;
    if xys
        set(ax,'DataAspectRatioMode','auto');
    end
else
    axis tight
    xr=get(ax,'xlim'); yr=get(ax,'ylim'); zr=get(ax,'zlim');
    set(ax,'nextplot','replace');
end
azel=get(ax,'view');
if abs(azel(2))==90
    renderer='ZBuffer';
else
    renderer='OpenGL'; c1=3;
end
set(fig,'Renderer',renderer);
set(ax,'LineStyleOrder',L,'ColorOrder',C,'DefaultLineTag',LT,...
    'DefaultSurfaceTag',ST,'DefaultSurfaceEdgeColor',EC,...
    'xlim',xr,'ylim',yr,'zlim',zr,'clim',get(ax,'CLim'));
if c1==3
    set(ax,'CameraViewAngle',get(ax,'CameraViewAngle'),...
        'PlotBoxAspectRatio',get(ax,'PlotBoxAspectRatio'));
end
if nargout
    hn=[H1(:);H2(:);H3(:)];
end


%-------------------------------------------------------------------------
% Local Functions
%-------------------------------------------------------------------------
% Update
function H=ea_LocalUpdate(fig,ax,ad,sf,flag)
global ColorOrder
p1=ad(:,1:3); p2=ad(:,4:6); c=ad(:,7:9); a=ad(:,13);
w=sf(1)*ad(:,10); h=sf(2)*ad(:,11); ip=sf(3)*ad(:,12);
H=get(ax,'children'); tag=get(H,'tag'); type=get(H,'type');
delete(H(strcmp(tag,'arrow3') & strcmp(type,'surface')));
set(fig,'nextplot','add'); set(ax,'nextplot','add'); H1=[];
if flag, map=get(fig,'colormap');                  % update colors
    M=(p1-p2); M=sqrt(sum(M.*M,2)); minM=min(M); maxM=max(M);
    H1=H(strcmp(tag,'arrow3') & strcmp(type,'line'));
    MagnitudeWarning=['Cannot perform magnitude coloring on lines ',...
        'that\nwere drawn with a single color, linestyle, and linewidth'];
    if length(H1)>1
        for i=1:length(H1)  % update line colors
            x=get(H1(i),'xdata'); y=get(H1(i),'ydata'); z=get(H1(i),'zdata');
            if length(x)>2                               % multiple lines
                warning('Arrow3:Magnitude',MagnitudeWarning), continue
            end
            m=sqrt((x(1)-x(2))^2+(y(1)-y(2))^2+(z(1)-z(2))^2);
            c=ea_LocalInterp(minM,maxM,map,m); set(H1(i),'color',c);
        end
    elseif length(H1)==1
        warning('Arrow3:Magnitude',MagnitudeWarning)
    end
    c=ea_LocalInterp(minM,maxM,map,M);
end
set(ax,'ColorOrder',c);                            % update surfaces
ColorOrder=[];
if strcmp(get(ax,'tag'),'Arrow3ConePlot')
    H=ea_arrow3(p1,p2,'o' ,w,h,'cone',a);            % update cones
else
    H=ea_arrow3(p1,p2,'o0',w,h,ip,a);
end
H=[H1(:);H(:)];
set(ax,'nextplot','replace');


%-------------------------------------------------------------------------
% SetSurface
function ea_LocalSetSurface(xys,xs,ys,dx,dy,xr,yr,x,y,z,a,c,H,n,m)
if xys
    x=x*dx+xr(1); y=y*dy+yr(1);
    if xs
        x=10.^x;
    end
    if ys
        y=10.^y;
    end
end
cd=zeros(n,m,3); cd(:,:,1)=c(1); cd(:,:,2)=c(2); cd(:,:,3)=c(3);
set(H,'XData',x,'YData',y,'ZData',z,'CData',cd,'FaceAlpha',a);


%-------------------------------------------------------------------------
% ColorTable
function [vc,cn]=ea_LocalColorTable(n,beta)
vc='kymcrgbadefhijlnpqstuvzw';                     % valid color codes
%                k               y               m               c
cn=[0.00,0.00,0.00; 1.00,1.00,0.00; 1.00,0.00,1.00; 0.00,1.00,1.00;
    %                r               g               b               a
    1.00,0.00,0.00; 0.00,1.00,0.00; 0.00,0.00,1.00; 0.42,0.59,0.24;
    %                d               e               f               h
    0.25,0.25,0.25; 0.00,0.50,0.00; 0.70,0.13,0.13; 1.00,0.41,0.71;
    %                i               j               l               n
    0.29,0.00,0.51; 0.00,0.66,0.42; 0.50,0.50,0.50; 0.50,0.20,0.00;
    %                p               q               s               t
    0.75,0.75,0.00; 1.00,0.50,0.00; 0.00,0.75,0.75; 0.80,0.34,0.00;
    %                u               v               z               w
    0.50,0.00,0.13; 0.75,0.00,0.75; 0.38,0.74,0.99; 1.00,1.00,1.00];

% Named Simulink Colors (zaql)
% LightBlue = 0.38  0.74  0.99 = aZure
% DarkGreen = 0.42  0.59  0.24 = Asparagus
% Orange    = 1.00  0.50  0.00 = kumQuat
% Gray      = 0.50  0.50  0.50 = Light gray
%
% Default ColorOrder Property Colors (bersvpd)
% Color1    = 0.00  0.00  1.00 = Blue
% Color2    = 0.00  0.50  0.00 = Evergreen
% Color3    = 1.00  0.00  0.00 = Red
% Color4    = 0.00  0.75  0.75 = Sky blue
% Color5    = 0.75  0.00  0.75 = Violet
% Color6    = 0.75  0.75  0.00 = Pear
% Color7    = 0.25  0.25  0.25 = Dark gray

if n, clf reset                                    % plot color table
    name={'blacK','Yellow','Magenta','Cyan',...
        'Red','Green','Blue','Asparagus',...
        'Dark gray','Evergreen','Firebrick','Hot pink',...
        'Indigo','Jade','Light gray','Nutbrown',...
        'Pear','kumQuat','Sky blue','Tawny',...
        'bUrgundy','Violet','aZure','White'};
    c=['yptn';'gjae';'czsb';'hmvi';'qrfu';'wldk'];
    set(gcf,'DefaultAxesXTick',[],'DefaultAxesYTick',[],...
        'DefaultAxesXTickLabel',[],'DefaultAxesYTickLabel',[],...
        'DefaultAxesXLim',[0,0.75],'DefaultAxesYLim',[0,0.75],...
        'DefaultRectangleEdgeColor','none');
    for i=1:24
        subplot(4,6,i)
        box on;
        j=find(vc==c(i)); title(name{j});
        dark=ea_LocalBrighten(cn(j,:),-beta);
        light=ea_LocalBrighten(cn(j,:),beta);
        rectangle('Position',[0,0.00,0.75,0.25],'FaceColor',dark);
        rectangle('Position',[0,0.25,0.75,0.25],'FaceColor',cn(j,:));
        rectangle('Position',[0,0.50,0.75,0.25],'FaceColor',light);
        rectangle('Position',[0,0.00,0.75,0.75],'EdgeColor','k');
        if rem(i,6)==1
            set(gca,'YTickLabel',{'dark','normal','light'},...
                'YTick',[0.125,0.375,0.625]);
            if i==19
                text(0,-0.25,['{\bf\itARROW3}  Named Color Table  ',...
                    '( \beta = ',num2str(beta),' )']);
            end
        end
    end
end


%-------------------------------------------------------------------------
% ColorMap
function [C,failed]=ea_LocalColorMap(c,vc,cn,beta)
n=length(c); failed=0; C=zeros(n,3); i=1; j=1;
while 1
    if ~sum([vc,'_^']==c(i))
        failed=1;
        break;
    end
    if sum('_^'==c(i))
        if i+1>n
            failed=1;
            break;
        end
        if ~sum(vc==c(i+1))
            failed=1;
            break;
        end
        cc=cn(vc==c(i+1),:); gamma=beta;
        if c(i)=='_'
            gamma=-beta;
        end
        C(j,:)=ea_LocalBrighten(cc,gamma); i=i+2;
    else
        C(j,:)=cn(vc==c(i),:); i=i+1;
    end
    if i>n
        break;
    end
    j=j+1;
end
if n>j
    C(j+1:n,:)=[];
end


%-------------------------------------------------------------------------
% Brighten
function C=ea_LocalBrighten(c,beta)
if sum([c==0,c==1])==3 && sum(c==0)<3 && sum(c==1)<3
    if beta<0
        C=(1+beta)*c;
    else
        C=c;  C(C==0)=beta;
    end
else
    C=c.^((1-min(1-sqrt(eps),abs(beta)))^sign(beta));
end


%-------------------------------------------------------------------------
% Repmat
function B=ea_LocalRepmat(A,siz)
if length(A)==1
    B(prod(siz))=A;
    B(:)=A;
    B=reshape(B,siz);
else
    [m,n]=size(A);
    mind=(1:m)';
    nind=(1:n)';
    mind=mind(:,ones(1,siz(1)));
    nind=nind(:,ones(1,siz(2)));
    B=A(mind,nind);
end


%-------------------------------------------------------------------------
% Interp
function v=ea_LocalInterp(xmin,xmax,y,u)
[m,n]=size(y); h=(xmax-xmin)/(m-1); p=length(u); v=zeros(p,n);
k=min(max(1+floor((u-xmin)/h),1),m-1); s=(u-xmin)/h-k+1;
for j=1:n
    v(:,j)=y(k,j)+s.*(y(k+1,j)-y(k,j));
end
v(v<0)=0; v(v>1)=1;


%-------------------------------------------------------------------------
% Check for supported log scales
function [xs,ys,xys]=ea_LocalLogCheck(ax)
xs=strcmp(get(ax,'xscale'),'log');
ys=strcmp(get(ax,'yscale'),'log');
zs=strcmp(get(ax,'zscale'),'log');
if zs
    error('Z log scale not supported');
end
xys=xs+ys;
if xys
    azel=get(ax,'view');
    if abs(azel(2))~=90
        error('3D log plot not supported');
    end
end


%-------------------------------------------------------------------------
% Generate valid value for color, linestyle and linewidth
function [c,ls,lw]=LocalValidateCLSW(s)
if nargin<1
    c='k';
    ls='-';
    lw=0.5;
else
    % identify linestyle
    if strfind(s,'--')
        ls='--';
        s=strrep(s,'--','');
    elseif strfind(s,'-.')
        ls='-.';
        s=strrep(s,'-.','');
    elseif strfind(s,'-')
        ls='-';
        s=strrep(s,'-','');
    elseif strfind(s,':')
        ls=':';
        s=strrep(s,':','');
    elseif strfind(s,'*')
        ls='*';
        s=strrep(s,'*','');
    else
        ls='-';
    end

    % identify linewidth
    tmp=double(s);
    tmp=find(tmp>45 & tmp<58);
    if ~isempty(tmp)
        if any(s(tmp)=='/')
            lw='/';
        else
            lw=str2double(s(tmp));
        end
        s(tmp)='';
    else
        lw=0.5;
    end

    % identify color
    if ~isempty(s)
        s=lower(s);
        if length(s)>1
            c=s(1:2);
        else
            c=s(1);
        end
    else
        c='k';
    end
end
