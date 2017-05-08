function [atlases,colorbuttons,atlassurfs]=ea_showatlas(varargin)
% This function shows atlas data in the 3D-Scene viewer. It
% reads in all atlases found in the eAuto_root/atlases folder, calculates a
% convex hull around the nonzero area and renders this area as 3D surfaces.
% For a small part of contact statistics, the function uses
% inhull.m which is covered by the BSD-license (see below).
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

maxcolor=64; % change to 45 to avoid red / 64 to use all colors


resultfig=varargin{1};
if nargin==2
    options=varargin{2};
end
if nargin>2
    elstruct=varargin{2};
    options=varargin{3};
end


nm=[1:2]; % native and mni
try
    nmind=[options.atl.can,options.atl.ptnative]; % which shall be performed?
catch
    nmind=[1 0];
end
nm=nm(logical(nmind)); % select which shall be performed.

mcr=ea_checkmacaque(options);

for nativemni=nm % switch between native and mni space atlases.
    
    switch nativemni
        case 1
            adir=[ea_space(options,'atlases'),options.atlasset,filesep];
            mifix='';
        case 2
            adir=[[options.root,options.patientname,filesep],'atlases',filesep,options.atlasset,filesep];
            mifix='';
    end
    
    
    atlascnt=1;
    set(0,'CurrentFigure',resultfig)
    
    if ~exist([adir,'atlas_index.mat'],'file')
        
        atlases=ea_genatlastable([],fileparts(fileparts(adir)),options,mifix);
    else
        load([adir,'atlas_index.mat']);
        atlases=ea_genatlastable(atlases,fileparts(adir),options,mifix);
    end
    
    
    if options.writeoutstats
        try
            load([options.root,options.patientname,filesep,'ea_stats']);
            prioratlasnames=ea_stats.atlases.names;
        end
    end
    
    
    
    
    if isfield(atlases,'colormap');
        
        try
            jetlist=eval(atlases.colormap);
            
        catch
            jetlist=atlases.colormap;
            
        end
        colormap(atlases.colormap);
        
        
        
    else
        try
            jetlist=options.colormap;
            atlases.colormap=jetlist;
            colormap(jetlist)
        catch
            atlases.colormap='jet';
            jetlist=jet;
        end
    end
    
    
    
    
    setinterpol=1;
    
   ht=getappdata(resultfig,'atlht');
   if ~isempty(ht) % sweep nonempty atlases toolbar
       delete(ht.Children(:));
   else
       ht=uitoolbar(resultfig);
   end
    atlcntbutton=uipushtool(ht,'CData',ea_get_icn('atlases',options),'TooltipString','Atlas Control Figure','ClickedCallback',{@ea_openatlascontrol,atlases,resultfig,options});

    % prepare stats fields
    if options.writeoutstats
        
        for el=1:length(elstruct)
            for side=1:length(elstruct(el).coords_mm)
                ea_stats.conmat{el,side}=nan(size(elstruct(el).coords_mm{side},1),length(atlases.names));
                ea_stats.conmat_inside_vox{el,side}=nan(size(elstruct(el).coords_mm{side},1),length(atlases.names));
                ea_stats.conmat_inside_hull{el,side}=nan(size(elstruct(el).coords_mm{side},1),length(atlases.names));
                ea_stats.patname{el,side}=elstruct(el).name;
            end
        end
        ea_stats.atlases=atlases;
        ea_stats.electrodes=elstruct;
    end
    
    % iterate through atlases, visualize them and write out stats.
    for atlas=1:length(atlases.names)
        
        [~,sidestr]=detsides(atlases.types(atlas));
        for side=detsides(atlases.types(atlas));
            
%             if ischar(atlases.pixdim{atlas,side}) % we are dealing with fibers
%                 
%                 %[~,alnm]=fileparts(atlases.names{atlas});
%                 %                     atlassurfs(atlascnt,fib)=surface([thisfib(:,1);thisfib(:,1)+randn(len,1)*0.1],...
%                 %                         [thisfib(:,2);thisfib(:,2)+randn(len,1)*0.1],...
%                 %                         [thisfib(:,3);thisfib(:,3)+randn(len,1)*0.1],...
%                 %                         [thisfib(:,4);thisfib(:,4)+randn(len,1)*0.1],'facecol','no','edgecol','interp','linew',5);
%                 
%                 
%             else
                
                
                fv=atlases.fv{atlas,side};
                
                
                if ischar(options.prefs.hullsimplify)   % for 'auto' hullsimplify
                    % get to 700 faces
                    simplify=700/length(fv.faces);
                    if simplify < 1 % skip volumes with fewer than 700 faces
                        fv=reducepatch(fv,simplify);
                    end
                    
                else
                    if options.prefs.hullsimplify<1 && options.prefs.hullsimplify>0
                        
                        fv=reducepatch(fv,options.prefs.hullsimplify);
                    elseif options.prefs.hullsimplify>1
                        simplify=options.prefs.hullsimplify/length(fv.faces);
                        fv=reducepatch(fv,simplify);
                    end
                end
                
                rndfactor=1;
                try
                switch atlases.names{atlas,side}(end-2:end)
                    case 'nii'
                        rndfactor=2;
                    case {'trk','mat'}
                        rndfactor=0.2;
                end
                end
                
                try
                    if ~options.prefs.d3.colorjitter
                       rndfactor=0; 
                    end
                end
                    
                cdat=repmat(atlases.colors(atlas),length(fv.vertices),1); % C-Data for surface
                
                if size(cdat,2)==1
                   cdat=atlases.colormap(round(cdat),:); 
                end
               
                
                % add color jitter
                cdat=cdat+(randn(size(cdat,1),3)*rndfactor);
        
                XYZ=atlases.XYZ{atlas,side};
                pixdim=atlases.pixdim{atlas,side};
                colorc=nan;
                
                %
                %             end
                
                
                % show atlas label
                
                if size(XYZ.mm,1)>1 % exception for single-coordinate atlases...
                    try
                    [~,centroid]=kmeans(XYZ.mm(:,1:3),1);
                    catch
                        centroid=mean(XYZ(:,1:3),1);
                    end
                else
                    
                    try
                    centroid=XYZ.mm(:,1:3);
                    catch
                        centroid=[0,0,0];
                        warning('No centroid found.')
                    end
                end
                try
                    centroid=centroid(1,:);
                catch % empty file..
                    break
                end
                
                
                set(0,'CurrentFigure',resultfig);

                atlassurfs(atlascnt,1)=patch(fv,'FaceVertexCData',cdat,'FaceColor','interp','facealpha',0.7,'EdgeColor','none','facelighting','phong');
                
           % end
            
            % export label and labelbutton
            
            [~,thislabel]=fileparts(atlases.names{atlas});
            try % use try here because filename might be shorter than .nii
                
                if strcmp(thislabel(end-3:end),'.nii') % if it was .nii.gz, fileparts will only remove .gz
                    [~,thislabel]=fileparts(thislabel);
                end
            end
            atlaslabels(atlas,side)=text(centroid(1),centroid(2),centroid(3),ea_sub2space(thislabel),'VerticalAlignment','Baseline','HorizontalAlignment','Center','Color','w');
            
           
            if ~exist('labelbutton','var')
                labelbutton=uitoggletool(ht,'CData',ea_get_icn('labels',options),'TooltipString','Labels');
                labelcolorbutton=uipushtool(ht,'CData',ea_get_icn('colors',options),'TooltipString','Label Color');
            end
            
            % make fv compatible for stats
            
            
            
            caxis([1 64]);
            
            % prepare colorbutton icon
            try
                atlasc=squeeze(jetlist(ceil(atlases.colors(atlas)),:));  % color for toggle button icon
            catch
                ea_error('Atlas color not found.');
            end

            colorbuttons(atlascnt)=uitoggletool(ht,'CData',ea_get_icn('atlas',atlasc),'TooltipString',atlases.names{atlas},'ClickedCallback',{@atlasvisible,resultfig,atlascnt},'State','on');
            
            % set Tags
            set(colorbuttons(atlascnt),'tag',[thislabel,'_',sidestr{side}])
            set(atlassurfs(atlascnt,1),'tag',[thislabel,'_',sidestr{side}])
            set(atlassurfs(atlascnt,1),'UserData',atlaslabels(atlas,side))
            
            % gather contact statistics
            if options.writeoutstats
                try
                    
                    if isfield(atlases.XYZ{atlas,side},'val') % volumetric atlas
                    thresh=ea_detthresh(atlases,atlas,atlases.XYZ{atlas,side}.val);
                    
                    atsearch=KDTreeSearcher(XYZ.mm(XYZ.val>thresh,:));
                    else % fibertract
                        atsearch=KDTreeSearcher(XYZ.mm(:,1:3));
                    end
                    for el=1:length(elstruct)
                        
                        [~,D]=knnsearch(atsearch,ea_stats.electrodes(el).coords_mm{side});
                        %s_ix=sideix(side,size(elstruct(el).coords_mm{side},1));
                        
                        ea_stats.conmat{el,side}(:,atlas)=D;
                        Dh=D;
                        
                        try
                            in=inhull(ea_stats.electrodes(el).coords_mm{side},fv.vertices,fv.faces,1.e-13*mean(abs(fv.vertices(:))));
                            Dh(in)=0;
                            
                        end
                        ea_stats.conmat_inside_hull{el,side}(:,atlas)=Dh;
                        
                        D(D<mean(pixdim))=0; % using mean here but assuming isotropic atlases in general..
                        ea_stats.conmat_inside_vox{el,side}(:,atlas)=D;
                        
                        
                        
                    end
                catch
                    warning('Statistics for tract atlas parts are not implemented yet.');
                end
            end
            
            %normals{atlas,side}=get(atlassurfs(atlascnt),'VertexNormals');
            
                ea_spec_atlas(atlassurfs(atlascnt,1),atlases.names{atlas},atlases.colormap,setinterpol);
            atlascnt=atlascnt+1;
            
            set(gcf,'Renderer','OpenGL')
            axis off
            % set(gcf,'color','w');
            axis equal
            
            if rand(1)>0.8 % we don't want to show every buildup step due to speed but want to show some buildup.
                drawnow
            end
            
            
        end
    end
    
    setappdata(resultfig,'atlassurfs',atlassurfs);
    setappdata(resultfig,'colorbuttons',colorbuttons);
    setappdata(resultfig,'atlht',ht);
    % configure label button to work properly and hide labels as default.
    
    atlabelsvisible([],[],atlaslabels(:),'off');
    set(labelbutton,'OnCallback',{@atlabelsvisible,atlaslabels(:),'on'},'OffCallback',{@atlabelsvisible,atlaslabels(:),'off'},'State','off');
    set(labelcolorbutton,'ClickedCallback',{@setlabelcolor,atlaslabels});
    
    
    
    
    %     % save table information that has been generated from nii files (on first run with this atlas set).
    %     try
    %         atlases.fv=ifv;
    %         atlases.cdat=icdat;
    %         atlases.XYZ=iXYZ;
    %         atlases.pixdim=ipixdim;
    %         atlases.colorc=icolorc;
    %         atlases.normals=normals;
    %     end
    
    
    try
        setappdata(gcf,'atlases',atlases);
        %        setappdata(gcf,'iXYZ',atlases.XYZ);
        %        setappdata(gcf,'ipixdim',atlases.pixdim);
    end
    try
        atlases.rebuild=0; % always reset rebuild flag.
        save([adir,options.atlasset,filesep,'atlas_index.mat'],'atlases','-v7.3');
    end
    
    if options.writeoutstats
        if exist('prioratlasnames','var')
            if ~isequal(ea_stats.atlases.names,prioratlasnames)
                warning('Other atlasset used as before. Deleting VAT and Fiberinfo. Saving backup copy.');
                save([options.root,options.patientname,filesep,'ea_stats_new'],'ea_stats');
                load([options.root,options.patientname,filesep,'ea_stats']);
                save([options.root,options.patientname,filesep,'ea_stats_backup'],'ea_stats');
                movefile([options.root,options.patientname,filesep,'ea_stats_new.mat'],[options.root,options.patientname,filesep,'ea_stats.mat']);
            else
                save([options.root,options.patientname,filesep,'ea_stats'],'ea_stats');
            end
        else
            save([options.root,options.patientname,filesep,'ea_stats'],'ea_stats');
            
        end
    end
    
end


% open up atlas control viewer




function setlabelcolor(hobj,ev,robject)

co=uisetcolor;
set(robject,'Color',co);




function C=rgb(C) % returns rgb values for the colors.

C = rem(floor((strfind('kbgcrmyw', C) - 1) * [0.25 0.5 1]), 2);



function atlasvisible(hobj,ev,resultfig,atlscnt,onoff)
if ~exist('onoff','var')
    onoff=hobj.State;
end

atls=getappdata(resultfig,'atlassurfs');

if(getappdata(resultfig,'altpressed'))
    
    cbutn=getappdata(resultfig,'colorbuttons');
    set(cbutn,'State',onoff);
    for el=1:length(atls)
        for side=1:2
            
            for atlshorz=1:size(atls,2)
                try
                    set(atls(atlscnt,atlshorz), 'Visible', onoff);
                end
            end
        end
        
    end
else
    for atlshorz=1:size(atls,2)
        try
            set(atls(atlscnt,atlshorz), 'Visible', onoff);
        end
    end
end



% check if new atlas select window is open:
figHandles = findobj('Type','figure');
atlspres=0;
for f=1:length(figHandles)
    if strcmp(figHandles(f).Tag,'atlasselect')
        atlspres=1;
        break
    end
end
if atlspres
   atfig=figHandles(f);
   clear figHandles
    handles=getappdata(atfig,'handles');
    ea_synctree(handles)
end



function atlabelsvisible(hobj,ev,obj,onoff)
for el=1:numel(obj)
    try
        set(obj(el),'Visible',onoff);
    end
end


function oldatlasvisible(hobj,ev,atls)
set(atls, 'Visible', 'on');
%disp([atls,'visible clicked']);
function atlasinvisible(hobj,ev,atls)
set(atls, 'Visible', 'off');
%disp([atls,'invisible clicked']);

function [sides,sidestr]=detsides(opt)

switch opt
    case 1 % right hemispheric atlas
        sides=1;
        sidestr={'right'};
    case 2 % left hemispheric atlas
        sides=2;
        sidestr={[''],'left'};
    case 3
        sides=1:2;
        sidestr={'right','left'};
    case 4
        sides=1:2;
        sidestr={'right','left'};
    case 5
        sides=1; % midline
        sidestr={'midline'};
        
end



function coords=map_coords_proxy(XYZ,V)

XYZ=[XYZ';ones(1,size(XYZ,1))];

coords=V.mat*XYZ;
coords=coords(1:3,:)';


% function six=sideix(side,howmany)
% howmany=howmany/2;
% six=side*howmany-(howmany-1):side*howmany;


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
    if size(aNr,2) ~= length(j),
        aNr = repmat(aN,1,length(j));
    end
    in(j) = all((nrmls*testpts(j,:)' - aNr) >= -tol,1)';
end



