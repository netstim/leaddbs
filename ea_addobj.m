
function ea_addobj(src, evt, resultfig, obj, options)

addht = getappdata(resultfig,'addht');
if isempty(addht)
    addht = uitoolbar(resultfig);
    setappdata(resultfig, 'addht', addht);
end

if iscell(obj) % dragndrop for tract and roi, 'obj' is a cell of the files
    if all(cellfun(@numel, regexp(obj, '(\.mat|\.trk)$', 'match', 'once')))
        for i=1:length(obj)
            addfibertract(obj{i}, resultfig, addht, [], 0, options);
        end
    elseif all(cellfun(@numel, regexp(obj, '(\.nii|\.nii\.gz)$', 'match', 'once')))
        pobj.plotFigureH = resultfig;
        pobj.htH = addht;
        for i=1:length(obj)
            % addroi(obj{i}, resultfig, addht, options);
            pobj.color = ea_uisetcolor;
            ea_roi(obj{i}, pobj);
        end
    else
        warndlg('Unsupported file(s) found!');
    end
else  % uigetfile, 'obj' is the type of the files to be selected
    switch obj
        case 'tract'
            % open dialog
            [fina,pana]=uigetfile({'*.mat;*.trk', 'Fiber Files (*.mat,*.trk)'},'Choose Fibertract to add to scene...',[options.root,options.patientname,filesep],'MultiSelect','on');
            if isempty(fina) % user pressed cancel.
                return
            end
            if iscell(fina)
                for fi=1:length(fina)
                    addfibertract([pana,fina{fi}],resultfig,addht,[],0,options);
                end
            else
                if ~fina % user pressed cancel
                    return
                end
                addfibertract([pana,fina],resultfig,addht,[],0,options);
            end
        case 'roi' % atlas
            % open dialog
            [fina,pana]=uigetfile({'*.nii';'*.nii.gz'},'Choose .nii image to add to scene...',[options.root,options.patientname,filesep],'MultiSelect','on');
            if iscell(fina) % multiple files
                for fi=1:length(fina)
                    addroi([pana,fina{fi}],resultfig,addht,options);
                end
            else
                if ~fina % user pressed cancel.
                    return
                end
                addroi([pana,fina],resultfig,addht,options);
            end
        case 'tractmap'
            [tfina,tpana]=uigetfile('*.mat','Choose Fibertract to add to scene...',[options.root,options.patientname,filesep],'MultiSelect','off');
            [rfina,rpana]=uigetfile({'*.nii';'*.nii.gz'},'Choose .nii image to colorcode tracts...',[options.root,options.patientname,filesep],'MultiSelect','off');
            addtractweighted([tpana,tfina],[rpana,rfina],resultfig,addht,options)
    end
end

axis fill


function addtractweighted(tract,weight,resultfig,addht,options)

disp('Loading fibertracts...');
[fibers,idx,voxmm,mat]=ea_loadfibertracts(tract);
disp('Done.');

disp('Loading weight...');
nii=ea_load_nii(weight);
disp('Done.');

nii.img(isnan(nii.img))=0;
nzeros=abs(nii.img)>(std(nii.img(:)));
nzeros=find(nzeros(:));
weights=nii.img(nzeros);
cmweights=weights;
cmweights=cmweights-min(cmweights);
cmweights=cmweights/max(cmweights);
cmweights=(cmweights*(length(gray)-1))+1; % normalize to colormap
cmweights=squeeze(ind2rgb(round(cmweights),jet));

[nzX,nzY,nzZ]=ind2sub(size(nii.img),nzeros);
nzXYZ=[nzX,nzY,nzZ]; clear nzX nzY nzZ

disp('Selecting fibertracts...')
[ix,d]=knnsearch(nzXYZ,fibers(:,1:3));
disp('Done.');
ea_dispercent(0,'Assigning colors to fibers');
ftractcols=zeros(length(idx),3);
fibno=length(idx);
idcnt=1;
for ftract=1:fibno
    thisfibentries=idcnt:idcnt+idx(ftract)-1;
    idcnt=idcnt+idx(ftract);
    ftdists=d(thisfibentries);
    [mindist,mindistix]=min(ftdists);

    smallftdists=ftdists<2;

    if any(smallftdists)
        minidentifier=thisfibentries(smallftdists);
        weights=weights(ix(minidentifier));
        [~,maxix]=max(weights);
        ftractcols(ftract,:)=cmweights(ix(minidentifier(maxix)),:);
    end
    ea_dispercent(ftract/fibno);
end
ea_dispercent(1,'end');

ea_dispercent(0,'Plotting fibers')
cnt=1;
coloredfibs=find(sum(ftractcols,2))';
numcoloredfibs=length(coloredfibs);
for fib=coloredfibs
    ea_dispercent(cnt/numcoloredfibs);
    thisfib=fibers(fibers(:,4)==fib,1:3);
    thisfib=[thisfib,repmat(ftractcols(fib),size(thisfib,1),1)];
    [~,fv(cnt)]=ea_plot_fiber(thisfib',6,0,options);
    cnt=cnt+1;
end
ea_dispercent(1,'end');
fv=ea_concatfv(fv);
addobjr=patch(fv,'Facecolor', 'interp', 'EdgeColor', 'none','FaceAlpha',0.3);

% add toggle button:
[~, tfina] = fileparts(tract);
[~, rfina] = fileparts(weight);
addbutn=uitoggletool(addht,'CData',ea_get_icn('fiber'),'TooltipString',[tfina,' weighted by ',rfina],'OnCallback',{@ea_atlasvisible,addobjr},'OffCallback',{@ea_atlasinvisible,addobjr},'State','on');
%storeinfigure(resultfig,addht,addbutn,addobjr,addobj,fina,'roi',XYZ,0,options); % store rendering in figure.
drawnow

disp('Done.');


function addroi(addobj,resultfig,addht,options)

% set cdata
c = ea_uisetcolor;

if numel(c)==1 && c==0
    return;
end

% load nifti
nii=ea_load_nii(addobj);
nii.img(isnan(nii.img))=0;
nii.img(isinf(nii.img))=0;
% if ~all(abs(nii.voxsize)<=1)
%     ea_reslice_nii(addobj,addobj,[0.5,0.5,0.5],0,[],3);
%     nii=ea_load_nii(addobj);
% end
% nii.img=round(nii.img);

[xx,yy,zz]=ind2sub(size(nii.img),find(nii.img>0)); %(mean(nii.img(nii.img~=0))/3))); % find 3D-points that have correct value.

if ~isempty(xx)
    XYZ=[xx,yy,zz]; % concatenate points to one matrix.
    XYZ=ea_vox2mm(XYZ,nii.mat); % map to mm-space
end

%% new approach with meshgrid based on voxels, then transformed to Nifti coordinatesystem -> should work for all Niftis
for dim=1:3
    gv{dim}=[1:size(nii.img,dim)];
end
[X,Y,Z]=meshgrid(gv{1},gv{2},gv{3});

if options.prefs.hullsmooth
    nii.img = smooth3(nii.img,'gaussian',options.prefs.hullsmooth);
end

fv=isosurface(X,Y,Z,permute(nii.img,[2,1,3]),max(nii.img(:))/2);
fvc=isocaps(X,Y,Z,permute(nii.img,[2,1,3]),max(nii.img(:))/2);
fv.faces=[fv.faces;fvc.faces+size(fv.vertices,1)];
fv.vertices=[fv.vertices;fvc.vertices];

tmp_vertices = nii.mat * [fv.vertices,ones(size(fv.vertices,1),1)]';
fv.vertices = tmp_vertices(1:3,:)';
%%

if ischar(options.prefs.hullsimplify)
    % get to 700 faces
    simplify=700/length(fv.faces);
    fv=reducepatch(fv,simplify);
else
    if options.prefs.hullsimplify<1 && options.prefs.hullsimplify>0
        fv=reducepatch(fv,options.prefs.hullsimplify);
    elseif options.prefs.hullsimplify>1
        simplify=options.prefs.hullsimplify/length(fv.faces);
        fv=reducepatch(fv,simplify);
    end
end

try
    fv=ea_smoothpatch(fv,1,35);
catch
    try
        cd([ea_getearoot,'ext_libs',filesep,'smoothpatch']);
        mex ea_smoothpatch_curvature_double.c -v
        mex ea_smoothpatch_inversedistance_double.c -v
        mex ea_vertex_neighbours_double.c -v
        fv=ea_smoothpatch(fv);
    catch
        warndlg('Patch could not be smoothed. Please supply a compatible Matlab compiler to smooth VTAs.');
    end
end

%?atlasc=59; %rand*64;
jetlist=jet;

co=ones(1,1,3);
co(1,1,:)=c;
atlasc=double(rgb2ind(co,jetlist));

cdat=abs(repmat(atlasc,length(fv.vertices),1) ... % C-Data for surface
    +randn(length(fv.vertices),1)*2)';

% show atlas.
set(0,'CurrentFigure',resultfig);
addobjr=patch(fv,'CData',cdat,'FaceColor',c,'facealpha',0.7,'EdgeColor','none','facelighting','phong');
%ea_spec_atlas(addobjr,'',jetlist,1);

% add toggle button:
[~, fina] = fileparts(addobj);
addbutn=uitoggletool(addht,'CData',ea_get_icn('atlas',c),'TooltipString',fina,'OnCallback',{@ea_atlasvisible,addobjr},'OffCallback',{@ea_atlasinvisible,addobjr},'State','on');
storeinfigure(resultfig,addht,addbutn,addobjr,addobj,fina,'roi',XYZ,0,options); % store rendering in figure.
drawnow


function addfibertract(addobj,resultfig,addht,connect,ft,options)
if ischar(addobj) % filename is given ? load fibertracts.
    if strfind(addobj,'.mat')
        load(addobj);
        if exist('fibsin', 'var')
            fibers = fibsin;
            clear fibsin
        end
        if exist('fibers', 'var')
            if size(fibers,1) < size(fibers,2)
                fibers = fibers';
            end
            if size(fibers,2) == 4
                thisset = fibers(:,1:3);
                [~,~,idx] = unique(fibers(:,4));
                fibidx = accumarray(idx,1);
            elseif size(fibers,2) == 3
                thisset = fibers;
                fibidx = idx;
            else
                error('Wrong input fiber tracts format!');
            end
            clear fibers idx
        else
            error('No fiber tracts found!');
        end
    elseif strfind(addobj,'.trk')
        fileOut = [addobj(1:end-3) 'mat'];
        disp('Converting .trk to ftr.')
        [thisset,fibidx] = ea_trk2ftr(addobj);
        thisset = thisset';
    else
        error('File is neither a .mat nor .trk!')
    end
else % fibers are already loaded and stored in figure.
    thisset=addobj.fibs;
    fibidx=addobj.idx;
end

fib_copy.fibs=thisset; % backup of whole original fiberset will be stored in figure.
fib_copy.idx=fibidx;

if ~isempty(connect) % select fibers based on connecting roi info (i.e. delete all other fibers).
    for roi=1:length(connect.rois) % check connectivities
        in=inhull(thisset,connect.xyz{roi})';
        idxv = repelem(1:numel(fibidx), fibidx)';
        selectedfibs{roi}=unique(idxv(in));
    end
    selectedfibs=unique(cell2mat(selectedfibs(:)));
    thisset=mat2cell(thisset,fibidx,3)';
    thisset=thisset(selectedfibs); % choose selected fibers.
end

%% new visualization part
c = ea_uisetcolor;
if c == 0
    c=NaN;
end
addobjr=ea_showfiber(thisset,fibidx,c);

axis fill

[~, fina] = fileparts(addobj);
addbutn=uitoggletool(addht,'CData',ea_get_icn('fibers'),'TooltipString',fina,'OnCallback',{@ea_atlasvisible,addobjr},'OffCallback',{@ea_atlasinvisible,addobjr},'State','on');
storeinfigure(resultfig,addht,addbutn,addobjr,addobj,fina,'tract',fib_copy,ft,options); % store rendering in figure.


function storeinfigure(resultfig,addht,addbutn,obj,path,name,type,data,replace,options)
AL=getappdata(resultfig,'AL');
% AL has fields of
% ROI
% ROINAMES
% FTS
% FTSNAMES
% MENU
% GUI
% store info..

if isempty(AL) % initialize AL
    AL.FTS=cell(0);
    AL.ROI=cell(0);
    AL.FTSNAMES=cell(0);
    AL.FTSFILES=cell(0);
    AL.FTSDATA=cell(0);
    AL.FTSBUTN=cell(0);
    AL.ROINAMES=cell(0);
    AL.ROIFILES=cell(0);
    AL.ROIDATA=cell(0);
    AL.MENU=struct;
    AL.GUI=struct;
end

switch type
    case 'tract'
        if replace % fibertract has been there before and has now been selected to be plotted connecting to a roi.
            AL.FTS{replace}=obj;
            AL.FTSNAMES{replace}=name;
            AL.FTSFILES{replace}=path;
            AL.FTSBUTN{replace}=addbutn;
            %         for roi=1:length(AL.ROI) % add connectivity data
            %             AL.GUI.FTS(length(AL.FTS)).ROI(roi)=0;
            %         end
            AL.FTSDATA{replace}=data;
        else
            AL.FTS{end+1}=obj;
            AL.FTSNAMES{end+1}=name;
            AL.FTSFILES{end+1}=path;
            AL.FTSBUTN{end+1}=addbutn;
            for roi=1:length(AL.ROI) % add connectivity data
                AL.GUI.FTS(length(AL.FTS)).ROI(roi)=0;
            end
            AL.FTSDATA{end+1}=data;

        end
    case 'roi'
        AL.ROI{end+1}=obj;
        AL.ROINAMES{end+1}=name;
        AL.ROIFILES{end+1}=path;
        AL.ROIDATA{end+1}=data;
        for ft=1:length(AL.FTS) % add connectivity data
            AL.GUI.FTS(ft).ROI(end+1)=0;
        end
end

% build menu for toggling
try
    delete(AL.MENU.MAINMENU)
end

if ~isempty(AL.FTS) % only build fibertracking menu if there is at least one fiberset.
    AL.MENU.MAINMENU=uimenu(resultfig,'Label','Fibertracking');
    for ft=1:length(AL.FTS)
        AL.MENU.FTMENU(ft) = uimenu(AL.MENU.MAINMENU,'Label',AL.FTSNAMES{ft});
        for roi=1:length(AL.ROI)
            AL.MENU.ROIMENU(ft,roi)=uimenu(AL.MENU.FTMENU(ft),'Label',AL.ROINAMES{roi},'Callback',{@dotracking,ft,roi,resultfig,addht,options});

            set(AL.MENU.ROIMENU(ft,roi),'Checked',binary2onoff(AL.GUI.FTS(ft).ROI(roi))) % set checks on menu.

        end
    end
end

axis fill
% store in figure.

setappdata(resultfig,'AL',AL);


function dotracking(hobj,ev,ft,roi,resultfig,addht,options)

AL=getappdata(resultfig,'AL'); %structure containing all the ft and roi handles.

set(hobj,'Checked',binary2onoff(~AL.GUI.FTS(ft).ROI(roi))) % check if this has been attached and uncheck if unattached.
AL.GUI.FTS(ft).ROI(roi)=~AL.GUI.FTS(ft).ROI(roi); % flip GUI info.

% delete and update fibertract.
delete(AL.FTS{ft}(:))
delete(AL.FTSBUTN{ft}); % delete togglebutton

%% define connect struct here:

connect.rois=find(AL.GUI.FTS(ft).ROI);
connect.xyz=AL.ROIDATA(connect.rois);

setappdata(resultfig,'AL',AL);

if isempty(connect.rois) % don't even send connect structure since this will save time.
    addfibertract(AL.FTSDATA{ft},resultfig,addht,AL.FTSNAMES{ft},[],ft,options)
else
    addfibertract(AL.FTSDATA{ft},resultfig,addht,AL.FTSNAMES{ft},connect,ft,options)
end


function str=binary2onoff(bin)
str='off';
if bin
    str='on';
end


function indcol=detcolor(mat) % determine color based on traversing direction.

xyz=abs(diff(mat,1,2));
rgb=xyz/max(xyz(:));

rgb=[rgb,rgb(:,end)];
rgbim=zeros(1,size(rgb,2),3);
rgbim(1,:,:)=rgb';
indcol=double(rgb2ind(rgbim,jet));


function  dispercent(varargin)
%
percent=round(varargin{1}*100);

if nargin==2
    if strcmp(varargin{2},'end')
        fprintf('\n')
        fprintf('\n')

        fprintf('\n')

    else
        fprintf(1,[varargin{2},':     ']);


    end
else
    fprintf(1,[repmat('\b',1,(length(num2str(percent))+1)),'%d','%%'],percent);
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
    if size(aNr,2) ~= length(j),
        aNr = repmat(aN,1,length(j));
    end
    in(j) = all((nrmls*testpts(j,:)' - aNr) >= -tol,1)';
end
