function ML=cv_build_graph_fig(resultfig,directory,selectedparc,handles,options)

% read in atlas structures:
aID = fopen([ea_space(options,'labeling'),selectedparc,'.txt']);
atlas_lgnd=textscan(aID,'%d %s');


D=length(atlas_lgnd{2});

parc=ea_load_nii([ea_space(options,'labeling'),selectedparc,'.nii']);

% build centroid mesh
parc.img=round(parc.img);
centr=zeros(D,3);
for reg=1:D
   [xx,yy,zz]=ind2sub(size(parc.img),find(parc.img==atlas_lgnd{1}(reg)));
    XYZ=[xx,yy,zz,ones(length(zz),1)]';

    XYZ=parc.mat*XYZ;
    centr(reg,:)=ea_nanmean(XYZ(1:3,:),2)';
end

% load CM into tCM variable
selcm=get(handles.wmatmodality,'String');
selcm=selcm{get(handles.wmatmodality,'Value')};
if strcmp(selcm,'Choose...')
    cmpth=get(handles.chosenwmatstr,'String');
else
    cmpth=[directory,'connectomics',filesep,selectedparc,filesep,selcm,'.mat'];
end
try
    CM=load(cmpth);
catch
    ea_error('Chosen connectivity matrix does not exist. File has vanished.');
end

fs=fieldnames(CM); % use first variable of file.
tCM=CM.(fs{1});
clear CM;
try
    tCM(logical(eye(D)))=nan;
catch
    ea_error('Chosen connectivity matrix does not match size of selected parcellation atlas. Please choose a different .mat file with only one variable matching the size of selected parcellation.');
end

% check if matrix has negative entries
isbipolar=any(tCM(:)<0);

ntCM=tCM; ntCM(isnan(ntCM))=0;
isdirected=~isequal(ntCM',ntCM);
clear ntCM

% normalize matrix entries:
ntCM=tCM;
ntCM(~isnan(tCM))=ea_normal(ntCM(~isnan(tCM)));


% threshold matrix
if strcmp(get(handles.wmatthresh,'String'),'auto')
    thresh=ea_nanmean(tCM(:))+3*ea_nanstd(tCM(:));
else
    thresh=str2double(get(handles.wmatthresh,'String'));
end
set(handles.wmatthreshis,'String',num2str(thresh));

try
tCM(abs(tCM)<thresh)=nan;
catch
    keyboard
end

if get(handles.wmatedges,'Value')
    % draw edges:
    stemfactor=0.5;
    if isdirected
        tipfactor=0.7;
    else
        tipfactor=0;
    end

    pcnt=1;
    ncnt=1;
    for edge=find(~isnan(tCM))'

        [ii,jj]=ind2sub(size(tCM),edge);
        if ii~=jj
            if isbipolar
                if tCM(edge)<0
                    cstr='r';
                    nfv(ncnt)=mArrow3(centr(ii,:),centr(jj,:),'color',cstr,'stemWidth',abs(ntCM(edge))*stemfactor,'tipWidth',abs(ntCM(edge))*tipfactor);
                    ncnt=ncnt+1;
                else
                    cstr='b';
                    pfv(pcnt)=mArrow3(centr(ii,:),centr(jj,:),'color',cstr,'stemWidth',ntCM(edge)*stemfactor,'tipWidth',ntCM(edge)*tipfactor);
                    pcnt=pcnt+1;
                end
            else
                cstr='b';

                pfv(pcnt)=mArrow3(centr(ii,:),centr(jj,:),'color',cstr,'stemWidth',ntCM(edge)*stemfactor,'tipWidth',ntCM(edge)*tipfactor);
                pcnt=pcnt+1;
            end
        end
    end

    if exist('pfv','var')
        pfv=ea_concatfv(pfv);
        set(0,'CurrentFigure',resultfig);
        hp=patch(pfv);

        set(hp,'Facecolor',[0.8,0.2,0.1]);
        set(hp,'EdgeColor','none');
        set(hp,'FaceAlpha',0.5);
    end
    if isbipolar
        if exist('nfv','var')
            nfv=ea_concatfv(nfv);
            set(0,'CurrentFigure',resultfig);
            hn=patch(nfv);
            set(hn,'Facecolor',[0.1,0.2,0.8]);
            set(hn,'EdgeColor','none');
            set(hn,'FaceAlpha',0.5);
        end
    end

end
if get(handles.wmatnodes,'Value')
    % show all areas which are connected above threshold


    tparc=parc;
    tnodes=ea_nansum(tCM,1);

    tnodeix=find(~isnan(tnodes));
    ptnodeix=tnodeix(tnodes(tnodeix)>0);
        ntnodeix=tnodeix(tnodes(tnodeix)<0);
    tparc.img(:)=0;
    for node=tnodeix
        tparc.img(parc.img==node)=tnodes(node);
    end


    pgraphsurf=ea_showconnectivitypatch(resultfig,parc,tparc.img,0,atlas_lgnd{2},ptnodeix,'autumn',get(handles.wmatnodes,'Value'),get(handles.wmatlabs,'Value'));

    ngraphsurf=ea_showconnectivitypatch(resultfig,parc,tparc.img*-1,0,atlas_lgnd{2},ntnodeix,'winter',get(handles.wmatnodes,'Value'),get(handles.wmatlabs,'Value'));

end
% export structures as ML
ML=struct;
try ML.pedges=hp; end
try ML.pnodes=pgraphsurf; end
if isbipolar
try ML.hn=hn; end
try ML.nnodes=ngraphsurf; end
end


function [fv] = mArrow3(p1,p2,varargin)
%mArrow3 - plot a 3D arrow as patch object (cylinder+cone)
%
% syntax:   h = mArrow3(p1,p2)
%           h = mArrow3(p1,p2,'propertyName',propertyValue,...)
%
% with:     p1:         starting point
%           p2:         end point
%           properties: 'color':      color according to MATLAB specification
%                                     (see MATLAB help item 'ColorSpec')
%                       'stemWidth':  width of the line
%                       'tipWidth':   width of the cone
%
%           Additionally, you can specify any patch object properties. (For
%           example, you can make the arrow semitransparent by using
%           'facealpha'.)
%
% example1: h = mArrow3([0 0 0],[1 1 1])
%           (Draws an arrow from [0 0 0] to [1 1 1] with default properties.)
%
% example2: h = mArrow3([0 0 0],[1 1 1],'color','red','stemWidth',0.02,'facealpha',0.5)
%           (Draws a red semitransparent arrow with a stem width of 0.02 units.)
%
% hint:     use light to achieve 3D impression
%

propertyNames = {'edgeColor'};
propertyValues = {'none'};

%% evaluate property specifications
for argno = 1:2:nargin-2
    switch varargin{argno}
        case 'color'
            propertyNames = {propertyNames{:},'facecolor'};
            propertyValues = {propertyValues{:},varargin{argno+1}};
        case 'stemWidth'
            if isreal(varargin{argno+1})
                stemWidth = varargin{argno+1};
            else
                warning('mArrow3:stemWidth','stemWidth must be a real number');
            end
        case 'tipWidth'
            if isreal(varargin{argno+1})
                tipWidth = varargin{argno+1};
            else
                warning('mArrow3:tipWidth','tipWidth must be a real number');
            end
        otherwise
            propertyNames = {propertyNames{:},varargin{argno}};
            propertyValues = {propertyValues{:},varargin{argno+1}};
    end
end

%% default parameters
if ~exist('stemWidth','var')
    ax = axis;
    if numel(ax)==4
        stemWidth = norm(ax([2 4])-ax([1 3]))/300;
    elseif numel(ax)==6
        stemWidth = norm(ax([2 4 6])-ax([1 3 5]))/300;
    end
end
if ~exist('tipWidth','var')
    tipWidth = 3*stemWidth;
end
tipAngle = 22.5/180*pi;
tipLength = tipWidth/tan(tipAngle/2);
ppsc = 50;  % (points per small circle)
ppbc = 250; % (points per big circle)

%% ensure column vectors
p1 = p1(:);
p2 = p2(:);

%% basic lengths and vectors
x = (p2-p1)/norm(p2-p1); % (unit vector in arrow direction)
y = cross(x,[0;0;1]);    % (y and z are unit vectors orthogonal to arrow)
if norm(y)<0.1
    y = cross(x,[0;1;0]);
end
y = y/norm(y);
z = cross(x,y);
z = z/norm(z);

%% basic angles
theta = 0:2*pi/ppsc:2*pi; % (list of angles from 0 to 2*pi for small circle)
sintheta = sin(theta);
costheta = cos(theta);
upsilon = 0:2*pi/ppbc:2*pi; % (list of angles from 0 to 2*pi for big circle)
sinupsilon = sin(upsilon);
cosupsilon = cos(upsilon);

%% initialize face matrix
f = NaN([ppsc+ppbc+2 ppbc+1]);

%% normal arrow
if norm(p2-p1)>tipLength
    % vertices of the first stem circle
    for idx = 1:ppsc+1
        v(idx,:) = p1 + stemWidth*(sintheta(idx)*y + costheta(idx)*z);
    end
    % vertices of the second stem circle
    p3 = p2-tipLength*x;
    for idx = 1:ppsc+1
        v(ppsc+1+idx,:) = p3 + stemWidth*(sintheta(idx)*y + costheta(idx)*z);
    end
    % vertices of the tip circle
    for idx = 1:ppbc+1
        v(2*ppsc+2+idx,:) = p3 + tipWidth*(sinupsilon(idx)*y + cosupsilon(idx)*z);
    end
    % vertex of the tiptip
    v(2*ppsc+ppbc+4,:) = p2;

    % face of the stem circle
    f(1,1:ppsc+1) = 1:ppsc+1;
    % faces of the stem cylinder
    for idx = 1:ppsc
        f(1+idx,1:4) = [idx idx+1 ppsc+1+idx+1 ppsc+1+idx];
    end
    % face of the tip circle
    f(ppsc+2,:) = 2*ppsc+3:(2*ppsc+3)+ppbc;
    % faces of the tip cone
    for idx = 1:ppbc
        f(ppsc+2+idx,1:3) = [2*ppsc+2+idx 2*ppsc+2+idx+1 2*ppsc+ppbc+4];
    end

%% only cone v
else
    tipWidth = 2*sin(tipAngle/2)*norm(p2-p1);
    % vertices of the tip circle
    for idx = 1:ppbc+1
        v(idx,:) = p1 + tipWidth*(sinupsilon(idx)*y + cosupsilon(idx)*z);
    end
    % vertex of the tiptip
    v(ppbc+2,:) = p2;
    % face of the tip circle
    f(1,:) = 1:ppbc+1;
    % faces of the tip cone
    for idx = 1:ppbc
        f(1+idx,1:3) = [idx idx+1 ppbc+2];
    end
end

%% draw
fv.faces = f;
fv.vertices = v;
%h = patch(fv);
% for propno = 1:numel(propertyNames)
%     try
%         set(h,propertyNames{propno},propertyValues{propno});
%     catch
%         disp(lasterr)
%     end
% end
