function [elrender,ellabel]=ea_showcorticalstrip(resultfig,elstruct,pt,options)
% This function renders the electrode as defined by options.elspec and
% coords_mm.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

coords_mm=elstruct.coords_mm;
trajectory=elstruct.trajectory;


if ~isfield(elstruct,'elmodel') % usually, elspec is defined by the GUI. In case of group analyses, for each patient, a different electrode model can be selected for rendering.

    elspec=options.elspec;
else % if elspec is defined for each electrode, overwrite options-struct settings here.
    o=ea_resolve_elspec(elstruct);
    elspec=o.elspec; clear o

end

if ~isfield(elstruct,'activecontacts')
    elstruct.activecontacts{1}=zeros(elspec.numel,1);
    elstruct.activecontacts{2}=zeros(elspec.numel,1);
end
        jetlist=jet;
     %   jetlist=jet;


for side=1:length(options.sides)
    trajvector=mean(diff(trajectory{side}));

    trajvector=trajvector/norm(trajvector);
try
    startpoint=trajectory{side}(1,:)-(2*(coords_mm{side}(1,:)-trajectory{side}(1,:)));
catch

end
    if options.d3.elrendering<3

        set(0,'CurrentFigure',resultfig);

        % draw patientname
        lstartpoint=startpoint-(0.03*(coords_mm{side}(1,:)-startpoint));
        ellabel(side)=text(lstartpoint(1),lstartpoint(2),lstartpoint(3),elstruct.name);


        % draw trajectory


        cnt=1;
        load([ea_getearoot,'templates',filesep,'electrode_models',filesep,elspec.matfname])
        A=[electrode.head_position,1;
            electrode.tail_position,1
            electrode.x_position,1
            electrode.y_position,1]; % points in model
        redomarkers=0;
        if ~isfield(elstruct,'markers') % backward compatibility to old electrode format
            redomarkers=1;

        else
            if isempty(elstruct.markers)

                redomarkers=1;
            end
        end
        if redomarkers
            for iside=options.sides
                elstruct.markers(iside).head=elstruct.coords_mm{iside}(1,:);
                elstruct.markers(iside).tail=elstruct.coords_mm{iside}(4,:);

                [xunitv, yunitv] = ea_calcxy(elstruct.markers(iside).head, elstruct.markers(iside).tail, 'null');
                elstruct.markers(iside).x = elstruct.coords_mm{iside}(1,:) + xunitv*(options.elspec.lead_diameter/2);
                elstruct.markers(iside).y = elstruct.coords_mm{iside}(1,:) + yunitv*(options.elspec.lead_diameter/2);
            end
        end

        B=[elstruct.markers(side).head,1;
           elstruct.markers(side).tail,1;
           elstruct.markers(side).x,1;
           elstruct.markers(side).y,1];
        setappdata(resultfig,'elstruct',elstruct);
        setappdata(resultfig,'elspec',elspec);

        X = mldivide(A,B); X=X';

        if options.d3.elrendering==2 % show a transparent electrode.
            aData=0.1;
        elseif options.d3.elrendering==1 % show a solid electrode.
            aData=1;
        end



        for ins=1:length(electrode.insulation)
            electrode.insulation(ins).vertices=X*[electrode.insulation(ins).vertices,ones(size(electrode.insulation(ins).vertices,1),1)]';
            electrode.insulation(ins).vertices=electrode.insulation(ins).vertices(1:3,:)';
            elrender{side}(cnt)=patch(electrode.insulation(ins));

        if isfield(elstruct,'group')
            usecolor=elstruct.groupcolors(elstruct.group,:);
        else
            usecolor=elspec.lead_color;
        end
        specsurf(elrender{side}(cnt),usecolor,aData);
            cnt=cnt+1;
        end
        for con=1:length(electrode.contacts)
            electrode.contacts(con).vertices=X*[electrode.contacts(con).vertices,ones(size(electrode.contacts(con).vertices,1),1)]';
            electrode.contacts(con).vertices=electrode.contacts(con).vertices(1:3,:)';
            elrender{side}(cnt)=patch(electrode.contacts(con));

            if options.d3.hlactivecontacts && ismember(con,find(elstruct.activecontacts{side})) % make active red contact without transparency
                specsurf(elrender{side}(cnt),[0.8,0.2,0.2],1);
            else
                specsurf(elrender{side}(cnt),elspec.contact_color,aData);
            end

            cnt=cnt+1;
        end

    else % simply draw pointcloud

        shifthalfup=0;
        % check if isomatrix needs to be expanded from single vector by using stimparams:
        try


        if size(options.d3.isomatrix{1},2)==elspec.numel-1 % 3 contact pairs
            shifthalfup=1;
        elseif size(options.d3.isomatrix{1},2)==elspec.numel % 4 contacts
            shifthalfup=0;
        else

            ea_error('Isomatrix has wrong size. Please specify a correct matrix.')
        end
        end



        if options.d3.prolong_electrode

            startpoint=trajectory{side}(1,:)-(options.d3.prolong_electrode*(coords_mm{side}(1,:)-trajectory{side}(1,:)));

        else
            startpoint=trajectory{side}(1,:);
        end
        set(0,'CurrentFigure',resultfig);

        % draw patientname
        lstartpoint=startpoint-(0.03*(coords_mm{side}(1,:)-startpoint));



        %ellabel(side)=text(lstartpoint(1),lstartpoint(2),lstartpoint(3),elstruct.name);
        ellabel=nan;
        pcnt=1;
        % draw contacts

     try
        minval=abs(min(options.d3.isomatrix{side}(:)));
        maxval=max(options.d3.isomatrix{side}(:));
     end
        for cntct=1:elspec.numel-shifthalfup

            if (options.d3.showactivecontacts && ismember(cntct,find(elstruct.activecontacts{side}))) || (options.d3.showpassivecontacts && ~ismember(cntct,elstruct.activecontacts{side}))
                if options.d3.hlactivecontacts && ismember(cntct,find(elstruct.activecontacts{side})) % make active red contact without transparency
                    useedgecolor=[0.8,0.5,0.5];
                    ms=10;
                elseif options.d3.hlactivecontacts && ~ismember(cntct,find(elstruct.activecontacts{side})) % make inactive grey and smaller contact without transparency
                    useedgecolor='none';
                    ms=5;
                else
                    useedgecolor='none';
                    ms=10;
                end
                % define color
                if options.d3.colorpointcloud
                    % draw contacts as colored cloud defined by isomatrix.

                    if ~isnan(options.d3.isomatrix{side}(pt,cntct))

                        usefacecolor=((options.d3.isomatrix{side}(pt,cntct)+minval)/(maxval+minval))*64;

%                         % ## add some contrast (remove these lines for linear
%                         % mapping)
%
%
%                         usefacecolor=usefacecolor-20;
%                         usefacecolor(usefacecolor<1)=1;
%                         usefacecolor=usefacecolor*2;
%                         usefacecolor(usefacecolor>64)=64;
%
%                         % ##

                        usefacecolor=ind2rgb(round(usefacecolor),jetlist);
                    else
                        usefacecolor=nan; % won't draw the point then.
                    end
                else
                    if isfield(elstruct,'group')
                        usefacecolor=elstruct.groupcolors(elstruct.group,:);
                    else
                        usefacecolor=elspec.contact_color;
                    end
                end

                if ~isnan(usefacecolor)
                    set(0,'CurrentFigure',resultfig);
                    if ~shifthalfup
                    elrender{side}(pcnt)=plot3(coords_mm{side}(cntct,1),coords_mm{side}(cntct,2),coords_mm{side}(cntct,3),'o','MarkerFaceColor',usefacecolor,'MarkerEdgeColor',useedgecolor,'MarkerSize',ms);
                    pcnt=pcnt+1;
                    else

                      elrender{side}(pcnt)=plot3(mean([coords_mm{side}(cntct,1),coords_mm{side}(cntct+1,1)]),...
                          mean([coords_mm{side}(cntct,2),coords_mm{side}(cntct+1,2)]),...
                      mean([coords_mm{side}(cntct,3),coords_mm{side}(cntct+1,3)]),...
                      'o','MarkerFaceColor',usefacecolor,'MarkerEdgeColor',useedgecolor,'MarkerSize',ms);
                    drawnow
                    pcnt=pcnt+1;
                    end
                else

                end
                hold on
            end
        end


    end

end


if ~exist('elrender','var')
    elrender=nan;
end


function m=maxiso(cellinp) % simply returns the highest entry of matrices in a cell.
m=0;
for c=1:length(cellinp)
    nm=max(cellinp{c}(:));
    if nm>m; m=nm; end
end

function m=miniso(cellinp)
m=inf;
for c=1:length(cellinp)
    nm=min(cellinp{c}(:));
    if nm<m; m=nm; end
end



function specsurf(varargin)

surfc=varargin{1};
color=varargin{2};
if nargin==3
    aData=varargin{3};
end

len=get(surfc,'ZData');

cd=zeros([size(len),3]);
cd(:,:,1)=color(1);
try % works if color is denoted as 1x3 array
    cd(:,:,2)=color(2);cd(:,:,3)=color(3);
catch % if color is denoted as gray value (1x1) only
    cd(:,:,2)=color(1);cd(:,:,3)=color(1);
end


cd=cd+0.01*randn(size(cd));

set(surfc,'FaceColor','interp');
set(surfc,'CData',cd);

try % for patches
    vertices=get(surfc,'Vertices');
    cd=zeros(size(vertices));
    cd(:)=color(1);
    cd=repmat(color,size(cd,1),size(cd,2)/size(color,2));
    set(surfc,'FaceVertexCData',cd);
end
set(surfc,'AlphaDataMapping','none');

set(surfc,'FaceLighting','phong');
set(surfc,'SpecularColorReflectance',0);
set(surfc,'SpecularExponent',10);
set(surfc,'EdgeColor','none')

if nargin==3
    set(surfc,'FaceAlpha',aData);
end

function C=rgb(C) % returns rgb values for the colors.

C = rem(floor((strfind('kbgcrmyw', C) - 1) * [0.25 0.5 1]), 2);
