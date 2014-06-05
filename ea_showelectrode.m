function elrender=ea_showelectrode(resultfig,elstruct,pt,options)
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


for side=options.sides
    trajvector=mean(diff(trajectory{side}));
    trajvector=trajvector/norm(trajvector);
    
    if options.d3.elrendering<3
    if options.d3.prolong_electrode
        startpoint=trajectory{side}(1,:)-(options.d3.prolong_electrode*(coords_mm{side}(1,:)-trajectory{side}(1,:)));

        else
        startpoint=trajectory{side}(1,:);
    end
            set(0,'CurrentFigure',resultfig); 

    
    % draw trajectory
    [elrender{side}(1),elrender{side}(2),elrender{side}(3)]=ea_cylinder(startpoint,coords_mm{side}(elspec.numel,:)-trajvector*(elspec.contact_length/2),elspec.lead_diameter/2,100,repmat(elspec.lead_color,1,3),1,0);
    
    
    if isfield(elstruct,'group')
        usecolor=elstruct.groupcolors(elstruct.group,:);
    else
        usecolor=elspec.lead_color;
    end
    
    
    if options.d3.elrendering==2 % show a transparent electrode.
        aData=0.1;
    elseif options.d3.elrendering==1 % show a solid electrode.
        aData=1;
    end
    
    
       specsurf(elrender{side}(1),usecolor,aData); specsurf(elrender{side}(2),usecolor,aData); specsurf(elrender{side}(3),usecolor,aData);

    cnt=4;
    
    
    % draw trajectory between contacts
    for cntct=1:elspec.numel-1
                set(0,'CurrentFigure',resultfig); 

        [elrender{side}(cnt),elrender{side}(cnt+1),elrender{side}(cnt+2)]=ea_cylinder(coords_mm{side}(cntct,:)-trajvector*(elspec.contact_length/2),coords_mm{side}(cntct+1,:)+trajvector*(elspec.contact_length/2),elspec.lead_diameter/2,100,repmat(elspec.lead_color,1,3),1,0);
        
        specsurf(elrender{side}(cnt),usecolor,aData); specsurf(elrender{side}(cnt+1),usecolor,aData); specsurf(elrender{side}(cnt+2),usecolor,aData);
        cnt=cnt+3;
    end
    
    
    % draw contacts
    for cntct=1:elspec.numel
                set(0,'CurrentFigure',resultfig); 

        
        [elrender{side}(cnt),elrender{side}(cnt+1),elrender{side}(cnt+2)]=ea_cylinder(coords_mm{side}(cntct,:)-trajvector*(elspec.contact_length/2),coords_mm{side}(cntct,:)+trajvector*(elspec.contact_length/2),elspec.contact_diameter/2,100,repmat(elspec.contact_color,1,3),1,0);
        if options.d3.hlactivecontacts && ismember(cntct,elstruct.activecontacts{side}) % make active red contact without transparency
        specsurf(elrender{side}(cnt),usecolor,1); specsurf(elrender{side}(cnt+1),usecolor,1); specsurf(elrender{side}(cnt+2),usecolor,1);
        else
        specsurf(elrender{side}(cnt),elspec.contact_color,aData); specsurf(elrender{side}(cnt+1),elspec.contact_color,aData); specsurf(elrender{side}(cnt+2),elspec.contact_color,aData);
        end
        cnt=cnt+3;
    end
    
    

    
    
    % draw tip
    
    if isfield(elstruct,'group')
        usecolor=elstruct.groupcolors(elstruct.group,:);
    else
        usecolor=elspec.tip_color;
    end
            set(0,'CurrentFigure',resultfig); 

    [cX,cY,cZ] = cylinder((repmat(elspec.tip_diameter/2,1,10)-([10:-1:1].^10/10^10)*(elspec.tip_diameter/2)));
    
    cZ=cZ.*(elspec.tip_length); % scale to fit tip-diameter
    
    % define two points to define cylinder.
    X1=coords_mm{side}(1,:)+trajvector*(elspec.contact_length/2);
    X2=X1+trajvector*elspec.tip_length;
    
    
    cX=cX+X1(1);
    cY=cY+X1(2);
    cZ=cZ-(2*elspec.tip_length)/2+X1(3);
    
    
    
    elrender{side}(cnt)=surf(cX,cY,cZ);
    
    % Calulating the angle between the x direction and the required direction
    % of cylinder through dot product
    angle_X1X2=acos( dot( [0 0 -1],(X2-X1) )/( norm([0 0 -1])*norm(X2-X1)) )*180/pi;
    
    % Finding the axis of rotation (single rotation) to roate the cylinder in
    % X-direction to the required arbitrary direction through cross product
    axis_rot=cross([0 0 -1],(X2-X1) );
    
    
    rotate(elrender{side}(cnt),axis_rot,angle_X1X2,X1)
    
    specsurf(elrender{side}(cnt),usecolor,aData);
    
    else % simply draw pointcloud
        pcnt=1;
                    jetlist=jet;
        % draw contacts
        for cntct=1:elspec.numel
            
            if (options.d3.showactivecontacts && ismember(cntct,elstruct.activecontacts{side})) || (options.d3.showpassivecontacts && ~ismember(cntct,elstruct.activecontacts{side}))
                if options.d3.hlactivecontacts && ismember(cntct,elstruct.activecontacts{side}) % make active red contact without transparency
                    useedgecolor=[0.8,0.5,0.5];
                    ms=10;
                elseif options.d3.hlactivecontacts && ~ismember(cntct,elstruct.activecontacts{side}) % make inactive grey and smaller contact without transparency
                    useedgecolor=[0.5,0.5,0.5];
                    ms=5;
                else
                    useedgecolor=[0.5,0.5,0.5];
                    ms=10;
                end
                % define color
                if options.d3.showisovolume && options.d3.isovscloud==1
                    % draw contacts as colored cloud defined by isomatrix.
                    if ~isnan(options.d3.isomatrix{side}(pt,cntct))
                        
                        usefacecolor=options.d3.isomatrix{side}(pt,cntct)*((64+miniso(options.d3.isomatrix(:)))/(maxiso(options.d3.isomatrix(:))+miniso(options.d3.isomatrix(:))));
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

                    elrender{side}(pcnt)=plot3(coords_mm{side}(cntct,1),coords_mm{side}(cntct,2),coords_mm{side}(cntct,3),'d','MarkerFaceColor',usefacecolor,'MarkerEdgeColor',useedgecolor,'MarkerSize',ms);
                pcnt=pcnt+1;
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