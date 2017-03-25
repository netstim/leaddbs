function varargout=ea_segmented_cylinder(varargin)
% This function creates a segmented cylinder as e.g. used by the boston
% vercise directed leads.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

h=figure;
hold on
for conins=1:2
    de=figure('visible','off');
    hold on
    switch conins
        case 1 % contacts are plotted
            ix=[15,105
                135,225
                255,345];
            fc='r'; usecolor=0.5;
        case 2 % insulations are plotted
            ix=[105,135
                225,255
                345,360]; % need to add  1,15 here..
            addix=[1,15];
            fc='b'; usecolor=1;
    end
    
    
    for out=[1,2]
        switch out
            case 1 % small circle inside electrode
                r=0.5;
            case 2 % outside electrode
                r=1;
        end
        [X,Y,Z]=cylinder(r,360); % starts at 3 o clock and goes counter-clockwise.
        
        for i=1:size(ix,1)
            
            if conins==2 && out==1 % small cylinder
                smcy=surf(X,Y,Z);
                smcy=surf2patch(smcy,'triangles');
                smcy.vertices=[0,0,0
                    0,0,1
                    smcy.vertices];
                for smci=3:length(smcy.vertices)-2
                    smcy.faces=[smcy.faces;
                        smci,smci+2,1+(~mod(smci,2)*1)];
                end
            end
            
            if i==3 && conins==2 % merge with other bit.
                
                
                
            segX{i}{out}=[X(:,ix(i,1):ix(i,2)),X(:,addix(1,1):addix(1,2))];
            segY{i}{out}=[Y(:,ix(i,1):ix(i,2)),Y(:,addix(1,1):addix(1,2))];
            segZ{i}{out}=[Z(:,ix(i,1):ix(i,2)),Z(:,addix(1,1):addix(1,2))];
                
            else
            
            segX{i}{out}=X(:,ix(i,1):ix(i,2));
            segY{i}{out}=Y(:,ix(i,1):ix(i,2));
            segZ{i}{out}=Z(:,ix(i,1):ix(i,2));
            end
            oute{i}{out}=surf(segX{i}{out},segY{i}{out},segZ{i}{out});
            oute{i}{out}=surf(segX{i}{out},segY{i}{out},segZ{i}{out});
            oute{i}{out}=surf(segX{i}{out},segY{i}{out},segZ{i}{out});
            fv{i}{out}=surf2patch(oute{i}{out},'triangles');
        end
    end
    close(de)
    
    
    for i=1:size(ix,1)
        fvx{conins,i}.faces=[fv{i}{1}.faces;fv{i}{2}.faces+length(fv{i}{1}.vertices)];
        
        fvx{conins,i}.vertices=[fv{i}{1}.vertices;fv{i}{2}.vertices];
        
        % add plates..
        
        % figure
        % patch(fvx{conins,i})
        fvx{conins,i}.faces=[fvx{conins,i}.faces
            2,1,length(fv{i}{1}.vertices)+2
            1,length(fv{i}{1}.vertices)+1,length(fv{i}{1}.vertices)+2];
        fvx{conins,i}.faces=[fvx{conins,i}.faces
            length(fv{i}{1}.vertices),length(fv{i}{1}.vertices)-1,length(fvx{conins,i}.vertices)
            length(fvx{conins,i}.vertices),length(fvx{conins,i}.vertices)-1,length(fv{i}{1}.vertices)-1];
        
        
        
        
        %close lids (top/bottom)..
        for face=1:length(fv{i}{1}.vertices)-2;
            fvx{conins,i}.faces=[fvx{conins,i}.faces;
                face+length(fv{i}{1}.vertices)+2,face+2,face;
                face+length(fv{i}{1}.vertices)+2,face+length(fv{i}{1}.vertices),face];
        end
        
        set(0,'CurrentFigure',h);
        p{conins,i}=patch(fvx{conins,i},'FaceColor',fc,'EdgeColor','none','CDataMapping','direct','CData',repmat([0.5,0.5,0.5],length(fvx{conins,i}.faces),1));
   
        specsurf(p{conins,i},usecolor,1)
        if conins==2 % plot ring of insulation
            fvx{conins,i+1}=smcy;
            p{conins,i+1}=patch(smcy,'FaceColor',fc,'EdgeColor','none','CDataMapping','direct','CData',repmat([0.5,0.5,0.5],length(smcy.faces),1));
                specsurf(p{conins,i+1},usecolor,1)
        end
        
    end
    clear segX segY segZ oute fv fvc
end

a=camlight('headlight');
set(gcf,'Renderer','OpenGL')
%[h]=ea_show_light(h);
%close(h)
varargout{1}=fvx;





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