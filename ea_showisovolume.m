function ea_showisovolume(resultfig,elstruct,options)
ht=uitoolbar(gcf);
        set(0,'CurrentFigure',resultfig); 

isom=options.d3.isomatrix;
jetlist=jet;
for side=options.sides
    

    cnt=1;
    for sub=1:length(elstruct)
        for cont=1:options.numcontacts
            if ~isnan(isom{side}(sub,cont));
            X{side}(cnt)=elstruct(sub).coords_mm{side}(cont,1);
            Y{side}(cnt)=elstruct(sub).coords_mm{side}(cont,2);
            Z{side}(cnt)=elstruct(sub).coords_mm{side}(cont,3);
            V{side}(cnt)=isom{side}(sub,cont);
            cnt=cnt+1;
            end
        end
    end
    
    X{side}=X{side}(:);        Y{side}=Y{side}(:);        Z{side}=Z{side}(:); V{side}=V{side}(:);
    assignin('base','X',X);
        assignin('base','Y',Y);
    assignin('base','Z',Z);

    
    bb(1,:)=[min(X{side}),max(X{side})];
    bb(2,:)=[min(Y{side}),max(Y{side})];
    bb(3,:)=[min(Z{side}),max(Z{side})];
    [XI,YI,ZI]=meshgrid(linspace(bb(1,1),bb(1,2),100),linspace(bb(2,1),bb(2,2),100),linspace(bb(3,1),bb(3,2),100));
    F = scatteredInterpolant(X{side},Y{side},Z{side},V{side});
    VI{side}=F(XI,YI,ZI);
    
    if options.d3.isovscloud==2 % show interpolated point mesh
        for xx=1:10:size(VI{side},1)
            for yy=1:10:size(VI{side},2)
                for zz=1:10:size(VI{side},3)
                    usefacecolor=VI{side}(xx,yy,zz)*((64+miniso(options.d3.isomatrix(:)))/(maxiso(options.d3.isomatrix(:))+miniso(options.d3.isomatrix(:))));
                    usefacecolor=ind2rgb(round(usefacecolor),jetlist);
                    
                    isopatch(side,xx,yy,zz)=plot3(XI(xx,yy,zz),YI(xx,yy,zz),ZI(xx,yy,zz),'.','Color',usefacecolor);
                end
            end
        end
        
    elseif options.d3.isovscloud==3 % show isovolume
        
        VI{side}=smooth3(VI{side},'gaussian',3);

        fv{side}=isosurface(XI,YI,ZI,VI{side},mean(VI{side}(:)));
        fv{side}=reducepatch(fv{side},0.5);
        cdat=repmat(60,(length(fv{side}.vertices)),1);
        isopatch(side,1)=patch(fv{side},'CData',cdat,'FaceColor',[0.8 0.8 1.0],'facealpha',0.7,'EdgeColor','none','facelighting','phong');
        jetlist=jet;
        ea_spec_atlas(isopatch(side,1),'isovolume',jet,1);
    
    end
    patchbutton(side)=uitoggletool(ht,'CData',ea_get_icn('isovolume',options),'TooltipString','Isovolume','OnCallback',{@isovisible,isopatch(side,:)},'OffCallback',{@isoinvisible,isopatch(side,:)},'State','on');

end



function isovisible(hobj,ev,atls)
set(atls, 'Visible', 'on');
%disp([atls,'visible clicked']);

function isoinvisible(hobj,ev,atls)
set(atls, 'Visible', 'off');
%disp([atls,'invisible clicked']);


function m=maxiso(cellinp) % simply returns the highest entry of matrices in a cell.
m=-inf;
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

