function ea_showisovolume(resultfig,elstruct,options)
ht=uitoolbar(gcf);
        set(0,'CurrentFigure',resultfig); 

isom=options.d3.isomatrix;
load([options.root,options.patientname,filesep,'LEAD_groupanalysis.mat']);

% check if isomatrix needs to be expanded from single vector by using stimparams:

if ~iscell(isom) % check if isomatrix is a cell ({[right_matrix]},{[left_matrix]}), if not convert to one.
    if min(size(isom))==1 && length(size(isom))==2 % single vector
        
        for side=1:2
            try
                stimmat{side}=cat(1,M.stimparams(:,1).U);
            catch
                warning('Stimulation parameters not set, using each electrode contact from lead.');
                stimmat{side}=ones(length(M.patient.list),4);
            end
            stimmat{side}=bsxfun(@times,stimmat{side}>0,isom);
        end
        
    end
    isom=stimmat;
end


%


jetlist=jet;
if size(isom{1},2)==size(M.elstruct(1).coords_mm{1},1)-1
    shifthalfup=1;
elseif size(isom{1},2)==size(M.elstruct(1).coords_mm{1},1)
    shifthalfup=0;
else
   error('Isomatrix has wrong size. Please specify a correct matrix.')
end
for side=options.sides
    cnt=1;
    for sub=1:length(elstruct)
        for cont=1:size(isom{1},2)
            if ~isnan(isom{side}(sub,cont));
                if ~shifthalfup
                    X{side}(cnt)=elstruct(sub).coords_mm{side}(cont,1);
                    Y{side}(cnt)=elstruct(sub).coords_mm{side}(cont,2);
                    Z{side}(cnt)=elstruct(sub).coords_mm{side}(cont,3);
                else % using pairs of electrode contacts (i.e. 3 pairs if there are 4 contacts)
                    X{side}(cnt)=mean([elstruct(sub).coords_mm{side}(cont,1),elstruct(sub).coords_mm{side}(cont+1,1)]);
                    Y{side}(cnt)=mean([elstruct(sub).coords_mm{side}(cont,2),elstruct(sub).coords_mm{side}(cont+1,2)]);
                    Z{side}(cnt)=mean([elstruct(sub).coords_mm{side}(cont,3),elstruct(sub).coords_mm{side}(cont+1,3)]);
                end
                V{side}(cnt)=isom{side}(sub,cont);
                
                cnt=cnt+1;
            end
        end
    end
    
    X{side}=X{side}(:);        Y{side}=Y{side}(:);        Z{side}=Z{side}(:); V{side}=V{side}(:);
    %assignin('base','X',X);
    %assignin('base','Y',Y);
    %assignin('base','Z',Z);

    
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
    
    if options.d3.exportisovolume % export to nifti volume
        
        Vol=spm_vol([options.earoot,'templates',filesep,'bb.nii']);
        nii{side}=spm_read_vols(Vol);
        nii{side}(:)=nan;
        XYZ=[X{side},Y{side},Z{side},ones(length(X{side}),1)]';
        XYZ=Vol.mat\XYZ; % to voxel space.
        XYZ=round(XYZ(1:3,:)');
        % repeat the above but in voxel space..
        clear bb
        bb(1,:)=[min(XYZ(:,1)),max(XYZ(:,1))];
        bb(2,:)=[min(XYZ(:,2)),max(XYZ(:,2))];
        bb(3,:)=[min(XYZ(:,3)),max(XYZ(:,3))];
        clear XI YI ZI
        [XI,YI,ZI]=meshgrid([bb(1,1):bb(1,2)],[bb(2,1):bb(2,2)],[bb(3,1):bb(3,2)]);
        
        F = scatteredInterpolant(XYZ(:,1),XYZ(:,2),XYZ(:,3),V{side});
        
        clear xix yix zix
        xix=bb(1,1):bb(1,2); yix=bb(2,1):bb(2,2); zix=bb(3,1):bb(3,2);
        
        nii{side}(xix,yix,zix)=F({xix,yix,zix});
        
        
        switch side
            case 1
                lr='right';
            case 2
                lr='left';
        end
        Vol.fname=[options.root,options.patientname,filesep,'iso_volume_',lr,'.nii'];
        Vol.dtype=[32 1];
        spm_write_vol(Vol,nii{side});
        if side==2; % write out combined volume.
           Vol.fname=[options.root,options.patientname,filesep,'iso_volume_combined.nii'];
           nii=nanmean(cat(4,nii{1},nii{2}),4);
            spm_write_vol(Vol,nii);
        end
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

