function ea_showisovolume(resultfig,elstruct,options)

% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn
set(0,'CurrentFigure',resultfig);

isobar=getappdata(resultfig,'isobar');
if isempty(isobar)
    isobar=uitoolbar(resultfig);
    setappdata(resultfig,'isobar',isobar);
end
hold on

if isfield(options.d3, 'regressorcolormap')
    cmap = options.d3.regressorcolormap;
else % default blue to red colormap
    cmap = ea_colorgradient(length(gray), [0,0,1], [1,1,1], [1,0,0]);
end

if size(options.d3.isomatrix{1},2)==get_maxNumContacts(elstruct)-1 % number of contact pairs
    shifthalfup=1;
elseif size(options.d3.isomatrix{1},2)==get_maxNumContacts(elstruct) % number of contacts
    shifthalfup=0;
else
    ea_cprintf('CmdWinErrors', 'Be careful! Isomatrix might have wrong size, or numbers of contacts are not consistent across patients.\n');
end

for iside=1:length(options.sides)
    side=options.sides(iside);
    
    cnt=1;
    for sub=1:length(elstruct)
        for cont=1:size(options.d3.isomatrix{1},2)
            if ~isnan(options.d3.isomatrix{side}(sub,cont))
                if ea_arenopoints4side(elstruct(sub).coords_mm,side)
                    %if there are no coordinates, don't parse anything, and skip to next
                    warning_printf=@(str_in) fprintf(['ATTENTION!! : ' str_in '\n']);
                    if side==1
                        warning_printf(['no isovolume will be exported for the right side of subj #' num2str(sub) ' as there is no lead in it.']);
                    elseif side==2
                        warning_printf(['no isovolume will be exported for the right side of subj #' num2str(sub) ' as there is no lead in it.']);
                    else
                        warning_printf(['no isovolume will be exported for side=' num2str(side) ' of subj #' num2str(sub) ' as there is no lead in it.']);
                    end
                else % ~isempty(elstruct(sub).coords_mm{side}) 
                    %if there are coordinates, parse them, otherwise skip to next
                    if ~shifthalfup
                        X{side}(cnt)=elstruct(sub).coords_mm{side}(cont,1);
                        Y{side}(cnt)=elstruct(sub).coords_mm{side}(cont,2);
                        Z{side}(cnt)=elstruct(sub).coords_mm{side}(cont,3);
                    else % using pairs of electrode contacts (i.e. 3 pairs if there are 4 contacts)
                        X{side}(cnt)=mean([elstruct(sub).coords_mm{side}(cont,1),elstruct(sub).coords_mm{side}(cont+1,1)]);
                        Y{side}(cnt)=mean([elstruct(sub).coords_mm{side}(cont,2),elstruct(sub).coords_mm{side}(cont+1,2)]);
                        Z{side}(cnt)=mean([elstruct(sub).coords_mm{side}(cont,3),elstruct(sub).coords_mm{side}(cont+1,3)]);
                    end
                    V{side}(cnt)=options.d3.isomatrix{side}(sub,cont);
                    cnt=cnt+1;
                end
            end
        end
    end

    if cnt==1
        %no elements/leads were present in this side, skip it
        continue
    end
    X{side}=X{side}(:);
    Y{side}=Y{side}(:);
    Z{side}=Z{side}(:);
    V{side}=V{side}(:);
    %assignin('base','X',X);
    %assignin('base','Y',Y);
    %assignin('base','Z',Z);

    bb(1,:)=[min(X{side}),max(X{side})];
    bb(2,:)=[min(Y{side}),max(Y{side})];
    bb(3,:)=[min(Z{side}),max(Z{side})];
    xvec=linspace(bb(1,1),bb(1,2),100);
    yvec=linspace(bb(2,1),bb(2,2),100);
    zvec=linspace(bb(3,1),bb(3,2),100);
    [XI,YI,ZI]=meshgrid(xvec,yvec,zvec);

    F = scatteredInterpolant(X{side},Y{side},Z{side},double(V{side}),'natural');
    F.ExtrapolationMethod='none';
    VI{side}=F(XI,YI,ZI);

    if options.d3.isovscloud==1 % show interpolated point mesh
        ipcnt=1;

        C = VI{side};
        C = (C-min(C(:)))./(max(C(:))-min(C(:))).*(length(cmap)-1);
        C(isnan(C)) = 0;
        C = C+1;

        for xx=1:10:size(VI{side},1)
            for yy=1:10:size(VI{side},2)
                for zz=1:10:size(VI{side},3)
                    if ~isnan(VI{side}(xx,yy,zz))
                        usefacecolor = C(xx,yy,zz);
                        usefacecolor = ind2rgb(round(usefacecolor),cmap);
                        isopatch(side,ipcnt)=plot3(XI(xx,yy,zz),YI(xx,yy,zz),ZI(xx,yy,zz),'o','MarkerFaceColor',usefacecolor,'MarkerEdgeColor',usefacecolor,'Color',usefacecolor);
                        ipcnt=ipcnt+1;
                    end
                end
            end
        end
    elseif options.d3.isovscloud==2 % show isovolume
        VI{side}=smooth3(VI{side},'gaussian',[15 15 15]);

        Vol=VI{side};
        Vol(isnan(Vol))=0;
        fv{side}=isosurface(XI,YI,ZI,Vol,0); % could use thresh instead of 0
        try fv{side}=ea_smoothpatch(fv{side},1,100); end

        C = VI{side};
        % thresh=ea_nanmean(VI{side}(:))-2*ea_nanstd(VI{side}(:));
        % thresh=nanmin(VI{side}(:));
        % C(C<thresh)=nan;
        C = (C-min(C(:)))./(max(C(:))-min(C(:))).*(length(cmap)-1);
        C(isnan(C)) = 0;
        C = C+1;

        nc=isocolors(XI,YI,ZI,C,fv{side}.vertices);
        nc=squeeze(ind2rgb(round(nc),cmap));
        isopatch(side,1)=patch(fv{side},'FaceVertexCData',nc,'FaceColor','interp','facealpha',0.7,'EdgeColor','none','facelighting','phong');

        ea_spec_atlas(isopatch(side,1),'isovolume',cmap,1);

        % export isovolume manually here:
        res=length(Vol);
        chun1=randperm(res); chun2=randperm(res); chun3=randperm(res);
        switch side
            case 1
                lr='right';
            case 2
                lr='left';
        end
        nii=ea_open_vol([ea_space,'bb.nii']);
        nii.mat=mldivide([(chun1);(chun2);(chun3);ones(1,res)]',[xvec(chun1);yvec(chun2);zvec(chun3);ones(1,res)]')';
        nii.img=permute(Vol,[2,1,3]);
        nii.voxsize=ea_detvoxsize(nii.mat);
        nii.dim=size(Vol);
        nii=rmfield(nii,'private');
        nii.fname=[options.root,options.patientname,filesep,options.d3.isomatrix_name,'_',lr,'.nii'];
        ea_write_nii(nii);
        patchbutton(side)=uitoggletool(isobar,'CData',ea_get_icn('isovolume'),'TooltipString',options.d3.isomatrix_name,'OnCallback',{@isovisible,isopatch(side,:)},'OffCallback',{@isoinvisible,isopatch(side,:)},'State','on');
    end
end


function isovisible(hobj,ev,atls)
set(atls, 'Visible', 'on');
%disp([atls,'visible clicked']);


function isoinvisible(hobj,ev,atls)
set(atls, 'Visible', 'off');
%disp([atls,'invisible clicked']);


function maxNumContacts = get_maxNumContacts(elstruct)
coords = {elstruct.coords_mm};
coords = horzcat(coords{:})';
maxNumContacts = max(cellfun(@(x) size(x,1), coords));
