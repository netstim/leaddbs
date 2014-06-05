function [coords,trajvector,fitline]=autocoord_side(patientname,options,side)

%% do auto_coordinate detection for one side (side==1 -> right electrode,
%% side==2 -> left electrode.

tra_nii=load_nii([options.root,patientname,filesep,patientname,'_tra_brain_A3_final.nii']);
cor_nii=load_nii([options.root,patientname,filesep,patientname,'_cor_brain_A3_final_wo_opt.nii']);


imat=zeros([size(tra_nii.img,1),size(tra_nii.img,2),size(tra_nii.img,3),2]);
imat(:,:,:,1)=tra_nii.img;
imat(:,:,:,2)=cor_nii.img;

showdis('Preparing contrasted volume...',options.verbose);

tra_nii.img=gencontrastimage(imat,options.axiscontrast);

fitline=[]; % empty initialization.
for refine=0:options.refinesteps
[fitline,trajvector,diams]=reconstruct_trajectory(fitline,tra_nii,patientname,side,refine,options);
end






%% determine height of last electrode

%detdiams=detrend(diams);


if options.verbose>1; di=figure('name','Finding local maxima in diameters...'); end
if options.verbose>1; plot(diams); end
if options.verbose>2; close(di); end





% %% display cor_nii on fitline
%
% display_cor_fit(trajvector,patientname,options);
%


%% find local maxima in diameters.


%% first, calculate distance between contacts.
% rename matfile to text
movefile([options.root,patientname,filesep,patientname,'_tra_brain_A3_mo_final.mat'],[options.root,patientname,filesep,patientname,'_tra_brain_A3_mo_final.txt']);
% load matrices
tramat=load([options.root,patientname,filesep,patientname,'_tra_brain_A3_mo_final.txt']);
% rename the file to .mat again
movefile([options.root,patientname,filesep,patientname,'_tra_brain_A3_mo_final.txt'],[options.root,patientname,filesep,patientname,'_tra_brain_A3_mo_final.mat']);



dist=calc_distance(options.eldist,trajvector,tramat(1:3,1:3),[options.root,patientname,filesep,patientname,'_tra_brain_A3_final.nii']);
% zdist is the distance between electrodes in z-direction.
zdist=dist/norm(trajvector);



diams(diams==0)=nan;

goodz=sample_cuboid(diams,fitline,trajvector,patientname,zdist,options);








% % [peaks,finalzs]=findpeaks(diams(imgsliceno+2:end)); % detdiams correctly assigns to the slices of the image. (detdiams(x)= slice x).
% % finalzs=finalzs+imgsliceno+1;
% %
% % diffstotarget=abs(diff([finalzs;repmat(33,1,length(finalzs))]));
% % goodz= diffstotarget==min(diffstotarget);



% diams(45:end)=0;
% diams(1:25)=0;
%
%
% [peaks,finalzs]=findpeaks(-diams);
%
%
%
% diffstotarget=abs(diff([finalzs;repmat(35,1,length(finalzs))]));
% goodz=find(diffstotarget==min(diffstotarget));
%
%
% goodz=finalzs(goodz);



%showdis(['Found good peak. Using: ',num2str(goodz(1)),'.'],options.verbose);





% determine coords by goodz, fitline and correction term.

correction=[0,0,0]; %[goodx,goody,0];
try
    coords=findcoords(goodz(1),fitline,trajvector,dist,correction,options);
catch
    showdis('Coords not found.',options.verbose);
    return
end








%
% %% Show figure of transformed image
%
%  The following function does work, but it isn't a nice way to display the images...
% [transnii,transcoords]=rotatenii(tra_nii,Vtra,normtrajvector,coords);

%% Show figure of image
if options.writeoutimages
    cuts=figure;
    wsize=150;
    
    
    
    %% write out axial images
    
    for el=1:4
        %subplot(2,2,el);
        colormap gray
        tra=permute(tra_nii.img,[2,1,3]);
        rcoords=round(coords);
        boundbox{1}=[zminus(rcoords(el,2),wsize):rcoords(el,2)+wsize];
        boundbox{2}=[zminus(rcoords(el,1),wsize):rcoords(el,1)+wsize];
        boundbox{3}=[rcoords(el,3)];
        imagesc(squeeze(flipud(tra(boundbox{1},boundbox{2},boundbox{3}))),[130 450]);
        axis square
        imcoords=coords(el,2:-1:1); % since the image is flipped (permute), flip the coords as well.
        imcoords=permute(((coords-rcoords)+wsize),[2,1]); % since the midpoint is the electrode, subtract wsize.
        imcoords(1)=wsize*2-imcoords(1); % correct for flipud.
        hold on
        
        for ol=1:length(options.overlays)
            if ~strcmp(options.ocolors{ol},'off')
                add_overlay(options.overlays{ol},options.ocolors{ol},options.oopacities(ol),boundbox,cuts,options);
            end
        end
        
        plot(imcoords(1),imcoords(2),'*w','MarkerSize',15);
        hold off
%        maximize(cuts);
        saveas(cuts,[options.root,patientname,filesep,'K',num2str(el*side-1),'_axial.png']);
        
    end
    
    %% write out coronar images
    
    for el=1:4
        %subplot(2,2,el);
        colormap gray
        cor=permute(cor_nii.img,[2,1,3]);
        rcoords=round(coords);
        boundbox{1}=[zminus(rcoords(el,2),wsize):rcoords(el,2)+wsize];
        boundbox{2}=[zminus(rcoords(el,1),wsize):rcoords(el,1)+wsize];
        boundbox{3}=[rcoords(el,3)];
        imagesc(squeeze(flipud(tra(boundbox{1},boundbox{2},boundbox{3}))),[130 450]);
        axis square
        imcoords=coords(el,2:-1:1); % since the image is flipped (permute), flip the coords as well.
        imcoords=permute(((coords-rcoords)+wsize),[2,1]); % since the midpoint is the electrode, subtract wsize.
        imcoords(1)=wsize*2-imcoords(1); % correct for flipud.
        hold on
        
        for ol=1:length(options.overlays)
            if ~strcmp(options.ocolors{ol},'off')
                add_overlay(options.overlays{ol},options.ocolors{ol},options.oopacities(ol),boundbox,cuts,options);
            end
        end
        
        plot(imcoords(1),imcoords(2),'*w','MarkerSize',15);
        hold off
        maximize(cuts);
        saveas(cuts,[options.root,patientname,filesep,'K',num2str(el*side-1),'_axial.png']);
        
    end
    
    
    
    
    % if options.verbose<4; close(cuts); end
end



function res=zminus(A,B)
res=A-B;
if res<0; res=0; end
