function cuts=ea_writeplanes(options)

% This function exports slice views of all electrode contacts reconstructed
% priorly. Images are written as .png image files. Bot transversal and
% coronar views are being exported. Additionally, overlays from atlas-data
% can be visualized via the function ea_add_overlay which uses all atlas
% files that are found in the eAuto_root/atlases directory.

% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn


% load prior results
try
    load([options.root,options.patientname,filesep,'ea_reconstruction']);
catch
    coords_mm=ea_read_fiducials([options.root,options.patientname,filesep,'ea_coords.fcsv'],options);
end

interpfactor=2;


cuts=figure('name',[options.patientname,': 2D cut views'],'numbertitle','off');
axis off
set(gcf,'color','w');
tracorpresent=zeros(3,1); % check if files are present.
switch options.modality
    case 1 % MR
        try
            Vtra=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.gtranii));
            tracorpresent(1)=1;
        catch
            try
                Vtra=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tranii));
                tracorpresent(1)=1;
            end
            
        end
        try
            Vcor=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.gcornii));
            tracorpresent(2)=1;
            
        catch
            try
                Vcor=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.cornii));
                tracorpresent(1)=1;
            end
        end
        try
            Vsag=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.gcornii));
            tracorpresent(3)=1;
            
        catch
            try
                Vsag=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.cornii));
                tracorpresent(3)=1;
            end
        end
    case 2 % CT
        Vtra=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tranii));
        Vcor=spm_vol(fullfile(options.root,options.prefs.patientdir,options.prefs.tranii));
end

for side=1:length(coords_mm)

coords{side}=Vtra.mat\[coords_mm{side},ones(size(coords_mm{side},1),1)]';
coords{side}=coords{side}(1:3,:)';
end
%XYZ_src_vx = src.mat \ XYZ_mm;


for side=options.sides
    %% write out axial images
    for tracor=find(tracorpresent)'
        
        for elcnt=1:options.elspec.numel
            el=elcnt+options.elspec.numel*(side-1);
            %subplot(2,2,el);
            
            % Show MR-volume
            colormap gray
            switch tracor
                
                case 1 % transversal images
                   
                    %title(['Electrode ',num2str(el-1),', transversal view.']);
                    
                    [slice,boundbox]=ea_sample_slice(Vtra,'tra',options.d2.bbsize,coords,el);
                    try
                    imagesc(slice,[nanmean(slice(slice>0))-3*nanstd(slice(slice>0)) nanmean(slice(slice>0))+3*nanstd(slice(slice>0))]);
                    catch
                        imagesc(slice);
                    end
                        axis square
                    hold on
                case 2 % coronar images
                   
                  [slice,boundbox]=ea_sample_slice(Vcor,'cor',options.d2.bbsize,coords,el);
                    try
                    imagesc(slice,[nanmean(slice(slice>0))-3*nanstd(slice(slice>0)) nanmean(slice(slice>0))+3*nanstd(slice(slice>0))]);
                    catch
                        imagesc(slice);
                    end
                        axis square
                    hold on
                 case 3 % coronar images
                   
                   [slice,boundbox]=ea_sample_slice(Vsag,'sag',options.d2.bbsize,coords,el);
                    try
                    imagesc(slice,[nanmean(slice(slice>0))-3*nanstd(slice(slice>0)) nanmean(slice(slice>0))+3*nanstd(slice(slice>0))]);
                    catch
                        imagesc(slice);
                    end
                        axis square
                    hold on
            end
            
            % Show overlays
            
            if options.d2.writeatlases
                cuts=ea_add_overlay(boundbox,cuts,side,tracor,options.patientname,options);
            end
            
            
            
            % Show coordinates
            
            plot((options.d2.bbsize+1)*interpfactor,(options.d2.bbsize+1)*interpfactor,'*','MarkerSize',15,'MarkerEdgeColor',[0.9 0.9 0.9],'MarkerFaceColor',[0.9 0.9 0.9],'LineWidth',2,'LineSmoothing','on');
            hold off
            set(gca,'LooseInset',get(gca,'TightInset'))
            % Save results
            
            switch tracor
                case 1
                    saveas(cuts,[options.root,options.patientname,filesep,options.elspec.contactnames{el},'_axial.png']);
                case 2
                    saveas(cuts,[options.root,options.patientname,filesep,options.elspec.contactnames{el},'_coronar.png']);
                case 3
                    saveas(cuts,[options.root,options.patientname,filesep,options.elspec.contactnames{el},'_saggital.png']);
            end
        end
    end
    
    
end

close(cuts)


function res=zminus(A,B)
res=A-B;
if res<0; res=0; end
