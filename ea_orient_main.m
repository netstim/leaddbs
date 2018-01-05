function roll_out = ea_orient_main(options,supervised)
%% Determine Orientation for BSCI directed leads from postoperative CT
% has an unsupervised and a supervised version

folder = [options.root options.patientname '\'];

if ~strcmp(options.elmodel,'Boston Scientific Vercise Directed') || options.modality == 1 % check for electrode type and postoperative imaging
    msg = sprintf(['Automatic rotation detection works only for Boston Scientific Vercise Directed leads and postoperative CT images.']);
    choice = questdlg(msg,'Invalid Electrode Type!','Abort','Abort');
    roll_out = [];
else
    %% import CTs and choose which CT to use
    if exist([folder 'rpostop_ct.nii']) == 2
        ct_reg = ea_load_nii([folder 'rpostop_ct.nii']);
        tmat_reg = ct_reg.mat;
    else
        error(['No rpostop_ct.nii found in folder: ' folder])
    end
    
    if exist([folder 'postop_ct.nii']) == 2
        ct_org = ea_load_nii([folder 'postop_ct.nii']);
        tmat_org = ct_org.mat;
        % use postop_ct as default but rpostop_ct if postop_ct has shearing
        if  ct_org.mat(2,3) ~= 0 || ct_org.mat(3,2) ~= 0 || ct_org.mat(1,1) < 0
            msg = sprintf(['Non-orthogonal rotation or shearing found inside the affine matrix of postop_ct.nii. You can either attempt:\n 1. Reslicing the postop_ct, which will be stored as postop_ct_resliced. \n 2. Use the rpostop_ct.nii instead (not recommended). \n 3. Abort and reslice the CT before coregistering and normalizing.']);
            choice = questdlg(msg,'Warning!','Reslice','rPostOpCT','Abort','Reslice');
            switch choice
                case 'Reslice'
                    if ~exist([folder 'postop_ct_resliced.nii'])
                        disp(['Reslicing postop_ct:'])
                        reslice_nii([folder 'postop_ct.nii'],[folder 'postop_ct_resliced.nii'])
                    end
                    ct_org = ea_load_nii([folder 'postop_ct_resliced.nii']);
                    tmat_org = ct_org.mat;
                    ct = ct_org;
                case 'rPostOpCT'
                    disp(['Using rpostop_ct.nii as reference image.'])
                    ct = ct_reg;
                case 'Abort'
                    error('Aborted due to non-orthogonal rotation or shearing found inside the affine matrix of postop_ct.nii')
            end
        else
            disp(['Using postop_ct.nii as reference image.'])
            ct = ct_org;
        end
        
    else
        msg = sprintf(['No postop_ct.nii found in folder: ' folder '\nScript will run using coregistered rpostop_ct.nii which may lead to inaccurate results.']);
        choice = questdlg(msg,'Warning!','Continue','Abort','Abort');
        switch choice
            case 'Continue'
                disp(['Using rpostop_ct.nii as reference image.'])
                ct = ct_reg;
            case 'Abort'
                error('Aborted by user')
        end
    end
    
    pixdim = ct.voxsize;
    
    %% import transformation matrices for CT coregistration
    try
        load([folder 'ea_coregctmethod_applied.mat']);
        
        if strcmp(coregct_method_applied{1},'ea_coregctmri_ants') || strcmp(coregct_method_applied{1},'ea_coregctmri_ants_refine')
            reg2org = load([folder 'anat_t12postop_ct_ants1.mat']);
        elseif strcmp(coregct_method_applied{1},'ea_coregctmri_brainsfit')
            reg2org.fixed = h5read([folder 'postop_ct2anat_t1_brainsfit_Inverse.h5','/TransformGroup/0/TranformFixedParameters']);
            reg2org.AffineTransform_float_3_3 = h5read([folder 'postop_ct2anat_t1_brainsfit_Inverse.h5','/TransformGroup/0/TranformParameters']);
        end
        tmat_reg2org =ea_antsmat2mat(reg2org.AffineTransform_float_3_3,reg2org.fixed);
    catch
        reg2org = load([folder 'Postop_CT_2_T1.mat']);
        tmat_reg2org =ea_antsmat2mat(reg2org.AffineTransform_double_3_3,reg2org.fixed);
        tmat_reg2org = inv(tmat_reg2org);
    end
    
    % colorscale for ct figures
    cscale = [-50 100];
    cscale2 = [-500 2000];
    
    % diameter of the slices shown in visualizations
    extractradius = 30;
    
    sides = {'right','left'};
    for side = options.elside
        disp(['Reconstructing rotation of ' sides{side} ' lead!'])
        
        % import lead information
        load([folder 'ea_reconstruction.mat']); % included in for-loop to make independent ea_save_reconstruction for both sides
        
        %% transform head/tail coordinates from native to image coordinates
        head_native = [reco.native.markers(side).head 1]';
        tail_native = [reco.native.markers(side).tail 1]';
        CTname = find(ct.fname=='\');
        CTname = ct.fname(CTname(end):end);
        if strcmp(CTname,'\postop_ct.nii') || strcmp(CTname,'\postop_ct_resliced.nii')
            % transform rpostop_ct -> postop_ct
            head_mm = (tmat_reg2org) * head_native;
            tail_mm = (tmat_reg2org) * tail_native;
            % transform postop_ct mm -> voxel
            head_vx = inv(tmat_org) * head_mm;
            tail_vx = inv(tmat_org) * tail_mm;
            tmat_vx2mm = tmat_org;
        elseif strcmp(CTname,'\rpostop_ct.nii')
            head_mm = head_native;
            tail_mm = tail_native;
            % transfrom rpostop_ct mm -> voxel
            head_vx = inv(tmat_reg) * head_mm;
            tail_vx = inv(tmat_reg) * tail_mm;
            tmat_vx2mm = tmat_reg;
        end
        
        %% determine location of the stereotactic marker and the directional
        % levels
        unitvector_mm = (tail_mm - head_mm)/norm(tail_mm - head_mm);
        marker_mm = round(head_mm + (10.25 .* unitvector_mm));
        dirlevel1_mm = round(head_mm + (2 .* unitvector_mm));
        dirlevel2_mm = round(head_mm + (4 .* unitvector_mm));
        
        % transform to vx
        marker_vx = round(inv(tmat_vx2mm) * marker_mm);
        dirlevel1_vx = round(inv(tmat_vx2mm) * dirlevel1_mm);
        dirlevel2_vx = round(inv(tmat_vx2mm) * dirlevel2_mm);        
                
        %% extract axial slices at the level of marker and directional electrodes
        if ~supervised
            artifact_marker = ct.img(marker_vx(1)-extractradius:marker_vx(1)+extractradius,marker_vx(2)-extractradius:marker_vx(2)+extractradius, marker_vx(3));
            artifact_dir1 = ct.img(dirlevel1_vx(1)-extractradius:dirlevel1_vx(1)+extractradius,dirlevel1_vx(2)-extractradius:dirlevel1_vx(2)+extractradius, dirlevel1_vx(3));
            artifact_dir2 = ct.img(dirlevel2_vx(1)-extractradius:dirlevel2_vx(1)+extractradius,dirlevel2_vx(2)-extractradius:dirlevel2_vx(2)+extractradius, dirlevel2_vx(3));
        elseif supervised
            %% Identify plane with optimal marker artifact
            tmp{1} = ct.img(marker_vx(1)-extractradius:marker_vx(1)+extractradius,marker_vx(2)-extractradius:marker_vx(2)+extractradius, marker_vx(3)-1);
            tmp{2} = ct.img(marker_vx(1)-extractradius:marker_vx(1)+extractradius,marker_vx(2)-extractradius:marker_vx(2)+extractradius, marker_vx(3));
            tmp{3} = ct.img(marker_vx(1)-extractradius:marker_vx(1)+extractradius,marker_vx(2)-extractradius:marker_vx(2)+extractradius, marker_vx(3)+1);
            
            h = figure('Name',['Lead ' sides{side}],'Position',[100 100 600 800],'Color','w');
            
            subplot(3,1,1)
            imagesc(tmp{1}')
            view(-180,90)
            hold on
            colormap gray
            caxis manual
            caxis(cscale)
            axis equal
            axis off
            title('Center -1')
            
            subplot(3,1,2)
            imagesc(tmp{2}')
            view(-180,90)
            hold on
            colormap gray
            caxis manual
            caxis(cscale)
            axis equal
            axis off
            title('Center')
            
            subplot(3,1,3)
            imagesc(tmp{3}')
            view(-180,90)
            hold on
            colormap gray
            caxis manual
            caxis(cscale)
            axis equal
            axis off
            title('Center +1')
            
            msg = sprintf(['Please select the slice on which the marker artifact is most clearly defined!']);
            choice = questdlg(msg,'Specify slices','Center -1','Center','Center +1', 'Center');
            switch choice
                case 'Center -1'
                    answer = 1;
                case 'Center'
                    answer = 2;
                case 'Center +1'
                    answer = 3;
            end
            
            artifact_marker = tmp{answer};
            marker_vx(3) = marker_vx(3) + answer - 3;
            marker_mm = tmat_vx2mm * marker_vx;
            clear tmp answer
            
            %% Identify plane with optimal dir level 1 artifact
            tmp{1} = ct.img(dirlevel1_vx(1)-extractradius:dirlevel1_vx(1)+extractradius,dirlevel1_vx(2)-extractradius:dirlevel1_vx(2)+extractradius, dirlevel1_vx(3)-1);
            tmp{2} = ct.img(dirlevel1_vx(1)-extractradius:dirlevel1_vx(1)+extractradius,dirlevel1_vx(2)-extractradius:dirlevel1_vx(2)+extractradius, dirlevel1_vx(3));
            tmp{3} = ct.img(dirlevel1_vx(1)-extractradius:dirlevel1_vx(1)+extractradius,dirlevel1_vx(2)-extractradius:dirlevel1_vx(2)+extractradius, dirlevel1_vx(3)+1);
            
            subplot(3,1,1)
            imagesc(tmp{1}')
            view(-180,90)
            hold on
            colormap gray
            caxis manual
            caxis(cscale)
            axis equal
            title('Center -1')
            
            subplot(3,1,2)
            imagesc(tmp{2}')
            view(-180,90)
            hold on
            colormap gray
            caxis manual
            caxis(cscale)
            axis equal
            title('Center')
            
            subplot(3,1,3)
            imagesc(tmp{3}')
            view(-180,90)
            hold on
            colormap gray
            caxis manual
            caxis(cscale)
            axis equal
            title('Center +1')
            
            msg = sprintf(['Please select the slice on which the directional level artifact is most clearly defined!']);
            choice = questdlg(msg,'Specify slices','Center -1','Center','Center +1', 'Center');
            switch choice
                case 'Center -1'
                    answer = 1;
                case 'Center'
                    answer = 2;
                case 'Center +1'
                    answer = 3;
            end
            
            artifact_dir1 = tmp{answer};
            dirlevel1_vx(3) = dirlevel1_vx(3) + answer - 2;
            dirlevel1_mm = tmat_vx2mm * dirlevel1_vx;
            clear tmp answer
            
            %% Identify plane with optimal dir level 2 artifact
            tmp{1} = ct.img(dirlevel2_vx(1)-extractradius:dirlevel2_vx(1)+extractradius,dirlevel2_vx(2)-extractradius:dirlevel2_vx(2)+extractradius, dirlevel2_vx(3)-1);
            tmp{2} = ct.img(dirlevel2_vx(1)-extractradius:dirlevel2_vx(1)+extractradius,dirlevel2_vx(2)-extractradius:dirlevel2_vx(2)+extractradius, dirlevel2_vx(3));
            tmp{3} = ct.img(dirlevel2_vx(1)-extractradius:dirlevel2_vx(1)+extractradius,dirlevel2_vx(2)-extractradius:dirlevel2_vx(2)+extractradius, dirlevel2_vx(3)+1);
            
            subplot(3,1,1)
            imagesc(tmp{1}')
            view(-180,90)
            hold on
            colormap gray
            caxis manual
            caxis(cscale)
            axis equal
            title('Center -1')
            
            subplot(3,1,2)
            imagesc(tmp{2}')
            view(-180,90)
            hold on
            colormap gray
            caxis manual
            caxis(cscale)
            axis equal
            title('Center')
            
            subplot(3,1,3)
            imagesc(tmp{3}')
            view(-180,90)
            hold on
            colormap gray
            caxis manual
            caxis(cscale)
            axis equal
            
            title('Center +1')
            
            msg = sprintf(['Please select the slice on which the directional level artifact is most clearly defined!']);
            choice = questdlg(msg,'Specify slices','Center -1','Center','Center +1', 'Center');
            switch choice
                case 'Center -1'
                    answer = 1;
                case 'Center'
                    answer = 2;
                case 'Center +1'
                    answer = 3;
            end
            
            artifact_dir2 = tmp{answer};
            dirlevel2_vx(3) = dirlevel2_vx(3) + answer - 2;
            dirlevel2_mm = tmat_vx2mm * dirlevel2_vx;
            clear tmp answer
            close(h)
        end
        %% define center of the artifacts
        
        center_marker = [(size(artifact_marker,1)+1)/2 (size(artifact_marker,1)+1)/2];
        center_dir1 = [(size(artifact_dir1,1)+1)/2 (size(artifact_dir1,1)+1)/2];
        center_dir2 = [(size(artifact_dir2,1)+1)/2 (size(artifact_dir2,1)+1)/2];
        
        %% allow for respecification of the artifact centers
        if supervised
            h = figure('Name',['Lead ' sides{side}],'Position',[100 100 600 800],'Color','w');
            imagesc(artifact_marker')
            axis equal
            axis off
            view(-180,90)
            hold on
            colormap gray
            caxis manual
            caxis(cscale2)
            scatter(center_marker(1),center_marker(2),'o','r');
            title('Please mark center of the artifact by doubleclicking')
            [a,b] = getpts;
            center_marker = [a(end) b(end)];
            marker_vx(1) = marker_vx(1)-extractradius + center_marker(1);
            marker_vx(2) = marker_vx(2)-extractradius + center_marker(2);
            marker_mm = tmat_vx2mm * marker_vx;
            clear a b
            close(h)

            h = figure('Name',['Lead ' sides{side}],'Position',[100 100 600 800],'Color','w');
            imagesc(artifact_dir1')
            axis equal
            axis off
            view(-180,90)
            hold on
            colormap gray
            caxis manual
            caxis(cscale2)
            scatter(center_dir1(1),center_dir1(2),'o','r');
            title('Please mark center of the artifact by doubleclicking')
            [a,b] = getpts;
            center_dir1 = [a(end) b(end)];
            dirlevel1_vx(1) = dirlevel1_vx(1)-extractradius + center_dir1(1);
            dirlevel1_vx(2) = dirlevel1_vx(2)-extractradius + center_dir1(2);
            dirlevel1_mm = tmat_vx2mm * dirlevel1_vx;
            clear a b
            close(h)

            h = figure('Name',['Lead ' sides{side}],'Position',[100 100 600 800],'Color','w');
            imagesc(artifact_dir2')
            view(-180,90)
            axis equal
            axis off
            hold on
            colormap gray
            caxis manual
            caxis(cscale2)
            scatter(center_dir2(2),center_dir2(2),'o','r');
            title('Please mark center of the artifact by doubleclicking')
            [a,b] = getpts;
            center_dir2 = [a(end) b(end)];
            dirlevel2_vx(1) = dirlevel2_vx(1)-extractradius + center_dir2(1);
            dirlevel2_vx(2) = dirlevel2_vx(2)-extractradius + center_dir2(2);
            dirlevel2_mm = tmat_vx2mm * dirlevel2_vx;
            clear a b
            close(h)
        end
        %% get intensity profiles at radius around the centers of marker and directional levels
        radius = 3;
        [angle, intensity,vector] = ea_orient_intensityprofile(artifact_marker,center_marker,pixdim,radius);
        radius = 8;
        [angle1, intensity1,vector1] = ea_orient_intensityprofile(artifact_dir1,center_dir1,pixdim,radius);
        [angle2, intensity2,vector2] = ea_orient_intensityprofile(artifact_dir2,center_dir2,pixdim,radius);
        
        %% detect peaks and valleys
        [peak,markerfft] = ea_orient_intensitypeaksFFT(intensity,2);
        valley = ea_orient_intensitypeaksFFT(-intensity,2);
        
        [peak1,~] = ea_orient_intensitypeaksFFT(intensity1,3);
        
        [peak2,~] = ea_orient_intensitypeaksFFT(intensity2,3);
        
        %% determine which of the two peaks is the marker by comparing marker peaks to the 3 dirlevel peaks
        % for this extent intensity from 0:360 to -360:+720 to exclude failures
        % due to detected peaks close to 0 or 360
        peak1tmp = [(peak1-360) peak1 (peak1+360)];
        peak2tmp = [(peak2-360) peak2 (peak2+360)];
        diff1 = min(abs(peak1tmp - peak(1))) + min(abs(peak2tmp - peak(1)));
        diff2 = min(abs(peak1tmp - peak(2))) + min(abs(peak2tmp - peak(2)));
        if ~supervised
            %% take anterior peak
            %         if peak(1) > 90 && peak(1) < 270
            %             finalpeak(side) = peak(2);
            %         else
            %             finalpeak(side) = peak(1);
            %         end
            %% take better peak
            if diff1 <= diff2
                finalpeak(side) = peak(1);
            else
                finalpeak(side) = peak(2);
            end
        elseif supervised
            
            h = figure('Name',['Lead ' sides{side}],'Position',[100 100 600 800],'Color','w');
            hold on
            title('Specify marker direction')
            imagesc(artifact_marker')
            view(-180,-90)
            axis equal
            axis off
            xlim([0 2*extractradius])
            ylim([0 2*extractradius])
            colormap gray
            caxis manual
            caxis(cscale)
            scatter(center_marker(1),center_marker(2),'o','g')
            scatter(vector(peak(1),1), vector(peak(1),2),50,'m','filled');
            scatter(vector(peak(2),1), vector(peak(2),2),50,'c','filled');
            if diff1 <= diff2
                suggestion = 'Peak 1';
            else
                suggestion = 'Peak 2';
            end
            
            msg = sprintf(['Two possible marker directions have been identified:\nPeak 1 = ' num2str(rad2deg(angle(peak(1)))) ' ° (magenta) \nPeak 2 = ' num2str(rad2deg(angle(peak(2)))) ' ° (cyan)\nAutomatic determination suggests ' suggestion '. Please select most likely direction.']);
            choice = questdlg(msg,'Specify marker direction','Peak 1','Peak 2',suggestion);
            switch choice
                case 'Peak 1'
                    finalpeak(side) = peak(1);
                case 'Peak 2'
                    finalpeak(side) = peak(2);
            end
            clear diff1 diff2 peak1tmp peak2tmp
            close(h)
        end
        %% calculate lead yaw and pitch angles for correction at the end
        yaw = asin(unitvector_mm(1));
        pitch = asin(unitvector_mm(2)/cos(yaw));
        
        if rad2deg(abs(pitch)) > 40
            disp(['Warning: Pitch > 40° - Determining orientation might be inaccurate!'])
        end
        if rad2deg(abs(yaw)) > 40
            disp(['Warning: Yaw > 40° - Determining orientation might be inaccurate!'])
        end
        
        %% correction for yaw and pitch to get rollangle for [0 0 1] lead 
        peakangle(side) = angle(finalpeak(side));        
        peakangle_corr(side) = (sin(peakangle(side)) * cos(pitch)) / ((cos(peakangle(side)) * cos(yaw)) - (sin(peakangle(side)) * sin(yaw) * sin(pitch)));  % see Sitz et al. 2017           
        peakangle_corr(side) = atan(peakangle_corr(side));
        
        if peakangle(side) > pi && peakangle(side) < 2 * pi && peakangle_corr(side) > 0
            peakangle_corr(side) = -peakangle_corr(side);
        end
        
        roll = peakangle_corr(side);
        
        %% determine angles of the 6-valley artifact ('dark star') artifact around the directional markers
        
        % approach creating multiple grids for neighboring roll angles
        count = 1;
        shifts = [];
        for k = -30:30
            shifts(count) = k;
            rolltemp = roll + deg2rad(k);
            dir1_angles = ea_orient_darkstar(rolltemp,pitch,yaw,dirlevel1_mm);
            [sumintensity1(count)] = ea_orient_intensitypeaksdirmarker(intensity1,dir1_angles);
            rollangles(count) = rolltemp;
            dir2_angles = ea_orient_darkstar(rolltemp,pitch,yaw,dirlevel2_mm);
            [sumintensity2(count)] = ea_orient_intensitypeaksdirmarker(intensity2,dir2_angles);
            count = count +1;
            clear rolltemp
        end
        clear count
        
        [~, temp] = min(sumintensity1);
        roll1 = rollangles(temp);
        dir1_angles = ea_orient_darkstar(roll1,pitch,yaw,dirlevel1_mm);
        dir1_valleys = round(rad2deg(dir1_angles) +1);
        dir1_valleys(dir1_valleys > 360) = dir1_valleys(dir1_valleys > 360) - 360;
        
        [~, temp] = min(sumintensity2);
        roll2 = rollangles(temp);
        dir2_angles = ea_orient_darkstar(roll2,pitch,yaw,dirlevel2_mm);
        dir2_valleys = round(rad2deg(dir2_angles) +1);
        dir2_valleys(dir2_valleys > 360) = dir2_valleys(dir2_valleys > 360) - 360;
        clear temp
        
        %% final figure
        fig(side).figure = figure('Name',['Lead ' sides{side}],'Position',[100 100 800 800],'Color','w');
        
        if peakangle(side) > pi
            tempangle = peakangle(side) - 2 * pi;
        else
            tempangle = peakangle(side);
        end
        fig(side).txt1 = uicontrol('style','text','units','pixels','Background','w',...
            'position',[650,650,150,75],'FontSize',12,'HorizontalAlignment','left',...
            'string',sprintf(['Artifact angle:\n' num2str(rad2deg(tempangle)) ' °\nPrimary roll angle:\n' num2str(rad2deg(roll)) ' °']));
        clear tempangle
        fig(side).txt2 = uicontrol('style','text','units','pixels','Background','w',...
            'position',[650,450,150,50],'FontSize',12,'HorizontalAlignment','left',...
            'string',sprintf(['Corrected roll angle:\n' num2str(rad2deg(roll1)) ' °']));
        
        fig(side).chk1 = uicontrol('style','checkbox','units','pixels',...
            'position',[650,425,150,25],'string','Accept','FontSize',12,'Background','w');
        
        fig(side).txt3 = uicontrol('style','text','units','pixels','Background','w',...
            'position',[650,225,150,50],'FontSize',12,'HorizontalAlignment','left',...
            'string',sprintf(['Corrected roll angle:\n' num2str(rad2deg(roll2)) ' °']));
        
        fig(side).chk2 = uicontrol('style','checkbox','units','pixels',...
            'position',[650,200,150,25],'string','Accept','FontSize',12,'Background','w');
        
        fig(side).txt4 = uicontrol('style','text','units','pixels','Background','w',...
            'position',[60,60,720,50],'FontSize',12,'HorizontalAlignment','left',...
            'string',sprintf(['Use the checkboxes if you want to correct the primary roll angle by the orientation angles of the directional levels and choose whether you want to save, manually refine, or discard the results.']));
        
        SaveButton = uicontrol('Style', 'pushbutton', 'String', 'Accept & Save',...
            'Position', [150 20 150 25],'FontSize',12,...
            'Callback', @savedirection);
        ManualButton = uicontrol('Style', 'pushbutton', 'String', 'Manual Refine',...
            'Position', [325 20 150 25],'FontSize',12,...
            'Callback', @manualretry);
        DiscardButton = uicontrol('Style', 'pushbutton', 'String', 'Discard',...
            'Position', [500 20 150 25],'FontSize',12,...
            'Callback', @discarddirection);
        
        %% marker
        ax1 = subplot(3,3,1,'Position',[0.05 0.75 0.2 0.2]);
        hold on
        title('Marker')
        imagesc(artifact_marker')
        view(-180,-90)
        axis equal
        axis off
        colormap gray
        caxis manual
        caxis(cscale)
        scatter(center_marker(1),center_marker(2),'o','g')
        plot(vector(:,1),vector(:,2),'g')
        scatter(vector(finalpeak(side),1), vector(finalpeak(side),2),'g','filled');
        
        ax2 = subplot(3,3,2,'Position',[0.3 0.75 0.2 0.2]);
        hold on
        title('Intensity Profile','FontWeight','normal')
        plot(rad2deg(angle),intensity)
        plot(rad2deg(angle),markerfft)
        scatter(rad2deg(angle(peak)), intensity(peak),'g');
        scatter(rad2deg(angle(valley)), intensity(valley),'r');
        set(ax2,'yticklabel',{[]})
        
        %% graphics dir level one
        ax4 = subplot(3,3,4,'Position',[0.05 0.475 0.2 0.2]);
        hold on
        title('Directional Level 1')
        imagesc(artifact_dir1')
        view(-180,-90)
        axis equal
        axis off
        colormap gray
        caxis manual
        caxis(cscale)
        scatter(center_dir1(1),center_dir1(2),'o','g')
        plot(vector1(:,1),vector1(:,2),'g')
        scatter(vector1(dir1_valleys,1),vector1(dir1_valleys,2),'r')
        
        ax5 = subplot(3,3,5,'Position',[0.3 0.475 0.2 0.2]);
        hold on
        title('Intensity Profile','FontWeight','normal')
        plot(rad2deg(angle1),intensity1)
        scatter(rad2deg(angle1(dir1_valleys)), intensity1(dir1_valleys),'r');
        set(ax5,'yticklabel',{[]})
        
        ax6 = subplot(3,3,6,'Position',[0.55 0.475 0.2 0.2]);
        hold on
        title('Correction','FontWeight','normal')
        plot(rad2deg(rollangles),sumintensity1)
        scatter(rad2deg(rollangles(rollangles == roll)), sumintensity1(rollangles == roll),'g','filled');
        scatter(rad2deg(rollangles(rollangles == roll1)), sumintensity1(rollangles == roll1),'r');
        text(rad2deg(rollangles(rollangles == roll1)), sumintensity1(rollangles == roll1),['\leftarrow HU = ' num2str(round(sumintensity1(rollangles == roll1)))]);
        
        set(ax6,'yticklabel',{[]})
        
        %% graphics dir level two
        ax7 = subplot(3,3,7,'Position',[0.05 0.2 0.2 0.2]);
        hold on
        title('Directional Level 2')
        imagesc(artifact_dir2')
        view(-180,-90)
        axis equal
        axis off
        colormap gray
        caxis manual
        caxis(cscale)
        scatter(center_dir2(1),center_dir2(2),'o','g')
        plot(vector2(:,1),vector2(:,2),'g')
        scatter(vector2(dir2_valleys,1),vector2(dir2_valleys,2),'r')
        
        ax8 = subplot(3,3,8,'Position',[0.3 0.2 0.2 0.2]);
        hold on
        title('Intensity Profile','FontWeight','normal')
        plot(rad2deg(angle2),intensity2)
        scatter(rad2deg(angle2(dir2_valleys)), intensity2(dir2_valleys),'r');
        set(ax8,'yticklabel',{[]})
        
        ax9 = subplot(3,3,9,'Position',[0.55 0.2 0.2 0.2]);
        hold on
        title('Correction','FontWeight','normal')
        plot(rad2deg(rollangles),sumintensity2)
        scatter(rad2deg(rollangles(rollangles == roll)), sumintensity2(rollangles == roll),'g','filled');
        scatter(rad2deg(rollangles(rollangles == roll2)), sumintensity2(rollangles == roll2),'r');
        text(rad2deg(rollangles(rollangles == roll2)), sumintensity2(rollangles == roll2),['\leftarrow HU = ' num2str(round(sumintensity2(rollangles == roll2)))]);
        set(ax9,'yticklabel',{[]})
        
        linkaxes([ax2 ax5 ax8],'xy');
        set(ax2,'Xlim',[0 360]);
        set(ax2,'Ylim',[min([intensity intensity1 intensity2])-50 max([intensity intensity1 intensity2])+50]);
        
        linkaxes([ax6 ax9],'xy');
        set(ax6,'Xlim',[rad2deg(rollangles(1)) rad2deg(rollangles(end))]);
        set(ax6,'Ylim',[-1000 1000]);
        
        uiwait
        
        if SaveButton.UserData == 1
            savestate = 1;
        elseif DiscardButton.UserData == 1
            savestate = 0;
            retrystate = 0;
        elseif ManualButton.UserData == 1
            savestate = 0;
            retrystate = 1;
            disp(['Retry with manual refinement!'])
            roll_out_retry = ea_orient_main(options,1);
        end
        
        %% saving results
        if savestate == 1
            clear x y
            checkbox1 = get(fig(side).chk1,'Value');
            checkbox2 = get(fig(side).chk2,'Value');
            
            if ~checkbox1 && ~checkbox2
                roll_y = roll;
                disp(['Using roll angle defined by stereotactic marker: ' num2str(rad2deg(roll)) ' °'])
            elseif checkbox1 && ~checkbox2
                roll_y = roll1;
                disp(['Using corrected roll angle defined by directional level 1: ' num2str(rad2deg(roll1)) ' °'])
            elseif ~checkbox1 && checkbox2
                roll_y = roll2;
                disp(['Using corrected roll angle defined by directional level 2: ' num2str(rad2deg(roll2)) ' °'])
            elseif checkbox1 && checkbox2
                roll_y = mean([roll1 roll2]);
                disp(['Using mean corrected roll angle defined by both directional levels: ' num2str(rad2deg(mean([roll1 roll2]))) ' °'])
            end
            
            %% calculate y
            [M,~,~,~] = ea_orient_rollpitchyaw(roll_y,pitch,yaw);
            y = M * [0;1;0];
            head = head_mm(1:3);
            
            %% transform y to native space and back
            y = head + y;
            y(4) = 1;
            if strcmp(CTname,'\postop_ct.nii') || strcmp(CTname,'\postop_ct_resliced.nii')
                % transform postop_ct_mm -> rpostop_ct_mm
                y = inv(tmat_reg2org) * y;
            elseif strcmp(CTname,'\rpostop_ct.nii')
                y = y;
            end
            y = y(1:3)';
            head = head_native(1:3)';
            tail = tail_native(1:3)';
            y = y - head;
            
            %% project y down on t and calculate x
            y = y/norm(y);
            t = diff([head;tail])/norm(diff([head;tail]));
            y = y - (dot(y,t) / (norm(t) ^2)) * t;
            y = y/norm(y);
            x = cross(y,t);
            
            y = (y / norm(y)) * (options.elspec.lead_diameter / 2);
            x = (x / norm(x)) * (options.elspec.lead_diameter / 2);
            head = head_native(1:3)';
            y_out = y;
            y = head + y;
            x = head + x;
            
            reco.native.markers(side).y = y;
            reco.native.markers(side).x = x;
            %% for direct saving into manual reconstruction
            [coords,trajectory,markers]=ea_resolvecoords(reco.native.markers,options);
            ea_save_reconstruction(coords,trajectory,markers,options.elmodel,1,options)
            
            % %% for transfering to ea_manualreconstruction
            y_out(3) = 0;
            y_out = y_out / norm(y_out);
            roll_out = rad2deg(atan2(norm(cross([0 1 0],y_out)),dot([0 1 0],y_out)));
            if markers(side).y(1) > markers(side).head(1) % negative 90 points to right, positive 90 points to left
                roll_out = - roll_out;
            end
            disp(['Corrected roll angle roll = ' num2str(rad2deg(roll_y)) ' °, has been converted to orientation angle = ' num2str(roll_out) ' for compatibility with ea_mancorupdatescene.'])
            %% methods dump:
ea_methods(options,...
            ['Rotation of directional DBS leads was determined using the algorithm published by Sitz et al. 2017 as implemented in LEAD-DBS software.'],...
            {'Sitz, A., Hoevels, M., Hellerbach, A., Gierich, A., Luyken, K., Dembek, T.A., Klehr, M., Wirths, J., Visser-Vandewalle, V., & Treuer, H. (2017). Determining the orientation angle of directional leads for deep brain stimulation using computed tomography and digital x-ray imaging: A phantom study. Medical Physics, 44(9):4463-4473. http://dx.doi.org/10.1002/mp.12424'});
        elseif retrystate == 0
            disp(['Changes to rotation not saved'])
            roll_out = [];
        elseif retrystate == 1
            disp(['Rotation saved after manual refinement'])
            roll_out = roll_out_retry;
        end
        close(fig(side).figure)
        clear reco
    end
end

end

function savedirection(hObject,eventdata)
hObject.UserData = 1;
uiresume
end
function discarddirection(hObject,eventdata)
hObject.UserData = 1;
uiresume
end
function manualretry(hObject,eventdata)
hObject.UserData = 1;
uiresume
end
