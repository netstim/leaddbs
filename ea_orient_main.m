function roll_out = ea_orient_main(options,supervised)
%% Determine Orientation for BSCI directed leads from postoperative CT
% has an unsupervised and a supervised version

folder = [options.root options.patientname filesep];

if  options.modality == 1 % check for electrode type and postoperative imaging
    msg = sprintf(['Automatic rotation detection works only for postoperative CT images.']);
    choice = questdlg(msg,'No postoperative CT!','Abort','Abort');
    roll_out = [];
elseif strcmp(options.elmodel,'Boston Scientific Vercise Directed') || strcmp(options.elmodel,'St. Jude Directed 6172 (short)') || strcmp(options.elmodel,'St. Jude Directed 6173 (long)')
    if ismember(options.elmodel,{'St. Jude Directed 6172 (short)','St. Jude Directed 6173 (long)'})
        disp(['Warning: DiODe algorithm not validated for ' options.elmodel '.'])
    end

    markerposition = options.elspec.markerpos;
    markerlength = options.elspec.markerlen;
    electrodespacing = options.elspec.contact_length+options.elspec.contact_spacing;
    contactlength = options.elspec.contact_length;
    tipInsulationlength = options.elspec.tip_length*~options.elspec.tipiscontact;

    level1center = tipInsulationlength + electrodespacing+contactlength/2;
    level2center = tipInsulationlength + electrodespacing*2+contactlength/2;
    markercenter = markerposition + markerlength/2;

    %%
    load(options.elspec.matfname)
    %% import CTs and choose which CT to use
    if exist([folder options.prefs.ctnii_coregistered],'file') == 2
        ct_reg = ea_load_nii([folder options.prefs.ctnii_coregistered]);
        tmat_reg = ct_reg.mat;
    else
        error(['No coregistered CT (',options.prefs.ctnii_coregistered,') found in folder: ' folder])
    end
    tol=0.0001; % tolerance of (rounding) difference in qform and sform matrices.
    if exist([folder options.prefs.rawctnii_unnormalized],'file') == 2
        ct_org = ea_load_nii([folder 'postop_ct.nii']);
        if sum(sum(abs(ct_org.mat-ct_org.private.mat0)))<tol
            tmat_org = ct_org.mat;
        else
            msg = sprintf(['Warning: Different sForm and qForm matrices in Nifti-object. Please select the matrix you want to use.']);
            choice = questdlg(msg,'Warning!','sForm','qForm','sForm');
            switch choice
                case 'sForm'
                    tmat_org = ct_org.mat;
                case 'qForm'
                    tmat_org = ct_org.private.mat0;
            end
        end
        ct = ct_org;

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
    tmat_reg2org=eye(4); % default.
    try
        if strcmp(options.prefs.reco.mancoruse,'postop')
        load([folder 'ea_coregctmethod_applied.mat']);
        switch coregct_method_applied{end}
            case 'ea_coregctmri_fsl'
                %             tmat_reg2org = dlmread([folder 'anat_t12postop_ct_flirt1.mat']));
                disp(['Warning: Temporary fix to use DiODe algorithm with FLIRT. rpostop_ct is used so results may be slightly less accurate.'])
                ct = ct_reg;
            otherwise
                [tmat_reg2org,ctfname] = ea_getrawct2preniimat(options,1);
                ct=ea_load_nii(ctfname);
        end
        else
            ct = ct_reg;
        end
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

    sides = {'right','left','3','4','5','6','7','8'};
    for side = options.elside
        disp(['Reconstructing rotation of ' sides{side} ' lead!'])

        % import lead information
        load([folder 'ea_reconstruction.mat']); % included in for-loop to make independent ea_save_reconstruction for both sides

        %% transform head/tail coordinates from native to image coordinates
        head_native = [reco.native.markers(side).head 1]';
        tail_native = [reco.native.markers(side).tail 1]';
        CTname = find(ct.fname==filesep);
        CTname = ct.fname(CTname(end):end);
        if strcmp(CTname,[filesep,options.prefs.rawctnii_unnormalized]) || strcmp(CTname,[filesep,'postop_ct_resliced.nii'])
            % transform rpostop_ct -> postop_ct
            head_mm = (tmat_reg2org) * head_native;
            tail_mm = (tmat_reg2org) * tail_native;
            % transform postop_ct mm -> voxel
            head_vx = inv(tmat_org) * head_mm;
            tail_vx = inv(tmat_org) * tail_mm;
            tmat_vx2mm = tmat_org;
        elseif strcmp(CTname,[filesep,options.prefs.ctnii_coregistered])
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
        marker_mm = round(head_mm + (markerposition .* unitvector_mm));
        dirlevel1_mm = round(head_mm + (electrodespacing .* unitvector_mm));
        dirlevel2_mm = round(head_mm + (2 * electrodespacing .* unitvector_mm));

        % transform to vx
        marker_vx = round(tmat_vx2mm\marker_mm);
        dirlevel1_vx = round(tmat_vx2mm\dirlevel1_mm);
        dirlevel2_vx = round(tmat_vx2mm\dirlevel2_mm);

        %% extract axial slices at the level of marker and directional electrodes
        if ~supervised
            artifact_marker=ea_sample_slice(ct,'tra',extractradius,'vox',{marker_vx(1:3)'},1)';
            if ct.mat(1,1) < 0
                artifact_marker = flip(artifact_marker,1);
            end
            if ct.mat(2,2) < 0
                artifact_marker = flip(artifact_marker,2);
            end

            artifact_dir1=ea_sample_slice(ct,'tra',extractradius,'vox',{dirlevel1_vx(1:3)'},1)';
            if ct.mat(1,1) < 0
                artifact_dir1 = flip(artifact_dir1,1);
            end
            if ct.mat(2,2) < 0
                artifact_dir1 = flip(artifact_dir1,2);
            end

            artifact_dir2=ea_sample_slice(ct,'tra',extractradius,'vox',{dirlevel2_vx(1:3)'},1)';
            if ct.mat(1,1) < 0
                artifact_dir2 = flip(artifact_dir2,1);
            end
            if ct.mat(2,2) < 0
                artifact_dir2 = flip(artifact_dir2,2);
            end
        elseif supervised
            h = figure('Name',['Respicify Slices for ' sides{side} ' Lead'],'Position',[100 100 600 800],'Color','w');
            txt1 = uicontrol('style','text','units','pixels','Background','w',...
                'position',[10,770,600,25],'FontSize',12,'HorizontalAlignment','center','FontWeight','bold',...
                'string',sprintf(['Please specify the slice with the most clearly defined artifact:']));
            plusoneButton = uicontrol('Style', 'pushbutton', 'String', '+ 1 slice',...
                'Position', [425 640 150 25],'FontSize',12,...
                'Callback', @buttonPress);
            centerButton = uicontrol('Style', 'pushbutton', 'String', 'Center',...
                'Position', [425 400 150 25],'FontSize',12,...
                'Callback', @buttonPress);
            minusoneButton = uicontrol('Style', 'pushbutton', 'String', '- 1 slice',...
                'Position', [425 160 150 25],'FontSize',12,...
                'Callback', @buttonPress);
            %% Identify plane with optimal marker artifact
            tmp{1}=ea_sample_slice(ct,'tra',extractradius,'vox',{marker_vx(1:3)' - [0 0 1]},1)';
            tmp{2}=ea_sample_slice(ct,'tra',extractradius,'vox',{marker_vx(1:3)'},1)';
            tmp{3}=ea_sample_slice(ct,'tra',extractradius,'vox',{marker_vx(1:3)' + [0 0 1]},1)';

            [tmp] = ea_orient_respecifyslices(tmp,ct,cscale,electrode,9);

            uiwait
            if plusoneButton.UserData == 1
                answer = 3;
                plusoneButton.UserData = 0;
            elseif centerButton.UserData == 1
                answer = 2;
                centerButton.UserData = 0;
            elseif minusoneButton.UserData == 1
                answer = 1;
                minusoneButton.UserData = 0;
            end

            artifact_marker = tmp{answer};
            marker_vx(3) = marker_vx(3) + answer - 2;
            marker_mm = tmat_vx2mm * marker_vx;
            clear tmp answer

            %% Identify plane with optimal dir level 1 artifact
            tmp{1}=ea_sample_slice(ct,'tra',extractradius,'vox',{dirlevel1_vx(1:3)' - [0 0 1]},1)';
            tmp{2}=ea_sample_slice(ct,'tra',extractradius,'vox',{dirlevel1_vx(1:3)'},1)';
            tmp{3}=ea_sample_slice(ct,'tra',extractradius,'vox',{dirlevel1_vx(1:3)' + [0 0 1]},1)';

            [tmp] = ea_orient_respecifyslices(tmp,ct,cscale,electrode,[2 3 4]);

            uiwait
            if plusoneButton.UserData == 1
                answer = 3;
                plusoneButton.UserData = 0;
            elseif centerButton.UserData == 1
                answer = 2;
                centerButton.UserData = 0;
            elseif minusoneButton.UserData == 1
                answer = 1;
                minusoneButton.UserData = 0;
            end

            artifact_dir1 = tmp{answer};
            dirlevel1_vx(3) = dirlevel1_vx(3) + answer - 2;
            dirlevel1_mm = tmat_vx2mm * dirlevel1_vx;
            clear tmp answer

            %% Identify plane with optimal dir level 2 artifact
            tmp{1}=ea_sample_slice(ct,'tra',extractradius,'vox',{dirlevel2_vx(1:3)' - [0 0 1]},1)';
            tmp{2}=ea_sample_slice(ct,'tra',extractradius,'vox',{dirlevel2_vx(1:3)'},1)';
            tmp{3}=ea_sample_slice(ct,'tra',extractradius,'vox',{dirlevel2_vx(1:3)' + [0 0 1]},1)';

            [tmp] = ea_orient_respecifyslices(tmp,ct,cscale,electrode,[5 6 7]);

            uiwait
            if plusoneButton.UserData == 1
                answer = 3;
                plusoneButton.UserData = 0;
            elseif centerButton.UserData == 1
                answer = 2;
                centerButton.UserData = 0;
            elseif minusoneButton.UserData == 1
                answer = 1;
                minusoneButton.UserData = 0;
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
            h = figure('Name',['Respicify Artifact Centers for ' sides{side} ' Lead'],'Position',[100 100 600 800],'Color','w');
            txt1 = uicontrol('style','text','units','pixels','Background','w',...
                'position',[50,770,550,25],'FontSize',12,'HorizontalAlignment','center','FontWeight','bold',...
                'string',sprintf(['Please mark the center of the artifact by doubleclicking:']));
            imagesc(artifact_marker');
            axis equal
            axis off
            view(-180,90)
            hold on
            colormap(gray)
            caxis manual
            caxis(cscale2)
            scatter(center_marker(1),center_marker(2),'o','r');
            [a,b] = getpts;
            center_marker = [a(end) b(end)];
            marker_vx(1) = marker_vx(1)-extractradius + center_marker(1);
            marker_vx(2) = marker_vx(2)-extractradius + center_marker(2);
            marker_mm = tmat_vx2mm * marker_vx;
            clear a b

            hold off
            imagesc(artifact_dir1')
            axis equal
            axis off
            view(-180,90)
            hold on
            colormap(gray)
            caxis manual
            caxis(cscale2)
            scatter(center_dir1(1),center_dir1(2),'o','r');
            [a,b] = getpts;
            center_dir1 = [a(end) b(end)];
            dirlevel1_vx(1) = dirlevel1_vx(1)-extractradius + center_dir1(1);
            dirlevel1_vx(2) = dirlevel1_vx(2)-extractradius + center_dir1(2);
            dirlevel1_mm = tmat_vx2mm * dirlevel1_vx;
            clear a b

            hold off
            imagesc(artifact_dir2')
            view(-180,90)
            axis equal
            axis off
            hold on
            colormap(gray)
            caxis manual
            caxis(cscale2)
            scatter(center_dir2(2),center_dir2(2),'o','r');
            [a,b] = getpts;
            center_dir2 = [a(end) b(end)];
            dirlevel2_vx(1) = dirlevel2_vx(1)-extractradius + center_dir2(1);
            dirlevel2_vx(2) = dirlevel2_vx(2)-extractradius + center_dir2(2);
            dirlevel2_mm = tmat_vx2mm * dirlevel2_vx;
            clear a b
            close(h)
        end
        %% get intensity profiles at radius around the centers of marker and directional levels
        radius = 4;
        radius = radius *2;
        [angle, intensity,vector] = ea_orient_intensityprofile(artifact_marker,center_marker,pixdim,radius);

        radius = 8;
        radius = radius *2;
        [angle1, intensity1,vector1] = ea_orient_intensityprofile(artifact_dir1,center_dir1,pixdim,radius);
        [angle2, intensity2,vector2] = ea_orient_intensityprofile(artifact_dir2,center_dir2,pixdim,radius);

        %% detect peaks and valleys
        [peak,markerfft] = ea_orient_intensitypeaksFFT(intensity,2);
        valley = ea_orient_intensitypeaksFFT(-intensity,2);

        %% take solution which is closer to the default direction - change defaultdirection if you do not rountinely implant in anterior direction
        defaultdirection = 'anterior';              % also use: 'posterior', 'medial', 'lateral'
        if ~supervised
            %% take solution which is closer to the default direction - change defaultdirection if you do not rountinely implant in anterior direction
            switch defaultdirection
                case 'anterior'
                    %% take anterior peak
                    if peak(1) > 90 && peak(1) < 270
                        finalpeak(side) = peak(2);
                    else
                        finalpeak(side) = peak(1);
                    end
                case 'posterior'
                    %% take posterior peak
                    if peak(1) > 90 && peak(1) < 270
                        finalpeak(side) = peak(1);
                    else
                        finalpeak(side) = peak(2);
                    end
                case 'medial'
                    if strcmp(sides{side},'right') % in case of more than 2 leads additional leads have to be specified here
                        if peak(1) > 180
                            finalpeak(side) = peak(2);
                        else
                            finalpeak(side) = peak(1);
                        end
                    elseif strcmp(sides{side},'left') % in case of more than 2 leads additional leads have to be specified here
                        if peak(1) > 180
                            finalpeak(side) = peak(1);
                        else
                            finalpeak(side) = peak(2);
                        end
                    end
                case 'lateral'
                    if strcmp(sides{side},'right')
                        if peak(1) > 180
                            finalpeak(side) = peak(1);
                        else
                            finalpeak(side) = peak(2);
                        end
                    elseif strcmp(sides{side},'left')
                        if peak(1) > 180
                            finalpeak(side) = peak(2);
                        else
                            finalpeak(side) = peak(1);
                        end
                    end
            end
            %% take "better" peak (not validated) to determine which of the two peaks is the marker
            %% by comparing marker peaks to the 3 dirlevel peaks
            % for this extent intensity from 0:360 to -360:+720 to exclude failures
            % due to detected peaks close to 0 or 360
            %         [peak1,~] = ea_orient_intensitypeaksFFT(intensity1,3);
            %         [peak2,~] = ea_orient_intensitypeaksFFT(intensity2,3);
            %             peak1tmp = [(peak1-360) peak1 (peak1+360)];
            %             peak2tmp = [(peak2-360) peak2 (peak2+360)];
            %             diff1 = min(abs(peak1tmp - peak(1))) + min(abs(peak2tmp - peak(1)));
            %             diff2 = min(abs(peak1tmp - peak(2))) + min(abs(peak2tmp - peak(2)));
            %             if diff1 <= diff2
            %                 finalpeak(side) = peak(1);
            %             else
            %                 finalpeak(side) = peak(2);
            %             end
            %             clear diff1 diff2 peak1tmp peak2tmp peak1 peak2
        elseif supervised
            peak1tmp = rad2deg(angle(peak(1)));
            peak2tmp = rad2deg(angle(peak(2)));
            if peak1tmp > 180
                peak1tmp = peak1tmp -360;
            end
            if peak2tmp > 180
                peak2tmp = peak2tmp -360;
            end

            h = figure('Name',['Lead ' sides{side}],'Position',[100 100 600 800],'Color','w');
            txt1 = uicontrol('style','text','units','pixels','Background','w',...
                'position',[50,770,500,25],'FontSize',12,'HorizontalAlignment','center','FontWeight','bold',...
                'string',sprintf(['Please select the marker direction:']));
            Solution1Button = uicontrol('Style', 'pushbutton', 'String', 'Solution 1',...
                'Position', [100 80 175 25],'FontSize',12,...
                'Callback', @buttonPress);
            Solution2Button = uicontrol('Style', 'pushbutton', 'String', 'Solution 2',...
                'Position', [350 80 175 25],'FontSize',12,...
                'Callback', @buttonPress);

            ax1 = axes('Position',[0.1 0.25 0.8 0.8]);
            hold on
            imagesc(ax1,artifact_marker');

            view(-180,-90)
            axis equal
            axis off
            colormap(gray)
            caxis manual
            caxis(cscale)
            scatter(center_marker(1),center_marker(2),'o','g')
            scatter(vector(peak(1),1), vector(peak(1),2),80,'g','filled');
            scatter(vector(peak(2),1), vector(peak(2),2),80,'g','filled');
            quiver(center_marker(1),center_marker(2),vector(peak(1),1)-center_marker(1) ,vector(peak(1),2)-center_marker(2),2,'LineWidth',2,'Color','g','MaxHeadSize',2)
            quiver(center_marker(1),center_marker(2),vector(peak(2),1)-center_marker(1) ,vector(peak(2),2)-center_marker(2),2,'LineWidth',2,'Color','g','MaxHeadSize',2)

            plot(vector(:,1),vector(:,2),'g','LineStyle',':','LineWidth',2)
            for k = 1:length(valley)
                plot([center_marker(1) (center_marker(1) + 3 * (vector(valley(k),1)-center_marker(1)))],...
                [center_marker(2) (center_marker(2) + 3 * (vector(valley(k),2)-center_marker(2)))],'r','LineStyle','--','LineWidth',2)
            end
            text(vector(peak(1),1)-5 ,vector(peak(1),2),['Solution 1'],'FontSize', 18, 'Color','g','HorizontalAlignment','left','VerticalAlignment','middle')
            text(vector(peak(2),1)+5 ,vector(peak(2),2),['Solution 2'],'FontSize', 18, 'Color','g', 'HorizontalAlignment','right','VerticalAlignment','middle')

            txt2 = uicontrol('style','text','units','pixels','Background','w',...
                'position',[100,110,500,150],'FontSize',12,'HorizontalAlignment','left',...
                'string',sprintf(['Two possible solutions have been identified:\n\nSolution 1 = ' num2str(peak1tmp) ' deg \nSolution 2 = ' num2str(peak2tmp) ' deg \n\nPlease select the most likely direction.']));
            uiwait

            if Solution1Button.UserData == 1
                answer = 1;
                Solution1Button.UserData = 0;
            elseif Solution2Button.UserData == 1
                answer = 2;
                Solution2Button.UserData = 0;
            end
            finalpeak(side) = peak(answer);
            clear answer peak1tmp peak2tmp
            close(h)
        end
        %% calculate lead yaw and pitch angles for correction at the end
        yaw = asin(unitvector_mm(1));
        pitch = asin(unitvector_mm(2)/cos(yaw));

        if rad2deg(abs(pitch)) > 40
            disp(['Warning: Pitch > 40 deg - Determining orientation might be inaccurate!'])
        end
        if rad2deg(abs(yaw)) > 40
            disp(['Warning: Yaw > 40 deg - Determining orientation might be inaccurate!'])
        end

        %% correction for yaw and pitch to get rollangle for [0 0 1] lead
        peakangle(side) = angle(finalpeak(side));
        peakangle_corr(side) = (sin(peakangle(side)) * cos(pitch)) / ((cos(peakangle(side)) * cos(yaw)) - (sin(peakangle(side)) * sin(yaw) * sin(pitch)));  % see Sitz et al. 2017
        peakangle_corr(side) = atan(peakangle_corr(side));

        if peakangle(side) < pi && peakangle_corr(side) < 0 && peakangle(side) - peakangle_corr(side) > pi/2
            peakangle_corr(side) = peakangle_corr(side) + pi;
        end
        if peakangle(side) > pi && peakangle_corr(side) > 0 && peakangle(side) - peakangle_corr(side) > pi/2
            peakangle_corr(side) = peakangle_corr(side) - pi;
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
        fig(side).figure = figure('Name',['Lead ' sides{side}],'Position',[100 100 800 800],'Color','w','Toolbar','none');

        if peakangle(side) > pi
            tempangle = peakangle(side) - 2 * pi;
        else
            tempangle = peakangle(side);
        end
        fig(side).txt1 = uicontrol('style','text','units','pixels','Background','w',...
            'position',[650,650,180,75],'FontSize',12,'HorizontalAlignment','left',...
            'string',sprintf(['Artifact Angle:\n' num2str(rad2deg(tempangle),'%.1f') ' deg\nMarker Angle:\n' num2str(rad2deg(roll),'%.1f') ' deg']));

        fig(side).txt2 = uicontrol('style','text','units','pixels','Background','w',...
            'position',[650,450,180,50],'FontSize',12,'HorizontalAlignment','left',...
            'string',sprintf(['Level 2 Angle:\n' num2str(rad2deg(roll2),'%.1f') ' deg']));

        fig(side).chk2 = uicontrol('style','checkbox','units','pixels',...
            'position',[650,425,180,25],'string','Accept','FontSize',12,'Background','w');

        fig(side).txt3 = uicontrol('style','text','units','pixels','Background','w',...
            'position',[650,225,180,50],'FontSize',12,'HorizontalAlignment','left',...
            'string',sprintf(['Level 1 Angle:\n' num2str(rad2deg(roll1),'%.1f') ' deg']));

        fig(side).chk1 = uicontrol('style','checkbox','units','pixels',...
            'position',[650,200,180,25],'string','Accept','FontSize',12,'Background','w');

        fig(side).txt4 = uicontrol('style','text','units','pixels','Background','w',...
            'position',[60,60,720,40],'FontSize',12,'HorizontalAlignment','left',...
            'string',sprintf(['Use the checkboxes if the algorithm accurately detected the artifacts of the directional levels and if you want to use them to correct the marker angle. Then accept, manually refine, or discard the results.']));

        if rad2deg(abs(pitch)) > 40 || rad2deg(abs(yaw)) > 40
            fig(side).txt5 = uicontrol('style','text','units','pixels','Background','w','ForegroundColor','r',...
                'position',[60,100,720,40],'FontSize',12,'HorizontalAlignment','left',...
                'string',sprintf(['WARNING: The polar angle of the lead is larger than 40 deg and results could be inaccurate.\nPlease inspect the results carefully and use manual refinement if necessary.']));
        elseif rad2deg(abs(roll)) > 60
            fig(side).txt5 = uicontrol('style','text','units','pixels','Background','w','ForegroundColor','r',...
                'position',[60,100,720,40],'FontSize',12,'HorizontalAlignment','left',...
                'string',sprintf(['WARNING: The orientation of the lead is far from ' defaultdirection '.\nPlease verify whether the correct marker orientation has been chosen.']));
        else
             fig(side).txt5 = uicontrol('style','text','units','pixels','Background','w','ForegroundColor','k',...
                'position',[60,100,720,40],'FontSize',12,'HorizontalAlignment','left',...
                'string',sprintf(['No warnings: The polar angle and lead orientation are within normal ranges.']));
        end

        SaveButton = uicontrol('Style', 'pushbutton', 'String', 'Accept & Save',...
            'Position', [150 20 150 25],'FontSize',12,...
            'Callback', @buttonPress);
        ManualButton = uicontrol('Style', 'pushbutton', 'String', 'Manual Refine',...
            'Position', [325 20 150 25],'FontSize',12,...
            'Callback', @buttonPress);
        DiscardButton = uicontrol('Style', 'pushbutton', 'String', 'Discard',...
            'Position', [500 20 150 25],'FontSize',12,...
            'Callback', @buttonPress);
        %% marker
        ax1 = subplot(3,3,1);
        hold on
        title(ax1,'Marker')
        imagesc(artifact_marker')
        view(-180,-90)
        axis equal
        axis off
        colormap(gray)
        caxis manual
        caxis(cscale)
        scatter(ax1,center_marker(1),center_marker(2),'o','g')
        plot(ax1,vector(:,1),vector(:,2),'g','LineStyle',':')
        scatter(ax1,vector(finalpeak(side),1), vector(finalpeak(side),2),'g','filled');
        quiver(center_marker(1),center_marker(2),vector(finalpeak(side),1)-center_marker(1) ,vector(finalpeak(side),2)-center_marker(2),2,'LineWidth',1,'Color','g','MaxHeadSize',2)
        for k = 1:length(valley)
            plot(ax1,[center_marker(1) (center_marker(1) + 3 * (vector(valley(k),1)-center_marker(1)))],...
                [center_marker(2) (center_marker(2) + 3 * (vector(valley(k),2)-center_marker(2)))],'r','LineStyle','--')
        end
        xlimit = get(ax1,'Xlim');
        ylimit = get(ax1,'Ylim');
        text(mean(xlimit),ylimit(2) - 0.15 * mean(ylimit),'A','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
        text(mean(xlimit),ylimit(1) + 0.15 * mean(ylimit),'P','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
        text(xlimit(1) + 0.1 * mean(xlimit),mean(ylimit),'L','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
        text(xlimit(2) - 0.1 * mean(xlimit),mean(ylimit),'R','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')

        ax2 = subplot(3,3,2);
        hold on
        title(ax2,'Intensity Profile','FontWeight','normal')
        plot(ax2,rad2deg(angle),intensity)
        plot(ax2,rad2deg(angle),markerfft)
        scatter(ax2,rad2deg(angle(peak)), intensity(peak),'g');
        scatter(ax2,rad2deg(angle(valley)), intensity(valley),'r');
        set(ax2,'yticklabel',{[]})

        %% graphics dir level one
        ax4 = subplot(3,3,4);
        hold on
        title(ax4,'Directional Level 1')
        imagesc(artifact_dir1')
        view(-180,-90)
        axis equal
        axis off
        colormap(gray)
        caxis manual
        caxis(cscale)
        scatter(ax4,center_dir1(1),center_dir1(2),'o','g')
        plot(ax4,vector1(:,1),vector1(:,2),'g','LineStyle',':')
        scatter(ax4,vector1(dir1_valleys,1),vector1(dir1_valleys,2),'r')
        for k = 1:length(dir1_valleys)
            plot(ax4,[center_dir1(1) (center_dir1(1) + 1.5 * (vector1(dir1_valleys(k),1)-center_dir1(1)))],...
                [center_dir1(2) (center_dir1(2) + 1.5 * (vector1(dir1_valleys(k),2)-center_dir1(2)))],'r','LineStyle','--')
        end
        xlimit = get(ax4,'Xlim');
        ylimit = get(ax4,'Ylim');
        text(mean(xlimit),ylimit(2) - 0.15 * mean(ylimit),'A','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
        text(mean(xlimit),ylimit(1) + 0.15 * mean(ylimit),'P','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
        text(xlimit(1) + 0.1 * mean(xlimit),mean(ylimit),'L','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
        text(xlimit(2) - 0.1 * mean(xlimit),mean(ylimit),'R','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')

        ax5 = subplot(3,3,5);
        hold on
        title(ax5,'Intensity Profile','FontWeight','normal')
        plot(ax5,rad2deg(angle1),intensity1)
        scatter(ax5,rad2deg(angle1(dir1_valleys)), intensity1(dir1_valleys),'r');
        set(ax5,'yticklabel',{[]})

        ax6 = subplot(3,3,6);
        hold on
        title(ax6,'Similarity Index','FontWeight','normal')
        plot(ax6,rad2deg(rollangles),sumintensity1)
        scatter(ax6,rad2deg(rollangles(rollangles == roll)), sumintensity1(rollangles == roll),'g','filled');
        scatter(ax6,rad2deg(rollangles(rollangles == roll1)), sumintensity1(rollangles == roll1),'r');
        text(rad2deg(rollangles(rollangles == roll1)), sumintensity1(rollangles == roll1),['\leftarrow HU = ' num2str(round(sumintensity1(rollangles == roll1)))]);
        set(ax6,'yticklabel',{[]})

        %% graphics dir level two
        ax7 = subplot(3,3,7);
        hold on
        title(ax7,'Directional Level 2')
        imagesc(artifact_dir2')
        view(-180,-90)
        axis equal
        axis off
        colormap(gray)
        caxis manual
        caxis(cscale)
        scatter(ax7,center_dir2(1),center_dir2(2),'o','g')
        plot(ax7,vector2(:,1),vector2(:,2),'g','LineStyle',':')
        scatter(ax7,vector2(dir2_valleys,1),vector2(dir2_valleys,2),'r')
        for k = 1:length(dir2_valleys)
            plot(ax7,[center_dir2(1) (center_dir2(1) + 1.5 * (vector2(dir2_valleys(k),1)-center_dir2(1)))],...
                [center_dir2(2) (center_dir2(2) + 1.5 * (vector2(dir2_valleys(k),2)-center_dir2(2)))],'r','LineStyle','--')
        end
        xlimit = get(ax7,'Xlim');
        ylimit = get(ax7,'Ylim');
        text(mean(xlimit),ylimit(2) - 0.15 * mean(ylimit),'A','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
        text(mean(xlimit),ylimit(1) + 0.15 * mean(ylimit),'P','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
        text(xlimit(1) + 0.1 * mean(xlimit),mean(ylimit),'L','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
        text(xlimit(2) - 0.1 * mean(xlimit),mean(ylimit),'R','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')

        ax8 = subplot(3,3,8);
        hold on
        title(ax8,'Intensity Profile','FontWeight','normal')
        plot(ax8,rad2deg(angle2),intensity2)
        scatter(ax8,rad2deg(angle2(dir2_valleys)), intensity2(dir2_valleys),'r');
        set(ax8,'yticklabel',{[]})

        ax9 = subplot(3,3,9);
        hold on
        title(ax9,'Similarity Index','FontWeight','normal')
        plot(ax9,rad2deg(rollangles),sumintensity2)
        scatter(ax9,rad2deg(rollangles(rollangles == roll)), sumintensity2(rollangles == roll),'g','filled');
        scatter(ax9,rad2deg(rollangles(rollangles == roll2)), sumintensity2(rollangles == roll2),'r');
        text(rad2deg(rollangles(rollangles == roll2)), sumintensity2(rollangles == roll2),['\leftarrow HU = ' num2str(round(sumintensity2(rollangles == roll2)))]);
        set(ax9,'yticklabel',{[]})

        linkaxes([ax2 ax5 ax8],'xy');
        set(ax2,'Xlim',[0 360]);
        set(ax2,'Ylim',[min([intensity intensity1 intensity2])-50 max([intensity intensity1 intensity2])+50]);

        linkaxes([ax6 ax9],'xy');
        set(ax6,'Xlim',[rad2deg(rollangles(1)) rad2deg(rollangles(end))]);
        set(ax6,'Ylim',[-1000 1000]);

        set(ax1,'Position',[0.13 0.75 0.2 0.2])
        set(ax2,'Position',[0.345 0.75 0.2 0.2])
        set(ax7,'Position',[0.13 0.475 0.2 0.2])
        set(ax8,'Position',[0.345 0.475 0.2 0.2])
        set(ax9,'Position',[0.56 0.475 0.2 0.2])
        set(ax4,'Position',[0.13 0.2 0.2 0.2])
        set(ax5,'Position',[0.345 0.2 0.2 0.2])
        set(ax6,'Position',[0.56 0.2 0.2 0.2])
        %% graphics lead
        ax_elec = axes('Position',[0 0.2 0.1 0.75]);
        axis vis3d
        hold on
        for k = 1:length(electrode.insulation)
            patch('Faces',electrode.insulation(k).faces,'Vertices',electrode.insulation(k).vertices,'Edgecolor','none','Facecolor',[electrode.lead_color electrode.lead_color electrode.lead_color]);
        end
        for k = 1:length(electrode.contacts)
            patch('Faces',electrode.contacts(k).faces,'Vertices',electrode.contacts(k).vertices,'Edgecolor','none','Facecolor',[electrode.contact_color electrode.contact_color electrode.contact_color]);
        end

        view(180,0)
        ylim([-1 1])
        xlim([-1 1])
        zlim([0 15])
        axis off
        axis equal

        camorbit(-rad2deg(tempangle),0)
        tempvec = [0; 1; 0];
        temp3x3 = ea_orient_rollpitchyaw(-tempangle,0,0);
        tempvec = temp3x3 * tempvec;
        text(tempvec(1),tempvec(2),markercenter,'M','FontSize',32,'HorizontalAlignment','center','VerticalAlignment','middle');
        text(tempvec(1),tempvec(2),level1center,'1','FontSize',32,'HorizontalAlignment','center','VerticalAlignment','middle');
        text(tempvec(1),tempvec(2),level2center,'2','FontSize',32,'HorizontalAlignment','center','VerticalAlignment','middle');
        clear tempangle

        set(ax_elec,'Position',[-0.16 0.21 0.43 0.73])

        %% get results
        if round(sumintensity1(rollangles == roll1)) <= -200
            checkbox1 = set(fig(side).chk1,'Value',1);
        end
        if round(sumintensity2(rollangles == roll2)) <= -200
            checkbox1 = set(fig(side).chk2,'Value',1);
        end

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
                disp(['Using roll angle defined by stereotactic marker: ' num2str(rad2deg(roll)) ' deg'])
            elseif checkbox1 && ~checkbox2
                roll_y = roll1;
                disp(['Using corrected roll angle defined by directional level 1: ' num2str(rad2deg(roll1)) ' deg'])
            elseif ~checkbox1 && checkbox2
                roll_y = roll2;
                disp(['Using corrected roll angle defined by directional level 2: ' num2str(rad2deg(roll2)) ' deg'])
            elseif checkbox1 && checkbox2
                roll_y = mean([roll1 roll2]);
                disp(['Using mean corrected roll angle defined by both directional levels: ' num2str(rad2deg(mean([roll1 roll2]))) ' deg'])
            end

            %% calculate y
            [M,~,~,~] = ea_orient_rollpitchyaw(roll_y,pitch,yaw);
            y = M * [0;1;0];
            head = head_mm(1:3);

            %% transform y to native space and back
            y = head + y;
            y(4) = 1;
            if strcmp(CTname,[filesep 'postop_ct.nii']) || strcmp(CTname,[filesep 'postop_ct_resliced.nii'])
                % transform postop_ct_mm -> rpostop_ct_mm
                y = inv(tmat_reg2org) * y;
            elseif strcmp(CTname,[filesep 'rpostop_ct.nii'])
                y = y;
            end

            head = head_native(1:3)';
            tail = tail_native(1:3)';
            y = y(1:3)' - head;

            %% Calculate direction of x and y markers
            [xunitv, yunitv] = ea_calcxy(head, tail, y);

            y = head + yunitv * (options.elspec.lead_diameter / 2);
            x = head + xunitv * (options.elspec.lead_diameter / 2);

            reco.native.markers(side).y = y;
            reco.native.markers(side).x = x;

            %% for direct saving into manual reconstruction
            [coords,trajectory,markers]=ea_resolvecoords(reco.native.markers,options);
            ea_save_reconstruction(coords,trajectory,markers,options.elmodel,1,options)

            % %% for transfering to ea_manualreconstruction
            yunitv(3) = 0;
            roll_out = rad2deg(atan2(norm(cross([0 1 0],yunitv)),dot([0 1 0],yunitv)));
            if markers(side).y(1) > markers(side).head(1) % negative 90 points to right, positive 90 points to left
                roll_out = - roll_out;
            end
            disp(['Corrected roll angle roll = ' num2str(rad2deg(roll_y)) ' deg, has been converted to orientation angle = ' num2str(roll_out) ' for compatibility with ea_mancorupdatescene.'])
            %% methods dump:
            ea_methods(options,...
                ['Orientation of directional DBS leads was determined using the algorithm published by Dembek et al. 2019 as implemented in Lead-DBS software.'],...
                {'T.A. Dembek, M. Hoevels, A. Hellerbach, A. Horn, J.N. Petry-Schmelzer, J. Borggrefe, J. Wirths, H.S. Dafsari, M.T. Barbe, V. Visser-Vandewalle & H. Treuer (2019). Directional DBS leads show large deviations from their intended implantation orientation. Parkinsonism Relat Disord. 2019 Oct;67:117-121. doi: 10.1016/j.parkreldis.2019.08.017.'});
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
else  % check for electrode type and postoperative imaging
    msg = sprintf(['No Valid Directional Lead Selected!']);
    choice = questdlg(msg,'No Directional Lead!','Abort','Abort');
    roll_out = [];
end


function buttonPress(hObject,eventdata)
hObject.UserData = 1;
uiresume
