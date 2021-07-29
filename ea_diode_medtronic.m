function [roll_y,y,solution] = ea_diode_mdt(side,ct,head_mm,unitvector_mm,tmat_vx2mm,elspec,supervised)
%% constant variables
% colorscale for ct figures
cscale = [-50 100];
cscale2 = [-500 2000];
% diameter of the slices shown in visualizations in vox
extractradius = 30;
% sidenames
sides = {'right','left','3','4','5','6','7','8'};

%% define electrode properties
% z position of the centers of directional levels 1 and 2 and the marker relative to
% head position
level1centerRelative = elspec.contact_length + elspec.contact_spacing;
level2centerRelative = (elspec.contact_length + elspec.contact_spacing) * 2;
markercenterRelative = elspec.markerpos - elspec.tip_length*~elspec.tipiscontact - elspec.contact_length/2;
load(elspec.matfname);

%% load CTs
%% Recalculate trajvector to optimize position at center of artifacts
% this part recalculates the exact lead position by calculating the
% center of mass of each CT slice from the head to 10mm above it and
% then performing a linear regression through these centers of mass

% sampling from head to 10mm above it in .5mm steps
samplelength = 20;
samplingvector_mm = vertcat([head_mm(1):unitvector_mm(1)./2:head_mm(1)+ samplelength*unitvector_mm(1)],...
    [head_mm(2):unitvector_mm(2)./2:head_mm(2)+ samplelength*unitvector_mm(2)],...
    [head_mm(3):unitvector_mm(3)./2:head_mm(3)+ samplelength*unitvector_mm(3)],...
    ones(1, samplelength*2+1));
samplingvector_vx = round(tmat_vx2mm\samplingvector_mm);

newcentervector_vx = nan(size(samplingvector_vx));
% for each slice calculate center of mass, if more than one
% center of mass is found, choose the one nearest to the
% original lead position
for k = 1:size(samplingvector_vx,2)
    [tmp,bb]=ea_sample_slice(ct,'tra',extractradius,'vox',{samplingvector_vx(1:3,k)'},1);
    tmp = tmp > 2000;
    tmpcent = regionprops(tmp,'Centroid');
    if numel(tmpcent) == 1
        tmpcent = round(tmpcent.Centroid);
        tmpcent = [bb{1}(tmpcent(1)),bb{2}(tmpcent(2)),samplingvector_vx(3,k),1]';
        newcentervector_vx(:,k) = tmpcent;
    elseif numel(tmpcent) > 1
        [~,tmpind] = min(sum(abs(vertcat(tmpcent.Centroid) - [(size(tmp,1)+1)/2 (size(tmp,1)+1)/2]),2));
        tmpcent = tmpcent(tmpind);
        clear tmpind
        tmpcent = round(tmpcent.Centroid);
        tmpcent = [bb{1}(tmpcent(1)),bb{2}(tmpcent(2)),samplingvector_vx(3,k),1]';
        newcentervector_vx(:,k) = tmpcent;
    elseif numel(tmpcent) == 0
        newcentervector_vx(:,k) = nan;
    end
end

newcentervector_mm = tmat_vx2mm * newcentervector_vx;
if numel(find(isnan(newcentervector_mm))) > 0.5 * numel(newcentervector_mm)
    error('Something went wrong with interpolation of the Lead - maybe wrong sForm/qForm was chosen?')
end
% fit linear model to the centers of mass and recalculate head
% and unitvector
new = [0:.5:samplelength];
xmdl = fitlm(new,newcentervector_mm(1,:));
ymdl = fitlm(new,newcentervector_mm(2,:));
zmdl = fitlm(new,newcentervector_mm(3,:));

head_mm = [predict(xmdl,0),predict(ymdl,0),predict(zmdl,0),1]';
other_mm = [predict(xmdl,10),predict(ymdl,10),predict(zmdl,10),1]';
unitvector_mm = (other_mm - head_mm)/norm(other_mm - head_mm);
clear tmpcent newcentervector_vx newcentervector_mm samplingvector_vx samplingvector_mm new k other_mm


% calculate locations of markers and directional levels
tail_mm = head_mm + (6 * unitvector_mm);
marker_mm = head_mm + (markercenterRelative .* unitvector_mm);
dirlevel1_mm = head_mm + (level1centerRelative .* unitvector_mm);
dirlevel2_mm = head_mm + (level2centerRelative .* unitvector_mm);

% transform to vx
marker_vx = tmat_vx2mm\marker_mm;
dirlevel1_vx = tmat_vx2mm\dirlevel1_mm;
dirlevel2_vx = tmat_vx2mm\dirlevel2_mm;

% in DiODe v2 only one directional level is used, starting at the
% middle between both directional levels and being optimized later
dirlevelnew_mm = mean([dirlevel1_mm,dirlevel2_mm]')';
dirlevelnew_vx = round(tmat_vx2mm\dirlevelnew_mm);

% calculate lead yaw and pitch angles for correction at the end
yaw = asin(unitvector_mm(1));
pitch = asin(unitvector_mm(2)/cos(yaw));
solution.polar1 = rad2deg(atan2(norm(cross([0;0;1],unitvector_mm(1:3))),dot([0;0;1],unitvector_mm(1:3))));
solution.polar2 = -rad2deg(atan2(unitvector_mm(2),unitvector_mm(1))) + 90;

if rad2deg(abs(pitch)) > 40
    disp(['Warning: Pitch > 40 deg - Determining orientation might be inaccurate!'])
end
if rad2deg(abs(yaw)) > 40
    disp(['Warning: Yaw > 40 deg - Determining orientation might be inaccurate!'])
end
if solution.polar1 > 40 || solution.polar2 > 40
    disp(['Warning: Polar > 40 deg - Determining orientation might be inaccurate!'])
end

%% extract axial slices at the level of marker and directional electrodes
if ~supervised
    %     %% select slice for marker where peak-valley-difference in FFT-intensity profile is largest
    %     count = 1;
    %     checkslices = [-1:0.5:1]; % check neighboring slices for marker +/- 1mm in .5mm steps
    %     for k = checkslices
    %         checklocation_mm = marker_mm + (unitvector_mm * k);
    %         checklocation_vx = round(tmat_vx2mm\checklocation_mm);
    %         tmp=ea_sample_slice(ct,'tra',extractradius,'vox',{checklocation_vx(1:3)' + [0 0 k]},1)';
    %         if ct.mat(1,1) < 0
    %             tmp = flip(tmp,1);
    %         end
    %         if ct.mat(2,2) < 0
    %             tmp = flip(tmp,2);
    %         end
    %         center_tmp = [(size(tmp,1)+1)/2 (size(tmp,1)+1)/2];
    %         radius = 4;
    %
    %         % calculate intensityprofile and its FFT for each slice
    %         [~, intensity,~] = ea_diode_intensityprofile(tmp,center_tmp,ct.voxsize,radius);
    %         [peak,tmpfft] = ea_diode_intensitypeaksFFT(intensity,2);
    %         valley = ea_diode_intensitypeaksFFT(-intensity,2);
    %         fftdiff(count) = mean(tmpfft(peak)) - mean(tmpfft(valley));
    %         count = count +1;
    %     end
    %
    %     % select slice with maximum difference in fft and respecify
    %     % marker accordingly
    %     [~,tmp_shift] = max(fftdiff);
    %     tmp_shift = checkslices(tmp_shift);
    %     marker_mm = marker_mm + (unitvector_mm * tmp_shift);
    %     marker_vx = tmat_vx2mm\marker_mm;
    %
    %     clear k count fftdiff valley peak radius intensity center_tmp tmp checkslices tmp_shift tmpfft checklocation_mm checklocation_vx
    %
    %     %% extract marker artifact from slice
    %     artifact_marker=ea_sample_slice(ct,'tra',extractradius,'vox',{round(marker_vx(1:3))'},1)';
    %     if ct.mat(1,1) < 0
    %         artifact_marker = flip(artifact_marker,1);
    %     end
    %     if ct.mat(2,2) < 0
    %         artifact_marker = flip(artifact_marker,2);
    %     end
    %     center_marker = [(size(artifact_marker,1)+1)/2 (size(artifact_marker,1)+1)/2];
    %
    %     %% extract intensity profile from marker artifact
    %     radius = 4;
    %     [angle, intensity,vector] = ea_diode_intensityprofile(artifact_marker,center_marker,ct.voxsize,radius);
    %
    %     %% detect peaks and valleys for marker artifact
    %     [peak,markerfft] = ea_diode_intensitypeaksFFT(intensity,2);
    %     valley = ea_diode_intensitypeaksFFT(-intensity,2);
    %
    %     %% Detect angles of the white streak of the marker (only for intensityprofile-based ambiguity features)
    %     valley_roll = ea_diode_angle2roll(angle(valley(1)),yaw,pitch);
    %     marker_angles = ea_diode_lightmarker(valley_roll,pitch,yaw,marker_mm);
    %
    %     solution.peaks = peak;
    %     solution.rolls_rad = [ea_diode_angle2roll(angle(peak(1)),yaw,pitch),ea_diode_angle2roll(angle(peak(2)),yaw,pitch)];
    %     solution.rolls_deg = rad2deg(solution.rolls_rad);
    %     solution.rolls_streak_deg = rad2deg(marker_angles);
    %
    %% Different methods to solve ambivalence for the marker
    % This is the main new feature of DiODe v2. A number of
    % different algorithms have been investigated. The first try to
    % resolve ambiguity due to different aspects of the intensity
    % profile of the marker - unfortunately these methods do not
    % function as well.
    % Next, two methods try to resolve ambiguity by calculating the
    % center of gravity (COG) of resampled, perpendicular slices exactly
    % through the marker and then investigate the direction of
    % shift from the expected marker center.
    % The last method is based on the asymmetry of the 3-line
    % "Darkstar" artifact generated by the directional level. It
    % resolves the Darkstar for both solutions and then compares
    % the goodness of fit.
    
    %% Intensity at FFT peak
    % compares the intensity at the location of the FFTpeak within
    % the intensity profile
    %             if intensity(peak(1)) > intensity(peak(2))
    %                 disp(['Intensity at FFT-Peak decides for peak 1'])
    %                 solution.FFTpeak = 1;
    %             else
    %                 disp(['Intensity at FFT-Peak decides for peak 2'])
    %                 solution.FFTpeak = 2;
    %             end
    
    %% Intensity at yaw/pitch corrected expected peak
    % compares the intensity at the location of the yaw/pitch
    % corrected peak within the intensity profile
    %             [~,corr_peak(1)] = min(abs(angle -marker_angles(1)));
    %             [~,corr_peak(2)] = min(abs(angle -marker_angles(2)));
    %             if intensity(corr_peak(1)) > intensity(corr_peak(2))
    %                 disp(['Intensity at corrected Peak decides for peak 1'])
    %                 solution.Correctedpeak = 1;
    %             else
    %                 disp(['Intensity at corrected Peak decides for peak 2'])
    %                 solution.Correctedpeak = 2;
    %             end
    %
    %% Maximum intensity in range
    % compares the maximum intensity between the valleys
    %             if max(intensity(valley(1):valley(2))) > max(intensity([1:valley(1),valley(2):length(intensity)]))
    %                 if peak(1) > valley(1) && peak(1) < valley(2)
    %                     disp(['Maximum Intensity decides for peak 1'])
    %                     solution.MAXpeak = 1;
    %                 else
    %                     disp(['Maximum Intensity decides for peak 2'])
    %                     solution.MAXpeak = 2;
    %                 end
    %             else
    %                 if peak(1) > valley(1) && peak(1) < valley(2)
    %                     disp(['Maximum Intensity decides for peak 2'])
    %                     solution.MAXpeak = 2;
    %                 else
    %                     disp(['Maximum Intensity decides for peak 1'])
    %                     solution.MAXpeak = 1;
    %                 end
    %             end
    
    %     %% ASM
    %     % compares the maximum intensity between the valleys in 3 radii
    %     ASMradii = [3,6,9];
    %     for k = 1:length(ASMradii)
    %         [~, ASMintensity(k,:),~] = ea_diode_intensityprofile(artifact_marker,center_marker,ct.voxsize,ASMradii(k));
    %     end
    %     ASMintensity = mean(ASMintensity);
    %     if max(ASMintensity(valley(1):valley(2))) > max(ASMintensity([1:valley(1),valley(2):length(ASMintensity)]))
    %         if peak(1) > valley(1) && peak(1) < valley(2)
    %             disp(['ASM decides for peak 1'])
    %             solution.ASM = 1;
    %         else
    %             disp(['ASM decides for peak 2'])
    %             solution.ASM = 2;
    %         end
    %     else
    %         if peak(1) > valley(1) && peak(1) < valley(2)
    %             disp(['ASM decides for peak 2'])
    %             solution.ASM = 2;
    %         else
    %             disp(['ASM decides for peak 1'])
    %             solution.ASM = 1;
    %         end
    %     end
    %
    %% Peak Shift from FFT
    % compares the angle-shift of the max intensity peak from the
    % peak in FFT and chooses the solution with the smaller shift
    %             tmp = [valley(1):valley(2)];
    %             [~,intpeak(1)] = max(intensity(tmp));
    %             intpeak(1) = tmp(intpeak(1));
    %             tmp = [1:valley(1),valley(2):length(intensity)];
    %             [~,intpeak(2)] = max(intensity(tmp));
    %             intpeak(2) = tmp(intpeak(2));
    %             intpeak = sort(intpeak,'ascend');
    %             clear tmp
    %             if max(abs(peak-intpeak)) > 90 && intpeak(1) < 90
    %                 intpeak(1) = intpeak(1) + 180;
    %             elseif max(abs(peak-intpeak)) > 90 && intpeak(2) > 270
    %                 intpeak(2) = intpeak(2) - 180;
    %             end
    %             intpeak = sort(intpeak,'ascend');
    %             if abs(peak(1) - intpeak(1)) < abs(peak(2) - intpeak(2))
    %                 disp(['Intensity peak shift from FFT decides for peak 1'])
    %                 solution.FFTshift = 1;
    %             else
    %                 disp(['Intensity peak shift from FFT decides for peak 2'])
    %                 solution.FFTshift = 2;
    %             end
    
    %% Peak Shift from Corrected Peak
    % compares the angle-shift of the max intensity peak from the
    % yaw- and shift-corrected peak and chooses the solution with
    % the smaller shift
    %             tmp = [valley(1):valley(2)];
    %             [~,intpeak(1)] = max(intensity(tmp));
    %             intpeak(1) = tmp(intpeak(1));
    %             tmp = [1:valley(1),valley(2):length(intensity)];
    %             [~,intpeak(2)] = max(intensity(tmp));
    %             intpeak(2) = tmp(intpeak(2));
    %             intpeak = sort(intpeak,'ascend');
    %             clear tmp
    %             if max(abs(corr_peak-intpeak)) > 90 && intpeak(1) < 90
    %                 intpeak(1) = intpeak(1) + 180;
    %             elseif max(abs(corr_peak-intpeak)) > 90 && intpeak(2) > 270
    %                 intpeak(2) = intpeak(2) - 180;
    %             end
    %             intpeak = sort(intpeak,'ascend');
    %             if abs(corr_peak(1) - intpeak(1)) < abs(corr_peak(2) - intpeak(2))
    %                 disp(['Intensity peak shift from corrected peak decides for peak 1'])
    %                 solution.Correctedshift = 1;
    %             else
    %                 disp(['Intensity peak shift from corrected decides for peak 2'])
    %                 solution.Correctedshift = 2;
    %             end
    
    %% Center of Mass method
    % this is where shit gets complicated
    
    % first, two orthogonal vectors, yvec which is the unitvector
    % pointing in the direction of 0 and x_vec, perpendicular
    % to it and unitvector are generated
%     rolltmp = ea_diode_angle2roll(0,yaw,pitch);
    
    rolltmp = 0;
    [M,~,~,~] = ea_diode_rollpitchyaw(rolltmp,pitch,yaw);
    yvec_mm = M * [0;1;0];
    xvec_mm = cross(unitvector_mm(1:3), yvec_mm);
    clear M
    
    %% create meshgrid for CT
    % coordinates in meshgrid format are created for the full
    % ct.img and a permuted ct is exported as Vnew (needed due to
    % the weird meshgrid format in Matlab)
    mincorner_mm = tmat_vx2mm * [1;1;1;1];
    maxcorner_mm = tmat_vx2mm * [size(ct.img)';1];
    
    [Xmm,Ymm,Zmm] = meshgrid([mincorner_mm(1):(maxcorner_mm(1)-mincorner_mm(1))./(size(ct.img,1)-1):maxcorner_mm(1)],...
        [mincorner_mm(2):(maxcorner_mm(2)-mincorner_mm(2))./(size(ct.img,2)-1):maxcorner_mm(2)],...
        [mincorner_mm(3):(maxcorner_mm(3)-mincorner_mm(3))./(size(ct.img,3)-1):maxcorner_mm(3)]);
    
    clear mincorner_mm maxcorner_mm
    
    Vnew = permute(ct.img,[2,1,3]);
    %% slice in CT
    % meshgrid based code to export a slice from the CT in Vnew
    % not needed any more
    
    %             extract_width = 10; % in mm
    %             samplingres = .1;
    %             [Xslice,Yslice] = meshgrid(...
    %                 [marker_mm(1)-extract_width:samplingres:marker_mm(1)+extract_width],...
    %                 [marker_mm(2)-extract_width:samplingres:marker_mm(2)+extract_width]);
    %             Zslice = repmat(marker_mm(3),size(Xslice));
    %
    %             figure
    %             newmarkerslice = slice(Xmm,Ymm,Zmm,Vnew,Xslice,Yslice,Zslice);
    %             set(newmarkerslice, 'EdgeColor','none')
    %             axis equal
    %             close
    
    %% slice perpendicular through the two markers
    % a 5mm slice with .1mm resolution is sampled perpendicular to
    % the lead over the range of the two markers, first the lower, than the
    % upper. COG is calculated for each slice at .2mm resolution. slices
    % with maximum COG distance to the geometric centers are calculated.
    % direction of deviation for these slices with respect to anterior
    % orientation are calculated.
    extract_width = 5; % in mm
    samplingres = .1;
    
    %% lower marker
    count = 1;
    for x = [-3:2*samplingres:0]
        center(count,:) = marker_mm + (x.*unitvector_mm);
        Xslice = ([-extract_width:samplingres:extract_width] .* xvec_mm(1)) + ([-extract_width:samplingres:extract_width] .* yvec_mm(1))' + center(count,1);
        Yslice = ([-extract_width:samplingres:extract_width] .* xvec_mm(2)) + ([-extract_width:samplingres:extract_width] .* yvec_mm(2))' + center(count,2);
        Zslice = ea_diode_perpendicularplane(unitvector_mm,center(count,:),Xslice,Yslice);
        myslice = interp3(Xmm,Ymm,Zmm,Vnew,Xslice,Yslice,Zslice);
        COG_mm(count,:) = ea_diode_calculateCOG((myslice >= 2000),Xslice,Yslice,Zslice);
        COG_dir(count,:) = (COG_mm(count,1:3)-center(count,1:3))/norm((COG_mm(count,1:3)-center(count,1:3)));
        distance(count) = norm(COG_mm(count,1:3)-center(count,1:3));        
        %% slice visualization if needed
%         figure
%         newmarkerslice = slice(Xmm,Ymm,Zmm,Vnew,Xslice,Yslice,Zslice);
%         set(newmarkerslice, 'EdgeColor','none')
%         hold on
%         scatter3(center(count,1),center(count,2),center(count,3),'k')
%         scatter3(center(count,1)+yvec_mm(1),center(count,2)+yvec_mm(2),center(count,3)+yvec_mm(3),'r')
%         plot3([center(count,1), center(count,1)+yvec_mm(1)],[center(count,2), center(count,2)+yvec_mm(2)],[center(count,3), center(count,3)+yvec_mm(3)],'r')
%         plot3([center(count,1), center(count,1)+unitvector_mm(1)],[center(count,2), center(count,2)+unitvector_mm(2)],[center(count,3), center(count,3)+unitvector_mm(3)],'g')
%         axis equal
%         scatter3(COG_mm(count,1),COG_mm(count,2),COG_mm(count,3),'g')
%         scatter3(marker_mm(1)+COG_dir(1),marker_mm(2)+COG_dir(2),marker_mm(3)+COG_dir(3),'b')
%         caxis([-500 3500])
%         view(0,90);
%         close
        %%
        count = count+1;
    end
    
    [~,bestslice] = max(distance);
    finalcenter{1} = center(bestslice,:);
    finalCOG{1} = COG_mm(bestslice,:);
    finalCOG_dir{1} = COG_dir(bestslice,:);
%     plot(distance)
%     close
    clear bestslice center COG_mm COG_dir distance count x
    %% upper marker
    count = 1;
    for x = [0:2*samplingres:3]
        center(count,:) = marker_mm + (x.*unitvector_mm);
        Xslice = ([-extract_width:samplingres:extract_width] .* xvec_mm(1)) + ([-extract_width:samplingres:extract_width] .* yvec_mm(1))' + center(count,1);
        Yslice = ([-extract_width:samplingres:extract_width] .* xvec_mm(2)) + ([-extract_width:samplingres:extract_width] .* yvec_mm(2))' + center(count,2);
        Zslice = ea_diode_perpendicularplane(unitvector_mm,center(count,:),Xslice,Yslice);
        myslice = interp3(Xmm,Ymm,Zmm,Vnew,Xslice,Yslice,Zslice);
        COG_mm(count,:) = ea_diode_calculateCOG((myslice >= 2000),Xslice,Yslice,Zslice);
        COG_dir(count,:) = (COG_mm(count,1:3)-center(count,1:3))/norm((COG_mm(count,1:3)-center(count,1:3)));
        distance(count) = norm(COG_mm(count,1:3)-center(count,1:3));
        
        %         COG_dir = (COG_mm-marker_mm(1:3))/norm((COG_mm-marker_mm(1:3)));
        %% slice visualization if needed
%         figure
%         newmarkerslice = slice(Xmm,Ymm,Zmm,Vnew,Xslice,Yslice,Zslice);
%         set(newmarkerslice, 'EdgeColor','none')
%         hold on
%         scatter3(center(count,1),center(count,2),center(count,3),'k')
%         scatter3(center(count,1)+yvec_mm(1),center(count,2)+yvec_mm(2),center(count,3)+yvec_mm(3),'r')
%         plot3([center(count,1), center(count,1)+yvec_mm(1)],[center(count,2), center(count,2)+yvec_mm(2)],[center(count,3), center(count,3)+yvec_mm(3)],'r')
%         plot3([center(count,1), center(count,1)+unitvector_mm(1)],[center(count,2), center(count,2)+unitvector_mm(2)],[center(count,3), center(count,3)+unitvector_mm(3)],'g')
%         axis equal
%         scatter3(COG_mm(count,1),COG_mm(count,2),COG_mm(count,3),'g')
%         scatter3(marker_mm(1)+COG_dir(1),marker_mm(2)+COG_dir(2),marker_mm(3)+COG_dir(3),'b')
%         caxis([-500 3500])
%         view(0,90);
%         close
        %%
        count = count+1;
    end
    
    [~,bestslice] = max(distance);
    finalcenter{2} = center(bestslice,:);
    finalCOG{2} = COG_mm(bestslice,:);
    finalCOG_dir{2} = COG_dir(bestslice,:);
%     plot(distance)
%     close
    clear bestslice center COG_mm COG_dir distance count x
    
    %%    
    finalangle(1) = rad2deg(atan2(norm(cross(finalCOG_dir{1},yvec_mm)),dot(finalCOG_dir{1},yvec_mm)));
    finalangle(2) = rad2deg(atan2(norm(cross(finalCOG_dir{2},yvec_mm)),dot(finalCOG_dir{2},yvec_mm)));
    diffangle = abs(diff(finalangle));
    
    
    if finalangle(2) <0 || finalangle(2) > 180
        keyboard
    end
    
    if finalCOG_dir{2}(1) > yvec_mm(1)
       finalangle(2) = -finalangle(2);
       finalangle(1) = finalangle(2) + diffangle;       
    end
    
    roll = deg2rad(mean(finalangle) -60);
%     if myroll < 0
%         myroll = 2*pi -myroll;
%     end
    
    
   
    
    
        
    %% slice perpendicular * 11
    % 10 times a 5mm slice with .1mm resolution is sampled perpendicular to
    % the lead from .5mm below to .5mm above the position of the
    % marker center oriented in the direction of x-vec and y-vec
    %             extract_width = 5; % in mm
    %             samplingres = .1;
    %             count = 1;
    %             COG_mm = [];
    %             COG_dir = [];
    %             for k = -0.5:0.1:0.5
    %                 markertmp = marker_mm +  (k .* unitvector_mm);
    %                 Xslice = ([-extract_width:samplingres:extract_width] .* xvec_mm(1)) + ([-extract_width:samplingres:extract_width] .* yvec_mm(1))' + markertmp(1);
    %                 Yslice = ([-extract_width:samplingres:extract_width] .* xvec_mm(2)) + ([-extract_width:samplingres:extract_width] .* yvec_mm(2))' + markertmp(2);
    %                 Zslice = ea_diode_perpendicularplane(unitvector_mm,marker_mm,Xslice,Yslice);
    %
    %                 myslice = interp3(Xmm,Ymm,Zmm,Vnew,Xslice,Yslice,Zslice);
    %                 COG_mm(:,count) = ea_diode_calculateCOG((myslice >= 2000),Xslice,Yslice,Zslice);
    %                 COG_dir(:,count) = (COG_mm(:,count)-markertmp(1:3))/norm((COG_mm(:,count)-markertmp(1:3)));
    %                 count = count +1;
    %             end
    %
    %             if sum(abs(yvec_mm-COG_dir),'all') < sum(abs(-yvec_mm-COG_dir),'all')
    %                 disp(['COGtrans11 decides for peak 1'])
    %                 solution.COGtrans11 = 1;
    %             else
    %                 disp(['COGtrans11 decides for peak 2'])
    %                 solution.COGtrans11 = 2;
    %             end
    %% slice parralel
    % a 1.5mm slice with .1mm resolution is sampled vertically
    % through the lead and through the marker center and oriented
    % in the direction of y-vec and unitvector
%     extract_width = 1.5; % in mm
%     samplingres = .1;
%     Xslice = ([-extract_width:samplingres:extract_width] .* unitvector_mm(1)) + ([-extract_width:samplingres:extract_width] .* yvec_mm(1))' + marker_mm(1);
%     Yslice = ([-extract_width:samplingres:extract_width] .* unitvector_mm(2)) + ([-extract_width:samplingres:extract_width] .* yvec_mm(2))' + marker_mm(2);
%     Zslice = ea_diode_perpendicularplane(xvec_mm,marker_mm,Xslice,Yslice);
%     
%     myslice = interp3(Xmm,Ymm,Zmm,Vnew,Xslice,Yslice,Zslice);
%     COG_mm = ea_diode_calculateCOG((myslice >= 2000),Xslice,Yslice,Zslice);
%     COG_dir = (COG_mm-marker_mm(1:3))/norm((COG_mm-marker_mm(1:3)));
%     
%     if sum(abs(yvec_mm-COG_dir)) < sum(abs(-yvec_mm-COG_dir))
%         disp(['COGsag decides for peak 1'])
%         solution.COGsag = 1;
%     else
%         disp(['COGsag decides for peak 2'])
%         solution.COGsag = 2;
%     end
    
%     %% slice visualization if needed
%                 figure
%                 newmarkerslice = slice(Xmm,Ymm,Zmm,Vnew,Xslice,Yslice,Zslice);
%                 set(newmarkerslice, 'EdgeColor','none')
%                 hold on
%                 scatter3(marker_mm(1),marker_mm(2),marker_mm(3),'k')
%                 scatter3(marker_mm(1)+yvec_mm(1),marker_mm(2)+yvec_mm(2),marker_mm(3)+yvec_mm(3),'r')
%                 plot3([marker_mm(1), marker_mm(1)+yvec_mm(1)],[marker_mm(2), marker_mm(2)+yvec_mm(2)],[marker_mm(3), marker_mm(3)+yvec_mm(3)],'r')
%                 plot3([marker_mm(1), marker_mm(1)+unitvector_mm(1)],[marker_mm(2), marker_mm(2)+unitvector_mm(2)],[marker_mm(3), marker_mm(3)+unitvector_mm(3)],'g')
%                 axis equal
%                 scatter3(COG_mm(1),COG_mm(2),COG_mm(3),'g')
%                 scatter3(marker_mm(1)+COG_dir(1),marker_mm(2)+COG_dir(2),marker_mm(3)+COG_dir(3),'b')
%                 caxis([-500 3500])
%                 close
    %% Slice parralel for visualization
    % a 10mm slice with .1mm resolution is sampled vertically
    % through the lead and through the marker center and oriented
    % in the direction of y-vec and unitvector for later
    % visualization    
    [M,~,~,~] = ea_diode_rollpitchyaw(roll,pitch,yaw);
    yvec_mm = M * [0;1;0];
    xvec_mm = cross(unitvector_mm(1:3), yvec_mm);
    clear M
    extract_width = 10; % in mm
    samplingres = .1;
    Xslice = ([-extract_width:samplingres:extract_width] .* unitvector_mm(1)) + ([-extract_width:samplingres:extract_width] .* yvec_mm(1))' + head_mm(1) + 7.5 * unitvector_mm(1);
    Yslice = ([-extract_width:samplingres:extract_width] .* unitvector_mm(2)) + ([-extract_width:samplingres:extract_width] .* yvec_mm(2))' + head_mm(2) + 7.5 * unitvector_mm(2);
    Zslice = ea_diode_perpendicularplane(xvec_mm,marker_mm,Xslice,Yslice);
    finalslice = interp3(Xmm,Ymm,Zmm,Vnew,Xslice,Yslice,Zslice);
    finalslice = finalslice';
    
    if rad2deg(roll) < 90 && rad2deg(roll) > -90
        finalslice = flipdim(finalslice,2);
    end
    %% darkstar method
    checkslices = [-2:0.5:2]; % check neighboring slices for marker
    
    % solution 1
    count = 1;
    myroll = roll;
    for x = checkslices
        checklocation_mm = dirlevelnew_mm + (unitvector_mm * x);
        checklocation_vx = round(tmat_vx2mm\checklocation_mm);
        artifact_tmp=ea_sample_slice(ct,'tra',extractradius,'vox',{checklocation_vx(1:3)'},1)';
        if ct.mat(1,1) < 0
            artifact_tmp = flip(artifact_tmp,1);
        end
        if ct.mat(2,2) < 0
            artifact_tmp = flip(artifact_tmp,2);
        end
        center_tmp = [(size(artifact_tmp,1)+1)/2 (size(artifact_tmp,1)+1)/2];
        radius = 8;
        
        [~, intensity_tmp,~] = ea_diode_intensityprofile(artifact_tmp,center_tmp,ct.voxsize,radius);
        %% determine angles of the 6-valley artifact ('dark star') artifact in each of the slices for +30:-30 deg
        for k = 1:61
            roll_shift = k-31;
            rolltemp = myroll + deg2rad(roll_shift);
            dirnew_angles = ea_diode_darkstar(rolltemp,pitch,yaw,checklocation_mm,radius);
            [sumintensitynew{1}(count,k)] = ea_diode_intensitypeaksdirmarker(intensity_tmp,dirnew_angles);
            %                     rollangles{1}(count,k) = rolltemp;
            rollangles{1}(count,k) = deg2rad(roll_shift);
        end
        count = count +1;
    end
    [~,darkstarangle(1)] = min(min(sumintensitynew{1},[],1));
    [~,darkstarslice(1)] = min(min(sumintensitynew{1},[],2));
    
    clear myroll roll_shift rolltemp dirnew_angles count
    
    % solution 2
    count = 1;
    myroll = roll + pi;
    for x = checkslices
        checklocation_mm = dirlevelnew_mm + (unitvector_mm * x);
        checklocation_vx = round(tmat_vx2mm\checklocation_mm);
        artifact_tmp=ea_sample_slice(ct,'tra',extractradius,'vox',{checklocation_vx(1:3)'},1)';
        if ct.mat(1,1) < 0
            artifact_tmp = flip(artifact_tmp,1);
        end
        if ct.mat(2,2) < 0
            artifact_tmp = flip(artifact_tmp,2);
        end
        center_tmp = [(size(artifact_tmp,1)+1)/2 (size(artifact_tmp,1)+1)/2];
        radius = 8;
        
        [~, intensity_tmp,~] = ea_diode_intensityprofile(artifact_tmp,center_tmp,ct.voxsize,radius);
        %% determine angles of the 6-valley artifact ('dark star') artifact in each of the slices for +30:-30 deg
        for k = 1:61
            roll_shift = k-31;
            rolltemp = myroll + pi + deg2rad(roll_shift);
            dirnew_angles = ea_diode_darkstar(rolltemp,pitch,yaw,checklocation_mm,radius);
            [sumintensitynew{2}(count,k)] = ea_diode_intensitypeaksdirmarker(intensity_tmp,dirnew_angles);
            %                     rollangles{2}(count,k) = rolltemp;
            rollangles{2}(count,k) = deg2rad(roll_shift);
        end
        count = count +1;
    end
    [~,darkstarangle(2)] = min(min(sumintensitynew{2},[],1));
    [~,darkstarslice(2)] = min(min(sumintensitynew{2},[],2));
    
    clear myroll roll_shift rolltemp dirnew_angles count
    
    % choose best slices and darkstar solution according to minimum
    % of sum intensity profile
    
    for k = 1:2
        sumintensitynew{k} = sumintensitynew{k}(darkstarslice(k),:);
        rollangles{k} = rollangles{k}(darkstarslice(k),:);
    end
    
    if min(sumintensitynew{1}(:)) < min(sumintensitynew{2}(:))
        disp(['Darkstar decides for peak 1'])
        solution.Darkstar = 1;
    else
        disp(['Darkstar decides for peak 2'])
        solution.Darkstar = 2;
    end
    
    %% Take anterior peak at the moment
%     if rad2deg(angle(peak(1))) <= 90 || rad2deg(angle(peak(1))) >= 270
%         finalpeak(side) = peak(1);
%         disp(['More anterior is peak 1'])
%         solution.Anterior = 1;
%     else
%         finalpeak(side) = peak(2);
%         disp(['More anterior is peak 2'])
%         solution.Anterior = 2;
%     end
%     
%     peakangle(side) = angle(finalpeak(side));
%     roll = ea_diode_angle2roll(peakangle(side),yaw,pitch);
% 
%     roll = myroll;
    peakangle(side) = rad2deg(roll);
    realsolution = 1;
    
    dirlevelnew_mm = dirlevelnew_mm + (unitvector_mm * checkslices(darkstarslice(realsolution)));
    dirlevelnew_vx = round(tmat_vx2mm\dirlevelnew_mm);
    
    artifact_dirnew = ea_sample_slice(ct,'tra',extractradius,'vox',{dirlevelnew_vx(1:3)'},1)';
    if ct.mat(1,1) < 0
        artifact_dirnew = flip(artifact_dirnew,1);
    end
    if ct.mat(2,2) < 0
        artifact_dirnew = flip(artifact_dirnew,2);
    end
    center_dirnew = [(size(artifact_dirnew,1)+1)/2 (size(artifact_dirnew,1)+1)/2];
    [anglenew, intensitynew,vectornew] = ea_diode_intensityprofile(artifact_dirnew,center_dirnew,ct.voxsize,radius);
    
    rollnew = roll + rollangles{realsolution}(darkstarangle(realsolution));
    dirnew_angles = ea_diode_darkstar(rollnew,pitch,yaw,dirlevelnew_mm,radius);
    dirnew_valleys = round(rad2deg(dirnew_angles) +1);
    dirnew_valleys(dirnew_valleys > 360) = dirnew_valleys(dirnew_valleys > 360) - 360;
    
    %             sumintensitynew = sumintensitynew(darkstarslice,:);
    %             rollangles = rollangles(darkstarslice,:);
end



%% final figure
fig(side).figure = figure('Name',['Lead ' sides{side}],'Position',[100 100 800 800],'Color','w','Toolbar','none');

if peakangle(side) > pi
    tempangle = peakangle(side) - 2 * pi;
else
    tempangle = peakangle(side);
end
fig(side).txt1 = uicontrol('style','text','units','normalized','Position',[.8,.8,.3,.1],...
    'Background','w', 'FontSize',12,'HorizontalAlignment','left',...
    'string',sprintf(['Artifact Angle:\n' num2str(rad2deg(tempangle),'%.1f') ' deg\nMarker Angle:\n' num2str(rad2deg(roll),'%.1f') ' deg']));

fig(side).txt3 = uicontrol('style','text','units','normalized','Position',[.8,.46,.3,.1],...
    'Background','w','FontSize',12,'HorizontalAlignment','left',...
    'string',sprintf(['Dir-Level Shift:\n' num2str(rad2deg(rollangles{[1 2] == realsolution}(darkstarangle([1 2] == realsolution))),'%.1f') ' deg\n' 'Corrected Angle:\n' num2str(rad2deg(rollnew),'%.1f') ' deg']));

fig(side).chk1 = uicontrol('style','checkbox','units','normalized','Position',[.8,.41,.3,.05],...
    'string','Accept','FontSize',12,'Background','w');

% fig(side).txt6 = uicontrol('style','text','units','pixels','Background','w',...
%     'position',[100,250,720,20],'FontSize',12,'HorizontalAlignment','left',...
%     'string',sprintf([...
%     'COM-Transversal Solution is: ' num2str(round(solution.rolls_deg(solution.COGtrans),1)) ' deg\n' ...
%     ]));
% if solution.COGsag ~= solution.COGtrans
%     txtcolor = [1 .5 .25];
% else
%     txtcolor = 'k';
% end
% fig(side).txt7 = uicontrol('style','text','units','pixels','Background','w',...
%     'position',[100,230,720,20],'FontSize',12,'HorizontalAlignment','left','ForegroundColor',txtcolor,...
%     'string',sprintf([...
%     'COM-Sagittal Solution is: ' num2str(round(solution.rolls_deg(solution.COGsag),1)) ' deg\n' ...
%     ]));
% if solution.Darkstar ~= solution.COGtrans
%     txtcolor = [1 .5 .25];
% else
%     txtcolor = 'k';
% end
% fig(side).txt8 = uicontrol('style','text','units','pixels','Background','w',...
%     'position',[100,210,720,20],'FontSize',12,'HorizontalAlignment','left','ForegroundColor',txtcolor,...
%     'string',sprintf([...
%     'STARS Solution is: ' num2str(round(solution.rolls_deg(solution.Darkstar),1)) ' deg\n' ...
%     ]));
% if solution.ASM ~= solution.COGtrans
%     txtcolor = [1 .5 .25];
% else
%     txtcolor = 'k';
% end
% fig(side).txt9 = uicontrol('style','text','units','pixels','Background','w',...
%     'position',[100,190,720,20],'FontSize',12,'HorizontalAlignment','left','ForegroundColor',txtcolor,...
%     'string',sprintf([...
%     'ASM Solution is: ' num2str(round(solution.rolls_deg(solution.ASM),1)) ' deg\n' ...
%     ]));


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
% hold on
% title(ax1,'Marker')
% imagesc(artifact_marker')
% view(-180,-90)
% axis equal
% axis off
% colormap(gray)
% caxis manual
% caxis(cscale)
% plot(ax1,vector(:,1),vector(:,2),'g','LineStyle',':')
% scatter(ax1,vector(peak,1), vector(peak,2),'g');
% scatter(ax1,vector(finalpeak(side),1), vector(finalpeak(side),2),'g','filled');
% quiver(center_marker(1),center_marker(2),vector(finalpeak(side),1)-center_marker(1) ,vector(finalpeak(side),2)-center_marker(2),2,'LineWidth',1,'Color','g','MaxHeadSize',2)
% scatter(ax1,center_marker(1),center_marker(2),[],[0 0.4470 0.7410],'filled')
% for k = 1:length(valley)
%     plot(ax1,[center_marker(1) (center_marker(1) + 3 * (vector(valley(k),1)-center_marker(1)))],...
%         [center_marker(2) (center_marker(2) + 3 * (vector(valley(k),2)-center_marker(2)))],'r','LineStyle','--')
% end
% xlimit = get(ax1,'Xlim');
% ylimit = get(ax1,'Ylim');
% text(mean(xlimit),ylimit(2) - 0.15 * mean(ylimit),'A','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
% text(mean(xlimit),ylimit(1) + 0.15 * mean(ylimit),'P','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
% text(xlimit(1) + 0.1 * mean(xlimit),mean(ylimit),'L','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
% text(xlimit(2) - 0.1 * mean(xlimit),mean(ylimit),'R','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')

ax2 = subplot(3,3,2);
% hold on
% title(ax2,'Intensity Profile','FontWeight','normal')
% plot(ax2,rad2deg(angle),intensity)
% plot(ax2,rad2deg(angle),markerfft)
% scatter(ax2,rad2deg(angle(peak)), intensity(peak),'g');
% scatter(ax2,rad2deg(angle(finalpeak(side))), intensity(finalpeak(side)),'g','filled');
% scatter(ax2,rad2deg(angle(valley)), intensity(valley),'r');
% set(ax2,'yticklabel',{[]})

ax3 = subplot(3,3,3);
hold on
title(ax3,'Sagittal View','FontWeight','normal')

imagesc(finalslice)
axis equal
axis off
caxis([1500 3000])
if peakangle(side) < 90 || peakangle(side) > 270
    quiver(round(size(finalslice,2)/2), round(size(finalslice,1)/2).*1.7, -round(size(finalslice,1)/8), 0, 2,'LineWidth',1.5,'Color','g','MaxHeadSize',2)
elseif peakangle(side) >= 90 && peakangle(side) <=270
    quiver(round(size(finalslice,2)/2), round(size(finalslice,1)/2).*1.7, +round(size(finalslice,1)/8), 0, 2,'LineWidth',1.5,'Color','g','MaxHeadSize',2)
end
scatter(ax3,round(size(finalslice,2)/2),round(size(finalslice,1)/2).*1.7,[],[0 0.4470 0.7410],'filled')
plot([round(size(finalslice,2)/2), round(size(finalslice,2)/2)], [round(size(finalslice,2)/2)-75, round(size(finalslice,2)/2)+100],'LineStyle','--','Color',[0 0.4470 0.7410])
xlimit = get(ax3,'Xlim');
ylimit = get(ax3,'Ylim');
text(xlimit(1) + 0.1 * mean(xlimit),mean(ylimit),'A','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
text(xlimit(2) - 0.1 * mean(xlimit),mean(ylimit),'P','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
%% graphics dir level one
ax4 = subplot(3,3,4);
hold on
title(ax4,'Directional Level')
imagesc(artifact_dirnew')
view(-180,-90)
axis equal
axis off
colormap(gray)
caxis manual
caxis(cscale)
plot(ax4,vectornew(:,1),vectornew(:,2),'g','LineStyle',':')
scatter(ax4,vectornew(dirnew_valleys,1),vectornew(dirnew_valleys,2),'r')
for k = 1:length(dirnew_valleys)
    plot(ax4,[center_dirnew(1) (center_dirnew(1) + 1.5 * (vectornew(dirnew_valleys(k),1)-center_dirnew(1)))],...
        [center_dirnew(2) (center_dirnew(2) + 1.5 * (vectornew(dirnew_valleys(k),2)-center_dirnew(2)))],'r','LineStyle','--')
end
scatter(ax4,center_dirnew(1),center_dirnew(2),[],[0 0.4470 0.7410],'filled')
xlimit = get(ax4,'Xlim');
ylimit = get(ax4,'Ylim');
text(mean(xlimit),ylimit(2) - 0.15 * mean(ylimit),'A','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
text(mean(xlimit),ylimit(1) + 0.15 * mean(ylimit),'P','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
text(xlimit(1) + 0.1 * mean(xlimit),mean(ylimit),'L','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
text(xlimit(2) - 0.1 * mean(xlimit),mean(ylimit),'R','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')

ax5 = subplot(3,3,5);
hold on
title(ax5,'Intensity Profile','FontWeight','normal')
plot(ax5,rad2deg(anglenew),intensitynew)
scatter(ax5,rad2deg(anglenew(dirnew_valleys)), intensitynew(dirnew_valleys),'r');
set(ax5,'yticklabel',{[]})

ax6 = subplot(3,3,6);
hold on
title(ax6,'Similarity Index','FontWeight','normal')
plot(ax6,rad2deg(rollangles{[1 2] == realsolution}),sumintensitynew{[1 2] == realsolution})
plot(ax6,rad2deg(rollangles{[1 2] ~= realsolution}),sumintensitynew{[1 2] ~= realsolution},'r')
scatter(ax6,...
    rad2deg(rollangles{[1 2] == realsolution}(rollangles{[1 2] == realsolution} == 0)),...
    sumintensitynew{[1 2] == realsolution}(rollangles{[1 2] == realsolution} == 0),...
    'g','filled');
scatter(ax6,...
    rad2deg(rollangles{[1 2] == realsolution}(darkstarangle([1 2] == realsolution))),...
    sumintensitynew{[1 2] == realsolution}(darkstarangle([1 2] == realsolution)),...
    'r');
text(rad2deg(rollangles{[1 2] == realsolution}(darkstarangle([1 2] == realsolution))),...
    sumintensitynew{[1 2] == realsolution}(darkstarangle([1 2] == realsolution)),...
    ['\leftarrow HU = ' num2str(round(sumintensitynew{[1 2] == realsolution}(darkstarangle([1 2] == realsolution))))]);
set(ax6,'yticklabel',{[]})

%% Positioning of different sublots
linkaxes([ax2 ax5],'xy');
set(ax2,'Xlim',[0 360]);
set(ax2,'Ylim',[min([intensitynew])-50 max([intensitynew])+50]);

set(ax6,'Xlim',[rad2deg(rollangles{1}(1)) rad2deg(rollangles{1}(end))]);
set(ax6,'Ylim',[-1000 1000]);

set(ax1,'Position',[0.13 0.75 0.2 0.2])
set(ax2,'Position',[0.345 0.75 0.2 0.2])
set(ax3,'Position',[0.56 0.75 0.2 0.2])
set(ax4,'Position',[0.13 0.395 0.2 0.2])
set(ax5,'Position',[0.345 0.395 0.2 0.2])
set(ax6,'Position',[0.56 0.395 0.2 0.2])

%% graphics lead
ax_elec = axes('Position',[0 0.3 0.1 0.75]);
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
zlim([0 max(electrode.insulation(end-1).vertices(:,3))])
axis off
axis equal

camorbit(-rad2deg(tempangle),0)
tempvec = [0; 1; 0];
temp3x3 = ea_diode_rollpitchyaw(-tempangle,0,0);
tempvec = temp3x3 * tempvec;
%         text(tempvec(1),tempvec(2),markercenter,'M','FontSize',32,'HorizontalAlignment','center','VerticalAlignment','middle');
%         text(tempvec(1),tempvec(2),level1center,'1','FontSize',32,'HorizontalAlignment','center','VerticalAlignment','middle');
%         text(tempvec(1),tempvec(2),level2center,'2','FontSize',32,'HorizontalAlignment','center','VerticalAlignment','middle');
clear tempangle

set(ax_elec,'Position',[-0.16 0.38 0.43 0.6])

%% get results
if round(sumintensitynew{[1 2] == realsolution}(darkstarangle([1 2] == realsolution))) <= -200
    checkbox1 = set(fig(side).chk1,'Value',1);
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
    [roll_y_retry,y_retry] = ea_diode_manual(side,ct,head_mm,unitvector_mm,tmat_vx2mm,elspec);
end

%% saving results
if savestate == 1
    clear x y
    checkbox1 = get(fig(side).chk1,'Value');
    if ~checkbox1
        roll_y = roll;
        disp(['Using roll angle defined by stereotactic marker: ' num2str(rad2deg(roll)) ' deg'])
    elseif checkbox1
        roll_y = rollnew;
        disp(['Using corrected roll angle defined by directional level 1: ' num2str(rad2deg(rollnew)) ' deg'])
    end
    %% calculate y
    [M,~,~,~] = ea_diode_rollpitchyaw(roll_y,pitch,yaw);
    y = M * [0;1;0];
    head = head_mm(1:3);
    y = head + y;
    y(4) = 1;
elseif retrystate == 0
    disp(['Changes to rotation not saved'])
    roll_y = [];
    y = [];
elseif retrystate == 1
    disp(['Rotation saved after manual refinement'])
    roll_y = roll_y_retry;
    y = y_retry;
end
close(fig(side).figure)
end

function buttonPress(hObject,eventdata)
hObject.UserData = 1;
uiresume
end