function [roll_y,y,solution] = ea_diode_medtronic(side,ct,head_mm,unitvector_mm,tmat_vx2mm,elspec)
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
disp('Inspecting lower marker');
count = 1;
for x = [-4:2*samplingres:0]
    center(count,:) = marker_mm + (x.*unitvector_mm);
    Xslice = ([-extract_width:samplingres:extract_width] .* xvec_mm(1)) + ([-extract_width:samplingres:extract_width] .* yvec_mm(1))' + center(count,1);
    Yslice = ([-extract_width:samplingres:extract_width] .* xvec_mm(2)) + ([-extract_width:samplingres:extract_width] .* yvec_mm(2))' + center(count,2);
    Zslice = ea_diode_perpendicularplane(unitvector_mm,center(count,:),Xslice,Yslice);
    myslice = interp3(Xmm,Ymm,Zmm,Vnew,Xslice,Yslice,Zslice);
    COG_mm(count,:) = ea_diode_calculateCOG((myslice >= 2000),Xslice,Yslice,Zslice);
    COG_dir(count,:) = (COG_mm(count,1:3)-center(count,1:3))/norm((COG_mm(count,1:3)-center(count,1:3)));
    distance(count) = norm(COG_mm(count,1:3)-center(count,1:3));
    %% slice visualization if needed
    %             figure
    %             newmarkerslice = slice(Xmm,Ymm,Zmm,Vnew,Xslice,Yslice,Zslice);
    %             set(newmarkerslice, 'EdgeColor','none')
    %             hold on
    %             scatter3(center(count,1),center(count,2),center(count,3),'k')
    %             scatter3(center(count,1)+yvec_mm(1),center(count,2)+yvec_mm(2),center(count,3)+yvec_mm(3),'r')
    %             plot3([center(count,1), center(count,1)+yvec_mm(1)],[center(count,2), center(count,2)+yvec_mm(2)],[center(count,3), center(count,3)+yvec_mm(3)],'r')
    %             plot3([center(count,1), center(count,1)+unitvector_mm(1)],[center(count,2), center(count,2)+unitvector_mm(2)],[center(count,3), center(count,3)+unitvector_mm(3)],'g')
    %             axis equal
    %             scatter3(COG_mm(count,1),COG_mm(count,2),COG_mm(count,3),'g')
    %             scatter3(marker_mm(1)+COG_dir(1),marker_mm(2)+COG_dir(2),marker_mm(3)+COG_dir(3),'b')
    %             caxis([-500 3500])
    %             view(0,90);
    %             close
    %%
    count = count+1;
end

[~,bestslice] = max(distance);
finalcenter{2} = center(bestslice,:);
finalcenter_vx{2} = round(tmat_vx2mm\finalcenter{2}');
finalCOG{2} = COG_mm(bestslice,:);
finalCOG_vx{2} = round(tmat_vx2mm\[finalCOG{2} 1]');
finalCOG_dir{2} = COG_dir(bestslice,:);
figure
plot([-4:2*samplingres:0],distance)
close
clear bestslice center COG_mm COG_dir distance count x
%% upper marker
disp('Inspecting upper marker');
count = 1;
for x = [0:2*samplingres:4]
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
    %             figure
    %             newmarkerslice = slice(Xmm,Ymm,Zmm,Vnew,Xslice,Yslice,Zslice);
    %             set(newmarkerslice, 'EdgeColor','none')
    %             hold on
    %             scatter3(center(count,1),center(count,2),center(count,3),'k')
    %             scatter3(center(count,1)+yvec_mm(1),center(count,2)+yvec_mm(2),center(count,3)+yvec_mm(3),'r')
    %             plot3([center(count,1), center(count,1)+yvec_mm(1)],[center(count,2), center(count,2)+yvec_mm(2)],[center(count,3), center(count,3)+yvec_mm(3)],'r')
    %             plot3([center(count,1), center(count,1)+unitvector_mm(1)],[center(count,2), center(count,2)+unitvector_mm(2)],[center(count,3), center(count,3)+unitvector_mm(3)],'g')
    %             axis equal
    %             scatter3(COG_mm(count,1),COG_mm(count,2),COG_mm(count,3),'g')
    %             scatter3(marker_mm(1)+COG_dir(1),marker_mm(2)+COG_dir(2),marker_mm(3)+COG_dir(3),'b')
    %             caxis([-500 3500])
    %             view(0,90);
    %             close
    %%
    count = count+1;
end

[~,bestslice] = max(distance);
finalcenter{1} = center(bestslice,:);
finalcenter_vx{1} = round(tmat_vx2mm\finalcenter{1}');
finalCOG{1} = COG_mm(bestslice,:);
finalCOG_vx{1} = round(tmat_vx2mm\[finalCOG{1} 1]');
finalCOG_dir{1} = COG_dir(bestslice,:);
hold on
plot([0:2*samplingres:4],distance)
% close
clear bestslice center COG_mm COG_dir distance count x

%% loop over markers
for x=1:2
    %% extract artifacts for upper and lower marker
    finalartifact{x} = ea_sample_slice(ct,'tra',extractradius,'vox',{finalcenter_vx{x}(1:3)'},1)';
    if ct.mat(1,1) < 0
        finalartifact{x} = flip(finalartifact{x},1);
    end
    if ct.mat(2,2) < 0
        finalartifact{x} = flip(finalartifact{x},2);
    end
    %% extract intensity profile from marker artifact
    radius = 4;
    %     [angle{x}, intensity{x},vector{x}] = ea_diode_intensityprofile(finalartifact{x},...
    %         [(finalCOG_dir{x}(1) ./ ct.voxsize(1))+((size(finalartifact{x},1)+1)./2),(finalCOG_dir{x}(2) ./ ct.voxsize(2)) + ((size(finalartifact{x},2)+1)./2)],...
    %         ct.voxsize,radius);
    [angle{x}, intensity{x},vector{x}] = ea_diode_intensityprofile(finalartifact{x},...
        [((size(finalartifact{x},1)+1)./2),((size(finalartifact{x},2)+1)./2)],...
        ct.voxsize,radius);
    %% detect peaks and valleys for marker artifact
    [~,markerfft{x}] = ea_diode_intensitypeaksFFT(intensity{x},2);
    [valleyfft{x},~] = ea_diode_intensitypeaksFFT(-intensity{x},2);
    windowwidth = 20;
    for k = 1:length(valleyfft{x})
        tmpind = (valleyfft{x}(k)-windowwidth:valleyfft{x}(k)+windowwidth);
        tmpind(tmpind<1) = tmpind(tmpind<1) +360;
        tmpind(tmpind>360) = tmpind(tmpind>360) -360;
        [~,loctemp] = min(intensity{x}(tmpind));
        valleyreal{x}(k) = valleyfft{x}(k)-windowwidth+loctemp;
        clear loctemp tmpind
    end
    valleyreal{x}(valleyreal{x}<1) = valleyreal{x}(valleyreal{x}<1) +360;
    valleyreal{x}(valleyreal{x}>360) = valleyreal{x}(valleyreal{x}>360) -360;

    %% Detect angles of the white streak of the marker (only for intensityprofile-based ambiguity features)
    figure
    plot(rad2deg(angle{x}),intensity{x})
    hold on
    plot(rad2deg(angle{x}),markerfft{x})
    scatter(rad2deg(angle{x}(valleyreal{x}-1)), intensity{x}(valleyreal{x}-1),'r')
    close

    valleycalc{x} = [round(mean(valleyfft{x}))-90, round(mean(valleyfft{x}))+90];
    valley_roll(x) = ea_diode_angle2roll(angle{x}(valleycalc{x}(1)),yaw,pitch);
    marker_angles{x} = ea_diode_lightmarker(valley_roll(x),pitch,yaw,marker_mm);
end

%% Angles from COG_dir
finalangle(1) = rad2deg(atan2(norm(cross(finalCOG_dir{1},yvec_mm)),dot(finalCOG_dir{1},yvec_mm)));
finalangle(2) = rad2deg(atan2(norm(cross(finalCOG_dir{2},yvec_mm)),dot(finalCOG_dir{2},yvec_mm)));

if finalangle(1) <0 || finalangle(1) > 180
    keyboard
end
%% Protection against left<->right switches
if finalCOG_dir{1}(1) > yvec_mm(1) && finalangle(1) > 0
    finalangle(1) = -finalangle(1)+360;
end
if finalCOG_dir{2}(1) > yvec_mm(1) && finalangle(2) > 0
    finalangle(2) = -finalangle(2)+360;
end

%% Angles from intensityprofile
[~,tmpind] = min(abs(rad2deg(marker_angles{1})-finalangle(1)));
finalintensityangle(1) = rad2deg(marker_angles{1}(tmpind));
[~,tmpind] = min(abs(rad2deg(marker_angles{2})-finalangle(2)));
finalintensityangle(2) = rad2deg(marker_angles{2}(tmpind));
clear tmpind

if finalangle(2) < finalangle(1)
    finalangle(2) = finalangle(2) + 360;
end
if finalintensityangle(2) < finalintensityangle(1)
    finalintensityangle(2) = finalintensityangle(2) + 360;
end
% finalangle(finalangle > 180) = finalangle(finalangle > 180) - 360;
% finalintensityangle(finalintensityangle > 180) = finalintensityangle(finalintensityangle > 180) - 360;

diffangle = abs(diff(finalangle));
diffintensityangle = abs(diff(finalintensityangle));
roll = deg2rad(mean(finalangle) -60);
rollintensity = deg2rad(mean(finalintensityangle) -60);

if roll > pi
    roll = roll - 2*pi;
end
if rollintensity > pi
    rollintensity = rollintensity - 2*pi;
end

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
finalslice = flipdim(finalslice,2);

% if rad2deg(roll) < 90 && rad2deg(roll) > -90
%     finalslice = flipdim(finalslice,2);
% end
%% darkstar method
% checkslices = [-2:0.5:2]; % check neighboring slices for marker
% 
% % solution 1
% count = 1;
% myroll = roll;
% for x = checkslices
%     checklocation_mm = dirlevelnew_mm + (unitvector_mm * x);
%     checklocation_vx = round(tmat_vx2mm\checklocation_mm);
%     artifact_tmp=ea_sample_slice(ct,'tra',extractradius,'vox',{checklocation_vx(1:3)'},1)';
%     if ct.mat(1,1) < 0
%         artifact_tmp = flip(artifact_tmp,1);
%     end
%     if ct.mat(2,2) < 0
%         artifact_tmp = flip(artifact_tmp,2);
%     end
%     center_tmp = [(size(artifact_tmp,1)+1)/2 (size(artifact_tmp,1)+1)/2];
%     radius = 8;
% 
%     [~, intensity_tmp,~] = ea_diode_intensityprofile(artifact_tmp,center_tmp,ct.voxsize,radius);
%     %% determine angles of the 6-valley artifact ('dark star') artifact in each of the slices for +30:-30 deg
%     for k = 1:61
%         roll_shift = k-31;
%         rolltemp = myroll + deg2rad(roll_shift);
%         dirnew_angles = ea_diode_darkstar(rolltemp,pitch,yaw,checklocation_mm,radius);
%         [sumintensitynew{1}(count,k)] = ea_diode_intensitypeaksdirmarker(intensity_tmp,dirnew_angles);
%         %                     rollangles{1}(count,k) = rolltemp;
%         rollangles{1}(count,k) = deg2rad(roll_shift);
%     end
%     count = count +1;
% end
% [~,darkstarangle(1)] = min(min(sumintensitynew{1},[],1));
% [~,darkstarslice(1)] = min(min(sumintensitynew{1},[],2));
% 
% clear myroll roll_shift rolltemp dirnew_angles count
% 
% % solution 2
% count = 1;
% myroll = roll + pi;
% for x = checkslices
%     checklocation_mm = dirlevelnew_mm + (unitvector_mm * x);
%     checklocation_vx = round(tmat_vx2mm\checklocation_mm);
%     artifact_tmp=ea_sample_slice(ct,'tra',extractradius,'vox',{checklocation_vx(1:3)'},1)';
%     if ct.mat(1,1) < 0
%         artifact_tmp = flip(artifact_tmp,1);
%     end
%     if ct.mat(2,2) < 0
%         artifact_tmp = flip(artifact_tmp,2);
%     end
%     center_tmp = [(size(artifact_tmp,1)+1)/2 (size(artifact_tmp,1)+1)/2];
%     radius = 8;
% 
%     [~, intensity_tmp,~] = ea_diode_intensityprofile(artifact_tmp,center_tmp,ct.voxsize,radius);
%     %% determine angles of the 6-valley artifact ('dark star') artifact in each of the slices for +30:-30 deg
%     for k = 1:61
%         roll_shift = k-31;
%         rolltemp = myroll + pi + deg2rad(roll_shift);
%         dirnew_angles = ea_diode_darkstar(rolltemp,pitch,yaw,checklocation_mm,radius);
%         [sumintensitynew{2}(count,k)] = ea_diode_intensitypeaksdirmarker(intensity_tmp,dirnew_angles);
%         %                     rollangles{2}(count,k) = rolltemp;
%         rollangles{2}(count,k) = deg2rad(roll_shift);
%     end
%     count = count +1;
% end
% [~,darkstarangle(2)] = min(min(sumintensitynew{2},[],1));
% [~,darkstarslice(2)] = min(min(sumintensitynew{2},[],2));
% 
% clear myroll roll_shift rolltemp dirnew_angles count
% 
% % choose best slices and darkstar solution according to minimum
% % of sum intensity profile
% 
% for k = 1:2
%     sumintensitynew{k} = sumintensitynew{k}(darkstarslice(k),:);
%     rollangles{k} = rollangles{k}(darkstarslice(k),:);
% end
% 
% if min(sumintensitynew{1}(:)) < min(sumintensitynew{2}(:))
%     disp(['Darkstar decides for peak 1'])
%     solution.Darkstar = 1;
% else
%     disp(['Darkstar decides for peak 2'])
%     solution.Darkstar = 2;
% end
% 
% %% Take peak
% peakangle(side) = rad2deg(roll);
% realsolution = 1;
% 
% dirlevelnew_mm = dirlevelnew_mm + (unitvector_mm * checkslices(darkstarslice(realsolution)));
% dirlevelnew_vx = round(tmat_vx2mm\dirlevelnew_mm);
% 
% artifact_dirnew = ea_sample_slice(ct,'tra',extractradius,'vox',{dirlevelnew_vx(1:3)'},1)';
% if ct.mat(1,1) < 0
%     artifact_dirnew = flip(artifact_dirnew,1);
% end
% if ct.mat(2,2) < 0
%     artifact_dirnew = flip(artifact_dirnew,2);
% end
% center_dirnew = [(size(artifact_dirnew,1)+1)/2 (size(artifact_dirnew,1)+1)/2];
% [anglenew, intensitynew,vectornew] = ea_diode_intensityprofile(artifact_dirnew,center_dirnew,ct.voxsize,radius);
% 
% rollnew = roll + rollangles{realsolution}(darkstarangle(realsolution));
% dirnew_angles = ea_diode_darkstar(rollnew,pitch,yaw,dirlevelnew_mm,radius);
% dirnew_valleys = round(rad2deg(dirnew_angles) +1);
% dirnew_valleys(dirnew_valleys > 360) = dirnew_valleys(dirnew_valleys > 360) - 360;

%% final figure
fig(side).figure = figure('Name',['Lead ' sides{side}],'Position',[100 100 800 800],'Color','w','Toolbar','none');

fig(side).txt1 = uicontrol('style','text','units','normalized','Position',[.675,.75,.3,.15],...
    'Background','w', 'FontSize',12,'HorizontalAlignment','left',...
    'string',sprintf(['Upper Marker angle based on...\n...COG:\n' num2str(convertAngle(finalangle(1)),'%.1f') ' deg\n...Intensity Profiles:\n' num2str(convertAngle(finalintensityangle(1)),'%.1f') ' deg']));

fig(side).txt2 = uicontrol('style','text','units','normalized','Position',[.675,.5,.3,.15],...
    'Background','w', 'FontSize',12,'HorizontalAlignment','left',...
    'string',sprintf(['Lower Marker angle based on...\n...COG:\n' num2str(convertAngle(finalangle(2)),'%.1f') ' deg\n...Intensity Profiles:\n' num2str(convertAngle(finalintensityangle(2)),'%.1f') ' deg']));

if diffangle > 100 && diffangle < 140
    txtcolor = [34 177 75]./255;
    fig(side).txt3 = uicontrol('style','text','units','pixels','Background','w','ForegroundColor',txtcolor,...
    'position',[80,340,700,40],'FontSize',12,'HorizontalAlignment','left',...
    'string',sprintf([...
    'Angular difference between Upper and Lower Marker for COG is: ' num2str(diffangle,'%.1f') ' deg\nThis seems plausible!']));
else
    txtcolor = 'r';
    fig(side).txt3 = uicontrol('style','text','units','pixels','Background','w','ForegroundColor',txtcolor,...
    'position',[80,340,700,40],'FontSize',12,'HorizontalAlignment','left',...
    'string',sprintf([...
    'Angular difference between Upper and Lower Marker for COG is: ' num2str(diffangle,'%.1f') ' deg\nThis seems implausible!']));
end

if diffintensityangle > 100 && diffintensityangle < 140
    txtcolor = [34 177 75]./255;
    fig(side).txt4 = uicontrol('style','text','units','pixels','Background','w','ForegroundColor',txtcolor,...
    'position',[80,300,700,40],'FontSize',12,'HorizontalAlignment','left',...
    'string',sprintf([...
    'Angular difference between Upper and Lower Marker for Intensity Profiles is: ' num2str(diffintensityangle,'%.1f') ' deg\nThis seems plausible!']));
else
    txtcolor = 'r';
    fig(side).txt4 = uicontrol('style','text','units','pixels','Background','w','ForegroundColor',txtcolor,...
    'position',[80,300,700,40],'FontSize',12,'HorizontalAlignment','left',...
    'string',sprintf([...
    'Angular difference between Upper and Lower Marker for Intensity Profiles is: ' num2str(diffintensityangle,'%.1f') ' deg\nThis seems implausible!']));
end


if abs(solution.polar1) <= 40
    txtcolor = [34 177 75]./255;
elseif abs(solution.polar1) > 40 && abs(solution.polar1) <= 55
    txtcolor = [255 128 0]./255;
elseif abs(solution.polar1) > 55
    txtcolor = 'r';
end

fig(side).txt5 = uicontrol('style','text','units','pixels','Background','w','ForegroundColor',txtcolor,...
    'position',[80,260,700,20],'FontSize',12,'HorizontalAlignment','left',...
    'string',sprintf([...
    'Polar Angle is: ' num2str(round(abs(solution.polar1))) ' deg\n' ...
    ]));


if abs(solution.polar1) > 40
    fig(side).txt6 = uicontrol('style','text','units','pixels','Background','w','ForegroundColor','r',...
        'position',[80,220,700,40],'FontSize',12,'HorizontalAlignment','left',...
        'string',sprintf(['WARNING: The polar angle of the lead is larger than 40 deg and results could be inaccurate.\nPlease inspect the results carefully and use manual refinement if necessary.']));
else
    fig(side).txt6 = uicontrol('style','text','units','pixels','Background','w','ForegroundColor',[34 177 75]./255,...
        'position',[80,220,700,40],'FontSize',12,'HorizontalAlignment','left',...
        'string',sprintf(['No warnings: The polar angle is within normal range.']));
end

if max(ct.voxsize) < 1
    txtcolor = [34 177 75]./255;
else
    txtcolor = 'r';
end
fig(side).txt7 = uicontrol('style','text','units','pixels','Background','w','ForegroundColor',txtcolor,...
    'position',[80,200,700,20],'FontSize',12,'HorizontalAlignment','left',...
    'string',sprintf([...
    'CT Resolution is: ' num2str(round(ct.voxsize(1),2)) ' mm; ' num2str(round(ct.voxsize(2),2)) ' mm; ' num2str(round(ct.voxsize(3),2)) ' mm; '  '\n' ...
    ]));

if max(ct.voxsize) > 1
    fig(side).txt8 = uicontrol('style','text','units','pixels','Background','w','ForegroundColor','r',...
        'position',[80,160,700,40],'FontSize',12,'HorizontalAlignment','left',...
        'string',sprintf(['WARNING: CT resolution is larger than 1 mm and results could be inaccurate.\nPlease inspect the results carefully and use manual refinement if necessary.']));
else
    fig(side).txt8 = uicontrol('style','text','units','pixels','Background','w','ForegroundColor',[34 177 75]./255,...
        'position',[80,160,700,40],'FontSize',12,'HorizontalAlignment','left',...
        'string',sprintf(['No warnings: CT resolution is within optimal range.']));
end

fig(side).txt9 = uicontrol('style','text','units','pixels','Background','w','ForegroundColor','r',...
    'position',[80,80,700,80],'FontSize',14,'HorizontalAlignment','left',...
    'string',sprintf(['WARNING: DiODe for Medtronic leads was never validated in phantom measurements! Results might be bogus. If you want to help with development contact till.dembek@uk-koeln.de!']));



COGButton = uicontrol('Style', 'pushbutton', 'String', ['COG solution: ' num2str(round(rad2deg(roll))) ' deg'],...
    'Position', [80 20 240 25],'FontSize',12,...
    'Callback', @buttonPress);

IntensityButton = uicontrol('Style', 'pushbutton', 'String', ['Intensity Profile solution: ' num2str(round(rad2deg(rollintensity))) ' deg'],...
    'Position', [340 20 240 25],'FontSize',12,...
    'Callback', @buttonPress);

DiscardButton = uicontrol('Style', 'pushbutton', 'String', 'Discard',...
    'Position', [600 20 140 25],'FontSize',12,...
    'Callback', @buttonPress);

%% upper marker
ax1 = subplot(3,3,1);
hold on
title(ax1,'Upper Marker')
imagesc(finalartifact{1}')

view(-180,-90)
axis equal
axis off
colormap(gray)
caxis manual
caxis(cscale)
% caxis([-50 600])

centertmp = [(size(finalartifact{1},1)+1)./2, (size(finalartifact{1},2)+1)./2];
scatter(ax1,centertmp(1),centertmp(2),'g');
plot(ax1,vector{1}(:,1),vector{1}(:,2),'g','LineStyle',':')
scatter(ax1,vector{1}(valleyreal{1},1), vector{1}(valleyreal{1},2),'r');
tmpind = round(finalintensityangle(1))-1;
if tmpind < 1
    tmpind = tmpind +360;
elseif tmpind > 360
    tmpind = tmpind -360;
end
scatter(ax1,vector{1}(tmpind,1), vector{1}(tmpind,2),'g');
clear tmpind
plot(ax1,...
    [vector{1}(valleyreal{1}(1),1),vector{1}(valleyreal{1}(2),1)],...
    [vector{1}(valleyreal{1}(1),2),vector{1}(valleyreal{1}(2),2)],'r','LineStyle','--');
quiver(centertmp(1),centertmp(2),(6./ct.voxsize(1))*finalCOG_dir{1}(1),(6./ct.voxsize(2))*finalCOG_dir{1}(2),2,'LineWidth',1,'Color','g','MaxHeadSize',2)


clear centertmp COGtmp COGdirtmp
zoom(2)
%% lower marker
ax2 = subplot(3,3,2);
hold on
title(ax2,'Lower Marker')
imagesc(finalartifact{2}')
view(-180,-90)
axis equal
axis off
colormap(gray)
caxis manual
caxis(cscale)
% caxis([-50 600])


centertmp = [(size(finalartifact{2},1)+1)./2, (size(finalartifact{2},2)+1)./2];
scatter(ax2,centertmp(1),centertmp(2),'g');
plot(ax2,vector{2}(:,1),vector{2}(:,2),'g','LineStyle',':')
scatter(ax2,vector{2}(valleyreal{2},1), vector{2}(valleyreal{2},2),'r');
tmpind = round(finalintensityangle(2))-1;
if tmpind < 1
    tmpind = tmpind +360;
elseif tmpind > 360
    tmpind = tmpind -360;
end
scatter(ax2,vector{2}(tmpind,1), vector{2}(tmpind,2),'g');
clear tmpind
plot(ax2,...
    [vector{2}(valleyreal{2}(1),1),vector{2}(valleyreal{2}(2),1)],...
    [vector{2}(valleyreal{2}(1),2),vector{2}(valleyreal{2}(2),2)],'r','LineStyle','--');
quiver(centertmp(1),centertmp(2),(6./ct.voxsize(1))*finalCOG_dir{2}(1),(6./ct.voxsize(2))*finalCOG_dir{2}(2),2,'LineWidth',1,'Color','g','MaxHeadSize',2)

clear centertmp COGtmp COGdirtmp

zoom(2)
%%
ax3 = subplot(3,3,3);
hold on
title(ax3,'Sagittal View','FontWeight','bold')

imagesc(finalslice(:,round(size(finalslice,2)./4):3*round(size(finalslice,2)./4)+1));
axis equal
axis off
caxis([1500 3000])

    quiver(round(size(finalslice,2)/4), round(size(finalslice,1)/2).*1.7, -round(size(finalslice,1)/16), 0, 2,'LineWidth',1.5,'Color','g','MaxHeadSize',2)

scatter(ax3,round(size(finalslice,2)/4),round(size(finalslice,1)/2).*1.7,[],[0 0.4470 0.7410],'filled')
plot(ax3,[round(size(finalslice,2)/4), round(size(finalslice,2)/4)], [round(size(finalslice,2)/2)-75, round(size(finalslice,2)/2)+100],'LineStyle','--','Color',[0 0.4470 0.7410])
xlimit = get(ax3,'Xlim');
ylimit = get(ax3,'Ylim');
% text(xlimit(1) + 0.1 * mean(xlimit),mean(ylimit),'A','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
% text(xlimit(2) - 0.1 * mean(xlimit),mean(ylimit),'P','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
%% graphics dir level one
% ax4 = subplot(3,3,4);
% hold on
% title(ax4,'Directional Level')
% imagesc(artifact_dirnew')
% view(-180,-90)
% axis equal
% axis off
% colormap(gray)
% caxis manual
% caxis(cscale)
% plot(ax4,vectornew(:,1),vectornew(:,2),'g','LineStyle',':')
% scatter(ax4,vectornew(dirnew_valleys,1),vectornew(dirnew_valleys,2),'r')
% for k = 1:length(dirnew_valleys)
%     plot(ax4,[center_dirnew(1) (center_dirnew(1) + 1.5 * (vectornew(dirnew_valleys(k),1)-center_dirnew(1)))],...
%         [center_dirnew(2) (center_dirnew(2) + 1.5 * (vectornew(dirnew_valleys(k),2)-center_dirnew(2)))],'r','LineStyle','--')
% end
% scatter(ax4,center_dirnew(1),center_dirnew(2),[],[0 0.4470 0.7410],'filled')
% xlimit = get(ax4,'Xlim');
% ylimit = get(ax4,'Ylim');
% text(mean(xlimit),ylimit(2) - 0.15 * mean(ylimit),'A','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
% text(mean(xlimit),ylimit(1) + 0.15 * mean(ylimit),'P','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
% text(xlimit(1) + 0.1 * mean(xlimit),mean(ylimit),'L','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
% text(xlimit(2) - 0.1 * mean(xlimit),mean(ylimit),'R','Color','w','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
% 
% ax5 = subplot(3,3,5);
% hold on
% title(ax5,'Intensity Profile','FontWeight','normal')
% plot(ax5,rad2deg(anglenew),intensitynew)
% scatter(ax5,rad2deg(anglenew(dirnew_valleys)), intensitynew(dirnew_valleys),'r');
% set(ax5,'yticklabel',{[]})
% 
% ax6 = subplot(3,3,6);
% hold on
% title(ax6,'Similarity Index','FontWeight','normal')
% plot(ax6,rad2deg(rollangles{[1 2] == realsolution}),sumintensitynew{[1 2] == realsolution})
% plot(ax6,rad2deg(rollangles{[1 2] ~= realsolution}),sumintensitynew{[1 2] ~= realsolution},'r')
% scatter(ax6,...
%     rad2deg(rollangles{[1 2] == realsolution}(rollangles{[1 2] == realsolution} == 0)),...
%     sumintensitynew{[1 2] == realsolution}(rollangles{[1 2] == realsolution} == 0),...
%     'g','filled');
% scatter(ax6,...
%     rad2deg(rollangles{[1 2] == realsolution}(darkstarangle([1 2] == realsolution))),...
%     sumintensitynew{[1 2] == realsolution}(darkstarangle([1 2] == realsolution)),...
%     'r');
% text(rad2deg(rollangles{[1 2] == realsolution}(darkstarangle([1 2] == realsolution))),...
%     sumintensitynew{[1 2] == realsolution}(darkstarangle([1 2] == realsolution)),...
%     ['\leftarrow HU = ' num2str(round(sumintensitynew{[1 2] == realsolution}(darkstarangle([1 2] == realsolution))))]);
% set(ax6,'yticklabel',{[]})

%% Positioning of different subplots
% linkaxes([ax2 ax5],'xy');
% set(ax5,'Xlim',[0 360]);
% set(ax5,'Ylim',[min([intensitynew])-50 max([intensitynew])+50]);

% set(ax6,'Xlim',[rad2deg(rollangles{1}(1)) rad2deg(rollangles{1}(end))]);
% set(ax6,'Ylim',[-1000 1000]);

set(ax1,'Position',[0.45 0.75 0.2 0.2])
set(ax2,'Position',[0.45 0.50 0.2 0.2])
set(ax3,'Position',[0.1 0.475 0.445 0.5])
% set(ax4,'Position',[0.13 0.395 0.2 0.2])
% set(ax5,'Position',[0.345 0.395 0.2 0.2])
% set(ax6,'Position',[0.56 0.395 0.2 0.2])

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

view(-90,0)
ylim([-1 1])
xlim([-1 1])
zlim([0 max(electrode.insulation(end-1).vertices(:,3))+2])
axis off
axis equal

% if peakangle(side) > 180
%     tempangle = peakangle(side) - 360;
% else
%     tempangle = peakangle(side);
% end

% camorbit(-rad2deg(tempangle),0)
% tempvec = [0; 1; 0];
% temp3x3 = ea_diode_rollpitchyaw(-tempangle,0,0);
% tempvec = temp3x3 * tempvec;
% clear tempangle

set(ax_elec,'Position',[0.05 0.485 0.25 0.5])

%% get results
% if round(sumintensitynew{[1 2] == realsolution}(darkstarangle([1 2] == realsolution))) <= -200
%     checkbox1 = set(fig(side).chk1,'Value',1);
% end

uiwait

if COGButton.UserData == 1
    savestate = 1;
    roll_y = roll;
elseif IntensityButton.UserData == 1
    savestate = 1;
    roll_y = rollintensity;
elseif DiscardButton.UserData == 1
    savestate = 0;
    retrystate = 0;
end

%% saving results
if savestate == 1
%     clear x y
%     checkbox1 = get(fig(side).chk1,'Value');
%     if ~checkbox1
%         roll_y = roll;
%         disp(['Using roll angle defined by stereotactic marker: ' num2str(rad2deg(roll)) ' deg'])
%     elseif checkbox1
%         roll_y = rollnew;
%         disp(['Using corrected roll angle defined by directional level 1: ' num2str(rad2deg(rollnew)) ' deg'])
%     end
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

function angleout = convertAngle(anglein)
    %% converts deg-angles to range from -180 to +180
    if anglein > 180
        angleout = anglein - 360;
    elseif anglein < -180
        angleout = anglein + 360;
    else
        angleout = anglein;
    end
end
