function [roll_y,y,solution] = ea_diode_manualGUI(side,ct,head_mm,unitvector_mm,tmat_vx2mm,elspec)
%% constant variables
% colorscale for ct figures
cscale = [-50 100];
cscale2 = [1500 3000];
% diameter of the slices shown in visualizations in vox
extractradius = 30;
% sidenames
sides = {'right','left','3','4','5','6','7','8'};

%% define electrode properties
% z position of the centers of directional levels 1 and 2 and the marker relative to
% head position
% level1centerRelative = elspec.contact_length + elspec.contact_spacing;
% level2centerRelative = (elspec.contact_length + elspec.contact_spacing) * 2;
% switch elspec.matfname
%     case 'medtronic_b33005'
%         markercenterRelative = 12.75;
%     case 'medtronic_b33015'
%         markercenterRelative = 15.75;
%     otherwise
%         keyboard
% end

load(elspec.matfname);

%% load CTs
%% Recalculate trajvector to optimize position at center of artifacts
% this part recalculates the exact lead position by calculating the
% center of mass of each CT slice from the head to 10mm above it and
% then performing a linear regression through these centers of mass

% sampling from head to 10mm above it in .5mm steps
disp('Refining Lead Trajectory!')
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
% marker_mm = head_mm + (markercenterRelative .* unitvector_mm);
% transform to vx
% marker_vx = tmat_vx2mm\marker_mm;

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
% mincorner_mm = tmat_vx2mm * [1;1;1;1];
% maxcorner_mm = tmat_vx2mm * [size(ct.img)';1];
%
% [Xmm,Ymm,Zmm] = meshgrid([mincorner_mm(1):(maxcorner_mm(1)-mincorner_mm(1))./(size(ct.img,1)-1):maxcorner_mm(1)],...
%     [mincorner_mm(2):(maxcorner_mm(2)-mincorner_mm(2))./(size(ct.img,2)-1):maxcorner_mm(2)],...
%     [mincorner_mm(3):(maxcorner_mm(3)-mincorner_mm(3))./(size(ct.img,3)-1):maxcorner_mm(3)]);
%
% clear mincorner_mm maxcorner_mm
%
% Vnew = permute(ct.img,[2,1,3]);

mincorner_mm = tail_mm - ([20; 20; 20; 0] .* sign(diag(tmat_vx2mm))); % tmat_vx2mm * [1;1;1;1];
maxcorner_mm = tail_mm + ([20; 20; 20; 0] .* sign(diag(tmat_vx2mm))); % tmat_vx2mm * [size(ct.img)';1];

mincorner_vx = round(tmat_vx2mm \ mincorner_mm);
maxcorner_vx = round(tmat_vx2mm \ maxcorner_mm);

ct.imgtmp = ct.img(mincorner_vx(1):maxcorner_vx(1),mincorner_vx(2):maxcorner_vx(2),mincorner_vx(3):maxcorner_vx(3));
[Xmm,Ymm,Zmm] = meshgrid([mincorner_mm(1):(maxcorner_mm(1)-mincorner_mm(1))./(size(ct.imgtmp,1)-1):maxcorner_mm(1)],...
    [mincorner_mm(2):(maxcorner_mm(2)-mincorner_mm(2))./(size(ct.imgtmp,2)-1):maxcorner_mm(2)],...
    [mincorner_mm(3):(maxcorner_mm(3)-mincorner_mm(3))./(size(ct.imgtmp,3)-1):maxcorner_mm(3)]);

clear mincorner_mm maxcorner_mm

Vnew = permute(ct.imgtmp,[2,1,3]);

%% slice perpendicular through the two markers
% a 5mm slice with .1mm resolution is sampled perpendicular to
% the lead over the range of the two markers, first the lower, than the
% upper. COG is calculated for each slice at .2mm resolution. slices
% with maximum COG distance to the geometric centers are calculated.
% direction of deviation for these slices with respect to anterior
% orientation are calculated.
extract_width = 5; % in mm
extract_height = 10; % in mm
samplingres = .1;
allrolls = deg2rad(1:1:360);
sagittal_shift = 7.5;
sagittal_limits = [sagittal_shift - extract_height, sagittal_shift + extract_height] ./ samplingres;
sagittal_headmodifier = (extract_height - sagittal_shift)  ./ samplingres;
sagittal_trajmodifier = (extract_height - sagittal_shift + 15)  ./ samplingres;
%% Slices parallel to the lead shaft are calculated for all orientations
for k = 1:length(allrolls)
    disp(['Calculating parallel Slices ' num2str(k) '/360'])
    roll_act = allrolls(k);
    [M,~,~,~] = ea_diode_rollpitchyaw(roll_act,pitch,yaw);
    yvec_mm = M * [0;1;0];
    xvec_mm = cross(unitvector_mm(1:3), yvec_mm);
    clear M
    Xslice = ([-extract_height:samplingres:extract_height] .* unitvector_mm(1)) + ([-extract_width:samplingres:extract_width] .* yvec_mm(1))' + head_mm(1) + sagittal_shift * unitvector_mm(1);
    Yslice = ([-extract_height:samplingres:extract_height] .* unitvector_mm(2)) + ([-extract_width:samplingres:extract_width] .* yvec_mm(2))' + head_mm(2) + sagittal_shift * unitvector_mm(2);
    Zslice = ea_diode_perpendicularplane(xvec_mm,tail_mm,Xslice,Yslice);
    parallelslices{k} = interp3(Xmm,Ymm,Zmm,Vnew,Xslice,Yslice,Zslice);
    parallelslices{k} = parallelslices{k}';
    parallelslices{k} = flipdim(parallelslices{k},2);
end

%% get parralel slices from head to 20mm above head;
extract_width = 10;
myslices = [sagittal_shift-extract_height:samplingres:sagittal_shift+extract_height-samplingres];
roll_act = 0;
[M,~,~,~] = ea_diode_rollpitchyaw(roll_act,pitch,yaw);
yvec_mm = M * [0;1;0];
xvec_mm = cross(unitvector_mm(1:3), yvec_mm);
for k = 1:length(myslices)
    disp(['Calculating transversal Slices ' num2str(k) '/' num2str(length(myslices))])
    center_act = head_mm + myslices(k) * unitvector_mm;
    myZvals(k) = center_act(3);
    Xslice = ([-extract_width:samplingres:extract_width] .* xvec_mm(1)) + ([-extract_width:samplingres:extract_width] .* yvec_mm(1))' + center_act(1);
    Yslice = ([-extract_width:samplingres:extract_width] .* xvec_mm(2)) + ([-extract_width:samplingres:extract_width] .* yvec_mm(2))' + center_act(2);
    Zslice = ea_diode_perpendicularplane(unitvector_mm,center_act,Xslice,Yslice);
    transversalslices{k} = interp3(Xmm,Ymm,Zmm,Vnew,Xslice,Yslice,Zslice);
    transversalslices{k} = flipdim(transversalslices{k},1);
end

%% final figure
fig(side).figure = uifigure('Name',['Lead ' sides{side}],'Position',[100 100 800 800],'Color','w','Toolbar','none');

%% get data to plot circle
tmpvector = [0 1] .* 50;
vector = [];
center = [round(size(transversalslices{1},2)/2),round(size(transversalslices{1},1)/2)];
for k = 1:360
    theta = (2*pi/360) * (k);    
    rotmat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
    vector =vertcat(vector,(tmpvector * rotmat) + center);
    angle(k) = theta;
end
clear theta rotmat tmpvector center

%% parralel Plot axis
orientation_act = 360;
level_act = round(length(transversalslices)/2);
ax1 = uiaxes(fig(side).figure,'Position',[90 130 300 600]);
hold(ax1,'on')
parralelslice_act = imagesc(ax1,parallelslices{orientation_act});
scatter(ax1,round(size(parallelslices{orientation_act},2)/2),sagittal_headmodifier,[],[0 0.4470 0.7410],'filled')
plot(ax1,[round(size(parallelslices{orientation_act},2)/2), round(size(parallelslices{orientation_act},2)/2)], [sagittal_headmodifier, sagittal_trajmodifier],'LineStyle','--','Color',[0 0.4470 0.7410])
level_marking = plot(ax1,[1,round(size(parallelslices{orientation_act},2))], [level_act, level_act],'Color',[0 0.4470 0.7410],'Tag','level_marking');
hold(ax1,'off')
axis(ax1,'equal')
axis(ax1,'off')
set(ax1,'colormap',gray)
caxis(ax1,cscale2)

%% transversal Plot axis
tmpcenter = [round(size(transversalslices{level_act},2)/2),round(size(transversalslices{level_act},1)/2)];
redline_angles = [orientation_act-90, orientation_act+90] + [((orientation_act-90)<1) .* 360 , ((orientation_act+90)>360) .* (-360)];

ax2 = uiaxes(fig(side).figure,'Position',[450 430 280 280]);
hold(ax2,'on')
transversalslice1_act = imagesc(ax2,transversalslices{level_act});
scatter(ax2,tmpcenter(1),tmpcenter(2),[],[0 0.4470 0.7410],'filled')
% plot(ax2,vector(:,1),vector(:,2),'g','LineStyle',':')
arrow1_act = quiver(ax2,tmpcenter(1),tmpcenter(2),(vector(orientation_act,1)-tmpcenter(1))/2,(vector(orientation_act,2)-tmpcenter(2))/2,2,'LineWidth',1,'Color','g','MaxHeadSize',2);
redline1 = plot(ax2,vector(redline_angles,1),vector(redline_angles,2),'r','LineStyle','--');
hold(ax2,'off')
axis(ax2,'equal')
axis(ax2,'off')
view(ax2,[-180 -90])
set(ax2,'colormap',gray)
caxis(ax2,cscale2)

ax3 = uiaxes(fig(side).figure,'Position',[450 150 280 280]);
hold(ax3,'on')
transversalslice2_act = imagesc(ax3,transversalslices{level_act});
scatter(ax3,round(size(transversalslices{level_act},2)/2),round(size(transversalslices{level_act},1)/2),[],[0 0.4470 0.7410],'filled')
% plot(ax3,vector(:,1),vector(:,2),'g','LineStyle',':')
arrow2_act = quiver(ax3,tmpcenter(1),tmpcenter(2),(vector(orientation_act,1)-tmpcenter(1))/2,(vector(orientation_act,2)-tmpcenter(2))/2,2,'LineWidth',1,'Color','g','MaxHeadSize',2);
redline2 = plot(ax3,vector(redline_angles,1),vector(redline_angles,2),'r','LineStyle','--');
hold(ax3,'on')
axis(ax3,'equal')
axis(ax3,'off')
view(ax3,[-180 -90])
set(ax3,'colormap',gray)
caxis(ax3,cscale)

%% Texts
orientation_label = uilabel(fig(side).figure, 'Position', [450 120 280 40],'HorizontalAlignment','center','VerticalAlignment','bottom','Fontsize',24,...
    'Text', ['Orientation: ' num2str(orientation_act - ((orientation_act > 180).*360)) '°']);
%% Controls
RotationSlider = uislider(fig(side).figure,'Value',0,'Limits',[-180 180],'MajorTicks',[-180:90:180],...
    'Position', [90 140 300 3],'FontSize',12,...
    'ValueChangedFcn',@(RotationSlider,eventdata) rotationChanged(eventdata,parallelslices,parralelslice_act,transversalslices,arrow1_act,arrow2_act,redline1,redline2,vector,orientation_label),'ValueChangingFcn',@(RotationSlider,eventdata) rotationChanged(eventdata,parallelslices,parralelslice_act,transversalslices,arrow1_act,arrow2_act,redline1,redline2,vector,orientation_label));
LevelSlider = uislider(fig(side).figure,'Value',level_act,'Limits',[1 length(transversalslices)],'Orientation','vertical',...
    'Position', [420 167 3 534],'FontSize',12,...
    'ValueChangedFcn',@(LevelSlider,eventdata) levelChanged(eventdata,transversalslices,parallelslices,orientation_act,level_marking,transversalslice1_act,transversalslice2_act),'ValueChangingFcn',@(LevelSlider,eventdata) levelChanged(eventdata,transversalslices,parallelslices,orientation_act,level_marking,transversalslice1_act,transversalslice2_act));
SaveButton = uibutton(fig(side).figure,'push', 'Text', 'Accept & Save',...
    'Position', [150 20 150 25],'FontSize',12,...
    'ButtonPushedFcn', @buttonPress);
DiscardButton =  uibutton(fig(side).figure,'push', 'Text', 'Discard',...
    'Position', [500 20 150 25],'FontSize',12,...
    'ButtonPushedFcn', @buttonPress);

%% graphics lead
ax_elec = uiaxes(fig(side).figure,'Position',[20 200 100 510]);
hold on
for k = 1:length(electrode.insulation)
    patch(ax_elec,'Faces',electrode.insulation(k).faces,'Vertices',electrode.insulation(k).vertices,'Edgecolor','none','Facecolor',[electrode.lead_color electrode.lead_color electrode.lead_color]);
end
for k = 1:length(electrode.contacts)
    patch(ax_elec,'Faces',electrode.contacts(k).faces,'Vertices',electrode.contacts(k).vertices,'Edgecolor','none','Facecolor',[electrode.contact_color electrode.contact_color electrode.contact_color]);
end


ylim(ax_elec,[-1 1])
xlim(ax_elec,[-1 1])
zlim(ax_elec,[0 20])
axis(ax_elec,'equal')
view(ax_elec,[-90,0])
axis(ax_elec,'off')

uiwait

if SaveButton.UserData == 1
    savestate = 1;
elseif DiscardButton.UserData == 1
    savestate = 0;
    retrystate = 0;
end

%% saving results
if savestate == 1
    clear x y
    roll_y = deg2rad(RotationSlider.Value);
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
end
close(fig(side).figure)
end

function buttonPress(hObject,eventdata)
hObject.UserData = 1;
uiresume
end

function [orientation_act] = rotationChanged(eventdata,parallelslices,parralelslice_act,transversalslices,arrow1_act,arrow2_act,redline1,redline2,vector,orientation_label)
orientation_act = round(eventdata.Value);
if orientation_act < 1
    orientation_act = orientation_act + 360;
end
set(parralelslice_act,'CData',parallelslices{orientation_act})
%% rotate arrow in transversalviews
tmpcenter = [round(size(transversalslices{1},2)/2),round(size(transversalslices{1},1)/2)];
set(arrow1_act,'UData',(vector(orientation_act,1)-tmpcenter(1))/2,'VData',(vector(orientation_act,2)-tmpcenter(2))/2);
set(arrow2_act,'UData',(vector(orientation_act,1)-tmpcenter(1))/2,'VData',(vector(orientation_act,2)-tmpcenter(2))/2);
redline_angles = [orientation_act-90, orientation_act+90] + [((orientation_act-90)<1) .* 360 , ((orientation_act+90)>360) .* (-360)];
set(redline1,'XData',vector(redline_angles,1),'YData',vector(redline_angles,2));
set(redline2,'XData',vector(redline_angles,1),'YData',vector(redline_angles,2));
set(orientation_label,'Text',['Orientation: ' num2str(orientation_act - ((orientation_act > 180).*360)) '°']);
end

function [level_act] = levelChanged(eventdata,transversalslices,parallelslices,orientation_act,level_marking,transversalslice1_act,transversalslice2_act)
level_act = round(eventdata.Value);
set(transversalslice1_act,'CData',transversalslices{level_act})
set(transversalslice2_act,'CData',transversalslices{level_act})

%% redraw ax1 to shift level markings
set(level_marking, 'XData',[1,round(size(parallelslices{orientation_act},2))],'YData',[level_act, level_act])
end