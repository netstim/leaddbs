function [tmp] = ea_orient_respecifyslices(tmp,ct,cscale,electrode,activecontact)
if ct.mat(1,1) < 0
    tmp{1} = flip(tmp{1},1);
    tmp{2} = flip(tmp{2},1);
    tmp{3} = flip(tmp{3},1);
end
if ct.mat(2,2) < 0
    tmp{1} = flip(tmp{1},2);
    tmp{2} = flip(tmp{2},2);
    tmp{3} = flip(tmp{3},2);
end


subplot(3,1,3)
imagesc(tmp{1}')
view(-180,90)
hold on
colormap(gray(64))
caxis manual
caxis(cscale)
axis equal
axis off
title('Center -1')

subplot(3,1,2)
imagesc(tmp{2}')
view(-180,90)
hold on
colormap(gray(64))
caxis manual
caxis(cscale)
axis equal
axis off
title('Center')

subplot(3,1,1)
imagesc(tmp{3}')
view(-180,90)
hold on
colormap(gray(64))
caxis manual
caxis(cscale)
axis equal
axis off
title('Center +1')
%% show lead
ax_elec = axes('Position',[0.125 0.125 0.1 0.75]);
hold on
for k = 1:length(electrode.insulation)
    patch('Faces',electrode.insulation(k).faces,'Vertices',electrode.insulation(k).vertices,'Edgecolor','none','Facecolor',[electrode.lead_color electrode.lead_color electrode.lead_color]);
end
for k = 1:length(electrode.contacts)
    patch('Faces',electrode.contacts(k).faces,'Vertices',electrode.contacts(k).vertices,'Edgecolor','none','Facecolor',[electrode.contact_color electrode.contact_color electrode.contact_color]);
end
for k = activecontact
    patch('Faces',electrode.contacts(k).faces,'Vertices',electrode.contacts(k).vertices,'Edgecolor','none','Facecolor','r');
end
view (90,0)
ylim([-1 1])
xlim([-1 1])
zlim([0 15])
axis equal
axis off
%%
% msg = sprintf(['Please select the slice on which the artifact is most clearly defined!']);
% choice = questdlg(msg,'Specify slices','Center -1','Center','Center +1', 'Center');
% switch choice
%     case 'Center -1'
%         answer = 1;
%     case 'Center'
%         answer = 2;
%     case 'Center +1'
%         answer = 3;
% end
end
