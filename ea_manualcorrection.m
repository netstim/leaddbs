function ea_manualcorrection(mcfig,coords_mm,trajectory,patientname,options)
% This is the function that enables the user to manually correct the
% electrode positions that have been reconstructed by the algorithm. A
% handle for the figure, the coordinates of the contacts in millimeter
% notation, optionally manually measured coordinates, the fitted line in
% form of a 1x2 cell each containing a nx3 matrix that describes the line,
% the full path to the coronar nifti file, the name of the patient and the
% usual options struct must be handed to the function as parameters.
%
% Output parameters are the figure handle and the corrected coordinates and
% will be returned once the user presses the spacebar.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
%
% Andreas Horn

set(gcf,'color','k');
axis off
setappdata(gcf,'c_lims',[1800,2800]);
setappdata(gcf,'selectrode',0);
setappdata(gcf,'planecset',0);

set(mcfig,'KeyPressFcn',@ea_keystr);
set(mcfig, 'BusyAction','cancel', 'Interruptible','off');


if options.modality==1 % MR
try
    Vcor=spm_vol([options.root,options.prefs.patientdir,filesep,options.prefs.cornii]);
    cornii=[options.root,options.prefs.patientdir,filesep,options.prefs.cornii];
    
catch % if not present
    
    Vcor=spm_vol([options.root,options.prefs.patientdir,filesep,options.prefs.tranii]);
    cornii=[options.root,options.prefs.patientdir,filesep,options.prefs.tranii];
end
else %CT
    Vcor=spm_vol([options.root,options.prefs.patientdir,filesep,options.prefs.tranii]);
    cornii=[options.root,options.prefs.patientdir,filesep,options.prefs.tranii];
end
Vtra=spm_vol([options.root,options.prefs.patientdir,filesep,options.prefs.tranii]);
tranii=[options.root,options.prefs.patientdir,filesep,options.prefs.tranii];


setappdata(mcfig,'patientname',patientname);
setappdata(mcfig,'coords_mm',coords_mm);
setappdata(mcfig,'trajectory',trajectory);
setappdata(mcfig,'origtrajectory',trajectory);

setappdata(mcfig,'options',options);
setappdata(mcfig,'Vcor',Vcor);
setappdata(mcfig,'cornii',cornii);
setappdata(mcfig,'Vtra',Vtra);
setappdata(mcfig,'tranii',tranii);

% initialize scene
updatescene;




% initialize toolbar
ht=uitoolbar(mcfig);
captions=getappdata(gcf,'captions');
c_step=2;
minuscontrast=uipushtool(ht,'CData',ea_get_icn('contrastminus',options),'TooltipString','Decrease Contrast [C]','ClickedCallback',{@setcontrast,'c',nan});
pluscontrast=uipushtool(ht,'CData',ea_get_icn('contrastplus',options),'TooltipString','Increase Contrast [V]','ClickedCallback',{@setcontrast,'v',nan});
minusoffset=uipushtool(ht,'CData',ea_get_icn('extleft',options),'TooltipString','Increase Offset [B]','ClickedCallback',{@setcontrast,'b',nan});
plusoffset=uipushtool(ht,'CData',ea_get_icn('extright',options),'TooltipString','Decrease Offset [N]','ClickedCallback',{@setcontrast,'n',nan});

eltog(1)=uitoggletool(ht,'CData',ea_get_icn('el0',options),'TooltipString','Select Electrode 0 [0]','State','off','OnCallback',{@selectelectrode},'OffCallback',{@deselectelectrode});
eltog(4)=uitoggletool(ht,'CData',ea_get_icn('el3',options),'TooltipString','Select Electrode 3 [3]','State','off','OnCallback',{@selectelectrode},'OffCallback',{@deselectelectrode});
eltog(5)=uitoggletool(ht,'CData',ea_get_icn('el4',options),'TooltipString','Select Electrode 4 [4]','State','off','OnCallback',{@selectelectrode},'OffCallback',{@deselectelectrode});
eltog(8)=uitoggletool(ht,'CData',ea_get_icn('el7',options),'TooltipString','Select Electrode 7 [7]','State','off','OnCallback',{@selectelectrode},'OffCallback',{@deselectelectrode});

rightview=uipushtool(ht,'CData',ea_get_icn('elR',options),'TooltipString','Set view from Right [R]','ClickedCallback',{@ea_view,'r'});
leftview=uipushtool(ht,'CData',ea_get_icn('elL',options),'TooltipString','Set view from Left [L]','ClickedCallback',{@ea_view,'l'});
antview=uipushtool(ht,'CData',ea_get_icn('elA',options),'TooltipString','Set view from Anterior [A]','ClickedCallback',{@ea_view,'a'});
postview=uipushtool(ht,'CData',ea_get_icn('elP',options),'TooltipString','Set view from Posterior [P]','ClickedCallback',{@ea_view,'p'});

xview=uipushtool(ht,'CData',ea_get_icn('elX',options),'TooltipString','Set view from X-Direction [X]','ClickedCallback',{@ea_view,'x'});
yview=uipushtool(ht,'CData',ea_get_icn('elY',options),'TooltipString','Set view from Y-Direction [Y]','ClickedCallback',{@ea_view,'y'});

finish_mc=uipushtool(ht,'CData',ea_get_icn('done',options),'TooltipString','Finish manual corrections [space]','ClickedCallback',{@robotSpace});

captoggle=uitoggletool(ht,'CData',ea_get_icn('labels',options),'TooltipString','Orientation','OnCallback',{@objvisible,captions},'OffCallback',{@objinvisible,captions},'State','on');



setappdata(mcfig,'eltog',eltog);


elplot=getappdata(mcfig,'elplot');

try
    realcoords_plot=getappdata(mcfig,'realcoords_plot');
end
trajectory_plot=getappdata(mcfig,'trajectory_plot');
planes=getappdata(mcfig,'planes');



%% Manual height correction here:
%set(mcfig,'Position',[10 400 700 500])

drawnow;

disp('Manual correction: Press arrows to adjust, space to end adjustment. For more shortcuts hover over buttons in the menu-bar.');


ea_draggable(elplot(1),@updatescene,'endfcn',@updatescene);
ea_draggable(elplot(4),@updatescene,'endfcn',@updatescene);
ea_draggable(elplot(5),@updatescene,'endfcn',@updatescene);
ea_draggable(elplot(8),@updatescene,'endfcn',@updatescene);


%% export variables to figure
setappdata(gcf,'eltog',eltog);
setappdata(gcf,'coords_mm',coords_mm);






function ea_endfcn
% This subfunction terminates the process of manual correction and saves
% results.
    disp('Saving results.');
coords_mm=getappdata(gcf,'coords_mm');
trajectory=getappdata(gcf,'trajectory');
options=getappdata(gcf,'options');

close(gcf)

% save results

try
ea_write_fiducials(coords_mm,fullfile(options.root,options.patientname,'ea_coords.fcsv'),options);
end
elmodel=options.elmodel;
save([options.root,options.patientname,filesep,'ea_reconstruction'],'trajectory','coords_mm','elmodel');
disp('Done.');

if options.autoimprove
    disp('Storing results in template.');
    ea_export_templates(coords_mm{1}(1:4,:),trajectory{1},options.patientname,options,'r')
    ea_export_templates(coords_mm{2}(1:4,:),trajectory{2},options.patientname,options,'l')
    disp('Done.');
end
% continue with rest of the program schedule..

ea_write(options);




function ea_keystr(mcfig,event)
%    pause
%commnd=get (gcf, 'CurrentKey');


%% get vars
eltog=getappdata(gcf,'eltog');
elplot=getappdata(gcf,'elplot');
coords_mm=getappdata(gcf,'coords_mm');




commnd=event.Character;
switch lower(commnd)
    case ' '
        coords_mm=getappdata(mcfig,'coords_mm');
        ea_endfcn;
        return
    case {'x','a','p','y','l','r'} % view angles.
        coords_mm=getappdata(mcfig,'coords_mm');
        ea_view(nan,nan,commnd);
        
    case {'0','3','4','7'}
        selectrode=str2double(commnd)+1; % +1 because of name convention 0->1.
        oselectrode=getappdata(gcf,'selectrode');
        if selectrode==oselectrode % toggle had already been clicked -> deselect all.
            % reset all toggletools
            for i=[1,4,5,8]
                set(eltog(i),'State','off');
            end
            setappdata(gcf,'selectrode',0);
            updatescene;
        else
            % clear all toggletools.
            for i=[1,4,5,8]
                set(eltog(i),'State','off');
            end
            
            % set the correct toggletool again.
            set(eltog(selectrode),'State','on');
            
            % store selected electrode in appdata.
            setappdata(gcf,'selectrode',selectrode);
            updatescene;
        end
        clear oselectrode
    case {'c','v','b','n'}
        setcontrast(nan,nan,event.Character,event.Modifier);
%     case 'v'
%         increasecontrast(nan,nan,event);
%     case 'b'
%         increaseoffset(nan,nan,event);
%     case 'n'
%         decreaseoffset(nan,nan,event);
        
        
    otherwise % arrow keys, plus, minus
        
        if ismember(event.Key,{'rightarrow','leftarrow','uparrow','downarrow'}) || ismember(event.Character,{'+','-','*','_'})
        selectrode=getappdata(gcf,'selectrode');
        if ~selectrode % no electrode is highlighted, move electrodes alongside trajectory or increase/decrease spacing.
            coords_mm=getappdata(mcfig,'coords_mm');
            trajectory=getappdata(mcfig,'trajectory');
            
            coords_mm=ea_correctcoords(coords_mm,trajectory,event);
            
            setappdata(mcfig,'coords_mm',coords_mm);
            updatescene;
            coords_mm=getappdata(mcfig,'coords_mm');
        else % electrode is highlighted. Move in xy dirs.
            
            coords_mm=getappdata(gcf,'coords_mm');
            trajectory=getappdata(gcf,'trajectory');
            movedcoords=moveonecoord(coords_mm,selectrode,event); % move the correct coord to the correct direction.
            cnt=1;
            for side=1:length(coords_mm)
                for el=1:4
                    set(elplot(cnt),'XData',movedcoords{side}(el,1),'YData',movedcoords{side}(el,2),'ZData',movedcoords{side}(el,3))
                    cnt=cnt+1;
                end
            end
            setappdata(mcfig,'coords_mm',coords_mm);
            updatescene;
            coords_mm=getappdata(mcfig,'coords_mm');
            %update_coords(elplot(selectrode),coords_mm,trajectory,movedcoords); % refresh scene view (including update for all other electrodes).
            %setappdata(gcf,'coords_mm',coords_mm);
            %setappdata(gcf,'trajectory',trajectory);
            
        end
        end
end
%end
cnt=1;
for side=1:length(coords_mm)
    for el=1:4
        set(elplot(cnt),'XData',coords_mm{side}(el,1),'YData',coords_mm{side}(el,2),'ZData',coords_mm{side}(el,3))
        cnt=cnt+1;
    end
end
refreshdata(elplot,'caller')
drawnow










function hdtrajectory=genhd_inside(trajectory)

resolution=20;

hdtrajectory(:,1)=interp1q([1:length(trajectory)]',trajectory(:,1),[1:1/resolution:length(trajectory)]');
hdtrajectory(:,2)=interp1q([1:length(trajectory)]',trajectory(:,2),[1:1/resolution:length(trajectory)]');
hdtrajectory(:,3)=interp1q([1:length(trajectory)]',trajectory(:,3),[1:1/resolution:length(trajectory)]');



function updatescene


%% inputs:

patientname=getappdata(gcf,'patientname');
coords_mm=getappdata(gcf,'coords_mm');
trajectory=getappdata(gcf,'trajectory');
cornii=getappdata(gcf,'cornii');
options=getappdata(gcf,'options');
movedel=getappdata(gcf,'movedel');
trajectory_plot=getappdata(gcf,'trajectory_plot');
spacetext=getappdata(gcf,'spacetext');
planes=getappdata(gcf,'planes');
c_lims=getappdata(gcf,'c_lims');

elplot = getappdata(gcf,'elplot');
selectrode=getappdata(gcf,'selectrode');

if ~isempty(selectrode) && selectrode>0
    coordhandle=elplot(selectrode);
end





xdata = cell2mat(get(elplot,'xdata'));
ydata = cell2mat(get(elplot,'ydata'));
zdata = cell2mat(get(elplot,'zdata'));

if selectrode
    [coords_mm,trajectory]=update_coords(coordhandle,coords_mm,trajectory,[xdata,ydata,zdata]);
end

















%% plot main figure

delete(trajectory_plot);
delete(spacetext);



mainax=subplot(4,4,[1:3,5:7,9:11,13:15]); % main plot
set(gca, 'LooseInset', [0,0,0,0]);
init=getappdata(gcf,'init');
if isempty(init)
    view(0,0);
    axis off
    setappdata(gcf,'init',1)
end

%% plot coords
hold on



if isempty(elplot) % first time plot electrode contacts
    clear elplot
    cnt=1;
    for side=1:length(coords_mm)
        for el=1:size(coords_mm{side},1)
            elplot(cnt)=plot3(coords_mm{side}(el,1),coords_mm{side}(el,2),coords_mm{side}(el,3),'O','MarkerEdgeColor',[0.9 0.9 0.9],'MarkerFaceColor','none','MarkerSize',25);
            cnt=cnt+1;
        end
    end
    
end

emp_eldist(1)=mean([pdist([coords_mm{1}(1,:);coords_mm{1}(4,:)]),pdist([coords_mm{2}(1,:);coords_mm{2}(4,:)])])/3;
spacetext=text(0,0,-17,['Electrode Spacing: ',num2str(emp_eldist),' mm.'],'Color','w','BackgroundColor','k','HorizontalAlignment','center'); 
set(gcf,'name',[options.patientname,', Electrode Spacing: ',num2str(emp_eldist),' mm.']);
setappdata(gcf,'spacetext',spacetext);

% %% plot lines

for side=options.sides
    
    try
        
        if ~isempty(trajectory{side})
            
            if options.verbose>1; trajectory_plot(side)=plot3(trajectory{side}(:,1),trajectory{side}(:,2),trajectory{side}(:,3),'color',[0.5,0.5,0.5],'linew',1.5); end
            
        end
    end
end




delete(planes);
clear planes
planecnt=1;
%% plot slices in x and y planes

for doxx=0:1
    for side=options.sides
        try
            sample_width=20-doxx*5; % a bit smaller sample size in x direction to avoid overlap.
            meantrajectory=genhd_inside(trajectory{side});
            clear imat
            %% sample plane left and right from meantrajectory
            
            
            Vcor=getappdata(gcf,'Vcor');
            
            imat=ea_resample_planes(Vcor,meantrajectory',sample_width,doxx,0.1);
            
            colormap gray
            
            
            if doxx % span surface in x direction
                spanvector=[sample_width,0,0];
            else % span surface in y direction
                spanvector=[0,sample_width,0];
            end
            
            boundingbox=[meantrajectory(1,:)-spanvector;...
                meantrajectory(1,:)+spanvector;...
                meantrajectory(end,:)-spanvector;...
                meantrajectory(end,:)+spanvector];
            
            
            xx=[boundingbox(1,1),boundingbox(2,1);boundingbox(3,1),boundingbox(4,1)];
            yy=[boundingbox(1,2),boundingbox(2,2);boundingbox(3,2),boundingbox(4,2)];
            zz=[boundingbox(1,3),boundingbox(2,3);boundingbox(3,3),boundingbox(4,3)];
            
            
            alphamap=imat;
            alphamap(:)=0.9;
            
            if ~getappdata(gcf,'planecset') % initially and once set contrast based on image data.
                
                if options.modality==1
                c_lims=[ea_nanmean(imat(:))-nanstd(imat(:))-3*nanstd(imat(:)),ea_nanmean(imat(:))-nanstd(imat(:))+3*nanstd(imat(:))];
                elseif options.modality==2
                        c_lims=[1800,2800]; % Initial guess, CT
                end
                caxis(c_lims);
                setappdata(gcf,'c_lims',c_lims);
                setappdata(gcf,'planecset',1);
            end
            
            
            planes(planecnt)=surface('XData',xx,'YData',yy,'ZData',zz,'CData',imat,'alphadata',alphamap,'FaceAlpha', 'texturemap','FaceColor','texturemap','EdgeColor','none','alphadatamapping','none');
            
            planecnt=planecnt+1;
            if ~doxx
                
                if side==1
                    captions(1)=text(0,... % x
                        ((min(boundingbox(:,2))+max(boundingbox(:,2)))/2)+20,... % y
                        0,... % z
                        'A','Color','w','BackgroundColor','k');
                    captions(2)=text(0,... % x
                        ((min(boundingbox(:,2))+max(boundingbox(:,2)))/2)-20,... % y
                        0,... % z
                        'P','Color','w','BackgroundColor','k');
                    % captions(1)=text(0,0,0,... % z
                    %  'C','Color','w');
                    
                    captions(3)=text(40,... % x
                        0,... % y
                        0,... % z
                        'R','Color','w','BackgroundColor','k');
                    
                    captions(4)=text(-40,... % x
                        0,... % y
                        0,... % z
                        'L','Color','w','BackgroundColor','k');
                    
                    
                    setappdata(gcf,'captions',captions);
                    
                end
            else
                
            end
        end
    end
end
caxis([c_lims(1) c_lims(2)]);



%% plot axial planes on the right hand side of the figure


Vtra=getappdata(gcf,'Vtra');
for side=1:length(coords_mm)
    coords{side}=Vtra.mat\[coords_mm{side},ones(size(coords_mm{side},1),1)]';
    coords{side}=coords{side}(1:3,:)';
end

%title(['Electrode ',num2str(el-1),', transversal view.']);
wsize=10;
cmap=[1,4,5,8];
for subpl=1:4
    subplot(4,4,subpl*4)
    
    slice=ea_sample_slice(Vtra,'tra',wsize,coords,cmap(subpl));
    try
        imagesc(slice,[ea_nanmean(slice(slice>0))-3*nanstd(slice(slice>0)) ea_nanmean(slice(slice>0))+3*nanstd(slice(slice>0))]);
    catch
        imagesc(slice);
    end
    
    hold on
    
    
    
    if selectrode && cmap(subpl)==selectrode
        
        fc='r';
    else
        fc='none';
    end
    
    plot((wsize+1)*2,(wsize+1)*2,'O','MarkerSize',15,'MarkerEdgeColor',[0.9 0.9 0.9],'MarkerFaceColor',fc,'LineWidth',2,'LineSmoothing','on');
    hold off
    axis square
    axis off
    caxis([c_lims(1) c_lims(2)]);
    
end




%% outputs

% try
%     setappdata(resultfig,'realcoords_plot',realcoords_plot);
% end

set(gcf,'CurrentAxes',mainax);
setappdata(gcf,'planes',planes);







%% outputs

setappdata(gcf,'elplot',elplot);
setappdata(gcf,'movedel',movedel);
setappdata(gcf,'coords_mm',coords_mm);
% try
%     setappdata(resultfig,'realcoords_plot',realcoords_plot);
% end
setappdata(gcf,'trajectory_plot',trajectory_plot);
setappdata(gcf,'planes',planes);
setappdata(gcf,'trajectory',trajectory);






function movedel=whichelmoved(coordhandle)

elplot = getappdata(gcf,'elplot');

for el=[1,4,5,8];
    
    if coordhandle==elplot(el)
        movedel=el;
    end
    
end

function   [nucoords,nutrajectory]=update_coords(coordhandle,coords,trajectory,movedcoords)

% movedcoords are NOT in cell format!

movel=whichelmoved(coordhandle);

elplot=getappdata(gcf,'elplot');
origtrajectory=getappdata(gcf,'origtrajectory');
%whichonemoved=1;

nucoords=coords;
nutrajectory=trajectory;

switch movel
    
    case 1 % lower right
        olddist=abs(pdist([coords{1}(3,:);coords{1}(4,:)]));
        refel=4;
        changecoord=1:3;
        cchangecoord=1;
        usefit=1;
        spin=-1;
        
    case 5 % lower left
        olddist=abs(pdist([coords{2}(3,:);coords{2}(4,:)]));
        refel=4;
        changecoord=1:3;
        cchangecoord=2;
        usefit=2;
        spin=-1;
    case 4 % upper right
        olddist=abs(pdist([coords{1}(1,:);coords{1}(2,:)]));
        refel=1;
        changecoord=2:4;
        cchangecoord=1;
        usefit=1;
        spin=1;
    case 8 % upper left
        olddist=abs(pdist([coords{2}(1,:);coords{2}(2,:)]));
        refel=1;
        changecoord=2:4;
        cchangecoord=2;
        usefit=2;
        spin=1;
end


refpt=coords{cchangecoord}(refel,:);
movedpt=movedcoords(movel,:); % movedcoords are NOT in cell format!
helppt=[movedpt(1:2),refpt(3)]; % this point is constructed as a projection of the movedpoint to the z-axis at height of refpt.
zdistfromhelppt=sqrt((3*olddist)^2-(pdist([helppt;refpt]))^2); % use Pythagorean theorem to calculate distance on z-axis.
% here, 3*olddist is the hypotenuse, pdist between the helper point and the
% reference point is one leg of the triangle, zdistfromhelppt is the
% leg that is calculated (the distance from the helper point to the new
% moved point (fulfilling the prerequisite that new moved point is at the
% same distance than the old point).
movedpt=helppt+spin*[0 0 zdistfromhelppt]; % set movedpt to be equidistant.

trajpts=[refpt;movedpt];

nutraj=diff(trajpts);
nutraj=nutraj/norm(nutraj);
%disp(num2str(nutraj));

% generate new coords
for ccord=changecoord
    nucoords{cchangecoord}(ccord,:)=refpt+(nutraj*(olddist*abs(refel-ccord)));
    %for xyz=1:2
    %    nucoords(ccord,xyz)=movedpt(xyz)-(movedpt(xyz)-refpt(xyz))*((oldmovedpt(xyz)-coords(ccord,xyz))/(oldmovedpt(xyz)-refpt(xyz)));
    %end
    
    set(elplot(ccord),'xdata',nucoords{cchangecoord}(ccord,1));
    set(elplot(ccord),'ydata',nucoords{cchangecoord}(ccord,2));
    set(elplot(ccord),'zdata',nucoords{cchangecoord}(ccord,3));
end

% generate new trajectory
nutrajectory{usefit}=[];
cnt=1;

for ix=-50:0.5:50
    
    thispoint=refpt+(nutraj*ix);
    if thispoint(3)<max(origtrajectory{usefit}(:,3)) && thispoint(3)>min(origtrajectory{usefit}(:,3))
        nutrajectory{usefit}(cnt,:)=thispoint;
        cnt=cnt+1;
    end
    
    
end
nutrajectory{usefit}=sortrows(nutrajectory{usefit},-3); % Make sure that trajectory is always listed in correct order.





function coords_mm=moveonecoord(coords_mm,selectrode,command)

grone=[0.1,1]; % step-sizes

movedcoords=[coords_mm{1};coords_mm{2}];

switch command.Key
    case 'leftarrow'
        movedcoords(selectrode,1)=movedcoords(selectrode,1)-grone(1+ismember('shift',command.Modifier));
    case 'rightarrow'
        movedcoords(selectrode,1)=movedcoords(selectrode,1)+grone(1+ismember('shift',command.Modifier));
    case 'uparrow'
        movedcoords(selectrode,2)=movedcoords(selectrode,2)+grone(1+ismember('shift',command.Modifier));
        
    case 'downarrow'
        movedcoords(selectrode,2)=movedcoords(selectrode,2)-grone(1+ismember('shift',command.Modifier));
        
        
        
end
coords_mm{1}=movedcoords(1:4,:);
coords_mm{2}=movedcoords(5:8,:);

function objvisible(hobj,ev,atls)
set(atls, 'Visible', 'on');
%disp([atls,'visible clicked']);

function objinvisible(hobj,ev,atls)
set(atls, 'Visible', 'off');
%disp([atls,'invisible clicked']);


function setcontrast(hobj,ev,key,modifier)
c_lims=getappdata(gcf,'c_lims');
comms={'v','c','b','n'}; % key commands
perfs=[1 -1 % actions to perform on key commands
      -1 1
       1 1
      -1 -1];
kern=100; % gain to correct contrast

doshift=any(ismember('shift',modifier));

c_lims=c_lims+(perfs(ismember(comms,lower(key)),:)*(kern*(doshift+1)));
setappdata(gcf,'c_lims',c_lims);
updatescene;


function selectelectrode(hobj,ev)
% reset all toggletools

warning ('off','all');

eltog=getappdata(gcf,'eltog');
elplot=getappdata(gcf,'elplot');

for i=[1,4,5,8]
    set(eltog(i),'State','off');
    set(elplot(i),'MarkerFaceColor','none');
    if eltog(i)==hobj;
        selectrode=i;
        
    end
end

% set the clicked toggletool again.
set(eltog(selectrode),'State','on');
set(elplot(selectrode),'MarkerFaceColor','r');

% store selected electrode in appdata.
setappdata(gcf,'selectrode',selectrode);
warning ('off','all');




function deselectelectrode(hobj,ev)
% reset all toggletools
eltog=getappdata(gcf,'eltog');
elplot=getappdata(gcf,'elplot');

for i=[1,4,5,8]
    set(eltog(i),'State','off');
    set(elplot(i),'MarkerFaceColor','none');
    
end
setappdata(gcf,'selectrode',0);

function ea_view(hobj,ev,commnd)
switch commnd
    case 'p'
        view(0,0);
    case {'x','a'}
        view(180,0);
    case {'y','r'}
        view(90,0);
    case 'l'
        view(270,0);
end

function ea_finish(hobj,ev)
disp('Manual correction done.');


function robotSpace(hobj,ev) % simulates key presses using Java.
import java.awt.Robot;
import java.awt.event.*;
SimKey=Robot;
SimKey.keyPress(KeyEvent.VK_SPACE)
SimKey.keyRelease(KeyEvent.VK_SPACE)

function y = ea_nanmean(varargin)
if nargin==2
    x=varargin{1};
    dim=varargin{2};
elseif nargin==1
x=varargin{1};
    dim=1;
end
    
N = sum(~isnan(x), dim);
y = nansum(x, dim) ./ N;
