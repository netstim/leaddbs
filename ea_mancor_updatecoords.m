function [numarkers,nutrajectory]=ea_mancor_updatecoords(coordhandle,markers,trajectory,movedheadtail,options,mcfig)

% movedcoords are NOT in cell format!
movel=whichelmoved(coordhandle,mcfig);
mplot=getappdata(gcf,'mplot');
origtrajectory=getappdata(gcf,'origtrajectory');
%whichonemoved=1;
numarkers=markers;
nutrajectory=trajectory;

switch movel
    case 1 % lower contact
        olddist=abs(ea_pdist([markers(options.elside).head;markers(options.elside).tail]));
        cchangecoord=1;
        usefit=options.elside;
        if markers(options.elside).head(3) < markers(options.elside).tail(3)
            spin = -1;
        else % In case tail is lower than head
            spin = 1;
        end

        refpt=markers(options.elside).tail;
        movedpt=movedheadtail(1,:);
        move='head';
    case 2 % upper contact
        olddist=abs(ea_pdist([markers(options.elside).head;markers(options.elside).tail]));
        cchangecoord=2;
        usefit=options.elside;
        if markers(options.elside).head(3) < markers(options.elside).tail(3)
            spin = 1;
        else % In case tail is lower than head
            spin = -1;
        end

        refpt=markers(options.elside).head;
        movedpt=movedheadtail(2,:);
        move='tail';
end

helppt=[movedpt(1:2),refpt(3)]; % this point is constructed as a projection of the movedpoint to the z-axis at height of refpt.
zdistfromhelppt=sqrt(abs((olddist)^2-(ea_pdist([helppt;refpt]))^2)); % use Pythagorean theorem to calculate distance on z-axis.
% here, 3*olddist is the hypotenuse, ea_pdist between the helper point and the
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

try
    xc=refpt+(nutraj*(olddist));
    numarkers(options.elside)=setfield(numarkers(options.elside),move,xc);
        set(mplot(cchangecoord,1),'xdata',xc(1));
        set(mplot(cchangecoord,1),'ydata',xc(2));
        set(mplot(cchangecoord,1),'zdata',xc(3));
catch
    keyboard
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


function movedel=whichelmoved(coordhandle,mcfig)
try
    mplot = getappdata(mcfig,'mplot');
    options = getappdata(mcfig,'options');
    for el=1:2
        if coordhandle==mplot(el,1)
            movedel=sub2ind([2,2],el,1);
        end
    end
catch
    keyboard
end
