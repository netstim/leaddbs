function [X,electrode,err]=ea_mapelmodel2reco(options,elspec,elstruct,side,resultfig)
err=0;

load([ea_getearoot,'templates',filesep,'electrode_models',filesep,elspec.matfname,'.mat'],'electrode')
A = [electrode.head_position,1;
     electrode.tail_position,1
     electrode.x_position,1
     electrode.y_position,1]; % points in model
redomarkers=0;

if ~isfield(elstruct,'markers') % backward compatibility to old electrode format
    redomarkers=1;
else
    if isempty(elstruct.markers)
        redomarkers=1;
    end
end

if redomarkers
    for iside=1:length(options.sides)
        side=options.sides(iside);

        elstruct.markers(side).head=elstruct.coords_mm{side}(1,:);

        switch options.elmodel
            case {'Medtronic B33005'
                  'Medtronic B33015'
                  'Boston Scientific Vercise Directed'
                  'Abbott Directed 6172 (short)'
                  'Abbott Directed 6173 (long)'}
                elstruct.markers(side).tail=elstruct.coords_mm{side}(8,:);
            case {'Boston Scientific Vercise Cartesia HX'
                  'Boston Scientific Vercise Cartesia X'}
                elstruct.markers(side).tail=mean(elstruct.coords_mm{side}(10:12,:));
            otherwise
                elstruct.markers(side).tail=elstruct.coords_mm{side}(4,:);
        end

        [xunitv, yunitv] = ea_calcxy(elstruct.markers(side).head, elstruct.markers(side).tail);
        elstruct.markers(side).x = elstruct.coords_mm{side}(1,:) + xunitv*(options.elspec.lead_diameter/2);
        elstruct.markers(side).y = elstruct.coords_mm{side}(1,:) + yunitv*(options.elspec.lead_diameter/2);
    end
end

B = [elstruct.markers(side).head,1;
     elstruct.markers(side).tail,1;
     elstruct.markers(side).x,1;
     elstruct.markers(side).y,1];
setappdata(resultfig,'elstruct',elstruct);
setappdata(resultfig,'elspec',elspec);

X = mldivide(A,B);

% Check precision using native coords instead of scrf coords
if options.native && exist([options.root,options.patientname,filesep,'scrf',filesep,'scrf_converted.mat'],'file')
    load([options.root,options.patientname,filesep,'scrf',filesep,'scrf_converted.mat'], 'mat');
    Btest = B/mat';
    Xtest = mldivide(A,Btest);
elseif options.native && exist([options.root,options.patientname,filesep,'scrf',filesep,'scrf.mat'],'file') % legacy
    mat = ea_getscrfmat(options);
    save([directory,'scrf',filesep,'scrf_converted.mat'],'mat');
    Btest = B/mat';
    Xtest = mldivide(A,Btest);
else
    Btest = B;
    Xtest = X;
end

% perform tests if A has been transformed to B correctly.
% First we will make sure that the projection from A to B (Ab)
% results in something very similar to B.
Ab=A*Xtest;
vizprec=0.01;
if sum(abs(Ab(:)-Btest(:)))>vizprec % visualization precision in sums of millimeters.
    err=1;
end

% Second we will make sure that projected markers structure is
% still close to orthogonal:
% figure, plot3(Ab(1,1),Ab(1,2),Ab(1,3),'r*');
% hold on
% plot3(Ab(2,1),Ab(2,2),Ab(2,3),'b*');
% plot3(Ab(3,1),Ab(3,2),Ab(3,3),'k*');
% plot3(Ab(4,1),Ab(4,2),Ab(4,3),'g*');
% axis square
angvizprec=0.5; % angular visualization precision tolerance in sums of degrees
if 90-rad2deg(acos(dot(...
        (Ab(1,1:3)-Ab(2,1:3))/...
        norm(Ab(1,1:3)-Ab(2,1:3)),...
        (Ab(1,1:3)-Ab(3,1:3))/...
        norm((Ab(1,1:3)-Ab(3,1:3)))...
        )+...
        dot(...
        (Ab(1,1:3)-Ab(2,1:3))/...
        norm(Ab(1,1:3)-Ab(2,1:3)),...
        (Ab(1,1:3)-Ab(4,1:3))/...
        norm(Ab(1,1:3)-Ab(4,1:3))) ...
        )) > angvizprec
    err=1;
end
% end tests

X=X';
