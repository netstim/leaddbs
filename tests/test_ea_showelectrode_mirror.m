function tests = test_ea_showelectrode_mirror
% Rendering must carry active/passive contact selection across hemispheres.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repo = fileparts(fileparts(mfilename('fullpath')));
oldPath = path;
addpath(repo,fullfile(repo,'helpers'),fullfile(repo,'templates','electrode_models'));
testCase.addTeardown(@() path(oldPath));
end

function testPlainActiveContactsUseContralateralSelection(testCase)
verifySelection(testCase,'2d',true);
end

function testSphereActiveContactsUseContralateralSelection(testCase)
verifySelection(testCase,'3d',true);
end

function testPlainPassiveContactsUseContralateralSelection(testCase)
verifySelection(testCase,'2d',false);
end

function testSpherePassiveContactsUseContralateralSelection(testCase)
verifySelection(testCase,'3d',false);
end

function verifySelection(testCase,style,showActive)
for mirrored = [false true]
    [el,options] = fixture(mirrored,style,showActive);
    for side = 1:2
        fig = figure('Visible','off');
        cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
        ax = axes('Parent',fig);
        hold(ax,'on');
        options.sides = side;
        obj = struct('plotFigureH',fig,'elstruct',el,'pt',1);
        rendered = ea_showelectrode(obj,'dbs',options);

        % The fixture mirrors the coordinate arrays but preserves the source
        % active-contact arrays, matching ea_mirrorsides' established API.
        sourceSide = side;
        if mirrored
            sourceSide = 3-side;
        end
        selected = logical(el.activecontacts{sourceSide});
        if ~showActive
            selected = ~selected;
        end
        expected = el.coords_mm{side}(selected,:);
        verifyEqual(testCase,numel(rendered),size(expected,1));
        actual = zeros(numel(rendered),3);
        for k = 1:numel(rendered)
            h = rendered(k);
            if strcmp(style,'3d')
                actual(k,:) = [midrange(h.XData),midrange(h.YData),midrange(h.ZData)];
                faceColor = h.FaceColor;
            else
                actual(k,:) = [h.XData,h.YData,h.ZData];
                faceColor = h.MarkerFaceColor;
            end
            if showActive
                verifyEqual(testCase,faceColor,[0.9 0.1 0.1],'AbsTol',1e-12);
            else
                verifyEqual(testCase,faceColor,el.groupcolors(1,:),'AbsTol',1e-12);
            end
        end
        verifyEqual(testCase,actual,expected,'AbsTol',1e-10);
        clear cleanup;
    end
end
end

function value = midrange(data)
value = (min(data,[],'all')+max(data,[],'all'))/2;
end

function [el,options] = fixture(mirrored,style,showActive)
el.elmodel = 'Medtronic 3389';
el.name = 'synthetic';
el.group = 1;
el.groupcolors = [0.2 0.4 0.6];
el.activecontacts = {[1;0;0;0],[0;1;0;1]};
contact = (1:4)';
el.coords_mm = {[12+contact/10,1+contact/5,3+2*contact], ...
               [-18-contact/10,2+contact/5,4+2*contact]};
if mirrored
    el.name = 'synthetic_mirrored';
    el.coords_mm = el.coords_mm([2 1]);
    for side = 1:2
        el.coords_mm{side}(:,1) = -el.coords_mm{side}(:,1);
    end
end
for side = 1:2
    el.trajectory{side} = el.coords_mm{side}([1 end],:);
end
options.d3.elrendering = 3;
options.d3.mirrorsides = double(mirrored);
options.d3.pntcmap = parula(64);
options.d3.prolong_electrode = 0;
options.d3.isomatrix = {};
options.d3.colorpointcloud = 0;
options.d3.hlactivecontacts = 1;
options.d3.showactivecontacts = double(showActive);
options.d3.showpassivecontacts = double(~showActive);
options.prefs.d3.pointcloudstyle = style;
end
