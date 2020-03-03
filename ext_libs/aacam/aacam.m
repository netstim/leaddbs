function [blurim] = aacam(focusmode,param1,param2,iparam,updown)
%AACAM, A Matlab Camera Model for Depth of Field
%   blurim = AACAM(focusmode, param1, param2, iparam, updown)
%   Simulate depth of field in the current figure and axis, return an RGB
%   image.
%
%   focusmode = 'LR'
%     param1 is L (focal length) and param2 is R (aperture radius)
%   focusmode = 'L1L2'
%     objects at perp. distances within [param1, param2] are sharp
%   focusmode = 'P1P2'
%     objects at perp. distances witin the 3D points param1 and param2 are
%     sharp, useful if you can select two points in the scene that you wish
%     to see clearly.
%   focusmode = 'P1R'
%     the focual length is determined by the 3D point param1 and the radius
%     of the aperture is param2, useful if you know a single point that
%     should be sharp
%
%   iparam = [N, stderr, maxerr], i.e. the simulation is terminated after N
%   iterations, or when standard deviation of the change of the simulated
%   image is less than stderr, or when the max change of the simulated
%   image is less than maxerr.
%
%   updown = [Nup, Ndown], e.g. [4 2] to make all calculations to be
%   upsampled 4 times and then reduced 2 times in the returned image
%   blurim. The fastest option is [1 1] but often [2 2]?or [2 1] are better
%   choices.
%
%   Example:
%     X = randn(1000,3)*1.5;
%     scatter3(X(:,1),X(:,2),X(:,3),4,[1 1 1],'filled'); hold on;
%     t = [0:0.0001:2]; 
%     h = plot3(1.4*cos(t*164),1.4*sin(2+t*131), 1.4*sin(154*t)); 
%     set(h,'Color',[0 1 1]);
%     camproj('perspective'); daspect([1 1 1]);
%     camtarget([0,0,0]); campos([9,9,9]); camva(10);
%     set(gcf,'renderer','opengl');
%     set(gca,'Color',0*[0.9 0.9 0.9]);
%     sphere(20);
%     lighting phong; light; drawnow;
%     [blurim] = aacam('P1R',[1 1 1]/sqrt(3),0.5,[100 0.2 1],[2 2]);
%     figure; image(blurim); set(gca,'Position',[0 0 1 1]); axis equal;
%
%   Author: Anders Brun
%           anders@cb.uu.se
%


K = updown(1);
M = updown(2);
N = iparam(1);

% Some things we have to assume or arrange for ....
camproj('perspective');
axis vis3d;

% Make a copy of the current figure
saveas(gcf,'temp.fig');
FIG = openfig('temp.fig');
delete('temp.fig');
set(FIG,'visible','off'); % Turn on to reveal magic
AXE = gca;

% Find out current camera position and orientation
da = daspect;
cpos = campos./da;
ctar = camtarget./da;
cva = camva;

% Get axis position in pixels (used to compensate later)
set(AXE,'Units','pixels');
set(gcf,'Units','pixels');
caxpos = get(AXE,'Position');
screen_dx = caxpos(3);
screen_dy = caxpos(4);

% Get focal length F for this axis
F = min(screen_dx,screen_dy)/2/tan(cva/360*2*pi/2);

% Calculate real up vector
cdir = (ctar-cpos)/norm(ctar-cpos);
cup = camup; cup = cross(cross(cdir,cup),cdir); cup = cup / norm(cup);
cright = cross(cdir,cup);

% Set L and R according to the selected "focus mode"
if isequal(focusmode,'LR')
    L = param1;
    R = param2;
elseif isequal(focusmode,'L1L2')
    L1 = param1;
    L2 = param2;
    L = 2*L2*L1/(L2+L1);
    R = L2*L1/(F*(L2-L1));
elseif isequal(focusmode,'P1P2')
    P1 = param1 ./ da;
    P2 = param2 ./ da;
    mcnear = (P1-cpos)*(ctar-cpos)' / norm(ctar-cpos)^2 * (ctar-cpos) + cpos;
    mctar = (P2-cpos)*(ctar-cpos)' / norm(ctar-cpos)^2 * (ctar-cpos) + cpos;
    L2 = norm(cpos - mctar);
    L1 = norm(cpos - mcnear);
    if L2 < L1
        p = L1; L1 = L2; L2 = p;
    end
    L = 2*L2*L1/(L2+L1);
    R = L2*L1/(F*(L2-L1));
elseif isequal(focusmode,'P1R')
    P1 = param1 ./ da;
    R = param2;
    mctar = (P1-cpos)*(ctar-cpos)' / norm(ctar-cpos)^2 * (ctar-cpos) + cpos;
    L = norm(cpos - mctar);
end

% Create MM, a Nmask x Nmask model of our aperture (circular for now)
% ii and jj are also used later
Nmask = 100;
[ii,jj] = ndgrid(-Nmask:Nmask,-Nmask:Nmask);
ii = ii/Nmask; jj = jj / Nmask;
MM = ii.^2 + jj.^2 <= 1;

% Sample N position (indices) from MM, our aperture model
idxs = mysamplemask(MM,N);

C = 0;
blurim = 0;
fprintf('Keep the camera steady during exposure... :-)\n');
for k = 1:N
    fprintf('%2.2f%% done\n',100*k/N); 
    di = ii(idxs(k));
    dj = jj(idxs(k));
    % Integrate over this aperture position
    C = C + 1;
    odx = di*R;
    ody = dj*R;
    odelta = cright*odx + cup*ody;
    % Move camera theta radians in alpha direction
    campos((cpos + odelta).*da);
    camtarget((ctar + odelta).*da);
    % Move axis to compensate for object focal length
    ax = F*odx/L;
    ay = F*ody/L;
    set(AXE,'Position',caxpos + [+ax +ay 0 0]);
    % Save image
    tempfile = 'myaa_temp_screendump.png';
    set(FIG,'PaperPositionMode','auto');
    set(FIG,'InvertHardcopy','off');
    screen_DPI = get(0,'ScreenPixelsPerInch');
    drawnow;
    print(FIG,['-r',num2str(screen_DPI*K)], '-dpng', tempfile);
    p = imread(tempfile);
    delete(tempfile);
    
    % HSV enhance bright stuff like stars ...
    % I commented this out, but add your own gamma function here if you
    % like to make things like white snow flakes or highlights to appear
    % brighter when they are out of focus.
    
    %[H,S,V] = rgb2hsv(p);
    %p = single(p).*(1 + repmat(single((V==1) & (H==0)),[1 1 3])*5);
    p = single(p);
    % Compute "errors", i.e. how much the image changed this iteration
    err = single(imresize((blurim + p)/(C)-blurim/(C-1),1/M));
    stderr = std(err(:));
    maxerr = max(abs(err(:)));
    
    blurim = blurim + p;
    
    % Now break if any of the error criteria are fulfilled
    if (stderr < iparam(2) || maxerr < iparam(3)) && C > 2
        break;        
    end    
    % Otherwise, just keep iterating up to N samples from the aperture
end

% Downsample and correct strange image values
blurim = single(blurim)/C/256;
blurim = imresize(blurim,1/M);
blurim(blurim>1) = 1;
blurim(blurim<0) = 0;

% POP Current camera / axis parameters
close(FIG);



function p = mysamplemask(MM,N)
% Simple farthest point sampling of a binary mask MM, return N indices to
% the bitmap. Should be better than regular sampling and random sampling,
% but after all it is a matter of taste...

p = zeros(0,1);
isize = size(MM,1);
jsize = size(MM,2);
M = prod(size(MM));

[ii,jj] = ndgrid(1:isize,1:jsize);

D = -inf(isize,jsize);
D(MM==1) = +inf;

for k = 1:N
    maxidx = sortrows([-D(:) (1:M)']);
    maxidx = maxidx(1,2);
    p(end+1) = maxidx;
    D = min(D,(ii-ii(maxidx)).^2 + (jj-jj(maxidx)).^2);
end






