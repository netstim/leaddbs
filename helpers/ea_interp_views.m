function ov=ea_interp_views(v,interstages,linear)

if ~exist('linear','var')
    linear=0;
end

if linear
cnt=1;
for stage=1:length(v)-1
    
    az=linspace(v(stage).az,v(stage+1).az,interstages);
    el=linspace(v(stage).el,v(stage+1).el,interstages);
    camva=linspace(v(stage).camva,v(stage+1).camva,interstages);
    camup=[linspace(v(stage).camup(1),v(stage+1).camup(1),interstages);...
        linspace(v(stage).camup(2),v(stage+1).camup(2),interstages);...
        linspace(v(stage).camup(3),v(stage+1).camup(3),interstages)];
    camtarget=[linspace(v(stage).camtarget(1),v(stage+1).camtarget(1),interstages);...
        linspace(v(stage).camtarget(2),v(stage+1).camtarget(2),interstages);...
        linspace(v(stage).camtarget(3),v(stage+1).camtarget(3),interstages)];
    campos=[linspace(v(stage).campos(1),v(stage+1).campos(1),interstages);...
        linspace(v(stage).campos(2),v(stage+1).campos(2),interstages);...
        linspace(v(stage).campos(3),v(stage+1).campos(3),interstages)];

    for interframe=1:interstages
        ov(cnt).az=az(interframe);
        ov(cnt).el=el(interframe);
        ov(cnt).camva=camva(interframe);
        ov(cnt).camup=camup(:,interframe)';
        ov(cnt).camtarget=camtarget(:,interframe)';
        ov(cnt).campos=campos(:,interframe)';
        ov(cnt).camproj=v(stage).camproj;
        cnt=cnt+1;
    end
end

else % log in and out
    cnt=1;
    for stage=1:length(v)-1
        
        az=nonLinspace(v(stage).az,v(stage+1).az,interstages,'cos');
        el=nonLinspace(v(stage).el,v(stage+1).el,interstages,'cos');
        camva=nonLinspace(v(stage).camva,v(stage+1).camva,interstages,'cos');
        camup=[nonLinspace(v(stage).camup(1),v(stage+1).camup(1),interstages,'cos');...
            nonLinspace(v(stage).camup(2),v(stage+1).camup(2),interstages,'cos');...
            nonLinspace(v(stage).camup(3),v(stage+1).camup(3),interstages,'cos')];
        camtarget=[nonLinspace(v(stage).camtarget(1),v(stage+1).camtarget(1),interstages,'cos');...
            nonLinspace(v(stage).camtarget(2),v(stage+1).camtarget(2),interstages,'cos');...
            nonLinspace(v(stage).camtarget(3),v(stage+1).camtarget(3),interstages,'cos')];
        campos=[nonLinspace(v(stage).campos(1),v(stage+1).campos(1),interstages,'cos');...
            nonLinspace(v(stage).campos(2),v(stage+1).campos(2),interstages,'cos');...
            nonLinspace(v(stage).campos(3),v(stage+1).campos(3),interstages,'cos')];
        
        for interframe=1:interstages
            ov(cnt).az=az(interframe);
            ov(cnt).el=el(interframe);
            ov(cnt).camva=camva(interframe);
            ov(cnt).camup=camup(:,interframe)';
            ov(cnt).camtarget=camtarget(:,interframe)';
            ov(cnt).campos=campos(:,interframe)';
            ov(cnt).camproj=v(stage).camproj;
            cnt=cnt+1;
        end
    end
    
    
end


% -------------------------------------------------------------------------
% nonLinspace(mn, mx, num, spacetype) returns a vector of non-linearly 
% spaced elements based on spacing specified by spacetype. 
%
% nonLinVec = nonLinspace(mn, mx, num, 'exp10') returns a vector of
% elements with smaller spacing at the beginning of the vector and greater
% spacing at the end of the vector based on the curve y = 10^x.
%
% nonLinVec = nonLinspace(mn, mx, num, 'cos') returns a vector of elements
% with smaller spacing at the beginning and end of the vector, and greater
% spacing in the middle based on the curve y = 1/2(1-cos(x)).
%
% nonLinVec = nonLinspace(mn, mx, num, 'log10') returns a vector of
% elements with greater spacing at the beginning of the vector and smaller
% spacing at the end of the vector. 
% 
%   Inputs: 
%       mn        - The minimum value in the vector. 
%       mx        - The maximum value in the vector.
%       num       - The number of elements in the vector. 
%       spacetype - Specifies the type of spacing needed. 
%
%   Outputs:
%       nonLinVec - A vector consisting of elements with spacing specified 
%                   by spacetype.
%
%
% Created: 10/12/17 - Connor Ott
% Last Modified: 10/23/17 - Connor Ott
% -------------------------------------------------------------------------
function [ nonLinVec ] = nonLinspace( mn, mx, num, spacetype )
if strcmpi(spacetype, 'exp10')
    % exponentialliness is the upper bound of the original 10^x curve
    % before it is scaled to fit the limits requested by the user. Since
    % the concavity of 10^x changes in different parts of its domain,
    % different spacing is seen when using different bounds. After some
    % basic qualitative analysis, an exponentialliness of 20 seemed to be a
    % good fit for my purposes. Increasing this value will increase the
    % spacing towards the end of the vector and decrease it towards the
    % beginning. 
    exponentialliness = 20;
    nonLinVec = (mx-mn)/exponentialliness * ...
                (10.^(linspace(0, log10(exponentialliness+1), num)) - 1)...
                + mn;
            
elseif strcmpi(spacetype, 'cos')
    nonLinVec = (mx - mn)*(0.5*(1-cos(linspace(0, pi, num)))) + mn;
    
elseif strcmpi(spacetype, 'log10')
    % As with exponentialliness, this defines the bounds on the log10(x)
    % curve. Increasing loginess will decreasing the spacing towards the
    % end of the vector and increase it towards the beginning. 
    loginess = 1.5;
    nonLinVec = (mx - mn)/loginess* ...
                log10((linspace(0, 10^(loginess) - 1, num)+ 1)) + mn;
            
end
   
