function [coords,trajectory,markers]=ea_resolvecoords(varargin)
% This function inflates the 2 markers reconstructed in
% ea_reconstruct_trajectory to the contacts as specified by the lead model.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn
markers=varargin{1};
options=varargin{2};

if nargin>2
    resize=varargin{3};
else
    resize=0;
end
if nargin==4
    rszfactor=varargin{4};
end

load([ea_getearoot,'templates',filesep,'electrode_models',filesep,options.elspec.matfname]);
for side=1:length(markers) %valid for unilateral support
    if resize
        can_dist=ea_pdist([electrode.head_position;electrode.tail_position]);
        %emp_dist=ea_pdist([markers(side).head;markers(side).tail]);
        %A=squareform(pdist(electrode.coords_mm));
        switch options.elmodel
            case {'Medtronic B33005'
                  'Medtronic B33015'
                  'Boston Scientific Vercise Directed'
                  'Abbott Directed 6172 (short)'
                  'Abbott Directed 6173 (long)'}
                coords_temp(1,:) = electrode.coords_mm(1,:);
                coords_temp(2,:) = mean(electrode.coords_mm(2:4,:));
                coords_temp(3,:) = mean(electrode.coords_mm(5:7,:));
                coords_temp(4,:) = electrode.coords_mm(8,:);
                A=sqrt(ea_sqdist(coords_temp',coords_temp'));
                can_eldist=sum(sum(tril(triu(A,1),1)))/(3);
                clear coords_temp
            case {'Boston Scientific Vercise Cartesia HX'
                  'Boston Scientific Vercise Cartesia X'}
                coords_temp(1,:) = mean(electrode.coords_mm(1:3,:));
                coords_temp(2,:) = mean(electrode.coords_mm(4:6,:));
                coords_temp(3,:) = mean(electrode.coords_mm(7:9,:));
                coords_temp(4,:) = mean(electrode.coords_mm(10:12,:));
                A=sqrt(ea_sqdist(coords_temp',coords_temp'));
                can_eldist=sum(sum(tril(triu(A,1),1)))/(3);
                clear coords_temp
            otherwise
                A=sqrt(ea_sqdist(electrode.coords_mm',electrode.coords_mm'));
                can_eldist=sum(sum(tril(triu(A,1),1)))/(options.elspec.numel-1);
        end
        vec=(markers(side).tail-markers(side).head)/norm(markers(side).tail-markers(side).head);
        if nargin>3
            stretch=can_dist*(rszfactor/can_eldist);
        else
            stretch=can_dist;
        end
        markers(side).tail=markers(side).head+vec*stretch;
    end

    if ~isempty(markers(side).head)
        M=[markers(side).head,1;markers(side).tail,1;markers(side).x,1;markers(side).y,1];
        E=[electrode.head_position,1;electrode.tail_position,1;electrode.x_position,1;electrode.y_position,1];
        X=mldivide(E,M);

        coords_mm=[electrode.coords_mm,ones(size(electrode.coords_mm,1),1)];
        coords{side}=X'*coords_mm';
        coords{side}=coords{side}(1:3,:)';

        trajvector{side}=(markers(side).tail-markers(side).head)/norm(markers(side).tail-markers(side).head);
        trajectory{side}=[markers(side).head-trajvector{side}*5;markers(side).head+trajvector{side}*25];
        trajectory{side}=[linspace(trajectory{side}(1,1),trajectory{side}(2,1),50)',...
            linspace(trajectory{side}(1,2),trajectory{side}(2,2),50)',...
            linspace(trajectory{side}(1,3),trajectory{side}(2,3),50)'];
    end
end
