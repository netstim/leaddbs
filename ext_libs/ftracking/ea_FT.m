function fibers=FT(FA,VectorF,Roi,parameters)
% This function will perform straight forward fiber tracking from the
% whole brain using linear interpolation
% Alper Yaman's contribution: Tracking initiates from ROI points -not from
% whole voxels. I changed a few lines and canceled textdisplay for
% performance purposes. 
%
% fibers=FT(FA,VectorF,Roi,parameters);
%
% FA: A 3D matrix with the fractional anistropy (FA) of a DTI dataset
% VectorF: A 3D matrix with the main fiber direction of each voxel
% Roi: A 3D matrix with the Roi. The fibers traced from the whole brain
%       which go through this Roi are kept.
% parameters: parameters.FiberLengthMax, parameters.FiberLengthMin
%       parameters.DeviationAngleMax, parameters.Step, parameters.Sampling,
%       parameters.FAmin and parameters.FAmax
%
% fibers: a cell array, with the coordinates of all fibers through the Roi
%
% Example:
%  see FT_example.m
%
% This function is written by D.Kroon University of Twente (August 2008)

if(parameters.textdisplay), disp('Start FT function'); pause(0.1); end

% Number of fibers found
fibern=0;

if(parameters.textdisplay), disp('Perform fiber tracking through the whole tissue'); end
if(parameters.Sampling>1)
    Roi2=zeros(size(Roi));Roi2(1:parameters.Sampling:end,1:parameters.Sampling:end,:)=1;
    Roibig=Roi.*Roi2;
    rbig=find(Roibig);
    [x,y,z]=ind2sub(size(Roibig),rbig);
else
    rbig=find(Roi);
    [x,y,z]=ind2sub(size(Roi),rbig);
end

% Loop through the voxel locations of the volume
% fiber is used to store the coordinates of one fiber
fiber=zeros(parameters.FiberLengthMax,3);
sz=size(FA);
for ix=1:length(x)%2:parameters.Sampling:(sz(1)-1),
    %     if(parameters.textdisplay),  disp(['fiber tracking process : ' num2str(round(100*x/(size(Roi,1)-1))) '%']); pause(0.1); end
    % First fiber coordinate is the current position
    fiber(1,:)=[x(ix) y(ix) z(ix)];

    % initialize variables
    check=true; f_length=1; f_Roi=false; old_gradient=[]; f_angle=0;
    while(check)
        % Linear interpolation of current location
        xBas=[0 0 0 0 1 1 1 1]+floor(fiber(f_length,1));
        yBas=[0 0 1 1 0 0 1 1]+floor(fiber(f_length,2));
        zBas=[0 1 0 1 0 1 0 1]+floor(fiber(f_length,3));
        xCom=fiber(f_length,1)-floor(fiber(f_length,1));
        yCom=fiber(f_length,2)-floor(fiber(f_length,2));
        zCom=fiber(f_length,3)-floor(fiber(f_length,3));
        % Linear interpolation percentages.
        perc=[(1-xCom) * (1-yCom) * (1-zCom); (1-xCom) * (1-yCom) * zCom;
            (1-xCom) * yCom * (1-zCom); (1-xCom) * yCom * zCom;
            xCom * (1-yCom) * (1-zCom); xCom * (1-yCom) * zCom;
            xCom * yCom * (1-zCom); xCom * yCom * zCom;];

        % The gradient / fiber direction is determined from a voxel
        % neighborhood of 8.
        gradient=[0 0 0]; f_FA=0;
        for i=1:8,
            gradient(1)=gradient(1)+VectorF(xBas(i),yBas(i),zBas(i),1)*perc(i);
            gradient(2)=gradient(2)+VectorF(xBas(i),yBas(i),zBas(i),2)*perc(i);
            gradient(3)=gradient(3)+VectorF(xBas(i),yBas(i),zBas(i),3)*perc(i);
            %gradient=gradient+(squeeze(VectorF(xBas(i),yBas(i),zBas(i),:))*perc(i))'; 
            % This line is canceled and 3 lines above are added for performance purposes. 
        end
        xx=(zBas-1)*sz(1)*sz(2)+(yBas-1)*sz(1)+xBas;
        f_FA=f_FA+FA(xx)*perc;%This line is taken from the for loop above for performance purposes. 
        gradient=gradient./(sqrt(sum(gradient.^2))+eps);

        % Set a step in the direction of the gradient.
        fiber(f_length+1,:)=fiber(f_length,:)+(parameters.Step) .* gradient;

        % Determine if this fiber crosses the Roi
        fiber_r=round(fiber(f_length,:));
        f_Roi=f_Roi||Roi(fiber_r(1),fiber_r(2),fiber_r(3));

        % Calculate angle between the current and last on the (fiber) position
        if(f_length>1),
            f_angle=abs(acos(sum(gradient.*old_gradient)));
        end

        % Do the Fiber stop checks:
        % Stop if the fiber becomes to long
        if(f_length>=parameters.FiberLengthMax), check=false; end
        % Stop if the fiber takes a hard turn.
        if(f_angle>parameters.DeviationAngleMax), check=false; end
        % Stop if the anistropy of the voxel is too low
        if(f_FA<parameters.FAmin | f_FA>parameters.FAmax), check=false; end
        % Check if the fiber will grow outside of the volume
        if((sum(fiber(f_length+1,:)>(size(Roi)-1))+sum(fiber(f_length+1,:)<2))>0), check=false; end

        % Update the fiber length
        if(check), f_length=f_length+1; end

        % Keep the gradient for next step angle change calculation
        old_gradient=gradient;
    end
    % Keep the fiber if it is long enough and crossed the Roi
    if((f_length>parameters.FiberLengthMin)&&f_Roi)
        fibern=fibern+1;
        fibers{fibern}=fiber(1:f_length,:);
    end
end
if(parameters.textdisplay), disp('FT function finished'); pause(0.1); end

