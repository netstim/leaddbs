function Sigmoid = ea_SigmoidFromEfield(varargin)
vizz=0;
% map E-field to sigmoid function
% by Till Dembek

Efield = varargin{1};
if nargin == 1
    minEfieldthresh = 0.061; % Astrom et al. 2014, 7.5µm, 1.5V
    maxEfieldthresh = 0.351; % Astrom et al. 2014, 2.5µm, 5V

    reestimatesteepness = 0;
elseif nargin == 3
    minEfieldthresh = varargin{2};
    maxEfieldthresh = varargin{3};
    reestimatesteepness = 1;
else
    error('Wrong number of input arguments!')
end


maxval = 1;
%% recalculate to Lead-DBS units (V/m)
multiplicationfactor = 1000;
minEfieldthresh = minEfieldthresh .* multiplicationfactor;
maxEfieldthresh = maxEfieldthresh .* multiplicationfactor;

%%
x0 = (minEfieldthresh + maxEfieldthresh) ./ 2;

if reestimatesteepness
    % reestimate steepness
    sn =  [1:0.1:100] ./ multiplicationfactor;
    for k = 1:length(sn)
        ymin(k) = maxval./(1+exp(-sn(k) .* (minEfieldthresh - x0)));
        ymax(k) = maxval./(1+exp(-sn(k) .* (maxEfieldthresh - x0)));
    end
    [~,tmpind1] = min(abs(ymin-0.05));
    [~,tmpind2] = min(abs(ymax-0.95));
    steepness = sn(round(mean([tmpind1 tmpind2])));
else

    steepness = 20.3 ./ multiplicationfactor;
end

if vizz
    figure
    disp(['Steepness: ' num2str(round(steepness.*multiplicationfactor,1))])
    x = [0:0.01:0.5] .* multiplicationfactor;
    y = maxval./(1+exp(-steepness .* (x - x0)));
    %     figure
    plot(x,y)
    hold on
    line([x(1),x(end)],[0.05,0.05])
    line([x(1),x(end)],[0.95,0.95])
    line([minEfieldthresh,minEfieldthresh],[0,1])
    line([maxEfieldthresh,maxEfieldthresh],[0,1])
    ylim([0 1])
end

Sigmoid = maxval./(1+exp(-steepness .* (Efield - x0)));
end
