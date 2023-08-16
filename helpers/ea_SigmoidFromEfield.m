function Sigmoid = ea_SigmoidFromEfield(varargin)

    % map E-field to sigmoid function
    % by Till Dembek

    Efield = varargin{1};
    if nargin == 1
        minEfieldthresh = 0.061; % Astrom et al. 2014, 7.5µm, 1.5V
        maxEfieldthresh = 0.351; % Astrom et al. 2014, 2.5µm, 5V
    elseif nargin == 3
        minEfieldthresh = varargin{2};
        maxEfieldthresh = varargin{3};
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
    steepness = 20.3 ./ multiplicationfactor;

    %% reestimate steepness
%     sn =  [1:0.1:100]; %* multiplicationfactor;
%     for k = 1:length(sn)
%        ymin(k) = maxval./(1+exp(-sn(k) .* (minEfieldthresh - x0)));
%        ymax(k) = maxval./(1+exp(-sn(k) .* (maxEfieldthresh - x0)));
%     end
%     [~,tmpind1] = min(abs(ymin-0.05));
%     [~,tmpind2] = min(abs(ymax-0.95));
%     steepness = sn(round(mean([tmpind1 tmpind2])));    
%     x = [0:0.01:0.5] .* multiplicationfactor;
%     y = maxval./(1+exp(-steepness .* (x - x0)));
%     plot(x,y)    
%     ylim([0 1])

    Sigmoid = maxval./(1+exp(-steepness .* (Efield - x0)));
end