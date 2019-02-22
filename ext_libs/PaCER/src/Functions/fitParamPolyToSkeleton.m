%% Fits a parameterized polynomial [0,1] -> R^3 to a skeleton point cloud
% Params: 
%              skeleton - 3 x N skeleton pointcloud
%       OPTIONAL degree - degree of fitted polynomial (DEFAULT: 3)
%
% Andreas Husch
% Centre Hospitalier de Luxembourg / Luxembourg Centre for Systems
% Biomedicine, University of Luxembourg
% 2014  - 2017
% mail@andreashusch.de, husch.andreas@chl.lu

function [r3polynomial, avgTperMm] = fitParamPolyToSkeleton(skeleton, degree)
if(nargin < 2)
    disp('No degree given, using default degree=3');
    degree=3;
end
%% Determine (approx.) t (i.e. LHS of parameterized regressor); (arg length to t) 
diffVecs = diff(skeleton);
apprTotalLengthMm = 0;
cumLengthMm = zeros(length(diffVecs),1);
deltas = zeros(length(diffVecs),1);
for i=1:length(diffVecs) % TODO express as matrix mult. (norms for every diffVec)
   deltas(i) = norm(diffVecs(i,:));
   cumLengthMm(i) = sum(deltas(1:i));
   apprTotalLengthMm = apprTotalLengthMm + deltas(i); 
end
avgTperMm = 1 / apprTotalLengthMm;
t = cumLengthMm ./ apprTotalLengthMm; % now range t=[0 1]
t = [0;t]; % same no of gaps than points 

%% Regress 
% Construct Design Matrix
% e.g. T = [t.^4 t.^3 t.^2 t ones(length(t),1)]'; %  (fliplr...)
T = [ones(length(t),1)]; %#ok<NBRAK>

for d=1:degree
    T = [t.^d T]; %#ok<AGROW>
end
T = T';

%display(['Optimal AIC at degree ' num2str(polydeg(t, skeleton(:,1))) ' for first (x) direction']);
%display(['Optimal AIC at degree ' num2str(polydeg(t, skeleton(:,2))) ' for second (y) direction']);
%display(['Optimal AIC at degree ' num2str(polydeg(t, skeleton(:,3))) ' for third (z) direction']);

r3polynomial =  T' \ skeleton; % = skeleton' * pinv(T))' =(skeleton'*T'*inv

% Asserts
fittingErrs = sqrt(sum((r3polynomial' * T - skeleton').^2));
meanFittingError = mean(fittingErrs);
stdFittingError = std(fittingErrs);
maxFittingError = max(fittingErrs);

fprintf(['Max off-model: ' num2str(maxFittingError) ' Mean off-model:  ' num2str(meanFittingError) '\n']);

if(maxFittingError > 0.3 && maxFittingError > (meanFittingError + 3*stdFittingError))
       fprintf('Check for outliers / make sure that the polynomial degree choosen is appropriate.\n In most cases this should be fine.\n');
end

%figure, plot(sum((r3polynomial' * T - skeleton').^2))
end