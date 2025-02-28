function fiberFiltered = ea_filterfiber_len(ftr, minFibLen)
% Filter fibers based on the axon length threshold chosen in the settings GUI
% By Butenko and Li, konstantinmgtu@gmail.com

arguments
    ftr         % Lead-DBS structure for fibers
    minFibLen   {mustBeNumeric}  % minimal length in mm
end

if ~iscell(ftr)
    ftr = {ftr};
end

fiberFiltered = cell(size(ftr));
for i=1:length(fiberFiltered)
    fiberFiltered{i}.fibers = [];
    fiberFiltered{i}.idx = [];
end

% Inline function to calculate the length of the fiber segments
sumlength = @(x) sum(sqrt(sum(diff(x).^2,2)));

disp('Removing short fibers...');
for i=1:length(ftr)
    if isempty(ftr{i}) || isempty(ftr{i}.idx)
        disp('No fibers found, skipping...');
    else
        % Prepare stard and end index of each fiber
        endIndex = cumsum(ftr{i}.idx);
        startIndex = [1;endIndex+1];
        startIndex = startIndex(1:end-1);
        fibIndex = arrayfun(@(x,y) (x:y)', startIndex, endIndex, 'Uni', 0);
        
        % Calculate numbers of points to keep
        fibLen = cellfun(@(x) sumlength(ftr{i}.fibers(x,1:3)), fibIndex);
        disp([num2str(sum(fibLen>=minFibLen)), ' out of ', num2str(length(ftr{i}.idx)),' fibers retained...'])
        fiberFiltered{i}.fibers = ftr{i}.fibers(cell2mat(fibIndex(fibLen>=minFibLen)), :);
        fiberFiltered{i}.idx = ftr{i}.idx(fibLen>=minFibLen);
    end
end

fprintf('\nFinished!\n\n');
