function ea_ants_applytransforms_to_points(input, output, transform)
% Wrapper for antsApplyTransformsToPoints

ea_libs_helper;

basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    applyTransformsToPoints = [basedir, 'antsApplyTransformsToPoints.exe'];
elseif isunix
    applyTransformsToPoints = [basedir, 'antsApplyTransformsToPoints.', computer];
end

cmd = [applyTransformsToPoints, ...
    ' --dimensionality 2' ...   % dimensionality
    ' --precision 1' ...    % double precision
    ' --input ', input ...  % input csv file with x,y,z,t (at least) as the column header
    ' --output ', output ...    % warped output csv file
    ' --transform [',transform, ',1]'];  % [transformFileName,useInverse]
    
if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end

end
