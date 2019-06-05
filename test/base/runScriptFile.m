function result = runScriptFile(fileName)
% This function runs the test in fileName
% It can distinguish between skipped and Failed tests. A test is considered
% to be skipped if it throws a COBRA:RequirementsNotMet error.
%
% OUTPUTS:
%
%    result:    A structure array with the following fields:
%
%                - `.passed`: true if the test passed otherwise false
%                - `.skipped`: true if the test was skipped otherwise false
%                - `.failed`: true if the test failed, or was skipped, otherwise false
%                - `.status`: a string representing the status of the test ('failed','skipped' or'passed')
%                - `.fileName`: the fileName of the test
%                - `.time`: the duration of the test (if passed otherwise NaN)
%                - `.statusMessage`: Informative string about potential problems.
%                - `.Error`: Error message received from a failed or skipped test
%
% .. Author:   - Thomas Pfau Jan 2018.
%
% NOTE:
%    this function has been adapted from
%    https://github.com/opencobra/cobratoolbox/blob/master/test/runScriptFile.m
%    [b58da0c]

global CBT_MISSING_REQUIREMENTS_ERROR_ID

LEADDBS_TESTSUITE_TESTFILE = fileName;

% get the timinig (and hope these values are not overwritten.
LEADDBS_TESTSUITE_STARTTIME = clock();

try
    % run the file
    executefile(fileName);
catch ME
    % vatch errors and interpret them
    clearvars -except ME LEADDBS_TESTSUITE_STARTTIME LEADDBS_TESTSUITE_TESTFILE CBT_MISSING_REQUIREMENTS_ERROR_ID
    scriptTime = etime(clock(), LEADDBS_TESTSUITE_STARTTIME);
    result = struct('status', 'failed', 'failed', true, 'passed', false, 'skipped', false, 'fileName', ...
                    LEADDBS_TESTSUITE_TESTFILE, 'time', scriptTime, 'statusMessage', 'fail', 'Error', ME);
    if strcmp(ME.identifier, CBT_MISSING_REQUIREMENTS_ERROR_ID)
        % requirement missing, so the test was skipped.
        result.status = 'skipped';
        result.skipped = true;
        result.failed = false;
        result.statusMessage = ME.message;
    else
        % actual error in the test.
        result.skipped = false;
        result.status = 'failed';
        result.statusMessage = ME.message;
    end
    return
end

% get the timinig.
scriptTime = etime(clock(), LEADDBS_TESTSUITE_STARTTIME);

result = struct('status', 'passed', 'failed', false, 'passed', true, 'skipped', false, 'fileName', ...
                LEADDBS_TESTSUITE_TESTFILE, 'time', scriptTime, 'statusMessage', 'success', 'Error', MException('', ''));

end

function executefile(fileName)
% runs a script file (used to separate workspaces)
    run(fileName)
end
