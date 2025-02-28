function ea_c3d_reorient(input, output, orient, type)
% Reorient image to standard orthogonal orientation using c3d

arguments
    input       {mustBeFile}                % Input image
    output      {mustBeTextScalar} = ''     % Output image
    orient      {mustBeTextScalar} = 'RAS'  % Output image orientation
    type        {mustBeTextScalar} = ''     % Output data type
end

if isempty(output)
    output = input;
end

% Create a orientation mapping
orientMap = containers.Map({'R','L','A','P','S','I'}, ...
                           {'L','R','P','A','I','S'});

% Flip the orientation since c3d has the opposite definition 
orient = arrayfun(@(x) orientMap(x), upper(orient));

% Get c3d executable
basedir = [fileparts(mfilename('fullpath')), filesep];
c3d = ea_getExec([basedir, 'c3d'], escapePath = 1);

% Get output data type
if isempty(type) % Keep input data type in case not specified
    cmd = {c3d, ea_path_helper(input), '-info-full'};
    [~, header] = ea_runcmd(strjoin(cmd, ' '));
    datatype = regexp(header, '(?<=datatype = )\d+', 'match', 'once');
    switch datatype
        case '2'
            type = 'uchar';
        case '4'
            type = 'short';
        case '8'
            type = 'int';
        case '16'
            type = 'float';
        case '64'
            type = 'double';
        case '256'
            type = 'char';
        case '512'
            type = 'ushort';
        case '768'
            type = 'uint';
        otherwise % Use float as default data type
            type = 'float';
    end
else
    if ~ismember(type, {'char', 'uchar', 'short', 'ushort', 'int', 'uint', 'float', 'double'})
        ea_error('Output data type must be {char | uchar | short | ushort | int | uint | float | double}', showdlg=0, simpleStack=1);
    end
end

cmd = {c3d, ea_path_helper(input), '-orient', orient, '-type', type, '-o', ea_path_helper(output)};
ea_runcmd(cmd);
