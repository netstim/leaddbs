function ispresent = ea_refreshscrf(options,handles)

ispresent = isfile(options.subj.brainshift.transform.instore);

if ispresent
    stack = dbstack;
    if ismember('ea_compute_scrf', {stack.name}) ...
            || ~all(isfile(struct2cell(options.subj.brainshift.checkreg)))
        % Always regenerate checkreg figure after recomputing
        ea_gencheckregfigs(options, 'brainshift');
    end
end

try
    standardslice = imread(options.subj.brainshift.checkreg.standard);
catch
    standardslice = imread([ea_getearoot,'helpers',filesep,'gui',filesep,'scrf_msg.png']);
end

try
    refineslice = imread(options.subj.brainshift.checkreg.scrf);
catch
    refineslice = imread([ea_getearoot,'helpers',filesep,'gui',filesep,'scrf_msg.png']);
end

set(0,'CurrentFigure',handles.scrf);
handles.scrf.CurrentAxes = handles.standardax;
imshow(standardslice);
handles.scrf.CurrentAxes = handles.scfax;
imshow(refineslice);

% calculate and display transform matrix:
if isfile(options.subj.brainshift.transform.instore)
    mat = ea_getscrfmat(options);
    handles.affmatrix.String = sprintf('% 0.2f  % 0.2f  % 0.2f  % 0.2f  \n% 0.2f  % 0.2f  % 0.2f  % 0.2f  \n% 0.2f  % 0.2f  % 0.2f  % 0.2f  \n% 0.2f  % 0.2f  % 0.2f  % 0.2f  ', mat');
    save(options.subj.brainshift.transform.converted, 'mat');
end
