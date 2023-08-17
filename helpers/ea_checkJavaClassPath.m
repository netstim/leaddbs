function status = ea_checkJavaClassPath
% Check if static Java class path has been set for LeadDBS

classpath = {fullfile(ea_getearoot, 'ext_libs', 'dragndrop')
    fullfile(ea_getearoot, 'helpers', 'gui', 'ListBoxRenderer')};

spath = javaclasspath('-static');

if all(ismember(classpath, spath))
    status = 1;
else
    classpathFile = fullfile(prefdir, 'javaclasspath.txt');
    f = fopen(classpathFile, 'a');
    fprintf(f, '%s\n', classpath{:});
    fclose(f);

    warning('off', 'backtrace');
    warndlg('Java class path updated, please restart MATLAB!');
    warning('on', 'backtrace');

    status = 0;
end
