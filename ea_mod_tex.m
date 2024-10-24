function ea_mod_tex(fname,options)

if exist([fname,'.aux'], 'file')
    delete([fname,'.aux']);
end

if exist([fname,'.tex'], 'file')
    delete([fname,'.tex']);
end

fID=fopen([fname,'.tex'],'w');

[pth,fn,ext]=fileparts(fname);
[~,patientname]=fileparts(pth);

% document header:

fprintf(fID,'%s\n','\documentclass{article}');
fprintf(fID,'%s\n','\usepackage[no-math]{fontspec}');
fprintf(fID,'%s\n','\usepackage{xunicode,xltxtra}');
fprintf(fID,'%s\n','\usepackage[dvipdfmx]{media9}');
fprintf(fID,'%s\n','\usepackage{graphicx}');
fprintf(fID,'%s\n','\usepackage{float}');
fprintf(fID,'%s\n\n','\usepackage{hyperref}');
fprintf(fID,'%s\n\n','\usepackage{bookmark}');

fprintf(fID,'%s\n','\hypersetup{');
fprintf(fID,'%s\n','    colorlinks,');
fprintf(fID,'%s\n','    citecolor=black,');
fprintf(fID,'%s\n','    filecolor=black,');
fprintf(fID,'%s\n','    linkcolor=black,');
fprintf(fID,'%s\n','    urlcolor=black');
fprintf(fID,'%s\n\n','}');

fprintf(fID,'%s\n',['\title{Lead-DBS Electrode localization}\author{',patientname,'}']);

fprintf(fID,'%s\n','\begin{document}');
fprintf(fID,'%s\n\n','	\pagenumbering{gobble}');
fprintf(fID,'%s\n','	\maketitle');
fprintf(fID,'%s\n\n','	\newpage');
fprintf(fID,'%s\n','	\tableofcontents');
fprintf(fID,'%s\n\n','	\newpage');
fprintf(fID,'%s\n','	\pagenumbering{arabic}');
fprintf(fID,'%s\n','');
fprintf(fID,'%s\n','	\section{3D view}');
fprintf(fID,'%s\n','	    \input{Lead-DBS_Electrode_Localization_small}');
fprintf(fID,'%s\n','	    ');
fprintf(fID,'%s\n','	\section{2D slices}');
fprintf(fID,'%s\n','	');

for side=1:2
    % begin right hemisphere

    switch side
        case 1
            cntnms=options.elspec.contactnames(1:options.elspec.numContacts);
            pf='Right';
        case 2
            cntnms=options.elspec.contactnames(1+options.elspec.numContacts:2*options.elspec.numContacts);
            pf='Left';
    end

    fprintf(fID,'%s\n',['	\subsection{',pf,' hemisphere}']);
    for cnt=1:options.elspec.numContacts
        fprintf(fID,'%s\n',['	\subsubsection{',cntnms{cnt},'}']);
        for view=1:3
            switch view
                case 1
                    views='Axial';
                case 2
                    views='Coronal';
                case 3
                    views='Sagittal';
            end
            fprintf(fID,'%s\n','	\begin{figure}[H]');
            fprintf(fID,'%s\n','		\centering');
            fprintf(fID,'%s\n',['		\includegraphics[keepaspectratio,width=\textwidth,height=0.75\textheight]{',cntnms{cnt},'_',lower(views),'}']);
            fprintf(fID,'%s\n',['		\caption{',views,' view of ',cntnms{cnt},'}']);
            fprintf(fID,'%s\n','	\end{figure}');
        end
    end
end

fprintf(fID,'%s\n','\end{document}');
fclose(fID);
