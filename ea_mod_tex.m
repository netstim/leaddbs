function ea_mod_tex(fname)

delete([fname,'.aux']);
delete([fname,'.tex']);
fID=fopen([fname,'.tex'],'w');

[pth,fn,ext]=fileparts(fname);
[~,patientname]=fileparts(pth);


% document header:
fprintf(fID,'%s\n','\documentclass[11pt, titlepage, twoside]{article}');
fprintf(fID,'%s\n','\usepackage[top=20mm, bottom=20mm, left=20mm, right=20mm]{geometry}');
fprintf(fID,'%s\n','\geometry{a4paper}');
fprintf(fID,'%s\n','\geometry{portrait}');
fprintf(fID,'%s\n','\usepackage[utf8]{inputenc}	% Comment out for xelatex; must use with pdflatex');
fprintf(fID,'%s\n','\usepackage[T1]{fontenc}');
fprintf(fID,'%s\n','\usepackage[parfill]{parskip} % Activate to begin paragraphs with an empty line rather than an indent');
fprintf(fID,'%s\n','\usepackage{graphicx}');
fprintf(fID,'%s\n','\usepackage{amsmath}');
fprintf(fID,'%s\n','\usepackage{amssymb}');
fprintf(fID,'%s\n','%\usepackage{fontspec}');
fprintf(fID,'%s\n','\usepackage{epstopdf}');
fprintf(fID,'%s\n','\usepackage{float}');
fprintf(fID,'%s\n','\usepackage[hang]{footmisc}');
fprintf(fID,'%s\n','\usepackage{endnotes}');
fprintf(fID,'%s\n','\usepackage{listingsutf8}');
fprintf(fID,'%s\n','%\restylefloat{table}');
fprintf(fID,'%s\n','\usepackage{xcolor}');
fprintf(fID,'%s\n','\usepackage{colortbl}');
fprintf(fID,'%s\n','%\usepackage{booktabs}');
fprintf(fID,'%s\n','\usepackage{tabu}');
fprintf(fID,'%s\n','\usepackage{caption}');
fprintf(fID,'%s\n','\usepackage{textcomp}');
fprintf(fID,'%s\n','\usepackage{subcaption}\usepackage[style=numeric,backend=bibtex8]{biblatex}');
fprintf(fID,'%s\n','\addbibresource{references.bib}');
fprintf(fID,'%s\n','\captionsetup[figure]{labelfont=sc}');
fprintf(fID,'%s\n','\captionsetup[table]{labelfont=sc}');
fprintf(fID,'%s\n','\captionsetup[tabu]{labelfont=sc}');

fprintf(fID,'%s\n','\DeclareCaptionType{MPEquation}[][List of equations]');
fprintf(fID,'%s\n','\captionsetup[MPEquation]{labelformat=empty}');

fprintf(fID,'%s\n','\DeclareCaptionType{MPListing}[][List of listings]');
fprintf(fID,'%s\n','\captionsetup[MPListing]{labelformat=empty}');
fprintf(fID,'%s\n','\usepackage{fancyhdr}');
fprintf(fID,'%s\n','\setlength{\headheight}{16pt}');
fprintf(fID,'%s\n','\pagestyle{fancyplain}');
fprintf(fID,'%s\n','\fancyhf{}');
fprintf(fID,'%s\n','%\fancypagestyle{plain}{}');
fprintf(fID,'%s\n','\cfoot[\thepage]{\thepage}');

fprintf(fID,'%s\n','\usepackage{stringenc}');
fprintf(fID,'%s\n','\usepackage{pdfescape}');


fprintf(fID,'%s\n','% Helper for debugging pdflatex errors induced by invalid UTF8 characters:');
fprintf(fID,'%s\n','% output the actual character code causing the issue.');

fprintf(fID,'%s\n','\makeatletter');
fprintf(fID,'%s\n','\renewcommand*{\UTFviii@defined}[1]{%');
fprintf(fID,'%s\n','\ifx#1\relax');
fprintf(fID,'%s\n','\begingroup');
fprintf(fID,'%s\n','% Remove prefix "\u8:"');
fprintf(fID,'%s\n','\def\x##1:{}%');
fprintf(fID,'%s\n','% Extract Unicode char from command name');
fprintf(fID,'%s\n','% (utf8.def does not support surrogates)');
fprintf(fID,'%s\n','\edef\x{\expandafter\x\string#1}%');
fprintf(fID,'%s\n','\StringEncodingConvert\x\x{utf8}{utf16be}% convert to UTF-16BE');
fprintf(fID,'%s\n','% Hexadecimal representation');
fprintf(fID,'%s\n','\EdefEscapeHex\x\x');
fprintf(fID,'%s\n','% Enhanced error message');
fprintf(fID,'%s\n','\PackageError{inputenc}{Unicode\space char\space \string#1\space');
fprintf(fID,'%s\n','(U+\x)\MessageBreak');
fprintf(fID,'%s\n','not\space set\space up\space');
fprintf(fID,'%s\n','for\space use\space with\space LaTeX}\@eha');
fprintf(fID,'%s\n','\endgroup');
fprintf(fID,'%s\n','\else\expandafter');
fprintf(fID,'%s\n','#1%');
fprintf(fID,'%s\n','\fi');
fprintf(fID,'%s\n','}');
fprintf(fID,'%s\n','\makeatother');
fprintf(fID,'%s\n','');
fprintf(fID,'%s\n','% Make multi-paragraph footnotes align nicely on the left');
fprintf(fID,'%s\n','\setlength\footnotemargin{10pt}');
fprintf(fID,'%s\n','');
fprintf(fID,'%s\n','% Do not display an automatic Notes heading for endnotes, as they are contained');
fprintf(fID,'%s\n','% by a Manuscripts section with a heading that has possibly been edited by the user');
fprintf(fID,'%s\n','\def\enoteheading{}');
fprintf(fID,'%s\n','');
fprintf(fID,'%s\n','% Apply footnote and endnote numbering schemes (decimal, lowercase roman etc.)');
fprintf(fID,'%s\n','\renewcommand{\thefootnote}{\arabic{footnote}}');
fprintf(fID,'%s\n','\renewcommand{\theendnote}{\arabic{endnote}}');
fprintf(fID,'%s\n','% Adjust endnotes to use normal-sized numbering and text');
fprintf(fID,'%s\n','\renewcommand{\enotesize}{\normalsize}');
fprintf(fID,'%s\n','\renewcommand\enoteformat{%');
fprintf(fID,'%s\n','  \raggedright');
fprintf(fID,'%s\n','  \leftskip=1.8em');
fprintf(fID,'%s\n','  \makebox[0pt][r]{\theenmark. \rule{0pt}{\dimexpr\ht\strutbox+\baselineskip}}%');
fprintf(fID,'%s\n','}');
fprintf(fID,'%s\n','');


fprintf(fID,'%s\n','\begin{document}');

% title of document:

fprintf(fID,'%s\n',['\title{Lead-DBS Electrode localization}\author{',patientname,'}']);

fprintf(fID,'%s\n','\maketitle');
guid=ea_generate_guid;
fprintf(fID,'%s\n',['\newtabulinestyle {MPFigureStyle_',guid,'_innerBorder=0pt}']);
fprintf(fID,'%s\n',['\newtabulinestyle {MPFigureStyle_',guid,'_outerBorder=0pt}']);

fprintf(fID,'%s\n',['\section{2D slices}\label{MPSection:',ea_generate_guid,'}']);

for side=1:2
    % begin right hemisphere
    
    switch side
        case 1
            cnts=0:3;
            pf='Right';
        case 2
            cnts=8:11;
            pf='Left';
    end
    
    fprintf(fID,'%s\n',['\subsection{',pf,' hemisphere}\label{MPSection:',ea_generate_guid,'}']);
    
    for cnt=0:3
        fprintf(fID,'%s\n',['\subsubsection{K',num2str(cnt),'}\label{MPSection:',ea_generate_guid,'} ']);
        for view=1:3
            switch view
                case 1
                    views='Axial';
                case 2
                    views='Coronar';
                case 3
                    views='Sagittal';
            end
            fprintf(fID,'%s\n','\begin{figure}[htbp]');
            fprintf(fID,'%s\n','\centering');
            fprintf(fID,'%s\n',['\includegraphics[keepaspectratio,width=\textwidth,height=0.75\textheight]{K',num2str(cnt),'_',lower(views),'.png}']);
            fprintf(fID,'%s\n',['\caption{',views,' view of K',num2str(cnt),'}\label{MPFigure:',ea_generate_guid,'}']);
            fprintf(fID,'%s\n','\end{figure}');
        end
    end
    
end

%enter 3D view

fprintf(fID,'%s\n',['\section{3D view}\label{MPSection:',ea_generate_guid,'}']);
fprintf(fID,'%s\n','    \input{3D_scene_small.tex}%');


fprintf(fID,'%s\n','\end{document}');


fclose(fID);