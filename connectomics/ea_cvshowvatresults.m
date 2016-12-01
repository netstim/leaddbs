function ea_cvshowvatresults(resultfig,pX,directory,filesare,handles,pV,selectedparc,options)

% determine if fMRI or dMRI
mods=get(handles.vatmodality,'String');
mod=mods{get(handles.vatmodality,'Value')};
if strfind(mod,'_tc')
        ea_cvshowvatfmri(resultfig,pX,directory,filesare,handles,pV,selectedparc,options);
else
        ea_cvshowvatdmri(resultfig,directory,handles,selectedparc,options);
end