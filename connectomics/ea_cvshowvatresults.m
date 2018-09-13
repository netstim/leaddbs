function ea_cvshowvatresults(resultfig,pX,directory,filesare,handles,pV,selectedparc,options)

% determine if fMRI or dMRI
mods=get(handles.vatmodality,'String');
mod=mods{get(handles.vatmodality,'Value')};
if regexp(mod, '^Patient''s fMRI - ', 'once')
	ea_cvshowvatfmri(resultfig,pX,directory,filesare,handles,pV,selectedparc,mod,options);
else
	ea_cvshowvatdmri(resultfig,directory,handles,selectedparc,options);
end