function ea_cleanlegacy(hobj,evt,handles)

% Add unnecessary or legacy files in the following list that the user may then 
% delete using the tools menu entry in selected patient folders:

delfis={'lanat.nii'
    'lpostop_tra.nii'
    'lpostop_cor.nii'
    'lpostop_sag.nii'
    'lpostop_ct.nii'
    'tmp'};


dirs=getappdata(handles.leadfigure,'uipatdir');
prefs=ea_prefs;
for pt=1:length(dirs)
    for f=1:length(delfis)
        ea_delete([dirs{pt},filesep,delfis{f}]);
    end
    
    whichnormmethod=ea_whichnormmethod(dirs{pt});
    
    switch whichnormmethod
        
        case ea_getantsnormfuns
            
            ea_cleanspmwarps(dirs{pt},prefs)
            ea_cleanfslwarps(dirs{pt},prefs)
            
        case ea_getfslnormfuns
            
            ea_cleanspmwarps(dirs{pt},prefs)
            ea_cleanantswarps(dirs{pt},prefs)
            
        otherwise % all SPM functions
            ea_cleanantswarps(dirs{pt},prefs)
            ea_cleanfslwarps(dirs{pt},prefs)
    end
    
end
disp('Deleted all unnecessary files.');

function ea_cleanantswarps(dir,prefs)

ea_delete([dir,filesep,'glanatComposite.h5']);
ea_delete([dir,filesep,'glanatInverseComposite.h5']);
ea_delete([dir,filesep,'glanatComposite.nii.gz']);
ea_delete([dir,filesep,'glanatInverseComposite.nii.gz']);
ea_delete([dir,filesep,'lanatComposite.h5']);
ea_delete([dir,filesep,'lanatInverseComposite.h5']);
ea_delete([dir,filesep,'glanat1Warp.nii.gz']);
ea_delete([dir,filesep,'glanat1InverseWarp.nii.gz']);
ea_delete([dir,filesep,'glanat0GenericAffine.mat']);

function ea_cleanfslwarps(dir,prefs)

ea_delete([dir,filesep,'glanatInverseWarpField.nii']);
ea_delete([dir,filesep,'glanatWarpField.nii']);
ea_delete([dir,filesep,'glanatWarpCoef.nii']);


function ea_cleanspmwarps(dir,prefs)
ea_delete([dir,filesep,'y_*.nii']);
ea_delete([dir,filesep,'iy_*.nii']);
ea_delete([dir,filesep,'u_*.nii']);
ea_delete([dir,filesep,'c*anat*.nii']);
