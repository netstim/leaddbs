function ea_aggregateVTA(~,~,handles,exportwhat)


suffx='';
if ismember(exportwhat,{'vat_seed_compound_dMRI','vat_seed_compound_fMRI'})
    prefs=ea_prefs;
    switch prefs.lcm.vatseed
        case 'binary'
            suffx='';
        case 'efield'
            suffx='_efield';
        case 'efield_gauss'
            suffx='_efield_gauss';
    end
end

if ~strcmp(handles.seeddefpopup.String{handles.seeddefpopup.Value}(1:9),'Use VATs:')
    ea_error('Please select a stimulation to export first.');
else
    stimname=handles.seeddefpopup.String{handles.seeddefpopup.Value}(11:end);
end
uipatdir=getappdata(handles.leadfigure,'uipatdir');
fname=uigetdir('','Select where to save files...');
fname=[fname,filesep];

tf=fopen([fname,'export.txt'],'w');
options.native=0;
for pt=1:length(uipatdir)
    [pth,ptname]=fileparts(uipatdir{pt});
    copyfile([uipatdir{pt},filesep,'stimulation',filesep,ea_nt(options),stimname,filesep,exportwhat,suffx,'.nii'],[fname,ptname,'.nii']);
    fprintf(tf,'%s\n',[fname,ptname,'.nii']);
end
fclose(tf);


msgbox('Files successfully exported.');
