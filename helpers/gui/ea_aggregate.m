function ea_aggregate(~,~,handles,exportwhat)


uipatdir=getappdata(handles.leadfigure,'uipatdir');
fname=uigetdir('','Select where to save files...');
fname=[fname,filesep];
ea_dispercent(0,'Exporting files');
for pt=1:length(uipatdir)
    [pth,ptname]=fileparts(uipatdir{pt});
    switch exportwhat
        case 'allcheckreg'
            infile=dir([uipatdir{pt},filesep,'checkreg',filesep,'*.png']);
            subpath=['checkreg',filesep];
        case 'normcheckreg'
            whichnormmethod=ea_whichnormmethod(uipatdir{pt});
            infile=dir([uipatdir{pt},filesep,'checkreg',filesep,'*',whichnormmethod,'.png']);
            subpath=['checkreg',filesep];
    end
    for inf=1:length(infile)
        
        copyfile([uipatdir{pt},filesep,subpath,infile(inf).name],[fname,ptname,'_',infile(inf).name]);
    end
    ea_dispercent(pt/length(uipatdir));
end
ea_dispercent(1,'end');

msgbox('Files successfully exported.');