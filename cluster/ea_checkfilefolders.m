function ficell=ea_checknofilefolders(ficell,fname)

todel=[];
for fi=1:length(ficell)
    if ~exist([ficell{fi},filesep,fname],'file')
        todel=[todel,fi];
    end
end
ficell(todel)=[];