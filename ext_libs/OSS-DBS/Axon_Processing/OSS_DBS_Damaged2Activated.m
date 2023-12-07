function fibers = OSS_DBS_Damaged2Activated(settings,fibers,idx,sideidx)

%% determine mean active contact
amplitudes = settings.Phi_vector(sideidx,:);
amplitudes = (abs(amplitudes)./sum(abs(amplitudes),"all",'omitnan'))';
contacts = settings.contactLocation{sideidx};
meanactive = sum(contacts.*amplitudes,1,'omitnan');

%% calculate minimal distance of each fiber to the mean active contact
fiberstmp = fibers;
fiberstmp(:,1:3)=fiberstmp(:,1:3)-meanactive;
distances = vecnorm(fiberstmp(:,1:3),2,2);
distances = mat2cell(distances, idx);
mindistances = cellfun(@min,distances);
%% calculate mean distance of activated fibers
fibtable=table;
[fibtable.idx,tmpidx]= unique(fibers(:,4));
fibtable.state = fibers(tmpidx,5);
fibtable.distance = mindistances;
fibtable.statenew = fibtable.state;

meanactivateddistance = mean(fibtable.distance(fibtable.state==1));

fibtable.statenew(fibtable.state == -1 & fibtable.distance < meanactivateddistance) = 1;
fibtable.statenew(fibtable.state == -1 & fibtable.distance >= meanactivateddistance) = 0;
%fibtable.statenew(fibtable.state == -2 & fibtable.distance < meanactivateddistance) = 1;
%fibtable.statenew(fibtable.state == -2 & fibtable.distance >= meanactivateddistance) = 0;
fibtable.statenew(fibtable.state == -3 & fibtable.distance < meanactivateddistance) = 1;
fibtable.statenew(fibtable.state == -3 & fibtable.distance >= meanactivateddistance) = 0;

for i = 1:height(fibtable)
    fibers(fibers(:,4)==fibtable.idx(i),5)=fibtable.statenew(i);
end

end
