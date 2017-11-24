function ea_checkoutliers(~,~,handles)


answ=inputdlg('Detect electrodes how many standard-deviations distant from robust average coordinate?','Set Threshold',1,{'2.5'});
sfact=str2double(answ);
M=getappdata(handles.leadfigure,'M');
ptidx=get(handles.patientlist,'Value');

els=M.elstruct(ptidx);
for side=1:2
    heads{side}=nan(length(els),3); tails{side}=nan(length(els),3);
    for pt=1:length(els)
        heads{side}(pt,:)=els(pt).markers(side).head;
        tails{side}(pt,:)=els(pt).markers(side).tail;
    end
    heads{side}(any(isnan(heads{side}),2),:)=[];
    tails{side}(any(isnan(tails{side}),2),:)=[];
meanhead=ea_robustmean(heads{side});
meantail=ea_robustmean(tails{side});

% calc for head
hdists=squareform(pdist([meanhead;heads{side}]));
hdists=hdists(1,2:end);
% calc for tail
tdists=squareform(pdist([meantail;tails{side}]));
tdists=tdists(1,2:end);

cdists{side}=hdists+tdists;

cstd=std(cdists{side});

outlieridx{side}=cdists{side}>(cstd*sfact);
    
end

outliers=find(logical(outlieridx{1}+outlieridx{2}));
dists=[cdists{1}+cdists{2}];
[~,ix]=sort(dists(outliers),'descend');
outliers=outliers(ix);
for ol=outliers
    [~,olname]=fileparts(M.patient.list{ptidx(ol)});
   disp(['Outlier (>',num2str(sfact),'STD): ',olname,' (',num2str(dists(ol)),' mm).']);    
end


