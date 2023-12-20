function Ihat=ea_ihat_databaselookup_netmap(obj,vals,connval,patientsel,test,training,usemask)


% set up the fingerprint to predict and the variable of interest 
fp_topredict=connval(patientsel(test),usemask)';
I_test=obj.responsevar(patientsel(test));

vizz=0; % visualisation of sharpening. 0 = turned off, 1 = turned on


sims=ea_corr(fp_topredict,vals{1}(:,usemask)',resolvecorrtype(obj.basepredictionon));
disp(['Average of similarities: ',num2str(mean(sims)),'.']);
sh_exp_sim=1; % exponent by which to sharpen similarities weight; 1 = no sharpening, can be set to larger number
sh_exp_cert=1; % exponent by which to sharpen certainties weight (if available); 1 = no sharpening, can be set to larger number
adaptivemode=1; % if adaptivemode is switched on, the standard sharpening exponents will be ignored
percenttarget_w1 = 0.1; % for adaptive mode, 10% of peaks are above mean, can be adjusted
percenttarget_w2 = 0.1; % for w2 adaptive mode, 10% of peaks are above mean, can be adjusted
% if we don't sharpen (e.g. sh_exp=1) then most patients will be taken into
% account similarly


% adaptive 10%
w1=sims-ea_nanmin(sims(:)); % make all positive
w1=ea_minmax(w1); % normalize

if vizz
    h=figure;
end


if adaptivemode==1
    for leftoutpt=1:size(w1,1)
        thisw1=w1(leftoutpt,:);

        for exponent=1:100 %testing exponents to get results, if results not found, largest exponent is used
            currentstate=thisw1.^exponent;
            if vizz
                set(0,"CurrentFigure",h);
                plot(currentstate);
                hold on
                plot(repmat(nanmean(currentstate),1,size(currentstate,2)));                
                drawnow
                hold off
                pause(0.1);
            end
            bincurrentstate=currentstate>ea_nanmean(currentstate,2);
            percentpeaks=ea_nansum(bincurrentstate)/length(bincurrentstate);
%             disp(percentpeaks);
            % thisstd=std(ea_minmax(thisw1.^exponent));
            if percentpeaks<=percenttarget_w1
                sh_exp_sim=exponent;
                break
            end
        end
        %% similarities
%         disp(percentpeaks);
        thisw1=thisw1.^sh_exp_sim; % sharpen
        thisw1=thisw1./repmat(sum(thisw1,2),1,sum(training)); % make sure they sum up to 1
        w1(leftoutpt,:)=thisw1;
    end
else
    w1=w1.^sh_exp_sim; % sharpen
    w1=w1./repmat(sum(w1,2),1,sum(training)); % make sure they sum up to 1
end

% figure, plot(w1')
% title('Sharpened similarities')

if ~strcmp(obj.certainvarlabel,'None') % certainty sidecar
    %% certainties
    w2=obj.certainvar(patientsel(training));
    w2=repmat(w2,1,ea_nansum(test))';
    w2=w2-ea_nanmin(w2(:)); % make all positive
    w2=ea_minmax(w2);
%     figure, plot(w2')
%     title('Raw absolute improvement')
    if vizz
        h1=figure;
    end

    if adaptivemode==1
        for leftoutpt=1:size(w2,1)
            thisw2=w2(leftoutpt,:);

            for exponent=1:100
                currentstate=thisw2.^exponent;
                if vizz
                    set(0,"CurrentFigure",h1);
                    plot(currentstate);
                    hold on
                    plot(repmat(ea_nanmean(currentstate,2),1,size(currentstate,2)));                
                    drawnow
                    hold off
                    pause(0.1);
                end
                bincurrentstate=currentstate>ea_nanmean(currentstate,2);
                percentpeaks=sum(bincurrentstate)/length(bincurrentstate);
%                 disp(percentpeaks);
                if percentpeaks<=percenttarget_w2
                    sh_exp_cert=exponent;
                    break
                end
             end
        %% certainty
%         disp(percentpeaks)
        thisw2=thisw2.^sh_exp_cert; % sharpen
        thisw2=thisw2./repmat(ea_nansum(thisw2,2),1,ea_nansum(training)); % make sure they sum up to 1
        w2(leftoutpt,:)=thisw2;

        end
    
    else
    w2=w2.^sh_exp_cert; % sharpen
    w2=w2./repmat(ea_nansum(w2,2),1,ea_nansum(training)); % make sure they sum up to 1
    end

    % multiplying the weights, the mean could also be used
    w=w1.*w2;
    w=w./repmat(sum(w,2),1,sum(training));
    %w = ea_robustmean([w1;w2],1);
    %w = mean([w1;w2],1);
    w=w./repmat(ea_nansum(w,2),1,ea_nansum(training));
   
else

    w=w1;


end

% remove nan values from the variable of interest
nnan_imp=obj.responsevar(patientsel(training));
ix=~isnan(nnan_imp);
nnan_imp=nnan_imp(ix);

% prediction
Ihat=(w(:,ix)*nnan_imp);

%print the results
disp(['Estimate of improvement: ',num2str(Ihat'),'.']);
disp(['Actual improvement: ',num2str(I_test'),'.']);
disp(['The sharpening exponent for w1 is: ',num2str(sh_exp_sim'),'.']);
disp(['The sharpening exponent for w2 is: ',num2str(sh_exp_cert'),'.']);


function c=resolvecorrtype(base)


switch lower(base)
    case 'spatial correlations (spearman)'
        c='Spearman';
    case 'spatial correlations (pearson)'
        c='Pearson';
    case 'spatial correlations (bend)'
        c='Bend';
end
