function AUC = ea_logit_regression(Ihat_train, Ihat, Improvement, training, test)

% Fit logit model, compute ROC and find the optimal threshold.
% Compute confustion matrix for the test set (can be the same as training)

% if Ihat_train was not provided, then we have in-sample analysis
if Ihat_train == 0
    Ihat_train = Ihat(training);
end

% first, we fit a logit function for our binary prediction
mdl = fitglm(Ihat_train,Improvement(training),'Distribution','binomial','Link','logit');



% second, we run ROC curve analysis
scores = mdl.Fitted.Probability;
[X,Y,T,AUC,OPTROCPT] = perfcurve(Improvement(training),scores,1);
figure, plot(X,Y, 'k', 'linew', 1.5)
set(gcf,'color','w');
hold on
plot(OPTROCPT(1),OPTROCPT(2),'ro', 'MarkerSize',10)
xlabel('False positive rate') 
ylabel('True positive rate')
txt = ['AUC: ' num2str(AUC)];
text(0.7,0.1,txt)
title('ROC for Classification by Logistic Regression')

% optimal threshold on the classifier
scores_thresh = T((X==OPTROCPT(1))&(Y==OPTROCPT(2)));

% plot logit fit for training
vec_val = min(Ihat_train):1:max(Ihat_train);
figure
set(gcf,'color','w');
lims=[min(Ihat_train)-0.1*(max(Ihat_train)-min(Ihat_train)), max(Ihat_train)+0.1*(max(Ihat_train)-min(Ihat_train))];
subplot(4,1,1)
subtitle('Response');
col=ea_color_wes('lifeaquatic');
g=ea_raincloud_plot(Ihat_train(Improvement(training)==1)','color',col(3,:),'box_on',1);
a1=gca;
set(a1,'ytick',[])
set(gca, 'xlim', lims)
a1.YLabel.String='Response';
a1.XLabel.String='Fiberscore';
a1.Box='off';
title('Logistic Regression for Training Cohort')

subplot(4,1,[2 3])
plot(vec_val', predict(mdl,vec_val'),'k', 'linew', 1.5)
xlabel('Fiberscore'), ylabel('Response')
set(gca, 'xlim', lims); box off
%plot(vec_val', predict(mdl,vec_val'),Ihat_train,Improvement(training),'s')
%plot(Ihat_av,Improvement,'s')

subplot(4,1,4)
subtitle('Control');
col=ea_color_wes('lifeaquatic');
g=ea_raincloud_plot(Ihat_train(Improvement(training)==0)','color',col(1,:),'box_on',1);
a1=gca;
set(a1,'ytick',[])
set(gca, 'xlim', lims)
a1.YLabel.String='Control';
a1.XLabel.String='Fiberscore';
a1.Box='off';


% prediction for test based on the logit model
scores_test = predict(mdl,Ihat(test));
Ihat_prediction = scores_test > scores_thresh;

% get the confussion matrix (this can be done on the test set now)
figure
cm = confusionchart(logical(Improvement(test)), Ihat_prediction);
set(gcf,'color','w');

tp = sum((Ihat_prediction == 1) & (Improvement(test) == 1));
fp = sum((Ihat_prediction == 1) & (Improvement(test) == 0));
tn = sum((Ihat_prediction == 0) & (Improvement(test) == 0));
fn = sum((Ihat_prediction == 0) & (Improvement(test) == 1));

sensitivity = tp/(tp + fn);  % TPR
specificity = tn/(tn + fp);  % TNR

cm.Title = ['Sensitivity: ', sprintf('%.2f',sensitivity), '; ', 'Specificity: ', sprintf('%.2f',specificity)];

end
