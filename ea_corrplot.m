function fig=ea_corrplot(X,description,labels)
% this simple function is a small wrapper for a corrplot figure.

dim=size(X,2);

disp(description);
try
[R,P]=corrplot(X,'testR','on','varNames',labels)

set(gcf,'name',description);
set(gcf,'NumberTitle','off');
end
% fig=figure;
% 
% for xx=1:dim
%     for yy=1:dim
%         
%     if xx<yy % below diagonal - show correlations
% subplot(dim,dim,sub2ind([dim,dim],xx,yy));
% scatter(X(:,xx),X(:,yy),'bo');
%         lsline();
%     elseif xx>yy % higher than diagonal - show significancies
%         
%     elseif xx==yy % diagonal - show histograms.
%         
%         
%     end
%     end
% end
% 
% 
% 
