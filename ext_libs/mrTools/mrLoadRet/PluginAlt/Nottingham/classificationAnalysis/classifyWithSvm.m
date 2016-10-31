function [outclass weights]=classifyWithSvm(test_patt,train_patt,train_lab)


% mm = mean([train_patt,test_patt],2);
% % %is this necessary?
% test_patt=100*(test_patt-repmat(mm,1,size(test_patt,2)))./repmat(mm,1,size(test_patt,2));
% train_patt=100*(train_patt-repmat(mm,1,size(train_patt,2)))./repmat(mm,1,size(train_patt,2));

warning('off','optim:quadprog:SwitchToMedScale');
[gindex,groups,glevels] = grp2idx(train_lab);
ngroups = length(groups);
ncat=ngroups;
svm = cell(ncat-1,ncat);
weights=nan(ncat,ncat,size(train_patt,1));
for j = 1:(ncat-1)
    for k = (j+1):ncat
      % get all instances
      x1 = train_patt(:,gindex==j);
      x2 = train_patt(:,gindex==k);
      % remove the one instance
%       trainset1 = setdiff(1:size(x1,1),test_idx);
%       trainset2 = setdiff(1:size(x2,1),test_idx);
      % create the classifer
%       svm{i_sphere}{test_idx,j,k} = getsvm_nottingham_copy(x1',x2');%,kernelfun,kernelargs,C);
      svm{j,k} = getsvm_nottingham_copy(x1',x2');%,kernelfun,kernelargs,C);
      % store the weights
      weights(j,k,:) = svm{j,k}.w;
      weights(k,j,:) = -svm{j,k}.w;
    end
end
weights(:,end+1,:)=nansum(weights(:,1:end,:),2);

classifierout = nan(ncat,ncat,size(test_patt,2));
% classifierout=cell(1,size(test_patt,2));
% % now test each one of the left out response
% for testStim = 1:size(test_patt,2)
    % get each one of the test stimuli in turn
    testvec =  test_patt;
    for j = 1:ncat
      for k = j:ncat
        if (j == k) 
          % classification against oneself we define as 0
          classifierout(j,k,:) = zeros(1,size(test_patt,2));
        else
          % use the svm, calculated above w/out the
          % test stimulus to classify the stimulus against
          % each one of the other stimuli categories
          classifierout(j,k,:) = getsvm_nottingham_copy(testvec',svm{j,k});
          % in the matrix of outputs, the classification of the
          % other stimulus categories vs. this one, is symmetric
          % i.e. if we know the classification for stimulus 1 vs
          % stimulus 2, then the classification for stimulus 2 vs
          % stimulus 1 is just the negative of that.
          classifierout(k,j,:) = -classifierout(j,k,:);
        end
      end
      % now calculate the sum of the classifier against all other
      % stimuli types
%         classifierout{testStim}(jc,ncat+1) = sum(classifierout{testStim}(j,1:ncat));
    end
    % make a tuning curve--this gives thet summed classifier outputs
    % for each possible stimulus types.

%     temp_tuningcurve = squeeze(sum(classifierout));%{testStim}(:,ncat+1));
        temp_tuningcurve = squeeze(sum(classifierout,2));%{testStim}(:,ncat+1));

%     tuningcurve = temp_tuningcurve;
    % the maximum point on this tuning curve is the best guess at
    % how to classify this instance
    [~, classification] = max(temp_tuningcurve);

%     svm_corr{i_sphere}(test_idx,teststim) = classification{i_sphere}(test_idx,teststim) == teststim;
% end
outclass = glevels(classification);
warning('on','optim:quadprog:SwitchToMedScale')