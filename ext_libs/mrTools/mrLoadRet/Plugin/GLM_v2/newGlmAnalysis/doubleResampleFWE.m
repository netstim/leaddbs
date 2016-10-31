% doubleResampleFWE.m
%
%        $Id$
%      usage: count = doubleResampleFWE(statistic,d.indexSortActual,params,d)
%         by: julien besle, 
%       date: 12/01/2011
%    purpose: compute Single-step and Step-down minP adjusted p-values
%             using algorithm 4 of Ge et al. Test 2003 (steps 2-5)
%             the nVoxel*nTests*nResamples matrix of p-values must be provided
%             as well as the actual bootstrap P-values and their ordering
%             the procedure is applied test by test

function countMinP = doubleResampleFWE(bootstrapP,actualP,d,nResamples)

bootstrapP = permute(bootstrapP,[1 3 2])/nResamples;
countMinP = NaN(size(bootstrapP,1),size(bootstrapP,3));
        
for iTest = 1:size(bootstrapP,3)
  switch(d.method)
    case 'Single-step' %this is a single-step version where there is not need 
                       %to sort the bootstrap p according to the actual p-value
      bootstrapP(:,:,iTest) = repmat(min(bootstrapP(:,:,iTest),[],1),[size(bootstrapP,1) 1]);
      countMinP(:,iTest) = sum(bootstrapP(:,:,iTest)<=repmat(actualP(:,iTest),[1 nResamples]),2);
      
    case 'Adaptive Single-step' %this is a single-step version where there we need 
                       %to sort the bootstrap p according to the actual p-value
                       %in order to use the estimated number of true null H0
      sortedBootstrapP = bootstrapP(d.indexSortActual{iTest},:,iTest);
      numberH0 = size(sortedBootstrapP,1);
      sortedActualP = actualP(d.indexSortActual{iTest},iTest);                           % DEBUG
      thisCountMinP = NaN(size(sortedActualP));
      for j=numberH0:-1:numberH0-d.numberFalseH0(iTest)+1
        thisCountMinP(j) = sum(min(sortedBootstrapP(j-d.numberTrueH0(iTest)+1:j,:),[],1)<=sortedActualP(j));
      end
      sortedBootstrapP = repmat(min(sortedBootstrapP,[],1),[d.numberTrueH0(iTest) 1]);
      thisCountMinP(1:j-1) = sum(sortedBootstrapP<=repmat(sortedActualP(1:j-1),[1 nResamples]),2);
      countMinP(d.actualIsNotNaN(:,iTest),iTest) = thisCountMinP(d.indexReorderActual{iTest});

  
    case {'Step-down','Adaptive Step-down'}
      sortedBootstrapP = bootstrapP(d.indexSortActual{iTest},:,iTest);
      numberH0 = size(sortedBootstrapP,1);
      sortedActualP = actualP(d.indexSortActual{iTest},iTest);                           
      thisCountMinP = NaN(size(sortedActualP));
      thisQ = ones(1,nResamples);
      %compute consecutive minima such that for decreasing p values (increasing significance)
      %the min is taken over more and more voxels (with a maximum of d.numberTrueH0 if this has been estimated)
      for j = numberH0:-1:numberH0-d.numberTrueH0(iTest)+1
        thisQ = min(thisQ,sortedBootstrapP(j,:));
        thisCountMinP(j) = sum(thisQ<=sortedActualP(j));
      end
      %from this point, take the min only over the estimated number of true H0
      for j = numberH0-d.numberTrueH0(iTest):-1:1
        thisCountMinP(j) = sum(min(sortedBootstrapP(j:j+d.numberTrueH0(iTest)-1,:),[],1)<=sortedActualP(j));
      end
      countMinP(d.actualIsNotNaN(:,iTest),iTest) = thisCountMinP(d.indexReorderActual{iTest});
  end
end