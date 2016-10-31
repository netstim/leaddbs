% testDesignMatrix.m
%
%        $Id:$
%      usage: retval = testDesignMatrix(d,params)
%         by: julien besle
%       date: 11/04/2013
%    purpose: test for any zero-column in designMatrix (taken out of makeDesignMatrix)
%

function nullComponents = testDesignMatrix(designMatrix,nHdr,nHrfComponents,EVnames)

if any(all(~designMatrix,1))
  nullComponents = find((all(~designMatrix,1)));
  disp('(testDesignMatrix) The following EV components cannot be estimated because the corresponding column in the design is null');
  for iEV = 1:nHdr;
    whichComponents = find(ismember((iEV-1)*nHrfComponents+1:iEV*nHrfComponents,nullComponents));
    if ~isempty(whichComponents)
      fprintf('EV ''%s'': Components %s\n',EVnames{iEV},mat2str(whichComponents));
    end
  end
else
  nullComponents = [];
end
