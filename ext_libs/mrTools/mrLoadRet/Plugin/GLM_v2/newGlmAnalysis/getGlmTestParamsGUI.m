% getGlmTestParamsGUI.m
%
%        $Id$
%      usage: params = getGlmTestParamsGUI(params,useDefault)
%         by: julien besle, 
%       date: 09/11/2010
%    purpose: return statistical test parameters for GLM analysis
%

function params = getGlmTestParamsGUI(thisView,params,useDefault)

keepAsking = 1;

while keepAsking
  %Default params
  if isfield(params, 'contrasts') && strcmp(params.contrasts,'all')
    params.numberContrasts = params.numberEVs;
  end
  if fieldIsNotDefined(params, 'contrasts') || ~isnumeric(params.contrasts) || ...
    ~isequal(size(params.contrasts,2),params.numberEVs)
      params.contrasts=eye(params.numberContrasts,params.numberEVs);
  end
  %add or remove lines if number of contrasts does not match what's in the contrast matrix 
  if params.numberContrasts>size(params.contrasts,1)
    params.contrasts=[params.contrasts;zeros(params.numberContrasts-size(params.contrasts,1),size(params.contrasts,2))];
  elseif params.numberContrasts<size(params.contrasts,1)
    params.contrasts=params.contrasts(1:params.numberContrasts,:);
  end
  if fieldIsNotDefined(params,'tTestSide')
    params.tTestSide = 'Both';
  end

  if fieldIsNotDefined(params, 'fTestNames') || ~isequal(length(params.fTestNames),params.numberFtests) 
    if params.numberFtests
      params.fTestNames=cellstr(reshape(sprintf('fTest%2d',1:params.numberFtests),7,params.numberFtests)');
    else
      params.fTestNames={};
    end
  end
  if fieldIsNotDefined(params, 'restrictions') || ...
    ~isequal(length(params.restrictions),params.numberFtests) || ~isequal(size(params.restrictions{1},2),params.numberEVs)
    params.restrictions = {};
    for iFtest = 1:params.numberFtests
        params.restrictions{iFtest}=zeros(params.numberEVs-1,params.numberEVs);
    end
  end

  %create model HRF to get the number of components per EV
  %here we assume that all scans in this group have the same sampling parameters 
  %(which are needed to determine the number of components in the deconvolution case) 
  framePeriod = viewGet(thisView,'framePeriod',params.scanNum(1),viewGet(thisView,'groupNum',params.groupName));
  [hrfParams,hrf] = feval(params.hrfModel, params.hrfParams,framePeriod/params.scanParams{params.scanNum(1)}.estimationSupersampling,[],1);
  nComponents = size(hrf,2);
  if fieldIsNotDefined(params, 'componentsToTest') || ~isequal(nComponents,length(params.componentsToTest));
    params.componentsToTest = ones(1,nComponents);
  end
  if fieldIsNotDefined(params, 'componentsCombination')
    params.componentsCombination = 'Or';
  end
  if fieldIsNotDefined(params, 'fweAdjustment')
    params.fweAdjustment = 0;
  end
  if fieldIsNotDefined(params, 'fdrAdjustment')
    params.fdrAdjustment = 1;
  end
  if fieldIsNotDefined(params, 'alphaContrastOverlay')
      params.alphaContrastOverlay = 'FDR';
  end
  if fieldIsNotDefined(params, 'statisticalThreshold')
      params.statisticalThreshold = .05;
  end
  if fieldIsNotDefined(params, 'showAdvancedStatisticMenu') 
    params.showAdvancedStatisticMenu = 0;
  end
  
  tTestSideMenu = putOnTopOfList(params.tTestSide,{'Both','Right','Left'});
  componentsCombinationMenu = putOnTopOfList(params.componentsCombination,{'Add','Or'});
  alphaContrastOverlayMenu = putOnTopOfList(params.alphaContrastOverlay,{'FDR','FWE','Uncorrected','None'});
  
  contrastOptionsVisible = 'visible=0';
  if params.numberContrasts
    contrastOptionsVisible = 'visible=1';
  end
  tTestOptionsVisible = 'visible=0';
  if params.computeTtests
    tTestOptionsVisible = 'visible=1';
  end
  fTestOptionsVisible = 'visible=0';
  if params.numberFtests
    fTestOptionsVisible = 'visible=1';
  end
  testOptionsVisible = 'visible=0';
  if params.computeTtests || params.numberFtests
    testOptionsVisible = 'visible=1';
  end
  componentOptionsVisible = 'visible=0';
  if (params.computeTtests || params.numberFtests) && nComponents>1
      componentOptionsVisible = 'visible=1';
  end
  
  paramsInfo = {...
      {'EVnames', params.EVnames, 'type=stringarray','editable=0','Name of the EVs'},...
      {'contrasts', params.contrasts,contrastOptionsVisible,...
            'incdec=[-1 1]','incdecType=plusMinus',...
            'Each row is a linear combination of EVs that defines a contrast of interest. Each contrast will be output as an overlay.'},...
      {'tTestSide', tTestSideMenu,tTestOptionsVisible,'type=popupmenu', 'Sidedness of contrast T-tests (Both = two-sided, Right = one-sided positive, Left = one-sided negative)'},...
       };
  for iFtest = 1:params.numberFtests
    paramsInfo{end+1} = {fixBadChars(sprintf('fTest%2d',iFtest)), params.fTestNames{iFtest},fTestOptionsVisible ,['Name of F-test ' num2str(iFtest)]};
    paramsInfo{end+1} = {fixBadChars(sprintf('restriction%2d',iFtest)), params.restrictions{iFtest},fTestOptionsVisible,...
            'incdec=[-1 1]','incdecType=plusMinus',...
            ['Restriction matrix defining F-test ' num2str(iFtest) '.']};
  end
  paramsInfo = [paramsInfo {...
      {'componentsToTest', params.componentsToTest,componentOptionsVisible,...
            'Vector defining which HRF component of each EV are tested for significance. Put zeros to exclude components or a non-zero weight to include them. '},...
      {'componentsCombination', componentsCombinationMenu,componentOptionsVisible,'type=popupmenu', 'How to combine EV components. ''Add'' adds the weighted components into a single EV for contrasts/F-test. ''Or'' ignores the weights and tests contrasts/F-tests at any component that is not 0. Note that ''Or'' is not compatible with one-sided T-tests'},...
      {'fdrAdjustment', params.fdrAdjustment,testOptionsVisible,'type=checkbox', 'Adjusts test outcomes for false discovery rate control (see advanced menu to change the method; default is the adaptive step-up method)'},...
      {'fweAdjustment', params.fweAdjustment,testOptionsVisible,'type=checkbox', 'Adjusts test outcomes for familywise error rate control by Bonferroni-type methods (see advanced menu to change the method; default is the adaptive step-down method)'},...
      {'alphaContrastOverlay', alphaContrastOverlayMenu,testOptionsVisible,'type=popupmenu', 'Sets alpha overlay of contrast overlay(s) to one of the computed significance values.'},...
      {'statisticalThreshold', params.statisticalThreshold,testOptionsVisible, 'minmax=[0 1]', 'Threshold (expressed as a p-value) that will be applied to the statistical maps. The threshold is not applied to the data themselves but to the range property of the overlay and can therefore be modified after running the analysis. The threshold p-value is automatically converted to the corresponding Z or -log10(P) value.'},...
      {'showAdvancedStatisticMenu', params.showAdvancedStatisticMenu,testOptionsVisible,'type=checkbox', 'Advanced settings for correlated noise correction, test outputs, FDR/FWE corrections, bootstrap/permutation tests and confidence intervals.'},...
       }];

  if useDefault
    tempParams = mrParamsDefault(paramsInfo);
  else
    tempParams = mrParamsDialog(paramsInfo,'Statistics Menu (Contrasts, parametricT-tests and F-tests)');
  end

  % user hit cancel
  if isempty(tempParams)
    params = tempParams;
    return;
  end
  
  params = mrParamsCopyFields(tempParams,params);
  %check that contrasts are not empty
  actualNumberContrasts=0;
  if isempty(params.contrasts) %mrParamsDialog returns an empty string instead of an empty array
    params.contrasts = [];
  else
    for iContrast = 1:params.numberContrasts
      if ~any(params.contrasts(iContrast,:))
        mrWarnDlg('(getGlmTestParamsGUI) Discarding empty contrast');
      else
        actualNumberContrasts = actualNumberContrasts+1;
        params.contrasts(actualNumberContrasts,:) = params.contrasts(iContrast,:);
      end
    end
    params.contrasts = params.contrasts(1:actualNumberContrasts,:);
  end
  params.numberContrasts = actualNumberContrasts;
  if ~params.numberContrasts %don't compute T-tests if there are no contrasts
    params.computeTtests=0;
  end

  allContrasts = params.contrasts;
  %check if contrasts are pairwise orthogonal and issue a warning if not
  dotProducts = allContrasts*allContrasts';
  if any(dotProducts( ~(diag(ones(length(dotProducts),1))) )) %if any off-diagonal dot-product is non-zero
    mrWarnDlg('(getGlmTestParamsGUI) Contrasts are not pairwise orthogonal.');
  end
  
  %check that F-tests are not empty
  restrictions = {};
  fTestNames={};
  actualNumberFtests=0;
  for iFtest = 1:params.numberFtests
    thisRestriction=params.(fixBadChars(sprintf('restriction%2d',iFtest)));
    params = mrParamsRemoveField(params,fixBadChars(sprintf('restriction%2d',iFtest)));
    if ~any(any(thisRestriction))
      mrWarnDlg('(getGlmTestParamsGUI) Discarding F-test with empty restriction matrix.');
    else
      actualNumberFtests = actualNumberFtests+1;
      fTestNames{actualNumberFtests} = params.(fixBadChars(sprintf('fTest%2d',iFtest)));
      restrictions{actualNumberFtests,1} = thisRestriction;
      allContrasts = [allContrasts;thisRestriction];
    end
    params = mrParamsRemoveField(params,fixBadChars(sprintf('fTest%2d',iFtest)));
  end
  params = mrParamsCopyFields(mrParamsDefault({...
  {'fTestNames',fTestNames,'Self-explanatory'},...
  {'restrictions',restrictions,'Restrictions matrix defining F-tests. Each line of each matrix defines a contrast of EVs.'},...
  }),params);
  params.restrictions = restrictions; %copy this again, as mrParamsDefault changes the value of these two parameters
  params.fTestNames = fTestNames; %but we use it to create a paramInfo entry for them, which gives acces to the hlep info later on
  
  %check that contrasts in restriction matrices are orthogonal
  orthogonal=1;
  for iR = 1:length(restrictions)
    dotProducts = params.restrictions{iR}*params.restrictions{iR}';
    if any(dotProducts( ~(diag(ones(length(dotProducts),1))) )) %if any off-diagonal dot-product is non-zero
      mrWarnDlg(sprintf('(getGlmTestParamsGUI) Contrasts for restriction matrix %i are not pairwise orthogonal.',iR));
      orthogonal=0;
    end
  end
  
  %check that user has selected the correction (FDR or FWE) corresponding to their choice of alpha overlay
  if params.computeTtests && strcmp(params.alphaContrastOverlay,'FDR') && ~params.fdrAdjustment
    params.fdrAdjustment=1;
    mrWarnDlg('(getTestParamsGUI) setting params.fdrAdjustment to Yes because it is the selected contrast overlay alpha','Yes');
  end
  if params.computeTtests && strcmp(params.alphaContrastOverlay,'FWE') && ~params.fweAdjustment
    params.fweAdjustment=1;
    mrWarnDlg('(getTestParamsGUI) setting params.fweAdjustment to Yes because it is the selected contrast overlay alpha','Yes');
  end
  
  %check consistency of parameters
  if params.numberContrasts && params.computeTtests  && ~strcmp(params.tTestSide,'Both') && ...
      nnz(params.componentsToTest)>1 && strcmp(params.componentsCombination,'Or')
    mrWarnDlg('(getTestParamsGUI) One-sided T-tests on several EV components with ''Or'' combination are not implemented','Yes');
  elseif ~orthogonal %if there is at least one restriction matrix with non-orthogonal contrasts
  else
    params.numberFtests = actualNumberFtests;
    params.numberContrasts = actualNumberContrasts;
    tempParams = getGlmAdvancedStatisticParamsGUI(thisView,params,useDefault);
    % if empty, user hit cancel, go back
    if isempty(tempParams)
      if size(tempParams,2) %if the top close button has been pressed
        params=tempParams;
        return
      else
        keepAsking = 1;
      end
    else
      params = tempParams;
      keepAsking = 0;
    end
  end
  if keepAsking && useDefault
    params = [];
    return;
  end

end



