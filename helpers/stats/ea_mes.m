function stats=ea_mes(X,Y,esm,varargin)
% ** function stats=mes(X,Y,esm,varargin)
% computes measures of effect size or related quantities between two
% samples (2-sample analyses) or one sample and a null value (1-sample
% analyses). The output consists of effect size measure(s), t statistics,
% and 95% confidence intervals where applicable. All input parameters
% except X, Y and esm are optional and must be specified as parameter/value
% pairs in any order, e.g. as in
%      mes(X,Y,'glassdelta','isDep',1,'nBoot',3000);
% 
%                           INPUT ARGUMENTS
%                           ---------------
% stats=mes(X,Y,esm)  computes the effect size measure(s) specified in 
%   input variable esm (see table at bottom) between samples X and Y or
%   between sample X and null value Y. X and Y may have multiple columns;
%   in this case the numbers of columns in X and Y must be equal as
%   comparisons will be made between matching columns. Furthermore, it is
%   assumed that each row corresponds to one subject/case (see treatment of
%   missing values below). No correction for multiple comparisons is
%   applied.
%   In 2-sample analyses it is assumed that the samples are independent
%   (not paired); if they are dependent (paired) use optional input
%   argument isDep (see below). Both the effect size measures available
%   and their computation are contingent on whether the samples are
%   independent or dependent. 
%   Input variable esm must be a character array (e.g. 'hedgesg') or a cell
%   array (e.g. {'hedgesg','glassdelta'}). 
% stats=mes(...,'isDep',1)  assumes that the samples in X and Y are
%   dependent (paired); accordingly, X and Y must have an equal number of 
%   rows. If this parameter is omitted the assumption is one of independent 
%   samples. 
% stats=mes(...,'missVal','listwise')  omits NaNs and Infs in X and Y, 
%   summarily termed here as missing values, in a within-variable listwise
%   fashion (which is the default): if there is a missing value anywhere in
%   X, the entire row will be omitted from analysis. If X has only one
%   column, only the single data point will be omitted. If the data are
%   unpaired the same applies to Y, independent of foul values in X. In
%   case of paired data deletion of rows is done in truly listwise fashion:
%   matching rows in X and Y will be omitted. If 'missVal' is set to
%   'pairwise' only the individual missing data points will be omitted.
%   Note that in case of paired data matching data points in X or Y will be
%   omitted.
% stats=mes(...,'nBoot',n)  computes confidence intervals of the statistic 
%   in question via bootstrapping n times. n should be on the order of
%   thousands, otherwise bootstrapping will lead to inaccurate results. By
%   default or if the value of nBoot is set to zero or nonsensical values
%   (<0, infinity, nan) confidence intervals will be computed analytically.
%   In cases where bootstrapping is not requested and computation of
%   analytical confidence intervals is not implemented the corresponding
%   fields of output struct stats will contain NaNs.
% stats=mes(...,'exactCi',true)  computes exact analytical confidence 
%   intervals for effect size measures for which both exact and approximate
%   CIs can be computed. Exact confidence intervals are based on iterative
%   determination of noncentrality parameters of noncentral Chi square, t
%   or F distributions. As this can be a very time-consuming process, and
%   as good bias corrections for small sample sizes are implemented for a
%   number of effect size measures, by default approximate analytical CIs
%   will be computed. See the documentation for details. Note that setting
%   this option is without effect if bootstrapping is requested (see input
%   variable 'nBoot' above).
% stats=mes(...,'confLevel',0.90)  computes 90 % confidence intervals (ci)
%   of the statistic in question (95 % ci are the default; any value may be
%   specified). If ci cannot be computed analytically and bootstrappig was
%   not requested, the corresponding fields of output struct stats (see
%   below) will contain NaNs.
% stats=mes(...,'ROCtBoot',1)  computes bootstrap confidence intervals for 
%   the area under the receiver-operating curve according to the 'bootstrap
%   t' method, which is more conservative than the 'bootstrap percentile'
%   method, the default (see the documentation for details). This option
%   will be ignored if bootstrapping is not requested
% stats=mes(...,'trCutoff',1.5)  sets 1.5 as the cutoff value expressed in 
%   standard deviations of the combined data set beyond the grand mean in
%   the computation of tail ratios. For positive values the right tail
%   ratio will be computed, for negative values the left tail ratio (see
%   documentation for further information). Default is 1.
% stats=mes(...,'trMeth','analytic')  determines that the tail ratios shall be computed
%   'analytically', assuming normal distributions, in the computation of
%   tail ratios. Default is 'count', which means that the ratios will be
%   determined by counting the actual data points beyond the cutoff
%   (relatively insensitive to deviations from normality)
% stats=mes(...,'doPlot',1)  will produce very simple plots of the results 
%   (one figure per effect size measure requested)
%
% -------------------------------------------------------------------------
% TABLE 1: ANALYSES TO BE SPECIFIED IN INPUT VARIABLE esm
% -------------------------------------------------------------------------
% esm           QUANTITIY COMPUTED
%        -- measures between one sample and a null value: --
% 'g1'          standardized difference (sample mean - comparison value)
% 'U3_1'        fraction of values below comparison value
%        -- measures between two samples: --
% 'md'          mean difference
% 'hedgesg'     Hedges' g (standardized mean difference)
% 'glassdelta'  Glass's delta (standardized mean difference)
% 'mdbysd'      mean difference divided by std of difference score
% 'requiv'      point-biserial correlation coefficient 
% 'cles'        common language effect size 
% 'U1'          Cohen's U1
% 'U3'          Cohen's U3
% 'auroc'       receiver-operating characteristic: area under curve
% 'tailratio'   tail ratios
% 'rbcorr'      rank-biserial correlation coefficient 
%
% 
%                         OUTPUT ARGUMENTS
%                         ----------------
%   The results of the computations are placed into fields of output
%   argument stats with names corresponding to the analyses requested. For 
%   example, stats=mes(X,Y,{'requiv','hedgesg'}) will result in 
%     stats.requiv 
%     stats.hedgesg
%   Where applicable, confidence intervals will be placed in fields 
%     .requivCi
%     .hedgesgCi
%   and the computations underlying confidence intervals (either of 'none',
%   'bootstrap', 'bootstrap t' (auroc only), 'approximate analytical', or
%   'exact analytical') will be placed in fields
%     .requivCiType
%     .hedgesgCiType
%   Additional fields:
%     .n (sample sizes)
%     .confLevel (confidence level)
%     .tstat (t test statistics) 
%   stats will also have fields .isDep and .nBoot with the same values as
%   the corresponding input arguments, thus providing potentially crucial
%   post-analysis information. 

% -------------------------------------------------------------------------
% Measures of Effect Size Toolbox Version 1.4, January 2015
% Code by Harald Hentschke (University of T?bingen) and 
% Maik St?ttgen (University of Bochum)
% For additional information see Hentschke and St?ttgen, 
% Eur J Neurosci 34:1887-1894, 2011
% -------------------------------------------------------------------------

% ----- default values & varargin -----
% standard values of optional input arguments (see description above)
isDep=false;
missVal='listwise';
nBoot=0;
ROCtBoot=false;
exactCi=false;
confLevel=.95;
trCutoff=1;
trMeth='count';
doPlot=false;

% check variable number of input args:
% - even number (parameter AND value specified)?
% - convert to proper upper/lowercase spelling as above
% - check whether valid parameter was given
% - overwrite defaults by values, if specified
nvarg=numel(varargin);
v=who;
if nvarg
  if nvarg/2==round(nvarg/2)
    for g=1:2:nvarg
      % check which optional input parameter is given, ignoring case, and
      % impose spelling used internally
      ix=find(strcmp(lower(varargin{g}),lower(v)));
      if ~isempty(ix)
        varargin{g}=v{ix};
      else
        error(['invalid optional input parameter (' varargin{g} ') was specified']);
      end
    end
    % finally, the assignment of values to parameters
    pvpmod(varargin);
  else
    error('optional input parameters must be specified as parameter/value pairs, e.g. as in ''par1'',1')
  end
end

% ----- a few 'constants':
alpha=1-confLevel;
% minimal number of bootstrapping runs below which a warning will be issued
% (and bootstrapping will be skipped)
minNBootstrap=1000;

% -------------------------------------------------------------------------
% ------- PART I: CHECKS OF INPUT ARGS & OTHER PREPARATORY WORKS ----------
% -------------------------------------------------------------------------
% Variable list_analyses below is a list of all possible analyses:
% - column 1 holds the shortcut as strings
% - column 2 contains a numerical tag coding for the kinds of samples 
%   (1=one-sample, 2=two-sample)
% - column 3 is a short description of the analyses
% This is the 'master list' of all analyses against which input argument
% esm will be checked (see) below
list_analysis={...
  'g1',          1, 'standardized difference (sample mean - comparison value)';...
  'U3_1',        1, 'fraction of values below comparison value';...
  'md',          2, 'mean difference (unstandardized)';...
  'hedgesg',     2, 'Hedges'' g (standardized mean difference)';...
  'glassdelta',  2, 'Glass''s delta (standardized mean difference)';...
  'mdbysd',      2, 'mean difference divided by std of difference score';...
  'requiv',      2, 'pointbiserial correlation coefficient';...
  'cles',        2, 'common language effect size';...
  'auroc',       2, 'receiver-operating characteristic';...
  'U1',          2, 'Cohen''s U1';...
  'U3',          2, 'Cohen''s U3';...
  'tailratio',   2, 'tail ratios';...
  'rbcorr',      2, 'rank-biserial correlation coefficient'...
  };
% in case esm is a char array convert it to a cell array because the code
% below relies on it
if ischar(esm)
  esm={esm};
end
% the tag column, extracted
listTag=cat(1,list_analysis{:,2});
% indices to currently requested analyses
listIx=ismember(list_analysis(:,1),esm);
% unique tags of these analyses: must be a scalar of value 1 (1-sample
% tests) or 2 (2-sample tests)
uTag=unique(listTag(listIx));
% catch a few silly errors:
% - typos or non-existent analysis type
if ~any(listIx)
  error('An illegal type of analysis was specified');
end
% - mixed one- and two-samples tests
if numel(uTag)>1
  error('A mix of one-sample and two-sample tests was requested')
end
% - programmer's error 
if isempty(intersect(uTag,[1 2]))
  error('internal: uTag different from 1 or 2');
end
% illegal values for missVal
if isempty(find(strcmpi(missVal,{'listwise','pairwise'})))
  error('illegal value for input parameter ''missVal'' (choose ''listwise'' or ''pairwise)''');
end

% --- check input X and Y
[nRowX nColX]=size(X);
[nRowY nColY]=size(Y);
% reshape X and Y in a few selected scenarios
% - if X is a single row array
if nRowX==1 && nColX>1
  warning('input variable X is a single row array - reshaping');
  X=X(:);
  [nRowX nColX]=size(X);
end
% - if X is a single column array and Y a single row array, reshape Y as 
% well
if nColX==1 && nRowY==1 && nColY>1
  warning('input variable Y is a single row-array - reshaping')
  Y=Y(:);
  [nRowY nColY]=size(Y);
end
% (if Y is a single row array while X is not the situation is ambiguous, so
% we'd better let the code produce an error further below)

% scalar expansion: if X has several columns and Y is a scalar, expand it
if nColX>1 && nRowY==1 && nColY==1
  Y=repmat(Y,1,nColX);
  nColY=nColX;
end
% now perform strict checks
if nColX~=nColY
  error('input variables X and Y must have the same number of columns');
end
if isDep
  if nRowX~=nRowY
    error('for paired data input variables X and Y must have the same number of rows');
  end
end
if uTag==1 && nRowY>1
  error('1-sample analyses require input variable y to be a scalar (or a row array with as many columns as x)');
end
  
% deal with foul values:
% - first, convert infs to nans
X(~isfinite(X))=nan;
Y(~isfinite(Y))=nan;
[nanRowX,nanColX]=find(isnan(X));
[nanRowY,nanColY]=find(isnan(Y));
% - depending on how missing values shall be treated, set selected
% values or entire rows to NaN. Differentiate between 1-sample and 2-sample
% analyses
if ~isempty(nanRowX) || ~isempty(nanRowY)
  if uTag==1 
    % 1-sample case (easy)
    if ~isempty(nanRowY)
      error('in 1-sample analyses Y is not allowed to contain NaN or Inf');
    end
    if strcmp(lower(missVal),'listwise')
      X(nanRowX,:)=nan;
    end
  else
    % 2-sample case
    switch lower(missVal)
      case 'listwise'
        if isDep
          nanRowX=union(nanRowX,nanRowY);
          nanRowY=nanRowX;
        end
        X(nanRowX,:)=nan;
        Y(nanRowY,:)=nan;
      case 'pairwise'
        if isDep
          badIx=union(sub2ind([nRowX nColX],nanRowX,nanColX),sub2ind([nRowY nColY],nanRowY,nanColY));
          X(badIx)=nan;
          Y(badIx)=nan;
        end
        % nothing to do in case of unpaired data and pairwise elemination
    end
  end
end
% temp variables not needed anymore
clear badIx nanRowX nanColX nanRowY nanColY;

% --- test-specific checks
% - mostly, we must check for the 'isDep' option: 
% a) for some analyses, data in x and y are usually not considered 'paired'
% (e.g. ROC), so by issuing an error or warning the user is encouraged to
% rethink his/her analysis
% b) the exact procedure/results of bootstrapping and t-test depend on
% whether the data are paired or not, hence a correct specification may be
% essential
if any(ismember(esm,'auroc'))
  if isDep
    warning('''auroc'' (receiver-operating characteristic) is implemented for independent samples; confidence intervals for dependent data are likely not correct (see parameter ''isDep'')');
  end
end
if any(ismember(esm,'glassdelta'))
  if isDep
    warning('''glassdelta'' does not make sense for dependent samples (see parameter ''isDep'')');
  end
end
if any(ismember(esm,'mdbysd'))
  if ~isDep
    error('''mdbysd'' (mean difference divided by std of difference score) is not defined for independent samples (see parameter ''isDep'')');
  end
end

if strcmpi(trMeth,'analytic')
  doAnalyticalTails=true;
elseif strcmp(lower(trMeth),'count')
  doAnalyticalTails=false;
else
  error('illegal value specified for input parameter ''trMeth''');
end

% --- check bootstrapping settings
doBoot=false;
if isfinite(nBoot)
  if nBoot>=minNBootstrap;
    doBoot=true;
  else
    if nBoot~=0
      % warn only if nBoot small but different from zero because zero may
      % by a deliberate input value
      warning('number of bootstrap repetitions is not adequate - not bootstrapping');
    end
    nBoot=0;
  end
end
if ROCtBoot && ~doBoot
  warning('option ''ROCtBoot'' is only effective if bootstrapping is requested - check input parameter ''nBoot''');
end

% --- check other input arguments
if confLevel<=0 || confLevel>=1
  error('input variable ''confLevel'' must be a scalar of value >0 and <1');
end

% -------------------------------------------------------------------------
% ------- PART II: THE WORKS ----------
% -------------------------------------------------------------------------
% The outermost loop cycles through all columns of X (and matching columns
% of Y), performs bootstrapping (if requested) and then runs the tests

% preallocate results arrays:
% - ttest associated parameters:
allocationNanny=repmat(nan,[1 nColX]);
% value of test statistic 
stats.t.tstat=allocationNanny;
% p values (which will be kept only for original data)
stats.t.p=allocationNanny;
% degrees of freedom
stats.t.df=allocationNanny;
% create a few standard fields
stats.isDep=isDep;
stats.nBoot=nBoot;
stats.confLevel=confLevel;
stats.n=[nRowX; nRowY];
% reorder field names: tstats moves from first place to last
stats=orderfields(stats,[2:numel(fieldnames(stats)) 1]);

% ** effect size and confidence intervals cannot be preallocated here
% because they depend on the exact tests required (will be done at first
% run of loop)

% waitbar only if X (and Y) has more than one column
if nColX>1
  wbH=waitbar(0,'Computing...','name','Please wait');
end

% loop over columns of X (and Y)
for g=1:nColX
  % pick non-nan values in columns
  x=X(~isnan(X(:,g)),g);
  y=Y(~isnan(Y(:,g)),g);
  n1=size(x,1);
  n2=size(y,1);

  % -----------------------------------------------------------------------
  % ------------------ BOOTSTRAPPING (IF REQUESTED) -----------------------
  % -----------------------------------------------------------------------
  if doBoot
    % *********************************************************************
    % Here, x, representing individual columns of variable X, is expanded
    % in the second dimension: the first column contains the original data,
    % the others will be filled up with sampled (with replacement) data. y
    % will be expanded in the same (or similar) manner in case of 2-sample
    % tests; in case of 1-sample tests its (scalar) value will be
    % replicated to nBoot+1 columns 
    % *********************************************************************
    % indices to be used for randomly sampling (with replacement) the
    % original data
    bootSamX=ceil(n1*rand([n1 nBoot]));
    % if data are paired use same indices for y; if they are unpaired
    % compute additional one
    if isequal(uTag,2) && ~isDep
      bootSamY=ceil(n2*rand([n2 nBoot]));
    end
    % from second column onwards fill with sampled data
    x(:,2:nBoot+1)=x(bootSamX);
    if isequal(uTag,2)
      % do the same with y in case of 2-sample tests
      if isDep
        y(:,2:nBoot+1)=y(bootSamX);
      else
        y(:,2:nBoot+1)=y(bootSamY);
      end
    elseif isequal(uTag,1)
      % expand y in case of 1-sample tests 
      y(:,2:nBoot+1)=y;
    end
    % delete index
    clear bootSam*
  end
    
  % -----------------------------------------------------------------------
  % ----------------- ES COMPUTATIONS -------------------------------------
  % -----------------------------------------------------------------------
  % In this section, first some some preparatory computations will be
  % performed, notably terms needed for t statistics and several mes. Then,
  % the FOR loop below will sequentially work on all types of analysis as 
  % listed in variable esm.

  % preparatory computations: 
  % - means, variances, number of samples
  m1=mean(x,1);
  s1=var(x,0,1); % var and std are by default divided by n-1
  m2=mean(y,1);
  s2=var(y,0,1);
  % exand n1 and n2
  n1=repmat(n1,[1 nBoot+1]);
  n2=repmat(n2,[1 nBoot+1]);

  % - degrees of freedom, t statistics, sd of diff score & pooled var
  if isDep || isequal(uTag,1)
    % df in dependent case: n=n1-1=n2-1
    df=n1-1;
    % standard deviation of difference score
    if uTag==1
      % in 1-sample analyses, stdD=sd(x)
      stdD=sqrt(s1);
    else
      stdD=std(x-y);      
    end
    % standard error
    seD=stdD./sqrt(n1);
    % t statistic
    tst=(m1-m2)./seD;
  else
    % independent case: n=[n1+n2-2]
    df=n1+n2-2;
    % pooled (within-groups) variance
    sP=((n1-1).*s1 + (n2-1).*s2)./(n1+n2-2);
    % standard error 
    seP=sqrt(sP.*((n1+n2)./(n1.*n2)));
    % t statistic
    tst=(m1-m2)./seP;
  end
  % p value
  p=2*(1-tcdf(abs(tst),df(1)));
  
  % in output struct stats.t keep only first element of t statistic 
  % computed from original data
  stats.t.tstat(g)=tst(1);
  % ditto for p and sd
  stats.t.p(g)=p(1);
  stats.t.df(g)=df(1);

  % - inverse cumulative t distribution:
  % n1, n2 and therefore df are identical for the original data and the
  % sampled (bootstrapped) data, so compute inverse cumulative t
  % distribution only from the former. However, as the computations below
  % rely on ciFac having dimensions [1 by size(x,2)] expand it here
  ciFac=repmat(-tinv(alpha/2,df(1)),[1 size(x,2)]);
  
  % ***********************************************************************
  %                  loop over all ES computations
  % ***********************************************************************
  % PROGRAMMER'S NOTES:
  % -> At this point in the code, we are within the outermost FOR loop 
  %    (mentioned above) which cycles through the columns of X, g being the 
  %    column index (see code line saying 'for g=1:nColX') 
  % -> Variable x represents data from an individual column of X. In case
  %    of 1-sample tests, y is the corresponding null value; in case of
  %    2-sample tests, it contains data from the corresponding column of Y.
  % -> If bootstrapping was requested, x has [nBoot+1] columns: the
  %    original data are in column 1 and the bootstrapped data are in 
  %    columns 2 and above. The same holds for y in case of 2-sample tests. 
  % -> In case of bootstrapping and 1-sample tests, x has [nBoot+1] columns
  %    as described above; y is a [1 by nBoot+1] array full of identical
  %    values (namely, the null value of the current column of interest in
  %    x): it has been expanded above to facilitate computations here.
  % -> To sum up, the columns in x and y always match and correspond to
  %    each other; thus any new code to be inserted below will work if it
  %    takes the potentially 2D layout of x and y, with the standard Matlab
  %    columnar organization of data, into account 
  % -> sample sizes, means and variances of x and y have been precomputed
  %    above (variables n, m, s, respectively) and need not (in fact should 
  %    not) be computed here
  % -> same holds for t statistics (variable tst, see case 'requiv')
  % -> if a formula for analytical confidence intervals is known, the code
  %    must be embedded (for each statistic of course) in 
  %       if ~doBoot
  %       end
  %    (see e.g. case 'hedgesg'). If they cannot be computed analytically,
  %    NaNs are inserted in the results variable ci. If bootstrapping was
  %    requested, the NaNs will be replaced by meaningful values
  % -> nRowX and nColX are the number or rows and columns, respectively, of 
  %    the original input variables. As the loop hops through the columns
  %    of input variable X, x is by definition either a single-column array
  %    or an array with nBoot columns - not with nColX columns! Moreover,
  %    as any NaNs will have been eliminated here from x, n1 is the number
  %    of rows in x, not nRowX. To make a long story short, n1 can be used
  %    for e.g. preallocation or the like, and make sure in your new code
  %    that you do not confuse nBoot and nColX. size(x,2) is the safest
  %    bet. Of course, the same applies to Y and y.
  
  nEs=numel(esm);
  for tti=1:nEs
    curEs=esm{tti};
    % If analytical confidence intervals for the statistic in question
    % cannot be computed and no bootstrapping is requested, ci must be set
    % to nan. Do this here for all statistics to avoid lots of redundant
    % code. ci will be overwritten with meaningful values if i) there is a
    % formula for analytical confidence intervals implemented within the
    % respective case in the switch statement, or ii) confidence intervals
    % based on bootstrapping are computed right after the switch statement
    ci=[nan; nan];
    % a similar argument applies to variable ciType, which indicates on
    % which method the computation of ci is based
    if g==1
      if doBoot
        ciType='bootstrap';
      else
        ciType='none';
      end
    end
    
    switch curEs
      case 'g1'
        es=(m1-y)./sqrt(s1);
        if ~doBoot
          % exact ci (Smithson 2003, p.34, formula 4.4)
          ci=ncpci(tst,'t',n1-1,'confLevel',confLevel)'/sqrt(n1);
          if g==1
            ciType='exact analytical';
          end          
        end

      case 'U3_1'
        % number of values in x below (scalar value) y
        es=sum(x<repmat(y,[n1(1) 1]))./n1;
        % count values on threshold half
        es=es+.5*sum(x==repmat(y,[n1(1) 1]))./n1;
        
      case 'md'
        % mean difference (unstandardized), the most basic mes
        % imaginable...
        es=m1-m2;
        if ~doBoot
          if isDep
            % se=standard error of difference scores
            ci=cat(1,es-ciFac.*seD,es+ciFac.*seD);
          else
            % se=standard error of mean difference
            ci=cat(1,es-ciFac.*seP,es+ciFac.*seP);
          end
          if g==1
            % confidence intervals of mean differences computed from
            % central t distributions are not approximate, hence we term
            % them exact here, although exactness does not, as in the case
            % of other mes, imply computation via noncentral distributions
            ciType='exact analytical';
          end
        end
        
      case 'mdbysd'
        % mean difference divided by std of difference score (defined only
        % for dependent data)
        es=(m1-m2)./std(x-y);
        if ~doBoot
          % exact ci (Smithson 2003, p.34-46, formula 4.4)
          ci=ncpci(tst,'t',n1-1,'confLevel',confLevel)'/sqrt(n1);
          if g==1
            ciType='exact analytical';
          end
        end
              
      case 'hedgesg'
        % Hedges' g
        if isDep
          % n=n1=n2
          es=tst.*sqrt(2*stdD.^2./(n1.*(s1+s2)));
        else
          es=(m1-m2)./sqrt(sP);
        end
        % correct for bias due to small n (both dependent and independent
        % data, Kline 2004 (p. 102, 106); Nakagawa & Cuthill 2007)
        biasFac=(1-(3./(4*n1+4*n2-9)));
        es=es.*biasFac;
        if ~doBoot
          if isDep
            % approximate ci for paired data (Nakagawa & Cuthill 2007) -
            % note that n=n1=n2 and that correlation coeff r is needed; no
            % bias correction here
            se=sqrt((2-2*corr(x,y))./n1 + es.^2./(2*n1-1));
            ci=cat(1,es-ciFac.*se,es+ciFac.*se);
            if g==1
              ciType='approximate analytical';
            end
          else
            if exactCi
              % exact ci (Smithson 2003, p. 37), including bias correction
              ci=biasFac*ncpci(tst,'t',n1+n2-2,'confLevel',confLevel)'*sqrt((n1+n2)/(n1*n2));
              if g==1
                ciType='exact analytical';
              end
            else
              % approximate ci (Nakagawa & Cuthill 2007, eq. 17 in table
              % 3), including bias correction
              se=sqrt((n1+n2)./(n1.*n2) + (es.^2./(2*n1+2*n2-4)));
              ci=biasFac*cat(1,es-ciFac.*se,es+ciFac.*se);
              if g==1
                ciType='approximate analytical';
              end
            end
          end
        end
        
      case 'glassdelta'
        es=(m1-m2)./sqrt(s1);
        % analytical confidence intervals: only approximate
        if ~doBoot
          se=sqrt(es^2/(2*n2-2)+(n1+n2)/(n1*n2));
          ci=cat(1,es-ciFac.*se,es+ciFac.*se);
          if g==1
            ciType='approximate analytical';
          end
        end
        
      case 'requiv'
        % note: an alternative computation would work as follows: assign
        % discrete numbers to the two different groups (say, 0 to x and 1
        % to y), concatenate the group indices and the real data and then
        % plug the resulting two variables into corr (Pearson's
        % correlation). This yields identical results, but is
        % computationally more demanding and also more complicated when x
        % and y are 2D (bootstrapped). 
        es=tst./sqrt(tst.^2 + n1+n2-2);
        if ~doBoot
          % if data are independent & exact analytical ci requested...
          if ~isDep && exactCi
            % exact analytical confidence intervals for partial eta2, of
            % which requiv is a special case (Smithson 2003, p.43 [5.6])
            ci=ncpci(tst^2,'F',[1 stats.t.df],'confLevel',confLevel)';
            % don't forget to sqrt and to take care of sign
            ci=sqrt(ci./(ci+1+stats.t.df+1));
            if es<0
              ci=-1*fliplr(ci);
            end
            if g==1
              ciType='exact analytical';
            end
          else
            % all other cases: approximate CI via Z-transform
            tmp=.5*log((1+es)/(1-es));
            tmp=tmp+[-1; 1]*ciFac./sqrt(n1+n2-3);
            % transform back
            ci=(exp(2*tmp)-1)./(exp(2*tmp)+1);
            if g==1
              ciType='approximate analytical';
            end
          end
        end

      case 'cles'
        % 'For continuous data, it is the probability that a score sampled
        % at random from one distribution will be greater than a score
        % sampled from some other distribution'
        es=normcdf((m1-m2)./sqrt(s1+s2));
        
      case 'U1'
        % line arrays of minimal and maximal values in x
        maxX=max(x);
        minX=min(x);        
        % ditto for y
        maxY=max(y);
        minY=min(y);        
        es=(sum(x>repmat(maxY,n1(1),1))+sum(y>repmat(maxX,n2(1),1))...
          +sum(x<repmat(minY,n1(1),1))+sum(y<repmat(minX,n2(1),1)))./(n1+n2);
        
      case 'U3'
        % we need the medians of both groups
        med1=median(x);
        med2=median(y);
        es=sum(x<repmat(med2,n1(1),1))./n1;
        % count identical values half
        es=es+.5*sum(x==repmat(med2,n1(1),1))./n1;
        % in case both medians are equal, U3 must by definition be 0.5, but
        % the code lines above may have resulted in divergent values if
        % the lower group (x) contains heavily aliased data (e.g. histogram
        % data with one dominant bin). The following two lines correct for
        % this
        tmpIx=med1==med2;
        es(tmpIx)=.5;
        
      case 'auroc'
        % compute area under curve using approach in Bamber 1975 (J Math
        % Psych 12:387-415), eq. 3, also presented in Hanley & McNeil
        % Radiology 143:29-36, 1982 (p. 31)
        es=zeros(1,nBoot+1);
        % loop over elements (rows) in x
        for xIx=1:n1
          tmp=repmat(x(xIx,:),n2(1),1);
          es=es+sum(tmp>y)+0.5*sum(tmp==y);
        end
        es=es./(n1.*n2);
        if ~doBoot
          % analytical confidence intervals: start by computing standard
          % error of AUROC according to Hanley & McNeil (Radiology
          % 143:29-36, 1982), who refer to the classic work by Bamber 1975
          % (J Math Psych 12:387-415)
          se=se_auroc(es,n1,n2);
          % note that normality is assumed in this step
          tmp=norminv(1-alpha/2).*se;
          ci=cat(1,es-tmp,es+tmp);
          if g==1
            ciType='approximate analytical';
          end
        else
          if ROCtBoot
            % this portion of code performs preparatory steps for the
            % 'bootstrap t' method of estimating CI (Obuchowski & Lieber
            % Acad Radiol 5:561-571, 1998), using the approach of Hanley &
            % McNeil (Radiology 143:29-36, 1982) for the estimation of each
            % bootstrapped data set's standard error
            % - compute standard error of original and bootstrapped sets
            se=se_auroc(es,n1,n2);
            % - overwrite the bootstrapped auroc values with the
            % studentized pivot statistic
            es(2:end)=(es(2:end)-es(1))./se(2:end);
            % - replace Infs resulting from bootstrapped cases with zero
            % se by NaNs because Matlab function prctile will eliminate
            % these
            es(isinf(es))=nan;
            % - finally, multiply by se of original data and add auroc of
            % original data. This way, the CI can be directly extracted
            % from es(2:end), as is done with all other effect size
            % measures, too
            es(2:end)=es(2:end)*se(1)+es(1);
            if g==1
              ciType='bootstrap t';
            end
          else
            % nothing to do here: CI will be extracted below, which
            % corresponds to the standard 'percentile bootstrap' method
          end
        end
        clear tmp;
        
      case 'tailratio'
        % total mean
        m_t=(n1.*m1 + n2.*m2)./(n1+n2);
        % total std (compute in two steps)
        std_t=...
          n1.*(m1-m_t).^2 + n2.*(m2-m_t).^2 + (n1-1).*s1 + (n2-1).*s2;
        std_t=sqrt(std_t./(n1+n2-1));
        % (that one would work as well: std(cat(1,x,y)))
        % the cutoff (threshold)
        cutoff=m_t+trCutoff*std_t;
        % now compute the proportions of samples...
        if doAnalyticalTails
          % ...the 'analytical' way (as in Kline p.126) - this has the
          % advantage that in cases where two samples are completely
          % separated we don't get Infinity as result. However, this really
          % requires normality
          if trCutoff>=0
            es=(1-normcdf((cutoff-m1)./sqrt(s1)))./(1-normcdf((cutoff-m2)./sqrt(s2)));
          else
            es=(normcdf((m1-cutoff)./sqrt(s1)))./(normcdf((m2-cutoff)./sqrt(s2)));
          end
        else
          % ...by counting
          if trCutoff>=0
            % right tail ratio
            es=(sum(x>ones(n1(1),1)*cutoff)./n1)./(sum(y>ones(n2(1),1)*cutoff)./n2);
          else
            % left tail ratio
            es=(sum(x<ones(n1(1),1)*cutoff)./n1)./(sum(y<ones(n2(1),1)*cutoff)./n2);
          end
        end
        % note that there may be division by zero, particularly in
        % bootstrapped data, which prevents proper computation of
        % confidence intervals, so replace the infs by nans
        es([false isinf(es(2:end))])=nan;
        
      case 'rbcorr'
        % rank-biserial correlation coefficient
        % ? compared to all other analyses this one takes an awfully long
        % time because of function corr and the functions it calls
        % (tiedrank and so on) - possibly this could be sped up
        es=corr(cat(1,x,y),cat(1,zeros(n1(1),1),ones(n2(1),1)),'type','Spearman');
                
      otherwise
        error(['internal: unrecognized test not caught by input checks: ' curEs]);
    end
    
    % *********************************************************************
    % If data were NOT bootstrapped, all computations are done at this
    % point and the results can be placed into appropriate fields of output
    % variable stats. If they were bootstrapped, confidence intervals must
    % be computed and the ES statistics extracted from the first element of
    % variable es.
    % *********************************************************************

    % start by preallocating major results fields of output variable stats
    if g==1
      stats.(curEs)=allocationNanny;
      stats.([curEs 'Ci'])=[allocationNanny; allocationNanny];
    end
    
    if doBoot
      % determine confidence intervals from array of effect size measures
      % generated from bootstrapped data
      ci=prctile(es(2:end),[alpha/2  1-alpha/2]'*100);
      % retain first element; this is the es computed from real data
      es=es(1);
    end
    
    % finally, use dynamic fields to store currently computed measures in
    % output variable stats
    stats.(curEs)(g)=es;
    stats.([curEs 'Ci'])(:,g)=ci;
    if g==1
      stats.([curEs 'CiType'])=ciType;
    end
  end
  
  % update waitbar
  if nColX>1
    waitbar(g/nColX,wbH);
  end
end

% kill waitbar window
if nColX>1
  delete(wbH);
end

% call simple plot routine if requested
if doPlot
  simpleEsPlot(stats,esm);
end

% ========================= LOCAL FUNCTIONS ===============================
% ========================= LOCAL FUNCTIONS ===============================

function se=se_auroc(es,n1,n2)
% ** function se_auroc(es,n1,n2) computes standard error of AUROC according
% to Hanley & McNeil Radiology 143:29-36, 1982
% let's break down the computations somewhat & stick to their terminology
q1=es./(2-es);
q2=2*es.^2./(1+es);
es2=es.^2;
se=sqrt((es-es2+(n1-1).*(q1-es2)+(n2-1).*(q2-es2))./(n1.*n2));


function simpleEsPlot(stats,esm)
% a simple routine for plotting results
nEs=numel(esm);
for tti=1:nEs
  curEs=esm{tti};
  % assign values to generic local variables es and ci
  curEsCi=[curEs 'Ci'];
  es=stats.(curEs);
  ci=stats.(curEsCi);
  % one figure per statistic
  figure(tti)
  clf
  hold on
  if numel(es)==1
    plot(1,es,'ko');
    plot([1 1],ci,'b+');
  else
    plot(es,'ko-');
    plot(ci','b-');
  end
end
  
  
% ======================= REPOSITORY =======================================

% The code below is a 'traditional' and inefficient way to compute the area
% under the curve of the ROC

%         % in the following lines, we generate an array of criterion values,
%         % equally spaced from the minimum to the maximum in [x and y
%         % combined], with at most rocNBin values. The array thus
%         % generated is used for both original and bootstrapped (if any)
%         % data.
%         xyso=sort([x(:,1); y(:,1)]);
%         xysod=unique(sort(diff(xyso)));
%         % zero is not an acceptable bin width
%         if ~xysod(1)
%           xysod(1)=[];
%         end
%         % number of bins: use the lesser of [rocNBin, max-min divided
%         % by smallest acceptable difference between values]
%         nBin=min(rocNBin,round((xyso(end)-xyso(1))/xysod(1)));
%         % with the bins defined like this and function histc, used below,
%         % the last bin contains only values falling exactly on the last
%         % value. Shift the border of the last bin by a tiny amount to the
%         % right to ensure that i) the last bin will contain a count of zero
%         % (and can be discarded), and ii) the last but one bin will contain
%         % the maximal value in the sample(s) without extending far beyond
%         % this value
%         critVal=linspace(xyso(1),xyso(end)+eps(xyso(end)),nBin+1);
%         % now compute the cum histograms - separately for x and y
%         % because their sizes (number of rows) may differ
%         xch=cumsum(histc(x,critVal))./repmat(n1,[nBin+1 1]);
%         ych=cumsum(histc(y,critVal))./repmat(n2,[nBin+1 1]);
%         % the first bin always contains a nonzero value. We need to have a
%         % first bin with zero count (needed for proper calculation of the
%         % area under the curve), so append these
%         xch=cat(1,zeros(1,size(xch,2)),xch);
%         ych=cat(1,zeros(1,size(ych,2)),ych);
%         % close the curves by replacing the useless values in the last bin
%         % by coordinate (defined in x-y space) [1 0]
%         xch(end,:)=ones(1,size(xch,2));
%         ych(end,:)=zeros(1,size(ych,2));
%         % compute area under curve using polyarea 
%         es=polyarea(xch,ych);

% here are the 95% CI based on smax: Bamber, J Math Psychol 1975
% % -> appear unrealistically narrow
% smax=es.*(1-es)./(min(n1,n2)-1);
% tmp=norminv(1-alpha/2).*smax;



function [stats,varargout]=mes1way(X,esm,varargin)
% ** function [stats,varargout]=mes1way(X,esm,varargin)
% computes measures of effect size between two or more samples and places
% the results in output structure stats. A summary table including the
% effect size measures and one-way ANOVA results will be displayed in the
% command window. All input parameters except X and esm are optional and
% must be specified as parameter/value pairs in any order, e.g. as in
%      mes1way(X,'eta2','group',g);
%
% For information on assumptions of the underlying model and related
% information please see the notes below TABLE 1 further down as well as
% the documentation accompanying this code.
%
%                           INPUT ARGUMENTS
%                           ---------------
% stats=mes1way(X,esm)  computes the effect size measure(s) specified in 
%   input variable esm (see table at bottom) from input array X. X must
%   have two or more columns, each representing a different group, or may
%   be a single column-array (see below). Input variable esm must be a
%   character array (e.g. 'eta2') or a cell array (e.g.
%   {'eta2','partialeta2'}).
% stats=mes1way(...,'group',g)  allows for an alternative structure of 
%   input data X. X must be a single-column vector and g a single-column
%   vector of arbitrary numbers coding for the different groups.
% stats=mes1way(...,'isDep',1)  assumes that the samples in X are dependent
%   (repeated measures); accordingly, the number of samples in each group
%   of X must be equal. If this parameter is omitted the assumption is one
%   of independent samples. Both the effect size measures available and
%   their computation are contingent on whether the samples are independent
%   or dependent. NOTE: if data is dependent, it is assumed that input data
%   X are sorted according to subjects. That is, if X is a multi-column-
%   array all data from each subject reside in one row. If X is a single-
%   column array, the equivalent assumption is made: that data for subjects
%   are listed in the order in which they occur within each group.
% stats=mes1way(...,'missVal','listwise')  omits NaNs and Infs in X,
%   summarily termed here as missing values, in listwise fashion: if X is a
%   multi-column array and there is a missing value anywhere in X, the
%   entire row of X will be omitted from analysis. If X is in single-column
%   format, positions corresponding to the row(s) in question will be
%   omitted. If set to 'pairwise', only the individual missing data will be
%   omitted. Note: 'listwise' is mandatory (and therefore the default) for
%   dependent data. Likewise, 'pairwise' is mandatory and default for
%   independent unbalanced data. It is also the default for independent,
%   balanced data.
% stats=mes1way(...,'nBoot',n)  computes confidence intervals of 
%   the statistic in question via bootstrapping n times. n should be on the
%   order of thousands, otherwise bootstrapping will lead to inaccurate
%   results. By default or if the value of nBoot is set to zero or
%   nonsensical values (<0, infinity, nan) confidence intervals will be 
%   computed analytically (where possible).
% stats=mes1way(...,'confLevel',0.90)  computes 90 % confidence intervals 
%   of the statistic in question (95 % ci are the default; any value may be 
%   specified)
% stats=mes1way(...,'cWeight',c)  allows specification of contrast weights
%   for the computation of effect size measures like standardized contrasts
%   and eta squared for focussed comparisons (see TABLE 1). Input array c
%   must contain as many columns as there are groups; you may specifiy as
%   many rows (=contrasts) as you wish.
% stats=mes1way(...,'tDenom','msw')  computes, in case of dependent data, F
%   and p values of contrasts, as well as confidence intervals of psi and
%   g_psi, from MS_[between x subject]. The default is 'sd', which means
%   that the values are based on the contrasts' difference score. Please
%   see the documentation for detailed explanation.
% -------------------------------------------------------------------------
% TABLE 1: ANALYSES TO BE SPECIFIED IN INPUT VARIABLE esm
% -------------------------------------------------------------------------
% esm             QUANTITIY COMPUTED ((*) require input of contrast weights)
% 'psi'           unstandardized contrast (*)
% 'g_psi'         standardized contrast (*)  
% 'psibysd'       standardized contrast for dependent data (*)
% 'eta2'          eta squared
% 'partialeta2'   partial eta squared 
% 'omega2',       omega squared
% 'partialomega2' partial omega squared
%
% Notes:
%   i. Parameters marked by (*) need contrast weights to be
%   computed; by definition they do not exist for an omnibus effect. The
%   other parameters will be computed for both the omnibus effect and the
%   specified contrasts.
%  ii. 'g_psi' is the oneway equivalent of Hedges' g, namely, contrast
%   divided by the square root of the pooled within-conditions variance.
%   Its value is identical for dependent and independent data, but for
%   dependent data confidence intervals will be smaller in case of
%   correlations between the groups.
% iii. 'psibysd' is a standardized contrast for dependent data only, 
%   namely, contrast divided by the standard deviation of the difference
%   score. It is the oneway equivalent of mdbysd in mes.m.
% iv. Fixed factors are assumed.
% v. Independent data may be unbalanced.
% vi. In the computation of omega2 and partial omega2 for dependent data,
%   an additive model (no subject x treatment interaction) is assumed
%
%                         OUTPUT ARGUMENTS
%                         ----------------
% The results of the computations are placed into fields of output
% structure stats with names corresponding to the analyses requested. For
% example, stats=mes1way(X,{'eta2','partialeta2'}) will result in fields
%   .eta2 
%   .partialeta2
% PLEASE NOTE: if contrast weights were specified, all output arguments
% will be single- or two-column arrays holding results in the following
% (row) order:
% 1. omnibus effect
% 2. first contrast
% 3. second contrast
% and so on. If a parameter is not defined for the omnibus effect or for
% contrasts, the corresponding entries will be NaN.
% Where applicable, confidence intervals will be placed in fields
%   .eta2Ci 
%   .partialeta2Ci
% and the computations underlying confidence intervals (either of 'none',
% 'bootstrap', 'approximate analytical', or 'exact analytical') will be
% placed in fields
%   .eta2CiType 
%   .partialeta2CiType
% Additional fields:
%   .n (sample sizes)
%   .confLevel (confidence level)
%   .cWeight (the contrast weights)
% stats will also have fields .isDep and .nBoot with the same values as the
% corresponding input arguments, thus providing potentially crucial post-
% analysis information. The second, optional output argument is the summary
% table of results.

% -------------------------------------------------------------------------
% Measures of Effect Size Toolbox Version 1.4, January 2015
% Code by Harald Hentschke (University of T?bingen) and 
% Maik St?ttgen (University of Bochum)
% For additional information see Hentschke and St?ttgen, 
% Eur J Neurosci 34:1887-1894, 2011
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%                       STRUCTURE OF CODE
% The code is composed of two major parts: In PART I input arguments are
% checked and pre-processed, 'constants' set up, and all kinds of other
% preparatory and householding jobs accomplished. PART II contains the
% computations proper. Right at the beginning, summed squares, F and p
% values are computed in local function prepcomp, which also accomplishes
% assembly of bootstrapped data if requested. Effect sizes and confidence
% intervals are then computed and the results placed in structure array
% 'stats' as described above. The results are additionally collected in
% cell array 'table' which is displayed in the command window (and may also
% be retrieved as an output argument)
% -------------------------------------------------------------------------

% ----- default values & varargin -----
% standard values of optional input arguments
group=[];
isDep=false;
missVal='pairwise';
nBoot=0;
cWeight=[];
confLevel=.95;
tDenom='sd'; 

% check variable number of input args:
% - even number (parameter AND value specified)?
% - convert to proper upper/lowercase spelling as above
% - check whether valid parameter was given
% - overwrite defaults by values, if specified
nvarg=numel(varargin);
v=who;
if nvarg
  if nvarg/2==round(nvarg/2)
    for g=1:2:nvarg
      % check which optional input parameter is given, ignoring case, and
      % impose spelling used internally
      ix=find(strcmpi(varargin{g},v));
      if ~isempty(ix)
        varargin{g}=v{ix};
      else
        error(['invalid optional input parameter (' varargin{g} ') was specified']);
      end
    end
    % finally, the assignment of values to parameters
    pvpmod(varargin);
  else
    error('optional input parameters must be specified as parameter/value pairs, e.g. as in ''par1'',1')
  end
end

% ----- a few 'constants':
alpha=1-confLevel;
% minimal number of bootstrapping runs below which a warning will be issued
% (and bootstrapping will be skipped)
minNBootstrap=1000;

% -------------------------------------------------------------------------
% ------- PART I: CHECKS OF INPUT ARGS & OTHER PREPARATORY WORKS ----------
% -------------------------------------------------------------------------
% Variable list_analyses below is a list of all possible analyses:
% - column 1 holds the shortcut as strings
% - column 2 holds 1 for quantities REQUIRING a contrast
% - column 3 is a short description of the analyses
% This is the 'master list' of all analyses against which input argument
% esm will be checked (see) below
list_analysis={...
'psi',           1, 'unstandardized contrast';...
'g_psi',         1, 'standardized contrast';...
'psibysd',       1, 'standardized contrast for dependent data';...
'eta2',          0, 'eta squared';...
'omega2',        0, 'omega squared';...
'partialeta2',   0, 'partial eta squared';...
'partialomega2', 0, 'partial omega squared';...
};

% in case esm is a char array convert it to a cell array because the code
% below relies on it
if ischar(esm)
  esm={esm};
end
tmp=find(strcmp('psibysd',esm));
if ~isempty(tmp) && ~isDep
  warning('parameter ''psibysd'' can only be computed for dependent data - eliminating it from list of parameters to be computed')
  esm(tmp)=[];
end
% the tag column, extracted
listTag=cat(1,list_analysis{:,2});
% indices to currently requested analyses
listIx=ismember(list_analysis(:,1),esm);
% unique tags of these analyses
uTag=unique(listTag(listIx));
% catch a few silly errors:
% - typos or non-existent analysis type
if any(~ismember(esm,list_analysis(:,1)))
  error('An illegal type of analysis was specified');
end

% check whether contrast weight(s) were specified, whether they are
% required, and set according flags (exact checks of contrast weights will
% be done further below)
if isempty(cWeight)
  doContrast=false;
  if any(uTag)
    error('part of the requested analyses require contrast weights as input args')
  end
else
  doContrast=true;
end

% illegal values for missVal
if isempty(find(strcmpi(missVal,{'listwise','pairwise'})))
  error('illegal value for input parameter ''missVal'' (choose ''listwise'' or ''pairwise)''');
end
if isDep && ~strcmpi(missVal,'listwise')
  warning('dependent data - setting input parameter ''missVal'' to ''listwise''');
  missVal='listwise';
end

% --- check and possibly transform input X: if it is a multiple column
% array, transform it into a single column array and set up variable group.
% This may appear counterintuitive as the computations are much easier and
% elegantly done on a matrix, in a vectorized version. However, this holds
% true only for balanced data. In real life more often than not we have to
% deal with unbalanced data, a vectorized computation of which would be
% cumbersome and, in conjunction with bootstrapping, inefficient.
% Therefore, the code is designed to deal with the most general case, and
% requires reshaping multiple column arrays into a single column array.
[nRowX nColX]=size(X);
if nColX>1
  % *** X is a multiple column-array:
  if ~isempty(group)
    % - group specified although not necessary
    warning('input variable ''group'' will be ignored with X having several columns');
    % instead, it will be generated internally
  end
  nGroup=nColX;
  group=repmat(1:nGroup,nRowX,1);
  goodIx=isfinite(X);
  nSample=sum(goodIx);
  % transformation of X and group into single column arrays and generation
  % of groupIx may need to be done after the switch construct, unless the
  % 'pairwise' case below is run
  doTransform=true;
  % missing values (nans and infs) in X
  if any(any(~goodIx))
    switch lower(missVal)
      case 'listwise'
        % omit bad row(s)
        badRowIx=any(~goodIx,2);
        X(badRowIx,:)=[];
        group(badRowIx,:)=[];
        nRowX=size(X,1);
        nSample=repmat(nRowX,1,nColX);
        % (value of doTransform unchanged)
      case 'pairwise'
        % omit bad values, at the same time transforming X and group into
        % single col-format
        X=X(goodIx);
        group=group(goodIx);
        % generate groupIx, the indexes coding for group assignment
        tmp=cumsum([0 nSample]);
        for gIx=nGroup:-1:1
          groupIx{gIx}=tmp(gIx)+1:tmp(gIx+1);
        end
        doTransform=false;
        % nSample unchanged
    end
  end   
  if doTransform
    % transform X and group into single col-format
    X=X(:);
    group=group(:);
    for gIx=nGroup:-1:1
      groupIx{gIx}=(gIx-1)*nSample+1:gIx*nSample;
    end
  end
else
  % *** X is a single column-array:
  if isempty(group)
    % - group not specified 
    error('input variable ''group'' must be specified when X is a single column-array');
  else
    [tmpNRow tmpNCol]=size(group);
    % - group not a single col-array
    if tmpNCol>1
      error('input variable group must be a single-column array');
    end
    % - numel not matching
    if tmpNRow~=nRowX
      error('input variables group and X must have the same number of elements');
    end
    % - undefined elements in group
    if any(~isfinite(group))
      error('input variable group contains NaNs or Infs');
    end
  end
  % determine how many groups there are
  [uGroup,aIx]=unique(group);
  nGroup=numel(uGroup);
  % revert sorting implicitly done by unique because the user probably
  % expects the original order to be maintained
  uGroup=group(sort(aIx));
  % if group is unsorted in the sense that samples from any group do not
  % form a contiguous block it appears possible that the data are messed
  % up, so better generate an error here and force the user to rethink (and
  % sort) the data
  isContiguousGroups=numel(find(diff(group)))==nGroup-1;
  if ~isContiguousGroups
    error([mfilename ' expects samples from each group to form contiguous blocks. Please sort your data accordingly']);
  end
  if isDep && mod(nRowX,nGroup)
    error('input variable X supposedly holds dependent data, but sample sizes are unequal');
  end
  if strcmpi(missVal,'listwise') && mod(nRowX,nGroup)
    warning('listwise elimination of missing data is undefined with unbalanced data - employing pairwise elimination (if any)');
    missVal='pairwise';
  end
  % logical index to missing values (nans and infs) in X
  badIx= ~isfinite(X);
  % if listwise elimination is requested and we have dependent or balanced data ...
  if strcmpi(missVal,'listwise') && ~mod(nRowX,nGroup)
    if ~isempty(any(badIx))
      badIx=reshape(badIx,nRowX/nGroup,nGroup);
      badIx(any(badIx,2),:)=true;
      badIx=badIx(:);
      X(badIx)=[];
      group(badIx)=[];
      nRowX=size(X,1);
    end
    % scalar expansion of variable nSample
    nSample=repmat(nRowX/nGroup,1,nGroup);
    for gIx=nGroup:-1:1
      groupIx{gIx}=(gIx-1)*nSample+1:gIx*nSample;
    end
    % variable group is not needed anymore
    group=[];
  else
    % pairwise elimination requested or independent/unbalanced data 
    if ~isempty(any(badIx))
      X(badIx)=[];
      group(badIx)=[];
    end
    for gIx=nGroup:-1:1
      groupIx{gIx}=find(group==uGroup(gIx));
      nSample(gIx)=numel(groupIx{gIx});
    end
  end
end
% At this point, variables nGroup and nSample hold the number of groups (a
% scalar) and the number of samples in each group (a row array),
% respectively, and X is a single column-array of data devoid of NaNs and
% Infs; group association is coded by variable group

% just to be 100% sure, this is an internal check that elements of groupIx
% are really integers incrementing in steps of 1
for g=1:nGroup
  tmp=unique(diff(groupIx{g}));
  if numel(groupIx{g})>1 && ~(numel(tmp)==1 && tmp==1)
    error('internal: groupIx messed up - tell the programmer');
  end
end

% --- check bootstrapping settings
doBoot=false;
if isfinite(nBoot)
  if nBoot>=minNBootstrap;
    doBoot=true;
  else
    if nBoot~=0
      % warn only if nBoot small but different from zero because zero may
      % be a deliberate input value
      warning('number of bootstrap repetitions is not adequate - not bootstrapping');
    end
    nBoot=0;
  end
end

% --- check other input arguments
if confLevel<=0 || confLevel>=1
  error('input variable ''confLevel'' must be a scalar of value >0 and <1');
end

[nCo nCoCol]=size(cWeight);
if doContrast
  if nCoCol~=nGroup
    error('number of contrast weights (per row) does not match number of groups in input data');
  end
  if any(any(~isfinite(cWeight)))
    error('contrast weights contain foul data (NaN or Inf)');
  end
  % check whether we're dealing with standard sets
  if any(abs(sum(abs(cWeight),2)-2)>eps(2))
    warning('at least one set of contrast weights is not a standard set: standardized mean differences (g_psi, psibysd) computed with this set do not represent the difference between the averages of two subsets of means')
  end
  % is there more than one contrast weight of zero? This may reflect the
  % user's intention to have one or more conditions 'set aside', in which
  % case the computations will not be valid
  if any(sum(cWeight==0,2)>1)
    warning('at least one set of contrast weights contains more than one zeroes - if you wish to exclude the corresponding groups from analysis you should eliminate them prior to computation')
  end
end
if ~any(ismember(tDenom,{'sd','msw'}))
  error('bad choice for input parameter ''tDenom'' (choose ''sd'' or ''msw'')');
end
% -------------------------------------------------------------------------
% ------- PART II: THE WORKS ----------
% -------------------------------------------------------------------------

% create a few standard fields
stats.isDep=isDep;
stats.nBoot=nBoot;
stats.confLevel=confLevel;
stats.n=nSample;
stats.cWeight=cWeight;

% preparatory computations of SS, MS etc.
r=prepcomp(X,groupIx,nSample,isDep,doBoot,nBoot,cWeight,tDenom);

% if eta2, partialeta2, omega2 or partialomega2 or any combination thereof
% is requested AND analytical CI shall be computed AND the data are
% independent, compute CI for partial eta2 here because those of most other
% ESM can be derived from them (and computing all independently of each
% other would be a waste).
if any(ismember({'eta2','partialeta2','omega2','partialomega2'},esm)) && ~doBoot && ~isDep
  % exact analytical confidence intervals for eta2 and/or partial eta2
  % (Smithson 2003, p.43 [5.6]), omnibus effects of which are identical in
  % oneway independent designs
  pEta2Ci=ncpci(r.F,'F',[r.dfGroup,r.dfErr],'confLevel',confLevel);
  pEta2Ci=pEta2Ci./(pEta2Ci+r.dfGroup+r.dfErr+1);
  for g=1:nCo
    % exact ci for contrasts: these apply only to PARTIAL eta2 because 
    % partial eta2 for a contrast = SSpsi/(SSpsi+SSerr) whereas eta2 for a
    % contrast = SSpsi/SStot = SSpsi/(SSa+SSerr), involving three terms
    pEta2Ci(g+1,:)=ncpci(r.FPsi(g),'F',[1,r.dfErr],'confLevel',confLevel);
    pEta2Ci(g+1,:)=pEta2Ci(g+1,:)./(pEta2Ci(g+1,:)+1+r.dfErr+1);
  end
  % information on kind of CI
  pEta2CiType=repmat('exact analytical',nCo+1,1);
end

% now computations of effect size parameters:
nEs=numel(esm);
% a temporary container to hold [es ci_lo ci_up] in the order in
% which they are computed (for the table of results to be displayed)
esStore=repmat(nan,[nCo+1 nEs*3]);

for tti=1:nEs
  curEs=esm{tti};
  es=repmat(nan,nCo+1,1);
  % If analytical confidence intervals for the statistic in question cannot
  % be computed and no bootstrapping is requested, ci must be set to nan.
  % Do this here for all statistics to avoid redundant assignments. ci will
  % be overwritten with meaningful values if i) there is a formula for
  % analytical confidence intervals implemented within the respective case
  % in the switch statement, or ii) confidence intervals based on
  % bootstrapping are computed right after the switch statement
  ci=repmat(nan,nCo+1,2);
  % a similar argument applies to variable ciType, which indicates on which
  % method the computation of ci is based (which may differ between the
  % omnibus effect and contrasts, hence we need more than one row
  if doBoot
    ciType=repmat('bootstrap',nCo+1,1);
  else
    ciType=repmat('none',nCo+1,1);
  end
  switch curEs
    case 'psi'
      % nan as first row because a contrast is not an omnibus effect; then
      % the contrast(s)
      es=cat(1,repmat(nan,1,nBoot+1),r.psi);
      if ~doBoot
        % note: confidence intervals of contrasts can be computed in an
        % exact fashion from central t distributions and are therefore
        % tagged 'exact analytical' here, although exactness does not, as
        % in the case of other mes below, imply computation via noncentral
        % distributions
        if isDep
          % dependent samples: 
          if strcmpi(tDenom,'sd')
            % exact CI of contrast
            ci(2:nCo+1,:)=repmat(r.psi,1,2) + r.sD/sqrt(nSample(1)) * tinv([alpha/2 1-alpha/2],r.dfSubj);
          elseif strcmpi(tDenom,'msw')
            % exact CI of contrast
            for g=1:nCo
              % SE of contrast
              fac=sqrt(r.msGroupSubj*sum((cWeight(g,:)).^2./nSample));
              ci(g+1,:)=repmat(r.psi(g),1,2) + fac * tinv([alpha/2 1-alpha/2],r.dfGroupSubj);
            end
          end
          ciType=char('none', repmat('exact analytical',nCo,1));
        else
          % independent samples:
          for g=1:nCo
            % SE of contrast
            fac=sqrt(r.msErr*sum((cWeight(g,:)).^2./nSample));
            ci(g+1,:)=repmat(r.psi(g),1,2) + fac * tinv([alpha/2 1-alpha/2],r.dfErr);
          end
          ciType=char('none', repmat('exact analytical',nCo,1));
        end
      end
    
    case 'g_psi'
      % nan as first row because this mes is not defined for an omnibus
      % effect; then the contrasts divided by standardizer proposed by
      % Kline 2004, the square root of the pooled within-conditions
      % variance, making this the equivalent of Hedges' g
      es=cat(1,repmat(nan,1,nBoot+1),r.psi./repmat(sqrt(r.msErr),nCo,1));
      if ~doBoot
        if isDep
          % dependent samples: 
          if strcmpi(tDenom,'sd')
            % exact CI of contrast...
            tmp=repmat(r.psi,1,2) + r.sD/sqrt(nSample(1)) * tinv([alpha/2 1-alpha/2],r.dfSubj);
            % ...leading to approximate CI of g_psi
            ci(2:nCo+1,:)=tmp/sqrt(r.msErr);
          elseif strcmpi(tDenom,'msw')
            % exact CI of contrast...
            for g=1:nCo
              % SE of contrast
              fac=sqrt(r.msGroupSubj*sum((cWeight(g,:)).^2./nSample));
              % ...leading to approximate CI of g_psi
              tmp=repmat(r.psi(g),1,2) + fac * tinv([alpha/2 1-alpha/2],r.dfGroupSubj);
              ci(g+1,:)=tmp/sqrt(r.msErr);
            end
          end
          ciType=char('none', repmat('approximate analytical',nCo,1));
        else
          % independent samples: exact confidence intervals 
          for g=1:nCo
            fac=sqrt(sum((cWeight(g,:)).^2./nSample));
            ci(g+1,:)=fac*ncpci(r.tPsi(g),'t',r.dfErr,'confLevel',confLevel);
          end
          ciType=char('none', repmat('exact analytical',nCo,1));
        end
      end
      
    case 'psibysd'
      % psibysd is defined only for dependent data; as the criterion of
      % dependence has been dealt with in the code above we do not heed it
      % here
      % - nan as first row, then the contrasts divided by standard
      % deviation of difference scores
      es=cat(1,repmat(nan,1,nBoot+1),r.psi./r.sD);
      if ~doBoot
        % as for g_psi above: exact CI of contrast, from which approximate
        % CI for es can be derived
        tmp=repmat(r.psi,1,2) + r.sD/sqrt(nSample(1)) * tinv([alpha/2 1-alpha/2],r.dfSubj);
        ci(2:nCo+1,:)=tmp./repmat(r.sD,1,2);
        ciType=char('none', repmat('approximate analytical',nCo,1));
      end
      
    case 'eta2'
      % eta2 does not differentiate between dependent and independent
      % samples
      es=r.ssGroup./r.ssTot;
      if doContrast
        es=cat(1,es,r.ssPsi./repmat(r.ssTot,nCo,1));
      end
      if ~doBoot
        if ~isDep
          % as explained above, ci for the omnibus effect are computable,
          % those for the contrasts are not
          ci(1,:)=pEta2Ci(1,:);
          ciType=char('exact analytical', repmat('none',nCo,1));
        end
      end
      
    case 'partialeta2'
      if isDep
        % dependent case: SSerror=r.ssGroupSubj
        es=r.ssGroup./(r.ssGroup+r.ssGroupSubj);
        if doContrast
          es=cat(1,es,r.ssPsi./(r.ssPsi+repmat(r.ssGroupSubj,nCo,1)));
        end
      else
        % independent samples: main effect same as eta squared because
        % r.ssTot=r.ssGroup+r.ssErr
        es=r.ssGroup./(r.ssTot);
        % contrast: denominator differs from dependent case
        if doContrast
          es=cat(1,es,r.ssPsi./(r.ssPsi+repmat(r.ssErr,nCo,1)));
        end
        if ~doBoot
          % all ci computed above apply to partial eta2
          ci=pEta2Ci;
          ciType=pEta2CiType;
        end
      end      
                
    case 'omega2'
      if isDep
        % formula according to Kline p. 188 & Table 6.8, setting
        % MSeffect=r.msGroup and MSerror=r.msGroupSubj
        % (note: in rare cases of partial omega squared being SMALLER than
        % omega squared, the term
        %         r.msSubj-r.msGroupSubj
        % in the denominator of omega squared must be negative)
        es=(r.dfGroup/sum(nSample)*(r.msGroup-r.msGroupSubj))./...
          (r.dfGroup/sum(nSample)*(r.msGroup-r.msGroupSubj) + ...
          1/nGroup*(r.msSubj-r.msGroupSubj) + r.msGroupSubj);
        % condensed formula from Grissom & Kim 2012, p. 198, formula 6.17
        % es=(r.dfGroup*(r.msGroup-r.msGroupSubj))./(r.ssTot+r.msSubj);
        if doContrast
          tmp=(1/sum(nSample)*(r.ssPsi-repmat(r.msGroupSubj,nCo,1)))./...
            repmat(r.dfGroup/sum(nSample)*(r.msGroup-r.msGroupSubj) + ...
            1/nGroup*(r.msSubj-r.msGroupSubj) + r.msGroupSubj,nCo,1);
          es=cat(1,es,tmp);      
        end
      else
        es=(r.ssGroup-(nGroup-1).*r.msErr)./(r.ssTot+r.msErr);
        if doContrast
          % contrasts: derived from formula 6.30, p. 186, Kline 2004, and
          % setting df(group)=1 and MS=SS
          es=cat(1,es,(r.ssPsi-repmat(r.msErr,nCo,1))./repmat(r.ssTot+r.msErr,nCo,1));
        end
        if ~doBoot
          % exact analytical confidence intervals according to Fidler &
          % Thompson 2001, p. 593, based on those for partial eta2
          tmp=(r.ssTot*(1-pEta2Ci))/r.dfErr;
          ci(1,:)=(r.ssTot*pEta2Ci(1,:)-r.dfGroup*tmp(1,:))./(r.ssTot+tmp(1,:));
          % as in the case of eta2 and contrasts, we cannot compute ci for
          % contrasts for omega2 for the same reasons stated above
          ciType=char('exact analytical', repmat('none',nCo,1));
        end
      end
      
    case 'partialomega2'
      if isDep
        % formula according to Kline p. 188 & Table 6.8, setting
        % MSeffect=r.msGroup and MSerror=r.msGroupSubj (note that the only
        % difference to the corresponding case of omega2 is one term less
        % in the denominator, namely the subjects' variance component. 
        es=(r.dfGroup/sum(nSample)*(r.msGroup-r.msGroupSubj))./...
          (r.dfGroup/sum(nSample)*(r.msGroup-r.msGroupSubj) + ...
          r.msGroupSubj);
        if doContrast
          tmp=(1/sum(nSample)*(r.ssPsi-repmat(r.msGroupSubj,nCo,1)))./...
            repmat(r.dfGroup/sum(nSample)*(r.msGroup-r.msGroupSubj) + ...
            r.msGroupSubj,nCo,1);
          es=cat(1,es,tmp);      
        end
      else
        % independent samples: same as omega squared
        es=(r.ssGroup-(nGroup-1).*r.msErr)./(r.ssTot+r.msErr);
        if doContrast
          % contrasts: assembled from Kline p. 186 & Table 6.7
          tmp=(1/sum(nSample)*(r.ssPsi-repmat(r.msErr,nCo,1)))./...
            (1/sum(nSample)*(r.ssPsi-repmat(r.msErr,nCo,1)) + repmat(r.msErr,nCo,1));
          es=cat(1,es,tmp);
        end
        if ~doBoot
          % exact analytical confidence intervals according to Fidler &
          % Thompson 2001, p. 593, based on those for partial eta2
          tmp=(r.ssTot*(1-pEta2Ci))/r.dfErr;
          ci=(r.ssTot*pEta2Ci(1,:)-r.dfGroup*tmp(1,:))./(r.ssTot+tmp(1,:));
          if doContrast
            % same formula applied to CI of contrasts, taking into account df(group)=1 
            ci(2:nCo+1,:)=(r.ssTot*pEta2Ci(2:nCo+1,:)-1*tmp(2:nCo+1,:))./(r.ssTot+tmp(2:nCo+1,:));
          end          
          ciType=pEta2CiType;
        end
      end
      
    otherwise
      error(['internal error: unrecognized test not caught by input checks: ' curEs]);
  end
  
  % *********************************************************************
  % If data were NOT bootstrapped, all computations are done at this
  % point and the results can be placed into appropriate fields of output
  % variable stats. If they were bootstrapped, confidence intervals must
  % be computed and the ES statistics extracted from the first element of
  % variable es.
  % *********************************************************************
  if doBoot
    % determine confidence intervals from array of effect size measures
    % generated from bootstrapped data
    ci=prctile(es(:,2:end)',[alpha/2  1-alpha/2]'*100)';
    % retain first column; this is the es computed from real data
    es=es(:,1);
  end
  % place es and ci in esStore
  esStore(:,(tti-1)*3+1:tti*3)=[es ci];
   
  % finally, use dynamic fields to store currently computed measures in
  % output variable stats
  stats.(curEs)=es;
  stats.([curEs 'Ci'])=ci;
  stats.([curEs 'CiType'])=ciType;
end

% assemble table for display and optional output:
% the concatenated es and ci as cell array
esStore=num2cell(esStore);
% 'fillers' of empty cells
filla=cell(1,nEs*3);
% insert ci titles
ciTi=cell(2,nEs);
[ciTi{1,:}]=deal('ci_lo');
[ciTi{2,:}]=deal('ci_up');
% chars to be added to headers of the F and p columns
if isDep && nCo>0 && strcmpi(tDenom,'sd')
  ac=' (**)';
else
  ac='';
end

esm=cat(1,esm,ciTi);
esm=esm(:)';
% column order:
% source | SS | df | MS | F | p | effect sizes (in the order requested)
summaryTable={...
  'source',                'SS',            'df',          'MS',            ['F',ac]   ['p',ac] esm{1:end};...
  'between-groups',        r.ssGroup(1),    r.dfGroup,     r.msGroup(1),    r.F(1),    r.p(1),  esStore{1,:};...
  };
for g=1:nCo
  summaryTable=cat(1,summaryTable,...
    {['- psi' int2str(g)], r.ssPsi(g),      1,             r.msPsi(g),      r.FPsi(g), r.pPsi(g), esStore{1+g,:}});
end
summaryTable=cat(1,summaryTable,...
  {'within-groups',        r.ssErr(1),      r.dfErr,       r.msErr(1),      [],      [],        filla{1:end}});
if isDep
  summaryTable=cat(1,summaryTable,...
    {'- subj',             r.ssSubj(1),     r.dfSubj,      r.msSubj(1),     [],      [],        filla{1:end};...
    '- group x subj',      r.ssGroupSubj(1),r.dfGroupSubj, r.msGroupSubj(1),[],      [],        filla{1:end}});
end
summaryTable=cat(1,summaryTable,...
  {'total',                r.ssTot(1),       r.dfTot,       r.msTot(1),     [],      [],        filla{1:end}});

summaryTable  
if isDep && nCo>0 && strcmpi(tDenom,'sd')
  disp(char({...
    '(**) NOTE: F and p values and the confidence intervals of psi and';...
    '     g_psi are derived from the standard deviation of the contrasts'' ';...
    '     difference scores. For an alternative computation set optional ';...
    '     input parameter ''tDenom'' to ''msw''.';...
    '     See the documentation for detailed information.'}));
end

if nargout>1
  varargout{1}=summaryTable;
end

% ========================= LOCAL FUNCTIONS ===============================
% ========================= LOCAL FUNCTIONS ===============================

function r=prepcomp(x,groupIx,nSample,isDep,doBoot,nBoot,cWeight,tDenom)
% performs preparatory computations for oneway analyses on single column
% array x with group assignment coded by input var group

nGroup=numel(groupIx);
nCo=size(cWeight,1);
% ------------------ START BY BOOTSTRAPPING (IF REQUESTED) ----------------
if doBoot
  % *********************************************************************
  % Here, individual, groupwise resampled instances of variable X are 
  % generated, expanding X in the second dimension: the first column
  % contains the original data, the others will be filled up with sampled
  % (with replacement) data.
  % *********************************************************************
  % variable bootSamX generated below shall contain indices to be used for
  % randomly sampling (with replacement) the original data (*** requires
  % groupIx to be integers incrementing in steps of 1! ***)
  if isDep
    % dependent data: generate random sampling index for first group and
    % replicate for all other groups (which are identical in size)
    bootSamX=repmat(ceil(nSample(1)*rand([nSample(1) nBoot])),[nGroup,1]);
    % add offset due to arrangement of data in a single-column array
    for g=1:nGroup
      bootSamX(groupIx{g},:)=bootSamX(groupIx{g},:)+groupIx{g}(1)-1;
    end
  else
    % independent data: resample data within each group independent of data
    % points picked in other groups
    bootSamX=repmat(nan,sum(nSample),nBoot);
    for g=1:nGroup
      bootSamX(groupIx{g},1:nBoot)=groupIx{g}(1)-1+ceil(nSample(g)*rand([nSample(g) nBoot]));
    end
  end
  % from second column onwards fill with sampled data
  x(:,2:nBoot+1)=x(bootSamX);
  % delete resampling index 
  clear bootSam* 
end

% between-groups df
r.dfGroup=nGroup-1;
% within-groups df
r.dfErr=sum(nSample-1);
% within-subjects df
r.dfSubj=nSample(1)-1;
% group x subjects df 
r.dfGroupSubj=r.dfGroup*r.dfSubj;

% grand mean
r.meanGrand=mean(x);

r.ssGroup=0;
r.ssErr=0;
for g=nGroup:-1:1
  % means of groups
  r.meanGroup(g,:)=sum(x(groupIx{g},:))/nSample(g);
  % SS_a, or SS_effect, between-groups SS
  r.ssGroup=r.ssGroup + nSample(g)*(r.meanGroup(g,:)-r.meanGrand).^2;
  % SS_error, within-groups SS
  r.ssErr=r.ssErr + sum((x(groupIx{g},:)-repmat(r.meanGroup(g,:),nSample(g),1)).^2);
end
% between-groups MS
r.msGroup=r.ssGroup/r.dfGroup;
% within-groups MS
r.msErr=r.ssErr/r.dfErr;

% SS_t, the sum of both
r.ssTot=r.ssGroup+r.ssErr;
% corresponding df
r.dfTot=r.dfGroup+r.dfErr;
r.msTot=r.ssTot/r.dfTot;

% parameters for dependent data
if isDep
  nSubj=nSample(1);
  r.ssSubj=zeros(1,nBoot+1);
  for g=nSubj:-1:1
    % means of subjects
    r.meanSubj(g,:)=mean(x(g:nSubj:end,:));
    % SS_subj, between-subjects SS
    r.ssSubj=r.ssSubj+nGroup*(r.meanSubj(g,:)-r.meanGrand).^2;
  end
  % SS group x subject
  r.ssGroupSubj=r.ssErr-r.ssSubj;
  % MS
  r.msSubj=r.ssSubj/r.dfSubj;
  r.msGroupSubj=r.ssGroupSubj/r.dfGroupSubj;
end
    
% compute contrasts and their SS
if nCo>0
  for g=nCo:-1:1
    % contrasts
    r.psi(g,:)=sum(repmat(cWeight(g,:)',1,nBoot+1).*r.meanGroup);
    % factor needed for computation of SS_psi and t for contrasts
    r.cFac(g,:)=sum(cWeight(g,:).^2./nSample);
    % SS_psi: see Kline 2004 (formula 6.8, p.168)
    r.ssPsi(g,:)=r.psi(g,:).^2./r.cFac(g,:);
    if isDep
      % std of the difference score
      tmp=repmat(cWeight(g,:),nSample(1),1);
      tmp=repmat(tmp(:),1,nBoot+1);
      tmp=reshape(x.*tmp,[nSample(1) nGroup nBoot+1]);
      r.sD(g,:)=permute(std(sum(tmp,2)),[1 3 2]);
    end
  end
  % as df is always 1 MS=SS, but create field nonetheless
  r.msPsi=r.ssPsi;
end

% finally, t, F, p
if isDep
  r.F=r.msGroup./r.msGroupSubj;
  r.p=1-fcdf(r.F,r.dfGroup,r.dfGroupSubj);
  if nCo>0
    if strcmpi(tDenom,'sd')
      % compute p values of contrasts from t statistic (Kline 2004, p.
      % 168): note that this t statistic is based on the std of the
      % contrast's difference score (in the denominator, see line below).
      % This means that i) t takes into account only the variability
      % between the groups/subjects compared in the contrast; specifically,
      % groups with a contrast weight of zero do not play a role, and ii) a
      % potential lack of sphericity among groups in the full data set is
      % not a problem.
      r.tPsi=r.psi./(r.sD/sqrt(nSample(1)));
      r.pPsi=2*(1-tcdf(abs(r.tPsi),r.dfSubj));
      r.FPsi=r.tPsi.^2;
    elseif strcmpi(tDenom,'msw')
      % F and p values are computed the 'classical' way, which takes
      % variability from all groups into account
      r.FPsi=r.msPsi./repmat(r.msGroupSubj,nCo,1);
      r.pPsi=1-fcdf(r.FPsi,1,r.dfGroupSubj);
      r.tPsi=sign(r.psi).*sqrt(r.FPsi);
    end
  end
else
  r.F=r.msGroup./r.msErr;
  r.p=1-fcdf(r.F,r.dfGroup,r.dfErr);
  if nCo>0
    % compute p values of contrasts from independent data on the basis of t
    % statistic (Kline 2004, formula 6.7, p. 168)
    r.tPsi=r.psi./sqrt(repmat(r.msErr,nCo,1).*repmat(r.cFac,1,nBoot+1));
    r.pPsi=2*(1-tcdf(abs(r.tPsi),r.dfErr));    
    % F is by definition the square of t
    r.FPsi=r.tPsi.^2;
    % ? note: for the example in Kline 2004, table 6.4, the F values are
    % all correct but the contrast-associated p values do not match
  end
end


function [stats,varargout]=mes2way(X,group,esm,varargin)
% ** function [stats,varargout]=mes2way(X,group,esm,varargin)
% computes measures of effect size for samples in a factorial two-way
% design and places the results in output structure stats. A summary table
% including the effect size measures and two-way ANOVA results will be
% displayed in the command window. All input parameters except X, group and
% esm are optional and must be specified as parameter/value pairs in any
% order, e.g. as in
%      mes2way(X,g,'eta2','confLevel',.9);
%
% For information on assumptions of the underlying model and related
% information please see the notes below TABLE 1 further down as well as
% the documentation accompanying this code.
%
%                           INPUT ARGUMENTS
%                           ---------------
% stats=mes2way(X,group,esm)  computes the effect size measure(s) 
%   specified in input variable esm (see TABLE 1) from input array X. X
%   must be a single-column vector. NaNs or Infs will be eliminated in
%   'pairwise' fashion, that is, without concurrent elimination of other
%   values (which would be the case in 'listwise' fashion). group must be a
%   two-column array with the same number of rows as X. The numbers in the
%   columns of group correspond to the levels of the two factors. 
%   PLEASE NOTE: the actual values in group may be arbitrary, but rank
%   order, not the order in which they appear in group, determines the
%   assignment to factor levels. For example, if group is
%     3.1  1
%     3.1  2
%     1.5  1
%     1.5  2
%   the first two data points of X are assigned to the SECOND level of the
%   first factor (because 3.1 > 1.5). This is particularly relevant for all
%   calculations based on contrasts (see below). Input variable esm must be
%   a character array (e.g. 'eta2') or a cell array (e.g.
%   {'eta2','partialeta2'}); see TABLE 1 below.
% stats=mes2way(...,'fName',{'genotype','treat'})  assigns names to the two 
%   factors which will be displayed in the summary table of results; by 
%   default they will be given generic names 'fac1' and 'fac2'
% stats=mes2way(...,'isDep',[1 0])  assumes that the samples in X are 
%   dependent (repeated measures) with respect to the first factor. Other
%   legal input values are [0 1] (dependence with respect to the second
%   factor), [1 1] (completely within-subjects design) and [0 0]
%   (completely between-subjects design, the default). If there is
%   within-subjects dependence along one or both factors, data must be
%   balanced. Furthermore, it is assumed that input data X are sorted
%   according to subjects; that is, within each group data for subjects are
%   listed in the same order.
% stats=mes2way(...,'nBoot',n)  computes confidence intervals of the 
%   statistic in question via bootstrapping n times. n should be on the
%   order of thousands, otherwise bootstrapping will lead to inaccurate
%   results. By default or if the value of nBoot is set to zero or
%   nonsensical values (<0, infinity, nan) confidence intervals will be 
%   computed analytically (where possible).
% stats=mes2way(...,'confLevel',0.9)  computes 90 % confidence  
%   intervals of the statistic in question (95 % ci are the default; any 
%   value may be specified)
% stats=mes2way(...,'cWeight',c)  allows specification of contrast weights
%   for the computation of standardized contrasts (see TABLE 1). The
%   required size of input array c depends on the kind of comparison (see
%   e.g. Kline 2004 for definitions and the documentation accompanying this
%   code for concrete examples):
%   MAIN COMPARISON CONTRAST: c must be a single-column array (comparison
%     of levels of the first factor) or a single-row array (ditto for
%     second factor)
%   SIMPLE COMPARISON CONTRAST: c must match the analysis design, i.e. in a
%     2x3 analysis c must have two rows and three columns, but all elements
%     except the row or the column of interest must be NaN
%  INTERACTION CONTRAST: c must match the analysis design and it should
%     also be doubly centered, that is, row sums and column sums must all
%     be zero
% stats=mes2way(...,'doDataPlot',true)  will produce a plot of the data 
%     according to the layout of the analyis: levels of the factors are in
%     subplot rows/columns, repeated measures data points are plotted with
%     identical colors, backgound colors of the plots reflect contrast
%     weights
%
% -------------------------------------------------------------------------
% TABLE 1: ANALYSES TO BE SPECIFIED IN INPUT VARIABLE esm
% -------------------------------------------------------------------------
% esm             QUANTITIY COMPUTED ((*) require input of contrast weights)
% 'psi'           unstandardized contrast (*)
% 'g_psi'         standardized contrast (*)  
% 'eta2'          eta squared
% 'partialeta2'   partial eta squared 
% 'omega2',       omega squared
% 'partialomega2' partial omega squared
%
% Notes:
% i. Parameters in TABLE 1 marked by (*) need contrast weights to be
%  computed. The other parameters will be computed for main and interaction
%  effects and contrasts (if specified).
% ii. 'g_psi' is one possible twoway equivalent of Hedges' g, namely,
%  contrast divided by the square root of the pooled within-conditions
%  variance. Its value is identical for dependent and independent data.
% iii. Fixed factors and interactions between the main factors are assumed.
% iv. Data may be unbalanced along the between-subjects factor(s)
% iv. If data is unbalanced it is assumed that the imbalance is due to
%  random loss of data, not that unequal cell size is part of the effect.
%  Accordingly, type III errors are computed and an unweighted means
%  analysis (using harmonic means) is performed
%
%                         OUTPUT ARGUMENTS
%                         ----------------
% The results of the computations are placed into fields of output
% structure stats with names corresponding to the analyses requested. For
% example, stats=effectsz_oneway(X,{'eta2','partialeta2'}) will result in
% fields
%   .eta2 
%   .partialeta2
% and so on. 
% PLEASE NOTE: if contrast weights were specified, all output arguments
% will be single- or two-column arrays holding results in the following
% (row) order:
% 1. main effect of factor 1
% 2. main effect of factor 2
% 3. interaction effect
% 4. contrast-related effect
% If a parameter is not defined for main/interaction effects or for
% contrasts, the corresponding entries will be NaN.
% Where applicable, confidence intervals will be placed in fields
%   .eta2Ci 
%   .partialeta2Ci
% and the computations underlying confidence intervals (either of 'none',
% 'bootstrap', 'approximate analytical', or 'exact analytical') will be
% placed in fields
%   .eta2CiType 
%   .partialeta2CiType
% Additional fields:
%   .n (sample sizes)
%   .confLevel (confidence level)
%   .contrast (containing fields .weight, the contrast weights and .type,
%   the type of comparison)
% stats will also have fields .isDep and .nBoot with the same values as the
% corresponding input arguments, thus providing potentially crucial post-
% analysis information. The second, optional output argument is the summary
% table of results.

% -------------------------------------------------------------------------
% Measures of Effect Size Toolbox Version 1.4, January 2015
% Code by Harald Hentschke (University of T?bingen) and 
% Maik St?ttgen (University of Bochum)
% For additional information see Hentschke and St?ttgen, 
% Eur J Neurosci 34:1887-1894, 2011
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%                       STRUCTURE OF CODE
% The code is composed of two major parts: In PART I input arguments are
% checked and pre-processed, 'constants' set up, and all kinds of other
% preparatory and householding jobs accomplished. PART II contains the
% computations proper. Right at the beginning, summed squares, F and p
% values are computed in local function prepcomp, which also accomplishes
% assembly of bootstrapped data if requested. Effect sizes and confidence
% intervals are then computed and the results placed in structure array
% 'stats' as described above. The results are additionally collected in
% cell array 'table' which is displayed in the command window (and may also
% be retrieved as an output argument)
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% ------- PART I: CHECKS OF INPUT ARGS & OTHER PREPARATORY WORKS ----------
% -------------------------------------------------------------------------
% ----- default values & varargin -----
% standard values of optional input arguments
fName=[];
isDep=[0 0];
nBoot=0;
cWeight=[];
confLevel=.95;
doDataPlot=false;

% check variable number of input args:
% - even number (parameter AND value specified)?
% - convert to proper upper/lowercase spelling as above
% - check whether valid parameter was given
% - overwrite defaults by values, if specified
nvarg=numel(varargin);
v=who;
if nvarg
  if nvarg/2==round(nvarg/2)
    for g=1:2:nvarg
      % check which optional input parameter is given, ignoring case, and
      % impose spelling used internally
      ix=find(strcmpi(varargin{g},v));
      if ~isempty(ix)
        varargin{g}=v{ix};
      else
        error(['invalid optional input parameter (' varargin{g} ') was specified']);
      end
    end
    % finally, the assignment of values to parameters
    pvpmod(varargin);
  else
    error('optional input parameters must be specified as parameter/value pairs, e.g. as in ''par1'',1')
  end
end

% ----- a few 'constants':
alpha=1-confLevel;
% number of factors
nFactor=2;
% minimal number of bootstrapping runs below which a warning will be issued
% (and bootstrapping will be skipped)
minNBootstrap=1000;
% this is an internal switch dictating the way summed squares are computed
% (relevant only to unbalanced designs); currently valid options are
% 'unweighted' and 'Type III'. If set to 'unweighted', an 'unweighted means
% analysis' is performed, that is, grand mean and marginal means are means
% of the corresponding cells NOT weighted by the number of samples in them,
% and 'overall' cell size is the harmonic mean of all cell sizes, which is
% plugged into the traditional formulae for balanced designs. It seems not
% to be used a lot in the literature, so by default the 'Method 1/Type III'
% way according to Overall & Spiegel 1969 is pursued, which is also the
% default in Matlab.
ssType='Type III';

if nargin<3
  error('input arguments X, group and esm are mandatory');
end
% Variable list_analyses below is a list of all possible analyses:
% - column 1 holds the shortcut as strings
% - column 2 holds 1 for quantities REQUIRING a contrast
% - column 3 is a short description of the analyses
% This is the 'master list' of all analyses against which input argument
% esm will be checked (see) below
list_analysis={...
'psi'            1, 'unstandardized contrast';...
'g_psi',         1, 'standardized contrast';...
'eta2',          0, 'eta squared';...
'omega2',        0, 'omega squared';...
'partialeta2',   0, 'partial eta squared';...
'partialomega2', 0, 'partial omega squared';...
};

% in case esm is a char array convert it to a cell array because the code
% below relies on it
if ischar(esm)
  esm={esm};
end

% the tag column, extracted
listTag=cat(1,list_analysis{:,2});
% indices to currently requested analyses
listIx=ismember(list_analysis(:,1),esm);
% unique tags of these analyses
uTag=unique(listTag(listIx));
% catch a few silly errors:
% - typos or non-existent analysis type
if any(~ismember(esm,list_analysis(:,1)))
  error('An illegal type of analysis was specified');
end

% check whether contrast weight(s) were specified, whether they are
% required, and set according flags (exact checks of contrast weights will
% be done further below). Combine all this information in a struct
contrast.weight=cWeight;
if isempty(contrast.weight)
  contrast.type='none';
  contrast.do=false;
  if any(uTag)
    error('part of the requested analyses require contrast weights as input args')
  end
else
  contrast.do=true;
  % kind of contrast unknown at this point (see further below)
  contrast.type='unspec';
end

% --- check input X and group
[nRowX nColX]=size(X);
[nRowG nColG]=size(group);
if nColX~=1 
  error('X must be a single column-array');
end
% - number of rows not matching
if nRowG~=nRowX
  error('input variables group and X must have the same number of rows');
end
if nColG~=nFactor
  error('input variable group must have two columns representing the levels of the two factors');
end
% - undefined elements in group
if any(~isfinite(group))
  error('input variable group contains NaNs or Infs');
end
% - nans and infs in X
badIx=find(~isfinite(X));
if ~isempty(badIx)
  warning('eliminating NaNs and Infs from input variable X');
  X(badIx)=[];
  group(badIx,:)=[];
  nRowX=nRowX-numel(badIx);
  nRowG=nRowG-numel(badIx);
end

% determine how many factors and levels there are, on the occasion
% replacing the values in variable group (i.e. the levels of the factors)
% by integers 1,2,... because function dummyvar used below requires this
for g=1:nFactor
  [factor(g).level,nada,intG]=unique(group(:,g));
  group(:,g)=intG;
  factor(g).nLevel=numel(factor(g).level);
  % factor(g).levelName=cellstr([repmat('level',factor(g).nLevel,1) int2str((1:factor(g).nLevel)')]);
end
% assign labels to factors
if iscell(fName) && numel(fName)==2
  [factor.name]=deal(fName{:});
else
  % if fName is anything but the default empty matrix, warn
  if ~isempty(fName)
    warning('check names for factors - setting generic names')
  end
  [factor.name]=deal('fac1','fac2');
end
% number of individual groups (cells) (nonredundant combinations of levels
% of all factors)
nGroup=prod([factor.nLevel]);
% ****************************** NOTE: ************************************
% the structure/layout of all intermediate and final results hinges on the
% order of groups established below by unique, namely, the rank order of
% the (numeric) group labels!
% *************************************************************************
[uGroup,aIx,bIx]=unique(group,'rows');
% any empty cell?
if nGroup~=size(uGroup,1)
  error('there is at least one empty group (cell), which is not allowed in factorial designs');
end
% if group is unsorted in the sense that samples from any group do not
% form a contiguous block it appears possible that the data are messed
% up, so better generate an error here and force the user to rethink (and
% sort) the data
isContiguousGroups=sum(any(diff(group),2))==nGroup-1;
if ~isContiguousGroups
  error([mfilename ' expects samples from each group to form contiguous blocks. Please sort your data accordingly']);
end
% *** for each group, determine indexes & count samples. Arrange the
% results in 2D cell array groupIx according to the analysis design: first
% factor = rows, second factor = columns ***
% (use output from unique, linear indexing and transposition to embed
% properly)
groupIx=cell([factor.nLevel])';
nSample=zeros([factor.nLevel])';
for gIx=1:nGroup
  tmpIx=find(bIx==gIx);
  groupIx{gIx}=tmpIx;
  nSample(gIx)=numel(tmpIx);
end
% don't forget to transpose
groupIx=groupIx';
nSample=nSample';

% in contrast to the 2D layout of groupIx and nSample most variables (i.e.
% summed squares, etc.) will have a 1D layout, generated by accessing
% elements of groupIx via linear indexing  (the 2nd dim is reserved for
% bootstrapped data). So, in order to compute e.g. marginal means we need
% linear indexes: 1D equivalents of 2D indexes into entire rows and columns
% in groupIx: with the code below, e.g. the first row of s2i2 contains
% indexes into e.g. r.meanGroup (computed in local function prepcomp),
% corresponding to the first row of groupIx, that is, all data of the first
% level of the first factor, and so on. Likewise, the first column will
% contain indexes to the first level of the second factor. (naming: s2i2 as
% a shortcut for 'output of sub2ind, placed in 2D array')
s2i2=reshape(1:nGroup,[factor.nLevel]);

% just to be 100% sure, this is an internal check that elements of groupIx
% are really integers incrementing in steps of 1
for g=1:nGroup
  tmp=unique(diff(groupIx{g}));
  if numel(groupIx{g})>1 && ~(numel(tmp)==1 && tmp==1)
    error('internal: groupIx messed up - tell the programmer');
  end
end

% --- check bootstrapping settings
doBoot=false;
if isfinite(nBoot)
  if nBoot>=minNBootstrap;
    doBoot=true;
  else
    if nBoot~=0
      % warn only if nBoot small but different from zero because zero may
      % be a deliberate input value
      warning('number of bootstrap repetitions is not adequate - not bootstrapping');
    end
    nBoot=0;
  end
end

% --- check other input arguments
if confLevel<=0 || confLevel>=1
  error('input variable ''confLevel'' must be a scalar of value >0 and <1');
end

% deal with contrast weights
if contrast.do
  [cwN1 cwN2]=size(contrast.weight);
  % should there be Infs, convert to NaNs
  contrast.weight(isinf(contrast.weight))=nan;
  if cwN1==factor(1).nLevel && cwN2==1
    % single-column array: main comparison of first level
    contrast.type='main1';
  elseif cwN1==1 && cwN2==factor(2).nLevel
    % single-row array: main comparison of second level
    contrast.type='main2';
  elseif isequal([cwN1 cwN2],[factor.nLevel])
    % layout of contrast.weight mirrors analysis design, but is not a 1 by
    % cwN2 or cwN1 by 1 design
    cwIsFin=isfinite(contrast.weight);
    if all(all(cwIsFin))
      % all values are finite: interaction contrast
      contrast.type='interaction';
    elseif numel(find(any(cwIsFin,1)))==1 && numel(find(any(cwIsFin,2)))==factor(1).nLevel
      % only one column contains finite values: simple comparison of first
      % factor at specific level of second factor
      contrast.type='simple1';
    elseif numel(find(any(cwIsFin,2)))==1 && numel(find(any(cwIsFin,1)))==factor(2).nLevel
      % the inverse
      contrast.type='simple2';
    else
      error('check matrix of contrast weights - it contains NaNs or Infs in improper places');
    end
  else
    error('check layout of contrast weights');
  end
end
% final checks: values
if contrast.do
  if any(strcmp(contrast.type,{'simple1','simple2','main1','main2'}))
    if nansum(nansum(abs(contrast.weight)))-2>eps(2)
      warning('contrast weights for simple or main comparison are not a standard set: the standardized mean difference g_psi computed with this set does not represent the difference between the averages of two subsets of means')
    end
  elseif strcmp(contrast.type,{'interaction'})
    if ~sum(sum(abs(contrast.weight)))
      error('contrast weights are all zero');
    elseif any(abs(sum(contrast.weight,1))>eps(1)) || any(abs(sum(contrast.weight,2))>eps(1))
      warning('array of contrast weights must be doubly centered, otherwise the resulting contrast is NOT independent of the main effects');
    elseif sum(sum(abs(contrast.weight)))-4>eps(4)
      warning('the sum of the absolute values of contrast weights is unequal 4: the standardized mean difference g_psi computed with this set does not represent the difference between a pair of simple comparisons');
    end
  end
end
% check design (isDep) and sample sizes
isDep=logical(isDep(:)');
if ~isequal(size(isDep), [1 2])
  error('check input variable ''isDep'' - it must be a two-element array'),
end
% in repeated measures designs check whether sample sizes match
if any(isDep)
  if isDep(1) && size(unique(nSample,'rows'),1)~=1
    error(['cell sizes across within-subjects factor (' factor(1).name ') must all be equal']);
  end
  % note the transpose
  if isDep(2) && size(unique(nSample','rows'),1)~=1
    error(['cell sizes across within-subjects factor (' factor(2).name ') must all be equal']);
  end
end

% if plot of data is requested do it here before the computations
if doDataPlot
  mesdplot(X,groupIx,nSample,factor,isDep,fName,contrast)
end

% -------------------------------------------------------------------------
% ------- PART II: THE WORKS ----------
% -------------------------------------------------------------------------

% create a few standard fields
stats.isDep=isDep;
stats.nBoot=nBoot;
stats.confLevel=confLevel;
stats.n=nSample;
stats.contrast.type=contrast.type;
stats.contrast.weight=contrast.weight;

% preparatory computations of SS, MS etc.
r=prepcomp(X,group,groupIx,factor,s2i2,nSample,isDep,doBoot,nBoot,contrast,ssType,alpha);

% if partialeta2 or partialomega2 or both are requested AND analytical CI
% shall be computed AND the data are independent, compute CI for partial
% eta2 here because those of partial omega2 can be derived from them (and
% computing them independently of each other would be a waste).
% Note on the side: in twoway designs, exact confidence intervals can be
% computed (via noncentral F) only for the partialled versions of eta2 and
% omega2. This is because the denominator in the formula for e.g. eta2 is
% SStot=SSerr+SS1+SS2+SS12, meaning that eta2 depends on several
% noncentralities. Thus, even if it is mentioned only in passing (or not at
% all) in papers, if analytical ci for what the authors term eta squared or
% omega squared are specified, in all likelihood these are the ci for the
% partialled MES
if any(ismember({'partialeta2','partialomega2'},esm)) && ~doBoot && ~any(isDep)
  % exact analytical confidence intervals for partial eta2 (Smithson 2003, p.43 [5.6])
  % - main effects
  for fi=1:2
    tmp=ncpci(r.factor(fi).F,'F',[r.factor(fi).df,r.dfErr],'confLevel',confLevel);
    pEta2Ci(fi,1:2)=tmp./(tmp+r.factor(fi).df+r.dfErr+1);
  end
  % - interaction
  tmp=ncpci(r.factor12.F,'F',[r.factor12.df,r.dfErr],'confLevel',confLevel);
  pEta2Ci(3,:)=tmp./(tmp+r.factor12.df+r.dfErr+1);
  % - contrast
  if contrast.do
    tmp=ncpci(r.FPsi,'F',[1,r.dfErr],'confLevel',confLevel);
    pEta2Ci(4,:)=tmp./(tmp+1+r.dfErr+1);
  end
  % information on kind of CI
  pEta2CiType=repmat('exact analytical',3+contrast.do,1);
end

% now computations of effect size parameters 
nEs=numel(esm);
% a temporary container to hold es and ci for the table of results to be
% displayed; column order is [es ci_lo ci_up] en block for each es in the
% order in which they are computed; row order is factor1, factor2,
% factor1*factor2, contrast
esStore=repmat(nan,[3+contrast.do nEs*3]);

for tti=1:nEs
  curEs=esm{tti};
  es=repmat(nan,3+contrast.do,1);
  % If analytical confidence intervals for the statistic in question cannot
  % be computed and no bootstrapping is requested, ci must be set to nan.
  % Do this here for all statistics to avoid redundant assignments. ci will
  % be overwritten with meaningful values if i) there is a formula for
  % analytical confidence intervals implemented within the respective case
  % in the switch statement, or ii) confidence intervals based on
  % bootstrapping are computed right after the switch statement
  ci=repmat(nan,3+contrast.do,2);
  % a similar argument applies to variable ciType, which indicates on which
  % method the computation of ci is based (which may differ between main
  % and interaction effects and contrasts, hence we need more than one row
  if doBoot
    ciType=repmat('bootstrap',3+contrast.do,1);
  else
    ciType=repmat('none',3+contrast.do,1);
  end
  switch curEs
    case 'psi'
      % rows 1-3 are nans because they are reserved for main and
      % interaction effects; then the contrast
      es=cat(1,repmat(nan,3,nBoot+1),r.psi);
      if ~doBoot
        ci(4,:)=r.ciPsi;
        ciType=char(repmat('none',3,1),'exact analytical');
      end
      
    case 'g_psi'
      % rows 1-3 are nans because this mes is only defined for contrasts;
      % then the contrast divided by standardizer, the square root of the
      % pooled within-conditions variance, making this one possible
      % equivalent of Hedges' g
      es=cat(1,repmat(nan,3,nBoot+1),r.psi./sqrt(r.msErr));
      if ~doBoot
        % independent samples: exact confidence intervals 
        if ~any(isDep)
          % relation of ncp to g_psi: Kline formula 6.14, p. 177
          ci(4,:)=sqrt(r.denomSsPsi)*ncpci(r.tPsi,'t',r.df4nc_Psi,'confLevel',confLevel);
          ciType=char(repmat('none',3,1),'exact analytical');
        end
      end
      
    case 'eta2'
      % main effects
      es(1,1:nBoot+1)=r.factor(1).ss./r.ssTot;
      es(2,:)=r.factor(2).ss./r.ssTot;
      % interaction
      es(3,:)=r.factor12.ss./r.ssTot;
      % contrast
      if contrast.do
        es(4,:)=r.ssPsi./r.ssTot;
      end
      % in contrast to oneway analyses exact CI are not computable for
      % eta2
      
    case 'partialeta2'
      % main effects
      % (note the second term in the denominator, which is the
      % design-dependent error of the effect)
      es(1,1:nBoot+1)=r.factor(1).ss./(r.factor(1).ss+r.factor(1).ssEffErr);
      es(2,:)=r.factor(2).ss./(r.factor(2).ss+r.factor(2).ssEffErr);
      % interaction
      es(3,:)=r.factor12.ss./(r.factor12.ss+r.factor12.ssEffErr);
      % contrast
      if contrast.do
        % (note second term in denominator)
        es(4,:)=r.ssPsi./(r.ssPsi+r.ssEffErr_Psi);
      end
      % exact analytical ci apply only to completely between subjects
      % design
      if ~doBoot && ~any(isDep)
        ci=pEta2Ci;
        ciType=pEta2CiType;
      end
      
    case 'omega2'
      if ~any(isDep)
        % main effects
        es(1,1:nBoot+1)=(r.factor(1).ss-r.factor(1).df.*r.msErr)./(r.ssTot+r.msErr);
        es(2,:)=(r.factor(2).ss-r.factor(2).df.*r.msErr)./(r.ssTot+r.msErr);
        % interaction
        es(3,:)=(r.factor12.ss-r.factor12.df.*r.msErr)./(r.ssTot+r.msErr);
        if contrast.do
          % contrasts: df(effect)=1 and MS=SS
          es(4,:)=(r.ssPsi-1*r.msErr)./(r.ssTot+r.msErr);
        end
        % ci: see corresponding comments for eta2
      end
      
    case 'partialomega2'
      if ~any(isDep)
        % (simplified formula 6.31, p. 186, Kline 2004)
        tmp=prod([factor.nLevel])*r.hCellSz;
        % main effects
        es(1,1:nBoot+1)=r.factor(1).df*(r.factor(1).F-1)./(r.factor(1).df*(r.factor(1).F-1)+tmp);
        es(2,:)=r.factor(2).df*(r.factor(2).F-1)./(r.factor(2).df*(r.factor(2).F-1)+tmp);
        % interaction
        es(3,:)=r.factor12.df*(r.factor12.F-1)./(r.factor12.df*(r.factor12.F-1)+tmp);
        if contrast.do
          % contrast: df(effect)=1
          es(4,:)=1*(r.FPsi-1)./(1*(r.FPsi-1)+tmp);
        end
        if ~doBoot
          % exact analytical confidence intervals according to Fidler &
          % Thompson 2001, p. 593, based on those for partial eta2
          tmp=(r.ssTot*(1-pEta2Ci))/r.dfErr;
          % main effects, interaction and contrast all in one step as the
          % only free variable is df
          tmp2=cat(1,r.factor.df,r.factor12.df,find(contrast.do))*[1 1];
          ci=(r.ssTot*pEta2Ci-tmp2.*tmp)./(r.ssTot+tmp);
          ciType=pEta2CiType;
        end
      end
      
    otherwise
      error(['internal error: unrecognized test not caught by input checks: ' curEs]);
  end
  
  % *********************************************************************
  % If data were NOT bootstrapped, all computations are done at this
  % point and the results can be placed into appropriate fields of output
  % variable stats. If they were bootstrapped, confidence intervals must
  % be computed and the ES statistics extracted from the first element of
  % variable es.
  % *********************************************************************
  if doBoot
    % determine confidence intervals from array of effect size measures
    % generated from bootstrapped data
    ci=prctile(es(:,2:end)',[alpha/2  1-alpha/2]'*100)';
    % retain first column; this is the es computed from real data
    es=es(:,1);
  end
  % place es and ci in esStore
  esStore(:,(tti-1)*3+1:tti*3)=[es ci];
  
  % finally, use dynamic fields to store currently computed measures in
  % output variable stats
  stats.(curEs)=es;
  stats.([curEs 'Ci'])=ci;
  stats.([curEs 'CiType'])=ciType;
end

% assemble table for display and optional output:
% the concatenated es and ci as cell array
esStore=num2cell(esStore);
% 'fillers' of empty cells
filla=cell(1,nEs*3);
% insert ci titles
ciTi=cell(2,nEs);
[ciTi{1,:}]=deal('ci_lo');
[ciTi{2,:}]=deal('ci_up');
esm=cat(1,esm,ciTi);
esm=esm(:)';
% column order:
% source | SS | df | MS | F | p | effect sizes (in the order requested)
summaryTable={...
    'SOURCE',                'SS',              'df',           'MS',              'F',              'p',              esm{1:end}};
caseStr={'between-subjects','within-subjects'};
caseCondit=[false true];
% loop over between/within subjects cases
for caseIx=1:2
  st=caseStr{caseIx};
  if any(isDep==caseCondit(caseIx))
    summaryTable=cat(1,summaryTable,...
      {caseStr{caseIx},        '---',             '---',          '---',             '---',            '---',          filla{1:end}});
    for g=1:2
      if isDep(g)==caseCondit(caseIx)
        tmp=['- ' factor(g).name];
        summaryTable=cat(1,summaryTable,...
          {tmp,              r.factor(g).ss(1), r.factor(g).df, r.factor(g).ms(1), r.factor(g).F(1), r.factor(g).p(1), esStore{g,:}});
      end
    end
    if ~any(isDep) && caseIx==1 || (any(isDep) && caseIx==2) 
      tmp=['- ' factor(1).name '*' factor(2).name];
      summaryTable=cat(1,summaryTable,...
        {tmp,                r.factor12.ss(1),  r.factor12.df, r.factor12.ms(1), r.factor12.F(1), r.factor12.p(1), esStore{3,:}});
    end
  end
end
if contrast.do
  summaryTable=cat(1,summaryTable,...
    {'contrast',               '---',             '---',          '---',             '---',            '---',          filla{1:end}},...
    {['- ' contrast.type ],  r.ssPsi(1),        1,              r.msPsi(1),        r.FPsi(1),        r.pPsi(1),        esStore{4,:}});
end

summaryTable=cat(1,summaryTable,...
  {'within-cells',           r.ssErr(1),        r.dfErr,        r.msErr(1),        [],               [],        filla{1:end}});

switch sum(isDep)
  case 1
    % mixed within-subjects
    % - index to between-subjects factor
    bi=find(~isDep);
    % - same for within-subjects 
    wi=find(isDep);
    tmp1=['- subj within ' factor(bi).name];
    tmp2=['- ' factor(wi).name '*subj within ' factor(bi).name];    
    summaryTable=cat(1,summaryTable,...
      {tmp1,              r.ssSubjWithinNonRM(1), r.dfSubjWithinNonRM, r.msSubjWithinNonRM(1), [], [], filla{1:end}},...
      {tmp2,              r.ssSubjRM(1),          r.dfSubjRM,          r.msSubjRM(1),          [], [], filla{1:end}});
  
  case 2
    % completely within-subjects 
    summaryTable=cat(1,summaryTable,...
      {'- subjects',         r.ssSubj(1),         r.dfSubj,        r.msSubj(1), [], [], filla{1:end}});
    for g=1:2
      tmp=['- ' factor(g).name '*subj'];
      summaryTable=cat(1,summaryTable,...                      
        {tmp,              r.factor(g).ssSubj(1), r.factor(g).dfSubj, r.factor(g).msSubj(1), [], [], filla{1:end}});
    end
    tmp=['- ' factor(1).name '*' factor(2).name  '*subj'];
    summaryTable=cat(1,summaryTable,...
      {tmp,              r.factor12.ssSubj(1), r.factor12.dfSubj, r.factor12.msSubj(1), [], [], filla{1:end}});
end
% finally, add line for total error
summaryTable=cat(1,summaryTable,...
  {'total',          r.ssTot(1), r.dfTot, r.msTot(1), [], [], filla{1:end}});

% display summaryTable
summaryTable  

if nargout>1
  varargout{1}=summaryTable;
end
  
% ========================= LOCAL FUNCTIONS ===============================
% ========================= LOCAL FUNCTIONS ===============================


% ------------------ START BY BOOTSTRAPPING (IF REQUESTED) ----------------
if doBoot
  % *********************************************************************
  % Here, individual, groupwise resampled instances of variable X are 
  % generated, expanding X in the second dimension: the first column
  % contains the original data, the others will be filled up with sampled
  % (with replacement) data.
  % *********************************************************************
  % variable bootSamX generated below shall contain indices to be used for
  % randomly sampling (with replacement) the original data (*** requires
  % groupIx to be integers incrementing in steps of 1! ***)
  switch sum(isDep)
    case 0
      % completely within-subjects: resample data within each group
      % independently of data points picked in other groups
      bootSamX=repmat(nan,size(x,1),nBoot);
      for g=1:nGroup
        bootSamX(groupIx{g},1:nBoot)=groupIx{g}(1)-1+ceil(nSample(g)*rand([nSample(g) nBoot]));
      end
    case 1
      % mixed within-subjects: the innermost loop applies a random sampling
      % index (offset-corrected) to all groups within a given column of
      % groupIx; this is done column by column of groupIx. This amounts to
      % a within-subjects design along the first factor; if the
      % within-subjects axis is along the second factor, we have to
      % temporarily transpose groupIx and nSample 
      if isequal(isDep,[false true])
        groupIx=groupIx';
        nSample=nSample';
      end
      bootSamX=repmat(nan,size(x,1),nBoot);
      for cIx=1:size(groupIx,2)
        % generic random sampling index for groups in current column
        partSamIx=ceil(nSample(1,cIx)*rand([nSample(1,cIx) nBoot]));
        for rIx=1:size(groupIx,1)
          bootSamX(groupIx{rIx,cIx},:)=partSamIx+groupIx{rIx,cIx}(1)-1;
        end          
      end
      % retranspose
      if isequal(isDep,[false true])
        groupIx=groupIx';
        nSample=nSample';
      end
      clear partSamIx;
    case 2
      % completely within-subjects: generate random sampling index for
      % first group and replicate for all other groups (which are identical
      % in size)
      bootSamX=repmat(ceil(nSample(1)*rand([nSample(1) nBoot])),[nGroup,1]);
      % add offset due to arrangement of data in a single-column array
      for g=1:nGroup
        bootSamX(groupIx{g},:)=bootSamX(groupIx{g},:)+groupIx{g}(1)-1;
      end
  end
  % from second column onwards fill with sampled data
  x(:,2:nBoot+1)=x(bootSamX);
  % delete resampling index 
  clear bootSam* 
end

%  --------------------- COMPUTATIONS PROPER ------------------------------
% within-groups SS, the individual groups' means and SS_error
r.ssErr=0;
for g=nGroup:-1:1
  % means of groups
  r.meanGroup(g,:)=sum(x(groupIx{g},:))/nSample(g);
  % within-groups SS
  r.ssErr=r.ssErr+sum((x(groupIx{g},:)-repmat(r.meanGroup(g,:),nSample(g),1)).^2);
end
% within-groups df
r.dfErr=sum(sum(nSample-1));
% within-groups MS
r.msErr=r.ssErr/r.dfErr;

% mean of all data points
r.mean=mean(x,1);
% unweighted arithmetic mean of all cell means
r.meanGrand=mean(r.meanGroup);
% cell size: harmonic mean 
r.hCellSz=nGroup/sum(1./nSample(:));
% cell size: arithmetic mean 
r.aCellSz=mean(nSample(:));

% for both factors, compute df, marginal means and SS_a, or SS_effect, the
% between-groups SS (or effect sum of squares)
switch ssType

  case 'unweighted'
    % ** unweighted means method using simple formulae for balanced data
    % and harmonic mean of cell size to accomodate unbalanced data
    for fi=1:nFactor
      % permute if fi==2
      s2i2=permute(s2i2,mod((1:nFactor)+fi,nFactor)+1);
      r.factor(fi).df=factor(fi).nLevel-1;
      r.factor(fi).ss=0;
      for lix=factor(fi).nLevel:-1:1
        % unweighted marginal means
        r.factor(fi).margMean(lix,:)=mean(r.meanGroup(s2i2(lix,:),:));
        r.factor(fi).ss=r.factor(fi).ss+...
          factor(mod(fi,nFactor)+1).nLevel*r.hCellSz*(r.factor(fi).margMean(lix,:)-r.meanGrand).^2;
      end
      r.factor(fi).ms=r.factor(fi).ss/r.factor(fi).df;
      % re-permute if fi==2
      s2i2=permute(s2i2,mod((1:nFactor)+fi,nFactor)+1);
    end
    % interaction terms
    r.factor12.df=r.factor(1).df*r.factor(2).df;
    r.factor12.ss=0;
    for lix1=factor(1).nLevel:-1:1
      for lix2=factor(2).nLevel:-1:1
        r.factor12.ss=r.factor12.ss+...
          r.hCellSz*(r.meanGroup(s2i2(lix1,lix2),:)-...
          r.factor(1).margMean(lix1,:)-...
          r.factor(2).margMean(lix2,:)+...
          r.meanGrand).^2;
      end
    end
    r.factor12.ms=r.factor12.ss/r.factor12.df;
  
  case 'Type III'
    % ** Method 1/Type III SS
    % i. create design matrix and bring it in proper shape:
    dm=dummyvar(group);
    % temporary helper var used for indexing
    tmp=[0 factor.nLevel];
    for fi=1:nFactor
      % within each factor, subtract last level's column from other columns
      % in dm
      dm(:,(1:tmp(fi+1))+tmp(fi))=dm(:,(1:tmp(fi+1))+tmp(fi))-...
        repmat(dm(:,sum(tmp(fi:fi+1))),1,tmp(fi+1));
      % df of main factor
      r.factor(fi).df=factor(fi).nLevel-1;
      % marginal mean: permute groupIx and/or s2i2 if fi==2
      groupIx=permute(groupIx,mod((1:nFactor)+fi,nFactor)+1);
      s2i2=permute(s2i2,mod((1:nFactor)+fi,nFactor)+1);
      for lix=factor(fi).nLevel:-1:1
        % unweighted marginal mean (=mean of cell means across rows or
        % columns regardless of cell size)
        r.factor(fi).margMean(lix,:)=mean(r.meanGroup(s2i2(lix,:),:));
        % % marginal mean from raw scores - not recommended
        % r.factor(fi).margMean(lix,:)=mean(x(cat(1,groupIx{lix,:}),:));
      end
      % re-permute if fi==2
      groupIx=permute(groupIx,mod((1:nFactor)+fi,nFactor)+1);
      s2i2=permute(s2i2,mod((1:nFactor)+fi,nFactor)+1);
    end
    % within each factor, delete columns corresponding to last level in dm
    dm(:,cumsum(tmp(2:end)))=[];
    % now add interaction columns:
    ct=prod([r.factor.df]);
    for lix1=r.factor(1).df:-1:1
      for lix2=r.factor(2).df:-1:1
        dm(:,ct+sum([r.factor.df]))=dm(:,lix1).*dm(:,lix2+r.factor(1).df);
        ct=ct-1;
      end
    end
    % df of interaction
    r.factor12.df=r.factor(1).df*r.factor(2).df;
    
    % ii. compute regression-related SS
    % the vector of ones which must be added to obtain a constant term from
    % the regressions
    vo=ones(size(x,1),1);
    % regressions:
    % - a helper index, pointing to the column in the design matrix dm
    % corresponding to the last entry for 'main' effects
    tmpIx=sum([r.factor.df]);
    % - a,b,ab
    ss_abab=bbregress(x,[vo dm]);
    % - a,b
    ss_ab=bbregress(x,[vo dm(:,1:tmpIx)]);
    % - a,ab
    ss_aab=bbregress(x,[vo dm(:,[1:r.factor(1).df  tmpIx+1:end])]);
    % - b,ab
    ss_bab=bbregress(x,[vo dm(:,[r.factor(1).df+1:end])]);
    
    % iii. finally, SS and MS associated with main and interaction effects
    r.factor(1).ss=ss_abab-ss_bab;
    r.factor(2).ss=ss_abab-ss_aab;
    r.factor12.ss=ss_abab-ss_ab;
    r.factor(1).ms=r.factor(1).ss/r.factor(1).df;
    r.factor(2).ms=r.factor(2).ss/r.factor(2).df;
    r.factor12.ms=r.factor12.ss/r.factor12.df;
    clear dm ct vo tmp* ss*
end

% SS_t, the sum of all ss, computed straight from the data
r.ssTot=sum((bsxfun(@minus,x,r.mean)).^2);

% corresponding df, ms
r.dfTot=r.factor(1).df+r.factor(2).df+r.factor12.df+r.dfErr;
r.msTot=r.ssTot/r.dfTot;

% compute contrasts and associated SS
if contrast.do
  switch contrast.type
    case {'simple1','simple2'}
      % linear index to all non-nan elements in contrast.weight...
      mIx=find(isfinite(contrast.weight(:)));
      % ...to be used for calculation of contrast
      r.psi=sum(repmat(contrast.weight(mIx),1,nBoot+1).*r.meanGroup(mIx,:));
      % denominator for SS, which will later also be used for computation
      % of confidence interval
      r.denomSsPsi=sum(contrast.weight(mIx).^2./nSample(mIx));
      % SS
      r.ssPsi=r.psi.^2./r.denomSsPsi;
      
    case {'main1','main2'}
      % index to factor under consideration (note usage of last letter of
      % contrast.type)
      ix=sscanf(contrast.type(end),'%i');
      % index to other factor
      oix=mod(ix,2)+1;
      % instead of computing contrast from marginal means replicate
      % contrast weights according to design and compute cellwise, then
      % average
      cw=repmat(contrast.weight,circshift([1 factor(oix).nLevel],[0 ix-1]));
      r.psi=sum(repmat(cw(:),1,nBoot+1).*r.meanGroup)/factor(oix).nLevel;
      % denominator for SS
      r.denomSsPsi=sum(contrast.weight.^2./sum(nSample,oix));
      % SS
      r.ssPsi=r.psi.^2./r.denomSsPsi;
      
    case 'interaction'
      r.psi=sum(repmat(contrast.weight(:),1,nBoot+1).*r.meanGroup);
      % denominator for SS
      r.denomSsPsi=sum(contrast.weight(:).^2./nSample(:));
      % SS
      r.ssPsi=r.psi.^2./r.denomSsPsi;
  end
  % as df is always 1 MS=SS, but create field nonetheless
  r.msPsi=r.ssPsi;
end

% compute design-dependent SS, F, p for main and interaction effects as
% well as contrast
switch sum(isDep)
  case 0
    % *****************************
    % completely between-subjects:
    % *****************************
    % effect error SS, F, p for main and interaction effects
    for fi=1:2
      % the effect error SS - generally it depends on the design; it is
      % needed for the computation of partialeta2
      r.factor(fi).ssEffErr=r.ssErr;
      r.factor(fi).F=r.factor(fi).ms./r.msErr;
      r.factor(fi).p=1-fcdf(r.factor(fi).F,r.factor(fi).df,r.dfErr);
    end
    r.factor12.ssEffErr=r.ssErr;
    r.factor12.F=r.factor12.ms./r.msErr;
    r.factor12.p=1-fcdf(r.factor12.F,r.factor12.df,r.dfErr);
      
    % effect error SS, df, ci, F, t, p for contrast
    if contrast.do
      r.ssEffErr_Psi=r.ssErr;
      % the df needed for the computation of the t noncentrality parameter
      % for the contrast (NOT the df of the contrast, which is by
      % definition always 1)
      r.df4nc_Psi=r.dfErr;
      % confidence intervals (original data only)
      r.ciPsi=repmat(r.psi(1),1,2) + sqrt(r.msErr(1)*r.denomSsPsi) * tinv([alpha/2 1-alpha/2],r.dfErr);
      % F (original & bootstrapped, because bootstrapped Fs are needed for
      % the computation of CI of partial eta2 and omega2)
      r.FPsi=r.msPsi./r.msErr;
      % compute t from F and sign of contrast (original data only)
      r.tPsi=sign(r.psi(1)).*sqrt(r.FPsi(1));
      r.pPsi=1-fcdf(r.FPsi(1),1,r.dfErr(1));
    end
    
  case 1
    % *****************************
    % mixed within-subjects:
    % *****************************
    % identify the repeated-measures (RM) factor
    rmFacIx=find(isDep);
    nonRmFacIx=find(~isDep);
    % the code further below assumes the second factor to be
    % repeated-measures (i.e. its levels span the rows of groupIx). If the
    % first factor is RM, we have to temporarily transpose all variables
    % whose 2D layout reflects the factors and levels.
    if rmFacIx==1
      % first factor is repeated-measures: 
      % - marginal mean
      margMean=r.factor(2).margMean;
      % - transpose groupIx, nSample and s2i2
      groupIx=groupIx';
      nSample=nSample';
      s2i2=s2i2';
    else
      % second factor is repeated-measures:
      % - marginal mean
      margMean=r.factor(1).margMean;
    end
    % (note on the side: the quantities are computed in a 'raw' fashion
    % here, not by subtraction, as this may be necessary in a future
    % version of the code allowing for unbalanced data across the
    % non-repeated measures factor)
    nSubj=nSample(:,1);
    % df
    r.dfSubjWithinNonRM=sum(nSample(:,1)-1);
    r.dfSubjRM=r.dfSubjWithinNonRM*(factor(rmFacIx).nLevel-1);
    % initialize
    r.ssSubjWithinNonRM=0;
    r.ssSubjRM=0;
    % loop over levels of non-RM factor
    for rIx=factor(nonRmFacIx).nLevel:-1:1
      tmpHelper=cat(2,groupIx{rIx,:});
      % loop over subjects
      for sIx=nSubj(rIx):-1:1
        % index to entries of current subject in current level of non-RM factor
        tmpIx=tmpHelper(sIx,:);
        % mean of those
        tmpMn=mean(x(tmpIx,:));
        % SS of subjects within levels of the non-RM factor (error between)
        r.ssSubjWithinNonRM=r.ssSubjWithinNonRM+factor(rmFacIx).nLevel*...
          (tmpMn-margMean(rIx,:)).^2;
        % loop over RM factor
        for cIx=factor(rmFacIx).nLevel:-1:1
          % interaction term: RM factor x subjects within levels of the
          % non-RM factor (error within)
          r.ssSubjRM=r.ssSubjRM+...
            (x(tmpIx(cIx),:)-tmpMn-r.meanGroup(s2i2(rIx,cIx),:)+margMean(rIx,:)).^2;
        end
      end
    end
    % MS
    r.msSubjWithinNonRM=r.ssSubjWithinNonRM/r.dfSubjWithinNonRM;
    r.msSubjRM=r.ssSubjRM/r.dfSubjRM;
    if rmFacIx==1
      % re-transpose variables transposed above
      groupIx=groupIx';
      nSample=nSample';
      s2i2=s2i2';
    end
    % effect error SS, F, t, p
    % - between-subjects factor
    r.factor(nonRmFacIx).ssEffErr=r.ssSubjWithinNonRM;
    r.factor(nonRmFacIx).F=r.factor(nonRmFacIx).ms./r.msSubjWithinNonRM;
    r.factor(nonRmFacIx).p=1-fcdf(r.factor(nonRmFacIx).F,r.factor(nonRmFacIx).df,r.dfSubjWithinNonRM);
    % - within-subjects factor
    r.factor(rmFacIx).ssEffErr=r.ssSubjRM;
    r.factor(rmFacIx).F=r.factor(rmFacIx).ms./r.msSubjRM;
    r.factor(rmFacIx).p=1-fcdf(r.factor(rmFacIx).F,r.factor(rmFacIx).df,r.dfSubjRM);
    % - interaction
    r.factor12.ssEffErr=r.ssSubjRM;    
    r.factor12.F=r.factor12.ms./r.msSubjRM;
    r.factor12.p=1-fcdf(r.factor12.F,r.factor12.df,r.dfSubjRM);
    % - ? contrast: 
    if contrast.do
      % if it is an interaction contrast or a single factor-contrast across
      % the within-subjects factor...
      if strcmp(contrast.type,'interaction') || rmFacIx==sscanf(contrast.type(end),'%i')
        r.ssEffErr_Psi=r.ssSubjRM;
        r.df4nc_Psi=r.dfSubjRM;
        % confidence intervals (original data only)
        r.ciPsi=repmat(r.psi(1),1,2) + sqrt(r.msSubjRM(1)*r.denomSsPsi) * tinv([alpha/2 1-alpha/2],r.dfSubjRM);
        % F (original & bootstrapped)
        r.FPsi=r.msPsi./r.msSubjRM;
        % p (original data only)
        r.pPsi=1-fcdf(r.FPsi(1),1,r.dfSubjRM);
      else
        r.ssEffErr_Psi=r.ssSubjWithinNonRM;
        r.df4nc_Psi=r.dfSubjWithinNonRM;
        % confidence intervals (original data only)
        r.ciPsi=repmat(r.psi(1),1,2) + sqrt(r.msSubjWithinNonRM(1)*r.denomSsPsi) * tinv([alpha/2 1-alpha/2],r.dfSubjWithinNonRM);
        % F (original & bootstrapped)
        r.FPsi=r.msPsi./r.msSubjWithinNonRM;
        % p (original data only)
        r.pPsi=1-fcdf(r.FPsi(1),1,r.dfSubjWithinNonRM);
      end
      % compute t from F and sign of contrast
      r.tPsi=sign(r.psi).*sqrt(r.FPsi);
    end

  case 2
    % *****************************
    % completely within-subjects:
    % *****************************
    nSubj=nSample(1);
    r.dfSubj=nSubj-1;
    r.ssSubj=zeros(1,nBoot+1);
    for sIx=nSubj:-1:1
      % means of subjects
      r.meanSubj(sIx,:)=mean(x(sIx:nSubj:end,:));
      % SS_subj, between-subjects SS
      r.ssSubj=r.ssSubj+nGroup*(r.meanSubj(sIx,:)-r.meanGrand).^2;
    end
    r.msSubj=r.ssSubj/r.dfSubj;
    % loop over factors
    for fi=1:nFactor
      % permute if fi==2
      groupIx=permute(groupIx,mod((1:nFactor)+fi,nFactor)+1);
      % df
      r.factor(fi).dfSubj=r.dfSubj*r.factor(fi).df;
      % initialize
      tmpSs=0;
      % loop over levels
      for lix=factor(fi).nLevel:-1:1
        % loop over subjects
        for sIx=nSubj:-1:1
          % index to entries of current subject in current level of current factor
          tmpIx=intersect(sIx:nSubj:size(x,1),cat(1,groupIx{lix,:}));
          % mean of those
          tmpMn=mean(x(tmpIx,:));
          % SS 
          tmpSs=tmpSs+factor(mod(fi,nFactor)+1).nLevel*(tmpMn-r.meanGrand).^2;
        end
      end
      % re-permute groupIx if fi==2
      groupIx=permute(groupIx,mod((1:nFactor)+fi,nFactor)+1);
      % obtain factor x subject term by subtraction
      r.factor(fi).ssSubj=tmpSs-r.factor(fi).ss-r.ssSubj;
      % ms
      r.factor(fi).msSubj=r.factor(fi).ssSubj/r.factor(fi).dfSubj;
    end
    % (interaction between factors) x subjects term:
    % - df
    r.factor12.dfSubj=r.dfSubj*r.factor12.df;
    % - obtain SS by subtraction
    r.factor12.ssSubj=r.ssErr-sum(cat(1,r.factor(:).ssSubj))-r.ssSubj;
    % - ms
    r.factor12.msSubj=r.factor12.ssSubj/r.factor12.dfSubj;
    % effect error SS, F, p
    for fi=1:2
      r.factor(fi).ssEffErr=r.factor(fi).ssSubj;
      r.factor(fi).F=r.factor(fi).ms./r.factor(fi).msSubj;
      r.factor(fi).p=1-fcdf(r.factor(fi).F,r.factor(fi).df,r.factor(fi).dfSubj);
    end
    r.factor12.ssEffErr=r.factor12.ssSubj;
    r.factor12.F=r.factor12.ms./r.factor12.msSubj;
    r.factor12.p=1-fcdf(r.factor12.F,r.factor12.df,r.factor12.dfSubj);
    
    % effect error SS, df, ci, F, t, p for contrast
    if contrast.do
      if strcmp(contrast.type,'interaction')
        r.ssEffErr_Psi=r.factor12.ssSubj;
        r.df4nc_Psi=r.factor12.dfSubj;
        % confidence intervals (original data only)
        r.ciPsi=repmat(r.psi(1),1,2) + sqrt(r.factor12.msSubj(1)*r.denomSsPsi) * tinv([alpha/2 1-alpha/2],r.factor12.dfSubj);
        % F (original & bootstrapped)
        r.FPsi=r.msPsi./r.factor12.msSubj;
        % p (original data only)
        r.pPsi=1-fcdf(r.FPsi(1),1,r.factor12.dfSubj);
      else
        fi=sscanf(contrast.type(end),'%i');
        r.ssEffErr_Psi=r.factor(fi).ssSubj;
        r.df4nc_Psi=r.factor(fi).dfSubj;
        % confidence intervals (original data only)
        r.ciPsi=repmat(r.psi(1),1,2) + sqrt(r.factor(fi).msSubj(1)*r.denomSsPsi) * tinv([alpha/2 1-alpha/2],r.factor(fi).dfSubj);
        % F (original & bootstrapped)
        r.FPsi=r.msPsi./r.factor(fi).msSubj;
        % p (original data only)
        r.pPsi=1-fcdf(r.FPsi(1),1,r.factor(fi).dfSubj);
      end
      % compute t from F and sign of contrast (original data only)
      r.tPsi=sign(r.psi(1)).*sqrt(r.FPsi(1));
    end
end

function ssReg=bbregress(dep,idep)
% 'bare bones' regression which performs only computations essential in the
% context of ANOVA computations cast as a linear regression, and returns
% the regression summed squares
% - coefficients
c=idep\dep;
% - expected values
e=idep*c;
% - regression SS
ssReg=sum((e-repmat(mean(dep),size(dep,1),1)).^2);

function mesdplot(x,groupIx,nSample,factor,isDep,fName,contrast)
% ** function mesdplot(x,groupIx,nSample,factor,isDep,fName,contrast)
% is an accessory function for mes2way.m which plots data for a twoway
% factorial analysis. This function is meant to be called from within
% mes2way.m, therefore there are neither explanations of input variables
% nor error checks here as they would be redundant. If you use it from
% outside mes2way make sure the input variables are properly shaped.
% Specifics:
% - each group of samples resides in one subplot
% - levels of factor 1 down the columns
% - levels of factor 2 along the rows
% - samples are circles, group mean is a horizontal line
% - dependent (repeated measures) data are depicted with individual colors
% - background color of the subplots reflects sign and value of contrast
%   weights (if they were specified)
% - levels of factors are written in upper left corner of each subplot

% -------------------------------------------------------------------------
% Measures of Effect Size Toolbox Version 1.4, January 2015
% Code by Harald Hentschke (University of T?bingen) and 
% Maik St?ttgen (University of Bochum)
% For additional information see Hentschke and St?ttgen, 
% Eur J Neurosci 34:1887-1894, 2011
% -------------------------------------------------------------------------

% ------- PART I: PREPARATORY WORKS ----------
% total number of groups
nGroup=numel(nSample);
% number of factors
nL1=factor(1).nLevel;
nL2=factor(2).nLevel;
%  max & min value in whole data set to determine axis limits
mima=[min(x(:)) max(x(:))];
% ...stretched a bit  so that all data points are fully contained in plot
mima=mima+diff(mima)*[-.1 .1];
% maximal group (cell) sample size within data set
maxNs=max(nSample(:));
% define colors for individual data points in case of dependent data
switch num2str(isDep')'
  case '11'
    % completely within-subjects: sample sizes are equal in all groups
    cm=colorcube(nSample(1));
  case '10'
    % repeated measures along first factor:
    % cumulative sum of sample sizes needed for... 
    cumSampleSize=[0 cumsum(nSample(1,:))];
    % ...indexing this colormap
    cm=colorcube(cumSampleSize(end));
  case '01'
    % repeated measures along second factor:
    cumSampleSize=[0 cumsum(nSample(:,1))'];
    cm=colorcube(cumSampleSize(end));
  otherwise
    % we don't need a colormap
    cm=[];
end
% define markersize depending on average number of samples per cell
% (empirical values, change according to gusto)
ms=3+7*exp(-mean(nSample(:))*.02);
% set random number generator to fixed state (old syntax for downward
% compatibility - use rng('default') for newer Matlab versions)
rand('seed',0);
% abscissa values: centered at zero; horizontal spread depending on max
% number of samples
absci=(rand(maxNs,1)-.5)*maxNs/(maxNs+10);
% expand contrast weights (makes plotting below easier)
if contrast.do
  switch contrast.type
    case 'main1'
      contrast.weight=repmat(contrast.weight,1,nL2);
    case 'main2'
      contrast.weight=repmat(contrast.weight,nL1,1);
    otherwise
      % do nothing: contrast weights are already handily shaped
  end
end
      
% ------- PART II: PLOT ----------
figure
% defaults for all axes on plot
set(gcf,'DefaultlineMarkerSize',ms)
set(gcf,'DefaultaxesFontSize',8)
set(gcf,'DefaulttextFontSize',8);
for r=1:nL1
  for c=1:nL2
    sph=subplot(nL1,nL2,(r-1)*nL2+c);
    hold on
    set(gca,'xtick',[],'ygrid','on','box','on'); 
    if any(isDep)
      % define which colors to use for individual samples depending on
      % which factor is repeated measures
      switch num2str(isDep')'
        case '11'
          % completely within-subjects (dependent) data: identical colors
          % of samples across all groups
          colorMapIndex=1:nSample(r,c);
        case '10'
          % repeated measures along first factor
          colorMapIndex=cumSampleSize(c)+1:cumSampleSize(c+1);
        case '01'
          % repeated measures along second factor
          colorMapIndex=cumSampleSize(r)+1:cumSampleSize(r+1);          
        otherwise
          error('something funny happened to isDep')
      end
      % plot individual samples
      for g=1:nSample(r,c)
        ph=plot(absci(g),x(groupIx{r,c}(g)),'o');
        % make outline color a darker version of fill color
        set(ph,'color',cm(colorMapIndex(g),:)*.5,'markerfacecolor',cm(colorMapIndex(g),:))
      end
    else
      % independent data: open black symbols
      ph=plot(absci(1:nSample(r,c)),x(groupIx{r,c}),'ko');
    end
    % plot mean as horizontal line
    line([-.3 .3],mean(x(groupIx{r,c})*[1 1]),'color','k','linewidth',2)
    % same scaling for all plots
    axis([-.55 .55 mima])
    if contrast.do
      % current contrast weight
      cw=contrast.weight(r,c);
      % axis color indicating sign and value of contrast weights (if
      % contrast is specified):
      if isfinite(cw)
        if cw>0
          % positive values in red (maximally reaching a somewhat blunted
          % hue)
          cCol=[1 1-cw 1-cw]*.7+.3;
        elseif cw<0
          % negative values in blue (ditto)
          cCol=[1+cw 1+cw 1]*.7+.3;
        else
          % zero contrast weights in light gray
          cCol=[.88 .88 .88];
        end
        set(gca,'color',cCol)
      end
    end
    % illustrate assignment of rows/columns to factors
    if r==1 && c==1
      title([factor(2).name ' --->'])
      ylabel(['<--- ' factor(1).name]); 
    end
    % finally, show levels of factors in upper left corner
    text(-.5, mima(2)-.1*diff(mima),...
      ['[' num2str(factor(1).level(r)) ',' num2str(factor(2).level(c)) ']']);
  end
end


function stats=mestab(table,varargin)
% ** function stats=mestab(table,varargin) computes basal parameters and
% effect size measures from N by M tables of categorical outcomes (mostly 
% 2 by 2 tables) and according confidence intervals where applicable and
% possible. The exact analyses/parameters need not be specified; as the
% computations involved are trivial all parameters will be computed and
% placed in output variable stats. All input parameters except table are
% optional and must be specified as parameter/value pairs in any order,
% e.g. as in
%      mestab(table,'confLevel',.90)
%
%                           INPUT ARGUMENTS
%                           ---------------
% stats=mestab(table)  computes all implemented parameters and effect size 
%   measure(s) listed below in TABLE 1. table must be minimally a 2 by 2 
%   table of frequencies of occurrence with the following row and column 
%   order:
%                          CHARACTERISTIC 1
%                        | - PRESENT         - NOT PRESENT
%       ------------------------------------------------
%       CHARACTERISTIC 2 |      
%       - PRESENT        |
%       - NOT PRESENT    |
%
%                        or (as a specific example)
% 
%                        | TRUE POSITIVE     TRUE NEGATIVE
%       ------------------------------------------------
%       TEST POSITIVE    |
%       TEST NEGATIVE    |
% 
%   To give another example, the rows are the dichotomous values of the
%   independent variable (e.g. control and treatment group) and the columns
%   those of the dependent variable (e.g. outcome):
%
%                        | OUTCOME NEGATIVE  OUTCOME POSITIVE
%       ------------------------------------------------
%       CONTROL          |
%       TREATMENT        |
%
%   If table is larger than 2 by 2 the only parameter to be computed is
%   Cramer's V.
% stats=mestab(...,'confLevel',0.90)  computes 90 % confidence intervals of 
%   the statistic in question (95 % ci are the default)
%
%                         OUTPUT ARGUMENTS
%                         ----------------
%   The results of the computations are placed into fields of output
%   argument stats with names corresponding to the analyses listed in
%   TABLE 1, e.g. 
%     .riskDifference
%     .riskRatio
%     ...
%   Where applicable, confidence intervals will be placed in fields 
%     .riskDifferenceCi
%     .riskRatioCi
%     ...
% Additional output fields are:
%     .confLevel (confidence level)
% -------------------------------------------------------------------------
% TABLE 1: PARAMETERS COMPUTED
% -------------------------------------------------------------------------
% FIELD NAME           QUANTITIY COMPUTED, SYNONYMS
%        -- measures for 2 by 2 tables --
% riskDifference       risk difference, proportion difference
% riskRatio            risk ratio, rate ratio
% oddsRatio            odds ratio
% phi                  degree of association
% baseRate             base rate
% sensitivity          sensitivity
% specificity          specificity
% posPredictValue      positive predictive value
% negPredictValue      negative predictive value
% successTreat         binomial effect size display: success rate (treated)
% successCtrl          binomial effect size display: success rate (control)
%        -- measures for M by N tables --
% cramerV              Cramer's V
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Measures of Effect Size Toolbox Version 1.4, January 2015
% Code by Harald Hentschke (University of T?bingen) and 
% Maik St?ttgen (University of Bochum)
% For additional information see Hentschke and St?ttgen, 
% Eur J Neurosci 34:1887-1894, 2011
% -------------------------------------------------------------------------

% ----- default values & varargin -----
% standard values of optional input arguments
confLevel=.95;
% check variable number of input args:
% - even number (parameter AND value specified)?
% - convert to proper upper/lowercase spelling as above
% - check whether valid parameter was given
% - overwrite defaults by values, if specified
nvarg=numel(varargin);
v=who;
if nvarg
  if nvarg/2==round(nvarg/2)
    for g=1:2:nvarg
      % check which optional input parameter is given, ignoring case, and
      % impose spelling used internally
      ix=find(strcmp(lower(varargin{g}),lower(v)));
      if ~isempty(ix)
        varargin{g}=v{ix};
      else
        error(['invalid optional input parameter (' varargin{g} ') was specified']);
      end
    end
    % finally, the assignment of values to parameters
    pvpmod(varargin);
  else
    error('optional input parameters must be specified as parameter/value pairs, e.g. as in ''par1'',1')
  end
end

% -------------------------------------------------------------------------
% ------- PART I: CHECKS OF INPUT ARGS & OTHER PREPARATORY WORKS ----------
% -------------------------------------------------------------------------
alpha=1-confLevel;
% --- check input table
[nRow nCol]=size(table);
if nRow<2 || nCol<2
  error('table size must at least be 2 by 2');
end
% what kinda table?
if nRow==2 && nCol==2
  is2by2=true;
else
  is2by2=false;
end
% reject foul values
if any(any(~isfinite(table)))
  error('input variable table contains foul data (nan or inf))');
end
% --- check other input arguments
if confLevel<=0 || confLevel>=1
  error('input variable ''confLevel'' must be a scalar of value >0 and <1');
end

% -------------------------------------------------------------------------
% ------- PART II: THE WORKS ----------
% -------------------------------------------------------------------------
if is2by2
  % first, for convenience, assign entries in table to alphabetic letters
  % according to Kline 2004 (p. 146) (all lowercase, though).
  % note linear indexing of a 2D-var
  a=table(1);
  b=table(3);
  c=table(2);
  d=table(4);
  
  % compute basal parameters 
  stats.baseRate=(a+c)/(a+b+c+d);
  stats.sensitivity=a/(a+c);
  stats.specificity=d/(b+d);
  stats.posPredictValue=a/(a+b);
  stats.negPredictValue=d/(c+d);
  
  % now compute values based on proportions. As all computations involved
  % are trivial and variable sizes miniscule let's not worry about their
  % redundancy and proliferation, respectively, and use the same
  % nomenclature as Kline 2004
  p_c=a/(a+b);
  p_t=c/(c+d);
  % will be needed for confidence intervals
  z2tailFactor=norminv(alpha/2)*[1 -1];
  
  % risk difference  
  stats.riskDifference=p_c-p_t;
  % confidence intervals
  stats.riskDifferenceCi=stats.riskDifference + ...
    sqrt(p_c*(1-p_c)/(a+b) + p_t*(1-p_t)/(c+d)) * z2tailFactor;
  % risk ratio
  stats.riskRatio=p_c/p_t;
  % confidence intervals of risk ratio: for better readibility compute in
  % two steps
  tmp=log(stats.riskRatio) + ...
    sqrt((1-p_c)/((a+b)*p_c) + (1-p_t)/((c+d)*p_t)) * z2tailFactor;
  stats.riskRatioCi=exp(tmp);
  % odds ratio
  stats.oddsRatio=a*d/(b*c);
  % ci: same story 
  tmp=log(stats.oddsRatio) + ...
    sqrt( 1/((a+b)*p_c*(1-p_c)) + 1/((c+d)*p_t*(1-p_t))) * z2tailFactor;
  stats.oddsRatioCi=exp(tmp);
  % the measure of association termed phi [equivalent:
  % stats.phi=(a*d-b*c)/sqrt((a+b)*(c+d)*(a+c)*(b+d));]
    
  % ** use code for Cramer's V for confidence intervals
  [stats.phi stats.phiCi chi2]=cramerv(table,nRow,nCol,confLevel,is2by2);
  % finally, binomial effect size display (to quote Randolph & Edmondson
  % 2005, 'What would the correlationally equivalent effect of the
  % treatment be if 50% of the participants had the occurrence and 50% did
  % not and 50% received treatment and 50% did not?')
  stats.successTreat=.5+stats.phi/2;
  stats.successCtrl=.5-stats.phi/2;
else
  % the only thing to be computed here is Cramer's V (and chi2)
  [stats.cramerV stats.cramerVCi chi2]=cramerv(table,nRow,nCol,confLevel,is2by2);
end

% for the sake of being complete, add chi square to stats
stats.chi2=chi2;


% ========================= LOCAL FUNCTIONS ===============================
% ========================= LOCAL FUNCTIONS ===============================

function [es,esCi,chi2]=cramerv(table,nRow,nCol,confLevel,is2by2)
% this function computes Cramer's V, including exact analytical CI
% ** NOTE: in the case of 2 by 2 tables Cramer's V is identical to phi
% except possibly for the sign), which will be taken care of in the last
% lines
colSum=sum(table);
rowSum=sum(table,2);
n=sum(sum(table));
k=min(nRow,nCol);
df=(nRow-1)*(nCol-1);
% expected frequency of occurrence in each cell: product of row and
% column totals divided by total N
ef=(rowSum*colSum)/n;
% chi square stats
chi2=(table-ef).^2./ef;
chi2=sum(chi2(:));
% Cramer's V
es=sqrt(chi2/(n*(k-1)));
% CI (Smithson 2003, p. 40)
ncp=ncpci(chi2,'X2',df,'confLevel',confLevel);
esCi=sqrt((ncp+df)/(n*(k-1)));
% in case we are dealing with 2 by 2 tables heed sign
if is2by2
  if det(table)<0
    es=es*-1;
    esCi=fliplr(esCi)*-1;
  end
end


function ci=ncpci(x,fType,df,varargin)
%  ** function ci=ncpci(x,fType,df,varargin)
% iteratively approaches two-sided confidence intervals for the
% noncentrality parameter of a noncentral Chi square (abbreviated X2), F or
% t distribution with degrees of freedom df, given an abscissa value (X2, F
% or t value). This is achieved by varying the X2, F or t noncentrality
% parameter of the corresponding probability distribution function (pdf)
% until the given abscissa value is, within a certain precision, at the
% percentile values required for the confidence interval (2.5th and 97.5th
% percentile for lower and upper 95 % confidence intervals, respectively).
% All input parameters listed below except x, fType and df are
% optional and must be specified as parameter/value pairs, e.g. as in
%          ncpci(x,'t',df,'confLevel',.9)
%
%               >>> INPUT VARIABLES >>>
% NAME          TYPE/DEFAULT      DESCRIPTION
% x             double scalar     X2, F or t value
% fType         char              'X2','F' or 't'
% df            scalar or array   degrees of freedom 
%                                 (F pdf: [numerator denominator])
% confLevel     double, 0.95      confidence level
% prec          double scalar,    precision: iteration will run until the
%                1e-6             estimated percentile is <=prec away from
%                                 the requested percentile
% doAnimate     logical,false     if true, the iteration process will be
%                                 graphically displayed in a figure window
%                     
%               <<< OUTPUT VARIABLES <<<
% NAME          TYPE/DEFAULT           DESCRIPTION
% ci            2 element array        confidence intervals


% -------------------------------------------------------------------------
% Measures of Effect Size Toolbox Version 1.4, January 2015
% Code by Harald Hentschke (University of T?bingen) and 
% Maik St?ttgen (University of Bochum)
% For additional information see Hentschke and St?ttgen, 
% Eur J Neurosci 34:1887-1894, 2011
% -------------------------------------------------------------------------

% defaults
prec=1e-6;
confLevel=.95;
doAnimate=false;
% replace defaults by input, if any
pvpmod(varargin);

% convert df to cell for automatic expansion of parameters
df=num2cell(df);
% convert confidence level to alpha
alpha=1-confLevel;
% target p values 
pTarget=[1-alpha/2  alpha/2];

% --- error checks, assignments of function handles, etc.
% are we dealing with pdf defined only for positive abscissa values?
isPosPdf=ismember(fType,{'X2','F'});
% if so...
if isPosPdf && x<0
  error('input arg ''x'' is negative but must be positive for X2 and F distributions')
end
% start index for outermost loop below, determining whether lower CI shall
% be computed or not
loopStartIx=1;
switch fType
  case 'X2'
    curPdf=@ncx2pdf;    
    curCdf=@ncx2cdf;
    curInv=@chi2inv;
    % abscissa limits for plots (if doAnimate==true): first row for lower
    % CI, second row for upper CI
    abscissLim=[0 2*x;0 5*x];
    % check: if cdf of x with noncentrality parameter 0 is less than
    % 1-alpha/2 don't even start on the lower CI because the iteration will
    % not converge (that is, there is no lower CI for given values of x and
    % df)
    if chi2cdf(x,df{:})<1-alpha/2
      % lower CI cannot be constructed as it is too close to zero - set to
      % NaN
      ci=nan;
      loopStartIx=2;
    end
    
  case 'F'
    curPdf=@ncfpdf;    
    curCdf=@ncfcdf;
    curInv=@finv;
    abscissLim=[0 2*x;0 5*x];
    % similar check as above
    if fcdf(x,df{:})<1-alpha/2
      % lower CI cannot be constructed as it is too close to zero - set to
      % NaN
      ci=nan;
      loopStartIx=2;
    end
    
  case 't'
    curPdf=@nctpdf;
    curCdf=@nctcdf;
    curInv=@tinv;
    abscissLim=x+[-4 2;-2 4]*sqrt(abs(x));
    
  otherwise
    error('illegal distribution function specified');
end

if prec>.001
  warning('results will be inaccurate - set input parameter ''prec'' to a lower value');
end

if doAnimate
  fh=figure;
  ph0=plot(x,0,'k^');
  hold on
  set(ph0,'markerfacecolor','k','markersize',6);
  ph=[];
  ti={'lower CI','upper CI'};
end

% loop twice: first lower ci (but see above), then upper ci
for iIx=loopStartIx:2
  % determine initial values: there are probably better ways of estimating
  % the limits of ncp for X2 and F pdfs than the guesses below (which work
  % best if the X2/F/t value is small)
  switch fType
    case 'X2'
      if iIx==1
        % lower CI
        ncp=x+curInv(pTarget(iIx),df{:});
      else
        % upper CI
        ncp=5*x;
      end
      
    case 'F'
      if iIx==1
        ncp=x+curInv(pTarget(iIx),df{:});
      else
        ncp=10*x;
      end
      
    case 't'
      % as a rough first approximation, assume that lower/upper limit of
      % ncp is close to corresponding percentiles of central pdfs
      if iIx==1
        ncp=x+curInv(pTarget(iIx),df{:});
      else
        ncp=x-curInv(pTarget(iIx),df{:});
      end
  end
  
  % interval of first estimates: guessed ncp enlarged by x/2 on either side
  ncp=ncp+abs(x)*[-.5 .5];
  % p values of current estimates
  p=curCdf(x,df{:},ncp);
  % deviations of p of current noncentral x pdfs from target p value
  deltaP=p-pTarget(iIx);
  nIter=1;
  if doAnimate
    ph=plotPdf(x,ncp,ph,curPdf,df,iIx,nIter,abscissLim,ti);
  end
  % while desired precision is not reached...
  while ~any(abs(deltaP)<=prec)
    if all(deltaP>0)
      % shift interval to the right by one interval length
      ncp=[ncp(2) ncp(2)+abs(diff(ncp))];
    elseif all(deltaP<0)
      % shift left by one interval length
      ncp=[ncp(1)-abs(diff(ncp)) ncp(1)];
    else
      % halve interval around mean
      ncp=mean(ncp)+.25*abs(diff(ncp))*[-1 1];
    end
    % X2 and F distributions need an extra check: the lower ncp must be >=0
    if isPosPdf
      if ncp(1)<0
        ncp(1)=0;
      end
      % if both values of ncp are zero here the upper CI is zero, too, so
      % stop here
      if ~any(ncp)
        break
      end
    end
    % p values of current estimates
    p=curCdf(x,df{:},ncp);
    % deviations of p of current nc x pdfs from target
    deltaP=p-pTarget(iIx);
    nIter=nIter+1;
    if doAnimate
      ph=plotPdf(x,ncp,ph,curPdf,df,iIx,nIter,abscissLim,ti);
    end
  end
  % pick border which is closer to the target value
  [nada,ix]=min(abs(deltaP));
  ci(iIx)=ncp(ix);
end
% close figure
if doAnimate
  pause(1)
  close(fh)
end

% ======================== LOCAL FUNCTION =================================
function ph=plotPdf(x,ncp,ph,pdfH,df,iIx,nIter,abscissLim,ti)
% ** function ph=plotPdf(x,ncp,ph,pdfH,df,iIx,nIter,abscissLim,ti)
% If doAnimate==true, plotPdf plots x (first input arg to ncpci) and
% noncentral pdfs with the noncentrality parameter estimates of each
% iteration step 
abscissVal=linspace(abscissLim(iIx,1),abscissLim(iIx,2),200);
if isempty(ph)
  ph(1)=plot(abscissVal,pdfH(abscissVal,df{:},ncp(1)),'-');
  ph(2)=plot(abscissVal,pdfH(abscissVal,df{:},ncp(2)),'-');
  set(ph(1),'color',[.9 .3 .3]);
  set(ph(2),'color',[.3 .3 .9]);
else
  set(ph(1),'xdata',abscissVal,'ydata',pdfH(abscissVal,df{:},ncp(1)));
  set(ph(2),'xdata',abscissVal,'ydata',pdfH(abscissVal,df{:},ncp(2)));  
end
title([ti{iIx} ', iteration # ' int2str(nIter)])
% supposedly, in animation mode we would like to be able to follow the
% iterative process with our eyes, so slow things down
drawnow
pause(.1)


function pvpmod(x)
% PVPMOD             - evaluate parameter/value pairs
% pvpmod(x) assigns the value x(i+1) to the parameter defined by the
% string x(i) in the calling workspace. This is useful to evaluate 
% <varargin> contents in an mfile, e.g. to change default settings 
% of any variable initialized before pvpmod(x) is called.
%
% (c) U. Egert 1998

%############################################
% this loop is assigns the parameter/value pairs in x to the calling
% workspace.

if ~isempty(x)
   for i = 1:2:size(x,2)
      assignin('caller', x{i}, x{i+1});
   end;
end;

%############################################

