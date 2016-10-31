% getsvm_nottingham_copy
%
%      usage: svm = getsvm_nottingham_copy(x1,x2,kernelfun,kernalargs,C)
%         by: justin gardner
%       date: 08/03/05
%    purpose: construct a support vector machine
%             to classify which stimulus responses 
%             came from.
%
%             based on Pattern Classification                                   
%             Duda, Hart & Stork p. Chapter 5                                       
%                    - and -
%             A Tutorial on Support Vector Machines for Pattern Recognition
%             Christopher J.C. Burges
%             Data Mining and Knowledge Discovery 2, 121-167 (1998)
%                   
%                    - and -
%            An introduction to kernel-Based learning algorithms
%            Muller, Mika Ratsch Tsuda and Scholkopf 
%            IEEE Transactions on neural networks 12.2 (2001)
%
%      e.g.: svm = getsvm_nottingham_copy(x1,x2);
%            where x1 and x2 hold arrays of exemplars from the wo
%            groups to be classified. Each row is a different
%            exmplar, and each column holds the value of the
%            exemplar along another dimension.
%            The default method is to use a linear classifier
%  
%            radial basis functions with C=specified penalty for errors
%            svm = getsvm_nottingham_copy(x1,x2,'radialbasis',0.5,10);
%            2nd order polynomial
%            svm = getsvm_nottingham_copy(x1,x2,'polynomial',2,10);
%            
%            To classify new points with a previously constructed
%            classifier do
%            classification = getsvm_nottingham_copy(xnew,svm);
%
%            try svmtest to play with this function
%            
function svm = getsvm_nottingham_copy(x1,x2,kernelfun,kernelargs,C)

if (nargin == 2)
  % if the second argument is a svm structure then
  % we just want to classify the points in x1
  if (isstruct(x2))
    for i = 1:size(x1,1)
      svm(i) = thissvmclassify(x1(i,:),x2);
    end
    return
  end
  % otherwise these are the defaults
  kernelfun = 'linear';
  kernelargs = [];
  C = 1;
elseif (nargin == 3)
  kernelargs = [];
  C = 1;
elseif (nargin == 4)
  C = 1;
elseif (nargin ~= 5)
  help getsvm_nottingham_copy;
  return
end

% save parameters
svm.kernelfun = kernelfun;
svm.kernelargs = kernelargs;
svm.C = C;

% Fisher's linear discriminant
if (strcmp(lower(kernelfun),'fisher'))
  % do fisher linear discriminant and return
  svm = fisherLinearDiscriminant(x1,x2,svm);
  return
end

% Fisher's linear discriminant
if (strcmp(lower(kernelfun),'fisherpillow'))
  % do fisher linear discriminant and return
  svm = fisherLinearDiscriminantPillow(x1,x2,svm);
  return
end

% get number of points
n1 = size(x1,1);n2 = size(x2,1);n = n1+n2;

% construct data matrix
X = [x1 ; x2];

% construct classification array
Y = [ones(n1,1) ; -ones(n2,1)];

% construct quadratic matrix
H = zeros(n,n);
for i = 1:n
  for j = 1:n
    H(i,j) = eval(sprintf('Y(i)*Y(j)*%s(X(i,:),X(j,:),kernelargs)',kernelfun));
  end
end

% construct linear component
f = -ones(n,1);

% construct equality constraint 
Aeq = Y';
beq = 0;

% get inequality constraints
A = -eye(n);
b = zeros(n,1);

% now get solution using optimization toolbox
% quadratic programming solver
[alpha svm.quadprog.fval svm.quadprog.exitflag svm.quadprog.output svm.quadprog.lambda] = quadprog(H,f,[],[],Aeq,beq,zeros(1,n),ones(1,n)*C,[],optimset('Display','off'));

% round to zero if very small
alpha(alpha < .001*max(alpha)) = 0;

% round to C if close to C
alpha(abs(alpha-C) < min(C/1000,.0001)) = C;

% find real support vectors
supportvectors = (alpha>0);

% get support vectors
svm.sv = X(find(supportvectors),:);

% just keep alpha > 0, weight by category variable
svm.alpha = Y(supportvectors).*alpha(supportvectors);

% calculate number of support vectors
svm.n = length(svm.alpha);

% get weight vectors
svm.w = svm.alpha'*svm.sv;

% calculate threshold, Scholkopf p. 186
% also see Burges, end of page 10
for i = 1:length(alpha)
  innersum = 0;
  for j = 1:length(alpha)
    innersum = innersum + Y(j)*alpha(j)*eval(sprintf('%s(X(i,:),X(j,:),kernelargs)',kernelfun));
  end
  b(i) = Y(i) - innersum;
end

% only use ones that are not 0 or not at C value
if ~isempty(find((alpha>0)&(alpha~=C)))
  svm.b = mean(b((alpha>0)&(alpha~=C)));
else
  svm.b = mean(b((alpha>0)));
end

% save parameters
svm.kernelfun = kernelfun;
svm.kernelargs = kernelargs;
svm.C = C;

%%%%%%%%%%%%%%%%%%%%
% linear kernel
%%%%%%%%%%%%%%%%%%%%
function y = linear(x1,x2,p);

y = dot(x1,x2);


%%%%%%%%%%%%%%%%%%%%
% polynomial kernel
%%%%%%%%%%%%%%%%%%%%
function y = polynomial(x1,x2,p);

y = (dot(x1,x2)+1)^p;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% radial basis function kernel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = radialbasis(x1,x2,p);

% default value
if (isempty(p)),p = 1;end

y = exp(-sum(((x1-x2).^2)/p^2));

%%%%%%%%%%%%%%%%%%%%
% sigmoidal kernel
%%%%%%%%%%%%%%%%%%%%
function y = sigmoid(x1,x2,p);

% default value
if (length(p) == 1),p(2) = 1;end

y = tanh(p(1)*dot(x1,x2)+p(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% classify a new point with svm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = thissvmclassify(x,svm)

% fisherLinearDiscriminant
if (strcmp(lower(svm.kernelfun),'fisher') || strcmp(lower(svm.kernelfun),'fisherpillow'))
  g = dot(svm.w,x) + svm.bias;
  return
end

g = svm.b;
for i = 1:svm.n
  g = g+eval(sprintf('svm.alpha(i)*%s(x,svm.sv(i,:),svm.kernelargs)',svm.kernelfun));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get Fisher linear discriminant
% see page 242 of Duda & Hart
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function svm = fisherLinearDiscriminantOld(x1,x2,svm)

% get number of points
n1 = size(x1,1);n2 = size(x2,1);n = n1+n2;
% construct data matrix
Y = [ones(n1,1) x1 ; -ones(n2,1) -x2];
% construct margin
b = [(n/n1)*ones(n1,1) ; (n/n2)*ones(n2,1)];
% get least squares solution
svm.w = pinv(Y)*b;
% get bias
svm.bias = svm.w(1);
svm.w = svm.w(2:length(svm.w));
% set support vector stuff to default values
svm.n = 0;
svm.alpha = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get Fisher linear discriminant
% uses standard formulation see
% p 120 Duda & Hart
% also added a "ridge coefficient"
% see program compFLD for an
% explanation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function svm = fisherLinearDiscriminant(x1,x2,svm)

% project onto basis constructed for range of x1 and x2
% this is a good thing to do if we have more dimensions
% than repeats
project = 0;
if project
  basis = orth([x1;x2]');
  x1old = x1;x2old = x2;
  x1 = x1*basis;
  x2 = x2*basis;
end
% get ridge coefficent
ridge = std([x1(:) ; x2(:)])/2;
% get number of points
n1 = size(x1,1);n2 = size(x2,1);n = n1+n2;
% get means
m1 = mean(x1);m2 = mean(x2);
% compute scatter matrix
S1 = cov(x1);S2=cov(x2);
% compute total scatter matrix
Sw = S1+S2+diag(ones(size(x1,2),1)*ridge);
% now compute weights
%svm.w = diag(1./diag(Sw))*(m1-m2)';svm.inversionType = 'diag';
svm.w = pinv(Sw)*(m1-m2)';svm.inversionType = 'pinv';

% rotate back if projected on to range basis
if project
  svm.w = basis*svm.w;x1 = x1old;x2 = x2old;
end
% compute the dot product of each
% training set against the weight vector
projection1 = sum(x1.*(svm.w*ones(1,size(x1,1)))',2);
projection2 = sum(x2.*(svm.w*ones(1,size(x2,1)))',2);
% now put the bias point in betweeen the two distributions
% since this is ideal if the two distributions have equal
% standard deviation
svm.bias = -(mean(projection1)+(mean(projection2)-mean(projection1))/2);
svm.n = 0;
svm.alpha = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get Fisher linear discriminant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function svm = fisherLinearDiscriminantPillow(x1,x2,svm)

% setting the ridge regression to auto does 
% not seem to produce good results on these data
% i think because in compFLD, the dimension
% reduction happens before the computation of the
% auto ridge regression parameter thus giving it
% a much higher value than what it should be.
%svm.w = compFLD(x1,x2,'auto',1);
% so we just pass in the ridge regression parameter
%svm.w = compFLD(x1,x2,std([x1(:) ; x2(:)])/2,1);
svm.w = compFLD(x1,x2,0,0);
% note that the direction of the weight vector that
% compFLD spits back is not well defined, so sometimes
% the classification goes in the opposite direction.
% this should probably be fixed if this classifier
% is to be used seriously.
svm.bias = 0;
svm.n = 0;
svm.alpha = [];

% compute the dot product of each
% training set against the weight vector
projection1 = sum(x1.*(svm.w*ones(1,size(x1,1)))',2);
projection2 = sum(x2.*(svm.w*ones(1,size(x2,1)))',2);
% note that the direction of the weight vector that
% compFLD spits back is not well defined, so sometimes
% the classification goes in the opposite direction.
% this fixes the direction to always be such that
% the projection of x1 is positive and x2 is negative
if (sum(projection1 > 0)+sum(projection2 < 0)) < (length(x1)+length(x2))/2
%  disp(sprintf('Flipping weights'));
  svm.w = -svm.w;
end
% now put the bias point in betweeen the two distributions
% since this is ideal if the two distributions have equal
% standard deviation
svm.bias = -(mean(projection1)+(mean(projection2)-mean(projection1))/2);
