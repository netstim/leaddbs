function [params dispParams] = GTdefaults

% regularizers
params.conpot_lambda = 0.5;
params.fiberSmooth = 0.1;
params.voxelSmooth_vf = 0.5; 
params.voxelSmooth_Ds = 0.1; 
params.voxelSmooth_Ds_range = 20;
params.penalty_SW = 0*0.01;
params.curv_hardthres = acos(45/180*pi);
params.trackingguide_strength = 50;

% tracking guide diffusion model
%params.fixedTissueMode = [1 0 1 -1 1];  % [ D_para D_perp vfi -1 weightscale]
params.fixedTissueMode = [2 0.5 0.5 -1 1];  % [ D_para D_perp vfi -1 weightscale]

% iteration parameters
params.Tstart = 3;
params.Tend = 0.001;
params.numstep = 100;
params.numiterations = -4;

% distribution of proposal types
params.prop_p_birth = 0.4;
params.prop_p_death = 0.1;		
params.prop_p_shift = 2;
params.prop_p_Dmod = 0.8;     
params.prop_p_vfmod = 0.4; 
params.prop_p_conprob = 0.2;


% particle related parameters
params.p_weight = 1.5; %1.2;                 % maximum weight
params.p_len = 2;                       % length of segment
params.p_wid = 2;                       % oversampling factor
params.p_chempot = 0.01;                   % cost of a particle (L0-regularizer)
params.p_expnump = 0.000002;                   % lambda of Poisson process

% connection related parameters
params.c_likli = 0.5;                   % likeliniess that a connection occurs
params.c_kappa = 0.5;                   % trade-off between angular and distance penalty (dist 0...1 ang)
params.c_bound_attract = 0.5;            % wm/gm-boundary attracts open fiber terminals

% diffusion model
params.maxbval = 4000;                  % all data beyond this bvalue is neglected
params.b_weighting = 0.05;              % relative weighting of b0-images
params.restrictions = 0*2^0 +0*2^1;        % if 0x01 is set D_para_ex == D_orth_ex
                                        % if 0x10 is set D_para_ex == D_para_int
params.trace_equality_constraint = 0.01;

% approximation model parameters
params.alpha = 0;                       % alpha==1: sw implicitly modeled, alpha==0: no sw
params.ordermax = [6 2 4 2];            % [order_MM kappa_MM order_MS kappa_MS]

% length thresholds of extracted fibers 
params.fibrange = [10 inf];

% still experimental
params.directional_propsal_distrib = 0;

% sphere discretization
params.sphericalDiscNumber = 128; % possible values 32,48,64,128,256
        % note that, if you change this you have recompile (or just delete /tmp/mesoFT/sinterp*)
        
        
        
        
        
% number of cores used
params.numcores = 1;  % if ==1, openMP free compile

% just a  string of additonal compiler parameters
params.compilerparams = '-fopenmp' ;
params.libs = '-lgomp' ;
        

% parameters shown in the GUI
cnt = 1;
dispParams(cnt).tag = 'Tstart';
dispParams(cnt).name = 'starting Temp';
cnt = cnt +1;
dispParams(cnt).tag = 'Tend';
dispParams(cnt).name = 'stopping Temp';
cnt = cnt +1;
dispParams(cnt).tag = 'numiterations';
dispParams(cnt).name = '#its';
cnt = cnt +1;
dispParams(cnt).tag = 'numstep';
dispParams(cnt).name = 'number of ir';
cnt = cnt +1;
dispParams(cnt).tag = 'p_len';
dispParams(cnt).name = 'Segment length';
cnt = cnt +1;
dispParams(cnt).tag = 'p_wid';
dispParams(cnt).name = 'oversampling';
cnt = cnt +1;
dispParams(cnt).tag = 'p_chempot';
dispParams(cnt).name = 'Segment Cost';
cnt = cnt +1;
dispParams(cnt).tag = 'c_likli';
dispParams(cnt).name = 'Con. Potential';
cnt = cnt +1;
dispParams(cnt).tag = 'c_kappa';
dispParams(cnt).name = 'Curv/Dist balance';
cnt = cnt +1;
dispParams(cnt).tag = 'numcores';
dispParams(cnt).name = '#Threads';
cnt = cnt +1;
