% This code is useful for search over small search space to reduce the
% sensitivity of optmisation  to inital guess for parameters and avoid
% local minia.
% you can perform greedy search first. Then around the winning model (i.e., greed boundary UB and LB )
% optimisation search can be used! You can alternatively use spm_argmax !
% Amirhossein Jafarian  aj631@cam.ac.uk
clear all;
clc;
%================ please add SPM12=========================================
spm('defaults','eeg');
%====================load GCM= group DCMs==================================
% load Group DCM as a column cell array of dim(GCM) ={Nx1}
load('GCM_DCM.mat');
GCM           = G;

%================find all subsets of connections in the model==============
% Example; self-inh of populations
EF_CMM_NMDA                          =  {'H(1,1)', 'H(2,2)','H(3,3)','H(4,4)'} ;
PS                                   = PowerSet(excitatory_CMM_NMDA);

%=================Loop over all fields i.e subsets of connections==========

evalute_free_energy                  = {}; 

for i  = 1:length(PS) % avoid using parfor !!
    try,
        field                    = PS{i,1}(1,:);
        fun                      = @(x)cost_function(x,GCM,field);
        x                        = fmincon(fun,theta,[],[],[],[],lb,ub); 
        evalute_free_energy{i,1} = feval(@eval_fun,x,GCM,field);
        evalute_free_energy{i,2} = x; 
    end
end

%==========================Find the model with the greatest evidence=======
st=[];
for i = 1:size(PS,1)
    st(i,:) = evalute_free_energy{i,1}.F;
end
ind_higher_evidence_model  = find(st == max(st));
%==========================================================================


GCM           = G;

evalute_free_energy                  = {};

for i  = 1:length(PS) % avoid using parfor !!
    try,
        field                    = PS{i,1}(1,:);
        fun                      = @(x)cost_function(x,GCM,field);
        x                        = fmincon(fun,theta,[],[],[],[],lb,ub); 
        evalute_free_energy{i,1} = feval(@eval_fun,x,GCM,field);
        evalute_free_energy{i,2} = x; .
    end
end

%==========================Find the model with the greatest evidence=======
st=[];
for i = 1:size(PS,1)
    st(i,:) = evalute_free_energy{i,1}.F;
end
ind_higher_evidence_model  = find(st == max(st));
%==========================================================================

function [ P ] = PowerSet( S )
% Modified to exclude empty subset. Orginal code was written by
% Paulo Abelha (2023). PowerSet( S )
% (https://www.mathworks.com/matlabcentral/fileexchange/62752-powerset-s),
    n      =  numel(S);
    x      =  1:n;
    P      =  cell(1,2^n);
    p_ix   =  2;
    for nn = 1:n

        a = combnk(x,nn);

        for j=1:size(a,1)

            P{p_ix} = S(a(j,:));
            p_ix    = p_ix + 1;

        end

    end
    P           = P';
    P           = P(2:end);
end
%==========================================================================
function F = cost_function(x,GCM,field)
%This is cost function measure the - model eivdence 
% Input vector x that define the unknown parameters of MRS/PET map
% field is a subset of connections in DCM
% Output
% -Free energy output of spm_dmc_peb
%==========================================================================
% load the vector of GLU or GABA here for  subjects
load('gaba.mat')
%=======================Transformed the GUL or GABA based on vector x======

g           = gaba./mean(gaba)-1; % optional: relative changes wrt population 
% Define log scale parameterisation  alpha here
TM(:,1)     = 1./(1+exp( alpha(1)*exp(x(1))*(g-alpha(2)*exp(x(2)))))-alpha(3)*exp(x(3)); 

%=========================SPM_DCM_PEB set up===============================

M              = struct();
M.Q            = 'all';
N2             = 12-1; % total number in your GCM
M.X            = [ones(N2,1),(TMRS- mean(TMRS))];
PEB            = spm_dcm_peb(GCM,M,field);
F              = -PEB.F; % negative of F because we are minizing the cost function using fmin
end
%==========================================================================
function D = eval_fun(x,GCM,field)
%This is cost function measure the - model eivdence 
% Input vector x that define the unknown parameters of MRS/PET map
% field is a subset of connections in DCM
% Output
% Free energy output of spm_dmc_peb
%==========================================================================
% load the vector of GLU or GABA 
load('gaba.mat')

%=====================Transformed the GUL or GABA based on vector=========

g           = gaba./mean(gaba)-1; % optional: relative changes wrt population 
% Define log scale parameterisation  alpha here
TM(:,1)     = 1./(1+exp( alpha(1)*exp(x(1))*(g-alpha(2)*exp(x(2)))))-alpha(3)*exp(x(3)); 

%=========================SPM_DCM_PEB set up===============================

M              = struct();
M.Q            = 'all';
N2             = 11; % total number in your GCM
M.X            = [ones(N2,1),(TM- mean(TM))]; 
PEB            = spm_dcm_peb(GCM,M,field); % you can use PEB = spm_dcm_peb_bmc(PEB) as well 
% but since we performed over all permustations we addressed this steps. In
% other words all possible local minia have been checked 
D.F            = PEB.F; % this is to evaluate the results so no need
D.field        = field ; 
end