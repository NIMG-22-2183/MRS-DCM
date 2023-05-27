close all; 
clear all;
clc;
% This code is useful for greedy search over parameter space of Sigmoid map!
% This kind of search (greedy) is useful when there is no prior inofrmation about
% model/data like emeprical MRS data. 
% Amirhossein Jafarian  aj631@cam.ac.uk
%================ please add SPM12=========================================
spm('defaults','eeg');
%====================load GCM= group DCMs==================================
load('GCM_DCM.mat');
GCM1           = G{1,1}';
GCM2           = G{1,2}';
T1 = {};
T2 = {};

%====search over fine resoultion in parameter space (LB lower bound, UB
%upper bound, Res resolution of parameter changes.
Lp1            = [LB(1) :  Res(1) :  UB(1)]';
Lp2            = [LB(2) :  Res(2) :  UB(2)]'; 
Lp3            = [0]';   % We only consider two parameters for greedy search
L1             = length(Lp1); 
L2             = length(Lp2); 
L3             = length(Lp3); 
 
% Model parameters
%================find all subsets of connections in the model==============
% Exampl: self inhibition 
CMM_NMDA_connections          =  {'H(1,1)', 'H(2,2)','H(3,3)','H(4,4)'} ;
PS                            =  PowerSet(CMM_NMDA_connections);
PR                            =  length(PS); 

%=================Loop over all fields i.e subsets of connections==========
e_even                  = {}  ; 
f_even                  = []  ;
e_odd                   = {}  ; 
f_odd                   = []  ;
diff                    = []  ;
Pf_even                 = []  ;
Pe_even                 = {}  ;
Pf_odd                  = []  ;
Pe_odd                  = {}  ;
Pdiff                   = []  ;

parfor i  = 1:PR 
    for p1= 1:L1
        for p2= 1:L2
            for p3= 1:L3
                field                 =  PS{i,1}(1,:);
                x                     =  [Lp1(p1), Lp2(p2), Lp3(p3)];
                [F1 , Results1 ]      =  cost_function_greedy(x,GCM1,field);
                [F2 , Results2 ]      =  cost_function_greedy(x,GCM2,field);
                f_even(i,p1,p2,p3,:)  =  F1;
                e_even{i,p1,p2,p3}    =  Results1;
                f_odd(i,p1,p2,p3,:)   =  F2;
                e_odd{i,p1,p2,p3}     =  Results2;
                diff(i,p1,p2,p3,:)    =  F1 - F2;
            end
        end
    end 
end

%==========================Find the model with the greatest evidence=======
ind_higher_evidence_even   = find(f_even == max(max(max(max(f_even)))));
[d1,d2,d3,d4]              = ind2sub(size(f_even), ind_higher_evidence_even);                   
ind_even                   = [d1,d2,d3,d4];


%==========================Find the model with the greatest evidence=======
ind_higher_evidence_odd    = find(f_odd == max(max(max(max(f_odd)))));
[d1,d2,d3,d4]              = ind2sub(size(f_odd), ind_higher_evidence_odd);                   
ind_odd                    = [d1,d2,d3,d4];
%==========================================================================
for i = 1:PR
    kke(i,:)               = spm_vec(f_even(i,:)');
    kko(i,:)               = spm_vec(f_odd(i,:)');
end
%==========================================================================
field                      =  PS{ind_even(1),1}(1,:);
x1                         =  [Lp1(ind_even(2)), Lp2(ind_even(3)), Lp3(ind_even(4))];
[F1 , Results1 ]           =  cost_function_greedy(x1,GCM1,field);
field                      =  PS{ind_odd(1),1}(1,:);
x2                         =  [Lp1(ind_odd(2)), Lp2(ind_odd(3)), Lp3(ind_odd(4))];
[F2 , Results2 ]           =  cost_function_greedy(x2,GCM2,field);
 Re                        =  {Results1,Results2,F1,F2,ind_even,ind_odd, Lp1,Lp2};
%=========================================================================
figure(1)
scatter (spm_vec(f_odd),spm_vec(f_even))
refline
r      = corrcoef(spm_vec(f_odd),spm_vec(f_even));
r      = full(r(1,2));
str = sprintf('corr = %-0.2f',r);
title(str,'FontSize',16)
xlabel('Free Energy (odd)','FontSize',12), ylabel('Free Energy (even)','FontSize',12)
title(str,'FontSize',16)
set(gca,'FontSize',12);
axis square
box off
%==========================================================================
function [ P ] = PowerSet( S )
% Modified to exclude emty subset. Orginal code was written by
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
function [F , Results] = cost_function_greedy(x,GCM,field)
% load the vector of GLU or GABA 
load('gaba.mat') 
%=======================Transformed the GUL/GABA based on vector x======
g               = gaba./mean(gaba)-1;  % deviation from mean 
% Define log scale parameterisation  alpha here- before doing the search
TM(:,1)         = 1./(1+exp( alpha(1)*exp(x(1))*(g-alpha(2)*exp(x(2)))))-alpha(3)*exp(x(3)); 

%=========================SPM_DCM_PEB set up===============================
M              = struct();
M.Q            = 'all';
N2             = 11; 
M.X            = [ones(N2,1),(TMRS- mean(TMRS))];
[PEB, P]       = spm_dcm_peb(GCM,M,field);
F              = PEB.F;
Results        = {PEB; P; TMRS};
end

