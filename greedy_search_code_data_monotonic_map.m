
% This code is useful for greedy search constraint by using monotonic map only!
% This kind of search (greedy) is useful when there is no prior inofrmation about
% model/data like emeprical MRS prior. 
% Amirhossein Jafarian aj631@cam.ac.uk

%================ please add SPM12=========================================
close all; 
clear all;
clc;
spm('defaults','eeg');
%====================load GCM= group DCMs==================================
load('GCM_DCM.mat');
GCM1           = G{1,1}';
GCM2           = G{1,2}';
T1             = {};
T2             = {};
%====search over fine resoultion in parameter space (LB lower bound, UB
%upper bound, Res resolution of parameter changes.
Lp1            = [LB(1) :  Res(1) :  UB(1)]';
Lp2            = [LB(2) :  Res(2) :  UB(2)]'; 
Lp3            = [0]';   % We only consider two parameters for greedy search
L1             = length(Lp1); 
L2             = length(Lp2); 
L3             = length(Lp3); 

%=Inclusion map and replace non-monotonic with trivial/original emeprical prior===  
%=======================Transformed the GUL or GABA based on vector x======
load('glu.mat') %  This can be GABA,GLU. 
g             = glu./mean(glu)-1; % define the relative deviation from mean
g             = sort(g);
y             = [];
T             = [];
Ty            = []; 
ind_fail      = [];

for p1= 1:L1
    for p2= 1:L2
        for p3= 1:L3
            x                        =  [Lp1(p1), Lp2(p2), Lp3(p3)];
            y(p1,p2,p3,:)            =   x(1)*g.^2 +x(2)*g+x(3);
            r                        = ismonotonic(y); 
            if r  == 1
               T                     = glu./mean(glu)-1; 
               Ty(p1,p2,p3,:)        = x(1)*T.^2 +x(2)*T+x(3);
               ind_fail(p1,p2,p3,:)  = 1;  
            else
               T                     = glu./mean(glu)-1; 
               Ty(p1,p2,p3,:)        = T;
               ind_fail(p1,p2,p3,:)  = 0;
            end   
        end
    end
end

%================find all subsets of connections in the model==============
% Exampl Exc connections
CMM_NMDA_connections                  = {'H(3,1)', 'H(4,2)','H(3,4)'};
PS                                    =  PowerSet(CMM_NMDA_connections);
PR                                    =  length(PS); 

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

for i  = 1:PR 
    for p1= 1:L1
        for p2= 1:L2
            for p3= 1:L3
                     field                = PS{i,1}(1,:);
                    [F1 , Results1 ]      =  cost_function_greedy_m(x,GCM1,field,squeeze(Ty(p1,p2,p3,:)));
                    [F2 , Results2 ]      =  cost_function_greedy_m(x,GCM2,field,squeeze(Ty(p1,p2,p3,:)));
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

ind_higher_evidence_even  = find(f_even == max(max(max(max(f_even)))));
[d1,d2,d3,d4]              = ind2sub(size(f_even), ind_higher_evidence_even);                   
ind_even                   = [d1,d2,d3,d4];


%==========================Find the model with the greatest evidence=======

ind_higher_evidence_odd  = find(f_odd == max(max(max(max(f_odd)))));
[d1,d2,d3,d4]              = ind2sub(size(f_odd), ind_higher_evidence_odd);                   
ind_odd                    = [d1,d2,d3,d4];

%==========================================================================

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
%===============================functions =================================
function [ P ] = PowerSet( S )
% Modified to exclude emty subset. Orginal code was written by
% Paulo Abelha (2023). PowerSet( S )
% (https://www.mathworks.com/matlabcentral/fileexchange/62752-powerset-s),
    n       = numel(S);
    x       = 1:n;
    P       = cell(1,2^n);
    p_ix    = 2;
    
    for nn = 1:n
        a = combnk(x,nn);
        for j=1:size(a,1)
            P{p_ix} = S(a(j,:));
            p_ix    = p_ix + 1;
        end
    end
    P      = P';
    P      = P(2:end);
end

function monotonic = ismonotonic(x, strict, direction, dim)
% ISMONOTONIC(X) returns a boolean value indicating whether or not a vector is monotonic.  
% By default, ISMONOTONIC returns true for non-strictly monotonic vectors,
% and both monotonic increasing and monotonic decreasing vectors. For
% matrices and N-D arrays, ISMONOTONIC returns a value for each column in
% X.
% 
% ISMONOTONIC(X, 1) works as above, but only returns true when X is
% strictly monotonically increasing, or strictly monotonically decreasing.
% 
% ISMONOTONIC(X, 0) works as ISMONOTONIC(X).
% 
% ISMONOTONIC(X, [], 'INCREASING') works as above, but returns true only
% when X is monotonically increasing.
% 
% ISMONOTONIC(X, [], 'DECREASING') works as above, but returns true only
% when X is monotonically decreasing.
% 
% ISMONOTONIC(X, [], 'EITHER') works as ISMONOTONIC(X, []).
% 
% ISMONOTONIC(X, [], [], DIM) works as above, but along dimension DIM.
% 
% NOTE: Third input variable is case insensitive, and partial matching is
% used, so 'd' would be recognised as 'DECREASING' etc..
% 
% EXAMPLE:
%     x = [1:4; 6:-2:2 3]
%     ismonotonic(x)
%     ismonotonic(x, [], 'i')
%     ismonotonic(x, [], [], 2)
% 
%     x =
%          1     2     3     4
%          6     4     2     3
%     ans = 
%          1     1     1     1
%     ans =
%          1     1     0     0
%     ans =
%          1
%          0
% 
% SEE ALSO: is*
% 
% $ Author: Richie Cotton $     $ Date: 2010/01/20 $    $ Version: 1.2 $
% Basic error checking & default setup
if ~isreal(x) || ~isnumeric(x)
   warning('ismonotonic:badXValue', ...
       'The array to be tested is not real and numeric.  Unexpected behaviour may occur.'); 
end
if nargin < 2 || isempty(strict)
    strict = false;
end
if nargin < 3 || isempty(direction)
    direction = 'either';
end
% Accept partial matching for direction
lendir = length(direction);
if strncmpi(direction, 'increasing', lendir)
    testIncreasing = true;
    testDecreasing = false;
elseif strncmpi(direction, 'decreasing', lendir)
    testIncreasing = false;
    testDecreasing = true;
elseif strncmpi(direction, 'either', lendir) 
    testIncreasing = true;
    testDecreasing = true;
else
    warning('ismonotonic:badDirection', ...
        'The string entered for direction has not been recognised, reverting to ''either''.');
    testIncreasing = true;
    testDecreasing = true;
end
if nargin < 4 || isempty(dim)
    dim = find(size(x) ~= 1, 1);
    if isempty(dim)
        dim = 1;
    end
end
% Test for monotonic increasing
if testIncreasing    
    if strict
        comparison = @gt;
    else
        comparison = @ge;
    end    
    monotonicAscending = all(comparison(diff(x, [], dim), 0), dim);
else
    monotonicAscending = false;
end
% Test for monotonic decreasing
if testDecreasing    
    if strict
        fhComparison = @lt;
    else
        fhComparison = @le;
    end
    monotonicDescending = all(fhComparison(diff(x, [], dim), 0), dim);
else
    monotonicDescending = false;
end
monotonic = monotonicAscending | monotonicDescending;
end

function [F , Results] = cost_function_greedy_m(x,GCM,field,TMRS)
M              = struct();
M.Q            = 'all';
N2             = 11; 
M.X            = [ones(N2,1),(TMRS)];
[PEB, P]       = spm_dcm_peb(GCM,M,field);
F              = PEB.F;
Results        = {PEB; P; TMRS};
end

       


