
% Please run this script to reproduce first-level DCM plots in the MRS-DCM
% paper. These plots include (i) correlation plots between the mean expectation
% of inferred synaptic physiology and free energies of inferred DCM of test 
% and re-test data and (ii) predicted versus observed test and re-test data.
% A jafarian aj631@cam.ac.uk
%============= Please add SPM12 path to Matlab=============================
clc
clear all
close all
spm('defaults','eeg')
%=====================Load first level DCMs================================
load('GCM_DCM.mat') % This file is the test and retest odd/even DCMs 

 T     = G{1,1}';
 R     = G{1,2}';
 
 %======plots predicted versus observed test and variance explained=========
for i = 1 :11
    Hz = T{1, 1}.Hz ;
    figure(i)
    DCM = T{i,1};
    subplot(2,1,1),
    plot(Hz,abs(DCM.Hc{1,1}(:,1,1)),'k',Hz,abs(DCM.Hc{1,1}(:,1,1) + DCM.Rc{1,1}(:,1,1)),':r')
    
    %  Varience explained even;
    pss     =sum(sum(sum(abs(spm_vec(DCM.Hc{1,1}(:,1,1))).^2)));
    rss     =sum(sum(sum(abs(spm_vec(DCM.Rc{1,1}(:,1,1))).^2)));
    VT(i,:) = 100*pss/(rss+pss);
    
    
    DCM = R{i,1};
    subplot(2,1,2), 
    plot(Hz,abs(DCM.Hc{1,1}(:,1,1)),'k',Hz,abs(DCM.Hc{1,1}(:,1,1) + DCM.Rc{1,1}(:,1,1)),':r')
 
    
    %  Varience explained odd;

    pss     =sum(sum(sum(abs(spm_vec(DCM.Hc{1,1}(:,1,1))).^2)));
    rss     =sum(sum(sum(abs(spm_vec(DCM.Rc{1,1}(:,1,1))).^2)));
    VR(i,:) = 100*pss/(rss+pss);
    
end

%====Scatter plots of  individual parameters and Free energies=============

 Ininh                 = [1:num_sub]; 
 InEx                  = [1:num_sub];
 for i = 1:11
     inh1(i,1)  =  T{Ininh(i), 1}.Ep.H(1,1);
     inh2(i,1)  =  T{Ininh(i), 1}.Ep.H(2,2);
     inh3(i,1)  =  T{Ininh(i), 1}.Ep.H(3,3);
     inh4(i,1)  =  T{Ininh(i), 1}.Ep.H(4,4);
     inh5(i,1)  =  T{Ininh(i), 1}.Ep.H(1,3);
     inh6(i,1)  =  T{Ininh(i), 1}.Ep.H(2,3);
     inh7(i,1)  =  T{Ininh(i), 1}.Ep.H(4,3);
     
     ex1(i,1)   =  T{InEx(i), 1}.Ep.H(2,1);
     ex2(i,1)   =  T{InEx(i), 1}.Ep.H(3,1);
     ex3(i,1)   =  T{InEx(i), 1}.Ep.H(3,4);
     ex4(i,1)   =  T{InEx(i), 1}.Ep.H(4,2);
     
     inh1(i,2)  =  R{Ininh(i), 1}.Ep.H(1,1);
     inh2(i,2)  =  R{Ininh(i), 1}.Ep.H(2,2);
     inh3(i,2)  =  R{Ininh(i), 1}.Ep.H(3,3);
     inh4(i,2)  =  R{Ininh(i), 1}.Ep.H(4,4);
     inh5(i,2)  =  R{Ininh(i), 1}.Ep.H(1,3);
     inh6(i,2)  =  R{Ininh(i), 1}.Ep.H(2,3);
     inh7(i,2)  =  R{Ininh(i), 1}.Ep.H(4,3);
     
     ex1(i,2)   =  R{InEx(i), 1}.Ep.H(2,1);
     ex2(i,2)   =  R{InEx(i), 1}.Ep.H(3,1);
     ex3(i,2)   =  R{InEx(i), 1}.Ep.H(3,4);
     ex4(i,2)   =  R{InEx(i), 1}.Ep.H(4,2);
 end
for i = 1:11
     F(i,1)    =  T{Ininh(i), 1}.F;
     F(i,2)    =  R{InEx(i), 1}.F;
end

figure
scatter (F(:,1),F(:,2))
refline
r      = corrcoef(F(:,1),F(:,2));
r      = full(r(1,2));
str = sprintf('corr = %-0.2f',r);
title(str,'FontSize',16)
xlabel('Free Energy (odd)','FontSize',12), ylabel('Free Energy (even)','FontSize',12)
title(str,'FontSize',16)
set(gca,'FontSize',12);


figure
scatter (inh1(:,1),inh1(:,2))
refline
r      = corrcoef(inh1(:,1),inh1(:,2));
r      = full(r(1,2));
str = sprintf('corr = %-0.2f',r);
title(str,'FontSize',16)
xlabel('H(1,1) (odd data)','FontSize',12), ylabel('H(1,1)(even data)','FontSize',12)
title(str,'FontSize',16)
set(gca,'FontSize',12);

%============================
figure
scatter (inh2(:,1),inh2(:,2))
refline
r      = corrcoef(inh2(:,1),inh2(:,2));
r      = full(r(1,2));
str = sprintf('corr = %-0.2f',r);
title(str,'FontSize',16)
xlabel('H(2,2) (odd data)','FontSize',12), ylabel('H(2,2)(even data)','FontSize',12)
title(str,'FontSize',16)
set(gca,'FontSize',12);


%============================
figure
scatter (inh3(:,1),inh3(:,2))
refline
r      = corrcoef(inh3(:,1),inh3(:,2));
r      = full(r(1,2));
str = sprintf('corr = %-0.2f',r);
title(str,'FontSize',16)
xlabel('H(3,3) (odd data)','FontSize',12), ylabel('H(3,3)(even data)','FontSize',12)
title(str,'FontSize',16)
set(gca,'FontSize',12);

%=============================
figure
scatter (inh4(:,1),inh4(:,2))
refline
r      = corrcoef(inh4(:,1),inh4(:,2));
r      = full(r(1,2));
str = sprintf('corr = %-0.2f',r);
xlabel('H(4,4) (odd data)','FontSize',12), ylabel('H(4,4)(even data)','FontSize',12)
title(str,'FontSize',16)
set(gca,'FontSize',12);

%==========================
figure
scatter (inh5(:,1),inh5(:,2))
refline
r      = corrcoef(inh5(:,1),inh5(:,2));
r      = full(r(1,2));
str = sprintf('corr = %-0.2f',r);
title(str,'FontSize',16)
xlabel('H(1,3) (odd data)','FontSize',12), ylabel('H(1,3)(even data)','FontSize',12)
title(str,'FontSize',16)
set(gca,'FontSize',12);

%==========================
figure
scatter (inh6(:,1),inh6(:,2))
refline
r      = corrcoef(inh6(:,1),inh6(:,2));
r      = full(r(1,2));
str = sprintf('corr = %-0.2f',r);
title(str,'FontSize',16)
xlabel('H(2,3) (odd data)','FontSize',12), ylabel('H(2,3)(even data)','FontSize',12)
title(str,'FontSize',16)
set(gca,'FontSize',12);

%==========================
figure
scatter (inh7(:,1),inh7(:,2))
refline
r      = corrcoef(inh7(:,1),inh7(:,2));
r      = full(r(1,2));
str = sprintf('corr = %-0.2f',r);
title(str,'FontSize',16)
xlabel('H(4,3) (odd data)','FontSize',12), ylabel('H(4,3)(even data)','FontSize',12)
title(str,'FontSize',16)
set(gca,'FontSize',12);

%=======================
figure
scatter (ex1(:,1),ex1(:,2))
refline
r      = corrcoef(ex1(:,1),ex1(:,2));
r      = full(r(1,2));
str = sprintf('corr = %-0.2f',r);
title(str,'FontSize',16)
xlabel('H(2,1) (odd)','FontSize',12), ylabel('H(2,1) (even)','FontSize',12)
title(str,'FontSize',16)
set(gca,'FontSize',12);

%===========================
figure
scatter (ex2(:,1),ex2(:,2))
refline
r      = corrcoef(ex2(:,1),ex2(:,2));
r      = full(r(1,2));
str = sprintf('corr = %-0.2f',r);
title(str,'FontSize',16)
xlabel('H(3,1) (odd)','FontSize',12), ylabel('H(3,1) (even)','FontSize',12)
title(str,'FontSize',16)
set(gca,'FontSize',12);

%====================
figure
scatter (ex3(:,1),ex3(:,2))
refline
r      = corrcoef(ex3(:,1),ex3(:,2));
r      = full(r(1,2));
str = sprintf('corr = %-0.2f',r);
title(str,'FontSize',16)
xlabel('H(3,4) (odd)','FontSize',12), ylabel('H(3,4) (even)','FontSize',12)
title(str,'FontSize',16)
set(gca,'FontSize',12);

%==========================================================
figure
scatter (ex4(:,1),ex4(:,2))
refline
r      = corrcoef(ex4(:,1),ex4(:,2));
r      = full(r(1,2));
str = sprintf('corr = %-0.2f',r);
title(str,'FontSize',16)
xlabel('H(4,2) (odd)','FontSize',12), ylabel('H(4,2) (even)','FontSize',12)
title(str,'FontSize',16)
set(gca,'FontSize',12);
