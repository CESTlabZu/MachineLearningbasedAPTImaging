%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The Multiple-Pool Model Lorentzian Fitting Code for ML NOE(-1.6) imaging
%
% Authors: Leqi Yin, Zhongliang Zu
%
% Correspondance: zhongliang.zu@vumc.org 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

% create frequency offset 
maxf=1000;
step=25;
offset= -maxf:step:maxf;
k_4p7T=[-2000, -1750, -1500, -1250, offset, 1250, 1500,1750,2000]';

% Load data containing Z spectra(clean), R1W and fm.
load("")


% multiple-pool Lorentzian fit.
for i=1:length(Zspectra_matrix)%length of trainingdata

    sig=(1-Zspectra_matrix(:,i));

    R1W_AREX=R1W_cal_matrix(i);
    fm_AREX=fm_cal_matrix(i);

    x =k_4p7T;
    beta0= [0.9, 0/2, 560/2,               0.025, -1400/2, 200/2,        0.01, -800/2, 600/2,        0.001, 600/2, 400/2,          0.02, 1400/2, 1200/2,          0.1, 0/2, 10000/2]; % initial test
    lb=[0.02, -400/2, 120/2,               0, -1600/2, 160/2,             0, -1200/2, 200/2,           0, 400/2, 0/2,                0, 1000/2, 400/2,            0, -1600/2, 4000/2]; % lower bound
    ub=[1, 400/2,4000/2,                    0.2, -1200/2, 1200/2,          0.2,-400/2, 2000/2,         0.2, 800/2, 600/2,            1, 1800/2, 2000/2,         1, 1600/2, 40000/2]; % upper bound




    Delta=[1];
    options=optimset('lsqcurvefit') ;
    options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',5e4*length(x),'MaxIter',2e5) ;

    [beta,resnorm,residual,exitflag,output,lambda,jacobian] = ...
        lsqcurvefit(@matsolv, beta0, x, sig, lb, ub, options, Delta) ;



    % amide
    beta_amide=beta;
    sig_simur_amide=matsolv(beta_amide,x,Delta);
    beta_amide(4)=0;
    sig_simur_ref_amide=matsolv(beta_amide,x,Delta);

    mor_AREX_amide(:,i)=(1./(1-sig_simur_amide)-1./(1-sig_simur_ref_amide))*R1W_AREX*(1+fm_AREX);

    % amine
    beta_amine=beta;
    sig_simur_amine=matsolv(beta_amine,x,Delta);
    beta_amine(7)=0;
    sig_simur_ref_amine=matsolv(beta_amine,x,Delta);

    mor_AREX_amine(:,i)=(1./(1-sig_simur_amine)-1./(1-sig_simur_ref_amine))*R1W_AREX*(1+fm_AREX);

    % NOE3p5
    beta_NOE3p5=beta;
    sig_simur_NOE3p5=matsolv(beta_NOE3p5,x,Delta);
    beta_NOE3p5(13)=0;
    sig_simur_ref_NOE3p5=matsolv(beta_NOE3p5,x,Delta);

    mor_AREX_NOE3p5(:,i)=(1./(1-sig_simur_NOE3p5)-1./(1-sig_simur_ref_NOE3p5))*R1W_AREX*(1+fm_AREX);


    % MT
    beta_MT=beta;
    sig_simur_MT=matsolv(beta_MT,x,Delta);
    beta_MT(16)=0;
    sig_simur_ref_MT=matsolv(beta_MT,x,Delta);

    mor_MT(:,i)=(sig_simur_MT-sig_simur_ref_MT)./(1-(sig_simur_MT-sig_simur_ref_MT))*R1W_AREX;


    % NOE1p6
    beta_NOE1p6=beta;
    sig_simur_NOE1p6=matsolv(beta_NOE1p6,x,Delta);
    beta_NOE1p6(10)=0;
    sig_simur_ref_NOE1p6=matsolv(beta_NOE1p6,x,Delta);

    mor_NOE1p6(:,i)=(1./(1-sig_simur_NOE1p6)-1./(1-sig_simur_ref_NOE1p6))*R1W_AREX*(1+fm_AREX);





    sprintf("----------------------- %d",i)

end

