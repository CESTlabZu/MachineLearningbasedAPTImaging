%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Codes for Deep Learning based APT imaging using partially
% synthetic data
%
% Authors: Malvika Viswanathan, Zhongliang Zu
%
% Please contact zhongliang.zu@vumc.org incase you have any doubts with the
% code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Tissue mimicking data; matrix should contain Z spectra(clean), R1W
% and fm.
load('');

% creating random 11 training samples and 1000 testing samples.
randidx = randperm(size(matrix_input_all_noise,2),11);
k = randperm(size(matrix_input_all_noise,2),size(matrix_input_all_noise,2));
for i = 1:1000
   if any(randidx == k(i))
       break
   else
       testidx(i) = k(i);
       groundtruth(:,i) = spectrum(matrix_output1(k(i),1),matrix_output1(k(i),2),-1440);
   end
end

testing_ip = matrix_input_noise(:,testidx);
testing_input = matrix_input_all_noise(:,testidx);
testing_inputclean =  matrix_input_all(:,testidx);

fm = fm_cal_matrix_output_all(randidx,1);
r1 = R1W_cal_matrix_output_all(randidx,1);

% tt mean saturation power (uT)
tt_9p4T=1;
% pulse duration is the total saturation time (s)
pulseduration=5;
maxf=2000;
step=50;
offset= -maxf:step:maxf;
k_9p4T=[-4000, -3500, -3000, -2500, offset, 2500, 3000,3500,4000];
k_9p4T=k_9p4T';
satangle=tt_9p4T*42.6*360*pulseduration;


for i=1:length(matrix_input_all)   
    sig=(1-matrix_input_all(:,i));
    
    R1W_AREX=R1W_cal_matrix_output_all(i);
    fm_AREX=fm_cal_matrix_output_all(i);
    
    x =k_9p4T;
    beta0= [0.9, 0, 560,               0.025, -1400, 200,        0.01, -800, 600,        0.001, 600, 400,          0.02, 1400, 1200,          0.1, 0, 10000]; % initial test
    lb=[0.07, -40, 200,               0, -1500, 200,             0, -920, 200,           0, 520, 200,                0, 1200, 400,            0, 0, 4000]; % lower bound
    ub=[1, 40,1200,                    0.2, -1300, 1200,          0.2,-680, 2000,         0.2, 720, 600,            0.3, 1600, 2000,         0.3, 1200, 40000]; % upper bound
    
    
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
    % Save for comparison with DL prediction
    mor_MTR_amide(:,i)=(sig_simur_amide-sig_simur_ref_amide);
    mor_AREX_amide(:,i)=(1./(1-sig_simur_amide)-1./(1-sig_simur_ref_amide))*R1W_AREX*(1+fm_AREX);


    % amine
    beta_amine=beta;
    sig_simur_amine=matsolv(beta_amine,x,Delta);
    beta_amine(7)=0;
    sig_simur_ref_amine=matsolv(beta_amine,x,Delta);
    % Save AREX fitting for partial synthetic data simulation
    mor_AREX_amine(:,i)=(1./(1-sig_simur_amine)-1./(1-sig_simur_ref_amine))*R1W_AREX*(1+fm_AREX);
    mor_MTR_amine(:,i)=(sig_simur_amine-sig_simur_ref_amine);
    
    % MT
    beta_MT=beta;
    sig_simur_MT=matsolv(beta_MT,x,Delta);
    beta_MT(16)=0;
    sig_simur_ref_MT=matsolv(beta_MT,x,Delta);
    % Save AREX fitting for partial synthetic data simulation
    mor_MT=(sig_simur_MT-sig_simur_ref_MT);
    mor_AREX_MT(:,i) = (mor_MT./(1-mor_MT))*R1W_AREX;
    mor_dir_MT(:,i)=(sig_simur_MT-sig_simur_ref_MT);
    sprintf("----------------------- %d",i)
        
end