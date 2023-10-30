%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Codes for Machine Learning based APT imaging using partially
% synthetic data
%
% Authors: Malvika Viswanathan, Zhongliang Zu
%
% Please contact zhongliang.zu@vumc.org incase you have any doubts with the
% code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
% Declaring global variables

rat_power=3 % 1uT
B0_sythetic=9.4

ttt=[0.25 0.5 1];
tt= ttt(rat_power);
delB0=B0_sythetic/9.4;  %3/9.4
delw1=tt/1;
offppm=400; %127
gauss=100;

betal_W=(9.4/B0_sythetic)^(0.4);
gamma_W=1;
betal=1;
gamma=1;

pulseduration=5;
i_SNR=5.0;

exciteduration=0.002;
exciteangle=90;
satangle=tt*42.6*360*pulseduration;

TR=0.02;
max=2000;
step=50;

sep1_9p4T=3.6*offppm;
sep2_9p4T=3*offppm;
sep3_9p4T=2*offppm
sep4_9p4T=-1.6*offppm;
sep5_9p4T=-3.3*offppm;

fs2=0.003;
fs3=0.0003;
fs4=0.003;

ksw2 =5000;
ksw3=500;
ksw4=50;
ksw5=20;
kmw=25;

R1S=1/1.5;
R2S1=1/0.002;
R2S2=1/0.01;
R2S3=1/0.01;
R2S4=1/0.001;
R2S5=1/0.0005;
R1M=1/1.5;
R2M=1/0.00005;


offset= -max:step:max;
k=[-4000, -3500, -3000, -2500, offset, 2500, 3000,3500,4000]*delB0;

%% Creating Partial Synthetic Data
%
% Please save variable matrix_MTR_output1/matrix_MTR_output for training
% target and matrix_input/noise for training inputs.

% Load fitted parameters for tissue mimicking data or fitted invivo data 
% in the line below.
load('');
idx = randidx(1);
fm=fm(1); % from tissue mimicking data

R1W=r1(1); % from tissue mimicking data

num_T1W=5;
num_T2W=5;
num_fs1=5;
num_fs2=5;
num_fs5=5;
num_fm=5;
num_ksw1=5;
num_T2S=4;

T1W_matrix=[1.6, 1.8, 2.0, 2.2, 2.4];
T2W_matrix=[40, 60, 80, 100, 120]*0.001;
T2S_matrix=[0.002, 0.003, 0.004, 0.005];


fs2_matrix=[0.5, 0.75, 1.0, 1.25, 1.5];
fm_matrix=[0.4, 0.7, 1.0, 1.3, 1.6];

i=1;




for ii_T1W=1:num_T1W
    ii_T1W
    R1W_cal=1./T1W_matrix(ii_T1W);
  for ii_T2W=1:num_T2W
      ii_T2W
      R2W_cal=1./T2W_matrix(ii_T2W);
      for ii_T2S=1:num_T2S
          R2S_cal = 1./T2S_matrix(ii_T2S);
              for ii_fs1=1:num_fs1
                fs1_cal=0.0006+0.0002*(ii_fs1-1);
                for ii_fs2=1:num_fs2
                    fs2_cal=fs2_matrix(ii_fs2);
                       for ii_fs5=1:num_fs5
                           fs5_cal=0.004+0.003*(ii_fs5-1);
                           for ii_fm=1:num_fm
                               fm_cal=fm_matrix(ii_fm);
                               for ii_ksw1=1:num_ksw1
                                    ksw1_cal=40+30*(ii_ksw1-1);
                                 
                             


 R1W_cal_obs=(R1W_cal+fm_cal*fm*R1M)./(1+fm_cal*fm);
 cal_Lorentzian1_cal=(fs1_cal.*ksw1_cal.*(tt.*42.6*2*pi).^2./((tt.*42.6*2*pi).^2+(R2S_cal+ksw1_cal)*ksw1_cal+ksw1_cal./(R2S_cal+ksw1_cal).*((k+sep1_9p4T)*2*pi).^2));
 cal_Lorentzian2_normal_cal = fs2_cal*mean(mor_AREX_amine(:,idx),2)';
 cal_Lorentzian5_cal=(fs5_cal.*ksw5.*(tt.*42.6*2*pi).^2./((tt.*42.6*2*pi).^2+(R2S5+ksw5)*ksw5+ksw5./(R2S5+ksw5).*((k+sep5_9p4T)*2*pi).^2));
 cal_Lorentzian6_normal_cal = fm_cal*mean(mor_AREX_MT(:,idx),2)';
 cal_eff_cal=R1W_cal_obs.*((k)*2*pi).^2./((tt.*42.6*2*pi).^2+((k)*2*pi).^2)+R2W_cal.*(tt.*42.6*2*pi).^2./((tt.*42.6*2*pi).^2+((k)*2*pi).^2);

%  cal_Lorentzian2_tumor_cal = fs2_cal*mean(Zspectra_mor_amine_tumor,2)';
%  cal_Lorentzian6_tumor_cal = fm_cal*mean(Zspectra_mor_MT_tumor,2)';


 sscal_normal = R1W_cal_obs./(cal_eff_cal+cal_Lorentzian1_cal./(1+fm_cal*fm)+cal_Lorentzian2_normal_cal./(1+fm_cal*fm)+cal_Lorentzian5_cal./(1+fm_cal*fm)+cal_Lorentzian6_normal_cal).*(((k)*2*pi).^2./((tt.*42.6*2*pi).^2+((k)*2*pi).^2));
 SS_cal_normal(:, ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1)=R1W_cal_obs./(cal_eff_cal+cal_Lorentzian1_cal./(1+fm_cal*fm)+cal_Lorentzian2_normal_cal./(1+fm_cal*fm)+cal_Lorentzian5_cal./(1+fm_cal*fm)+cal_Lorentzian6_normal_cal).*(((k)*2*pi).^2./((tt.*42.6*2*pi).^2+((k)*2*pi).^2));
 SS_cal_value_ref_normal=R1W_cal_obs./(cal_eff_cal+0./(1+fm_cal*fm)+cal_Lorentzian2_normal_cal./(1+fm_cal*fm)+cal_Lorentzian5_cal./(1+fm_cal*fm)+cal_Lorentzian6_normal_cal).*(((k)*2*pi).^2./((tt.*42.6*2*pi).^2+((k)*2*pi).^2));
 
%  sscal_tumor = R1W_cal_obs./(cal_eff_cal+cal_Lorentzian1_cal./(1+fm_cal*fm)+cal_Lorentzian2_tumor_cal./(1+fm_cal*fm)+cal_Lorentzian5_cal./(1+fm_cal*fm)+cal_Lorentzian6_tumor_cal).*(((k)*2*pi).^2./((tt.*42.6*2*pi).^2+((k)*2*pi).^2));
%  SS_cal_tumor(:, ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1)=R1W_cal_obs./(cal_eff_cal+cal_Lorentzian1_cal./(1+fm_cal*fm)+cal_Lorentzian2_tumor_cal./(1+fm_cal*fm)+cal_Lorentzian5_cal./(1+fm_cal*fm)+cal_Lorentzian6_tumor_cal).*(((k)*2*pi).^2./((tt.*42.6*2*pi).^2+((k)*2*pi).^2));
%  SS_cal_value_ref_tumor=R1W_cal_obs./(cal_eff_cal+0./(1+fm_cal*fm)+cal_Lorentzian2_tumor_cal./(1+fm_cal*fm)+cal_Lorentzian5_cal./(1+fm_cal*fm)+cal_Lorentzian6_tumor_cal).*(((k)*2*pi).^2./((tt.*42.6*2*pi).^2+((k)*2*pi).^2));

 R1W_cal_matrix(ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1)=R1W_cal_obs;
 fm_cal_matrix(ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1)=fm_cal*fm;
 params_matrix(:, ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1) = [R1W_cal_obs; R2W_cal;R2S_cal; fs1_cal; fs2_cal; fs5_cal; fm_cal; ksw1_cal];
 
 % Direct subtraction of label and reference signals
 mtr_normal =(1-sscal_normal) - (1-SS_cal_value_ref_normal);
 cestr_normal(:, ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1) =(1-sscal_normal) - (1-SS_cal_value_ref_normal);
 mtr_amp_normal(1,ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1) = mtr_normal(1,16);
 mtr_width_normal(1,ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1)= fwhm(mtr_normal,k);
%  mtr_tumor = (1-sscal_tumor)-(1- SS_cal_value_ref_tumor);
%  cestr_tumor(:, ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1)=(1-sscal_tumor)-(1- SS_cal_value_ref_tumor);
%  mtr_amp_tumor(1,ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1) = mtr_tumor(1,16);
%  mtr_width_tumor(1,ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1)= fwhm(mtr_tumor,k);

 i = i+1;

                            
                               end
                           end
                       end
                end
          end
      end
  end
end

k_cut=[1:25 41:49 85:89];

matrix_input(:,:)=reshape(SS_cal_normal(k_cut,:,:,:,:,:,:,:),  [length(k_cut') num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1]);   
matrix_input_all(:,:)=reshape(SS_cal_normal,  [length(k') num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1]); 
% matrix_input2(:,:)=reshape(SS_cal_tumor(k_cut,:,:,:,:,:,:,:),  [length(k_cut') num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1]);   
% matrix_input_all2(:,:)=reshape(SS_cal_tumor,  [length(k') num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1]);
% matrix_input = [matrix_input1, matrix_input2];


sz = size(matrix_input_all(:,1));
matrix_input_all_noise = matrix_input_all + (0.002+0.002*(i_SNR-1))*(randn(sz));
sz = size(matrix_input(:,1));
matrix_input_noise = matrix_input + (0.002+0.002*(i_SNR-1))*(randn(sz));


R1W_cal_matrix_output_all=reshape(R1W_cal_matrix,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1 1]);
fm_cal_matrix_output_all=reshape(fm_cal_matrix,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1 1]);
matrix_MTR_output1(:,1)=reshape(mtr_amp_normal,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1 1]);   
matrix_MTR_output1(:,2)=reshape(mtr_width_normal,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1 1]);
% matrix_MTR_output2(:,1)=reshape(mtr_amp_tumor,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1 1]);   
% matrix_MTR_output2(:,2)=reshape(mtr_width_tumor,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1 1]);
% matrix_MTR_output(:,1)=[matrix_MTR_output1(:,1); matrix_MTR_output2(:,1)]; 
% matrix_MTR_output(:,2)=[matrix_MTR_output1(:,2); matrix_MTR_output2(:,2)];
matrix_CESTR_output1 = reshape(cestr_normal, [89 num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1]);
%matrix_CESTR_output2 = reshape(cestr_tumor, [89 num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1]);

