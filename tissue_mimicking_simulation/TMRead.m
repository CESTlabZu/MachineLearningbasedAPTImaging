clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Codes for Machine Learning based APT imaging using partially
% synthetic data
%
% Authors: Malvika Viswanathan, Zhongliang Zu
%
% Please contact zhongliang.zu@vumc.org incase you have any doubts with the
% code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tt mean saturation power (uT)
tt_9p4T=1;
% pulse duration is the total saturation time (s)
pulseduration=5;
i_SNR=5;  % add noise

gauss=100;

% create frequency offset 
maxf=2000;
step=50;

% frequency offset of each pool
sep1_9p4T=3.6*400;  %amide
sep2_9p4T=3*400; %fast amine
sep3_9p4T=2*400;  % creatine amine
sep4_9p4T=-1.6*400; % NOE at -1.6
sep5_9p4T=-3.3*400; % NOE at -1.6



% relaxations
R1S=1/1.5;
R2S1=1/0.002;
R2S2=1/0.01;
R2S3=1/0.01;
R2S4=1/0.001;
R2S5=1/0.0005;
R1M=1/1.5;
R2M=1/0.00005;

offset= -maxf:step:maxf;
k_9p4T=[-4000, -3500, -3000, -2500, offset, 2500, 3000,3500,4000];
k_9p4T=k_9p4T';
satangle=tt_9p4T*42.6*360*pulseduration;


fs3=0.0003;
fs4=0.003;
fs5=0.01;

ksw2 =5000;
ksw3=500;
ksw4=50;
ksw5=20;
kmw=25;


num_T1W=3;
num_T2W=3;
num_T2S=3;
num_T2M=3;
num_fs1=3;
num_fs2=3;
num_fs3=3;
num_fs5=3;
num_fm=3;
num_ksw1=3;
num_ksw2=3;
num_ksw3=3;

T1W_matrix=[1.5, 1.9, 2.3];
T2W_matrix=[30, 70, 110]*0.001;
T2S_matrix=[0.001, 0.0025, 0.004];
T2M_matrix = [30, 50, 70]*0.000001;

fs2_matrix=[0.5 1.0 1.5]*0.003;
fs3_matrix = [0.0001, 0.0003, 0.0005];

fm_matrix=[0.04 0.08 0.12];
ksw2_matrix = [3000, 5000, 7000];
ksw3_matrix = [300, 500, 700];



i=1;

for ii_T1W=1:num_T1W
    ii_T1W
    R1W=1./T1W_matrix(ii_T1W);
    for ii_T2W=1:num_T2W
      ii_T2W
      R2W=1./T2W_matrix(ii_T2W);
      for ii_T2S=1:num_T2S
          R2S_cal= 1./T2S_matrix(ii_T2S);
          for ii_T2M = 1:num_T2M
              R2M = 1./T2M_matrix(ii_T2M);
              for ii_fs1=1:num_fs1
                  fs1=0.0005+0.0004*(ii_fs1-1);
                  for ii_fs2=1:num_fs2
                      fs2=fs2_matrix(ii_fs2);
                      for ii_fs3=1:num_fs3
                          fs3=fs3_matrix(ii_fs3);
                          for ii_fs5=1:num_fs5
                              fs5=0.002+0.006*(ii_fs5-1);
                              for ii_fm=1:num_fm
                                  fm=fm_matrix(ii_fm);
                                  for ii_ksw1=1:num_ksw1
                                      ksw1=20+60*(ii_ksw1-1);
                                      for ii_ksw2=1:num_ksw2
                                          ksw2=ksw2_matrix(ii_ksw2);
                                          for ii_ksw3 = 1:num_ksw3
                                              ksw3=ksw3_matrix(ii_ksw3);
                                 
a25mspulse = runsteadysimgauss(ksw1, ksw2, ksw3, ksw4, ksw5, kmw, fs1, fs2, fs3,fs4, fs5, 1, fm, R1S, R2S_cal, R2S2, R2S3,R2S4, R2S5, R1W, R2W, R1M, R2M,sep1_9p4T*2*pi,sep2_9p4T*2*pi,sep3_9p4T*2*pi,sep4_9p4T*2*pi, sep5_9p4T*2*pi, pulseduration, gauss, satangle, 1, 2, 1, .00, 1, 1, k_9p4T*2*pi, 1);
a25mspulse_ref = runsteadysimgauss(ksw1,ksw2, ksw3, ksw4, ksw5, kmw, 0, fs2, fs3,fs4, fs5, 1, fm, R1S, R2S_cal, R2S2, R2S3,R2S4, R2S5, R1W, R2W, R1M, R2M,sep1_9p4T*2*pi,sep2_9p4T*2*pi,sep3_9p4T*2*pi,sep4_9p4T*2*pi, sep5_9p4T*2*pi, pulseduration, gauss, satangle, 1, 2, 1, .00, 1, 1, k_9p4T*2*pi, 1);
disp("step1")
aa_9p4T(:,ii_T1W,ii_T2W,ii_T2S,ii_T2M, ii_fs1, ii_fs2, ii_fs3,ii_fs5, ii_fm, ii_ksw1,ii_ksw2,ii_ksw3)=a25mspulse(:,6);
aa_9p4T_ref(:,ii_T1W,ii_T2W,ii_T2S,ii_T2M, ii_fs1, ii_fs2, ii_fs3,ii_fs5, ii_fm, ii_ksw1,ii_ksw2,ii_ksw3)=a25mspulse_ref(:,6);
disp("step2")
disp("Zspec")
Slab = a25mspulse(:,6);
Sref = a25mspulse_ref(:,6);
S0 = 1;
R1W_cal_obs=(R1W+(fm*R1M))./(1+fm); 
R1W_cal_matrix(ii_T1W,ii_T2W,ii_T2S,ii_T2M, ii_fs1, ii_fs2, ii_fs3,ii_fs5, ii_fm, ii_ksw1,ii_ksw2,ii_ksw3)=R1W_cal_obs;
fm_cal_matrix(ii_T1W,ii_T2W,ii_T2S,ii_T2M, ii_fs1, ii_fs2, ii_fs3,ii_fs5, ii_fm, ii_ksw1,ii_ksw2,ii_ksw3)=fm;
mtr = reshape(((1-Slab)-(1-Sref)), [89 1]);
params_matrix(:, ii_T1W,ii_T2W,ii_T2S,ii_T2M, ii_fs1, ii_fs2,ii_fs5, ii_fs3, ii_fm, ii_ksw1,ii_ksw2,ii_ksw3) = [R1W_cal_obs; R2W;R2S_cal;R2M; fs1; fs2; fs3;fs5; fm; ksw1; ksw2; ksw3];
mtr_amp(ii_T1W,ii_T2W,ii_T2S,ii_T2M, ii_fs1, ii_fs2, ii_fs3,ii_fs5, ii_fm, ii_ksw1,ii_ksw2,ii_ksw3) = mtr(16,:);
mtr_width(ii_T1W,ii_T2W,ii_T2S,ii_T2M, ii_fs1, ii_fs2, ii_fs3,ii_fs5, ii_fm, ii_ksw1,ii_ksw2,ii_ksw3)= fwhm(mtr,k_9p4T);
arex_amp(ii_T1W,ii_T2W,ii_T2S,ii_T2M, ii_fs1, ii_fs2, ii_fs3,ii_fs5, ii_fm, ii_ksw1,ii_ksw2,ii_ksw3) = ((1./Slab(16)) - (1./Sref(16))).*R1W_cal_obs*(1+fm);
arex_width(ii_T1W,ii_T2W,ii_T2S,ii_T2M, ii_fs1, ii_fs2, ii_fs3,ii_fs5, ii_fm, ii_ksw1,ii_ksw2,ii_ksw3)= fwhm((((1./Slab(1:30)) - (1./Sref(1:30))).*R1W_cal_obs*(1+fm)),k_9p4T);
X = sprintf("-------------------------------------%d",i);
disp(X)
i = i+1;
% 
                                            end
                                        end
                                    end
                               end
                          end
                        end
                    end
                end
           end
       end
  end
end


% Save outputs matrix_input_all, R1W_cal_matrix_output_all,
% fm_cal_matrix_output_all for lorentzian fitting

k_cut=[1:25 41:49 85:89];

matrix_input(:,:)=reshape(aa_9p4T(k_cut,:,:,:,:,:,:,:),  [length(k_cut.') num_T1W*num_T2W*num_T2S*num_T2M*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_ksw2*num_ksw3]);   
sz = size(matrix_input(:,1));
matrix_input_noise = matrix_input + (0.002+0.002*(i_SNR-1))*(randn(sz));
matrix_output1(:,1)=reshape(mtr_amp,  [num_T1W*num_T2W*num_T2S*num_T2M*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_ksw2*num_ksw3 1]);     
matrix_output1(:,2)=reshape(mtr_width,  [num_T1W*num_T2W*num_T2S*num_T2M*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_ksw2*num_ksw3 1]);


matrix_output2(:,1)=reshape(arex_amp,  [num_T1W*num_T2W*num_T2S*num_T2M*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_ksw2*num_ksw3 1]);     
matrix_output2(:,2)=reshape(arex_width,  [num_T1W*num_T2W*num_T2S*num_T2M*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_ksw2*num_ksw3 1]);


matrix_input_all(:,:)=reshape(aa_9p4T,  [length(k_9p4T) num_T1W*num_T2W*num_T2S*num_T2M*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_ksw2*num_ksw3]);   
sz = size(matrix_input_all(:,1));
matrix_input_all_noise = matrix_input_all + (0.002+0.002*(i_SNR-1))*(randn(sz));
toc

R1W_cal_matrix_output_all=reshape(R1W_cal_matrix,  [num_T1W*num_T2W*num_T2S*num_T2M*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_ksw2*num_ksw3 1]);
fm_cal_matrix_output_all=reshape(fm_cal_matrix,  [num_T1W*num_T2W*num_T2S*num_T2M*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_ksw2*num_ksw3 1]);
