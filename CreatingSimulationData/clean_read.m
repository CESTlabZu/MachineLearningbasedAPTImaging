clear;

tic
rat_power=3
B0_sythetic=9.4

ttt=[0.25 0.5 1];
tt= ttt(rat_power);
delB0=B0_sythetic/9.4;  %3/9.4
delw1=tt/1;
offppm=400; %127

betal_W=(9.4/B0_sythetic)^(0.4);
gamma_W=1;
betal=1;
gamma=1;


load('matlab_amine')
fs1_array(1,:)=[0.9726    1.6976    1.0851    1.2638    1.4219    1.8203    1.2727    1.1059]*10^(-3);
fs1_array(2,:)=[1.2047    1.0391    0.8226    1.1658    0.9194    1.1223    1.1429    0.8765]*10^(-3);

ksw1_array(1,:)=[81.4153   51.8419   55.5867   43.9763   49.0173   85.0796   88.1235  197.8926];
ksw1_array(2,:)=[63.9803   93.5020   75.5704   67.4344   87.1886  305.7937  117.1014   67.9304];

fs2_array(1,:)=[8.8134    1.1773    1.2529    3.5029    3.6070    1.0256    0.8471    0.7456]*10^(-3);
fs2_array(2,:)=[14.7529    4.6015    5.8083    2.9135    6.4230    3.3807    5.2381   11.2357]*10^(-3);

ksw2_array(1,:)=[1.0000    0.2631    0.3601    0.9326    0.9543    0.2518    0.2363    0.2152]*10^(4);
ksw2_array(2,:)=[1.0000    0.9698    1.0000    0.7914    1.0000    0.7264    1.0000    1.0000]*10^(4);


T2s1_array(1,:)=[1.5051    1.1567    1.2421    1.6341    1.6854    0.8884    0.8144    0.9684]*10^(-3);
T2s1_array(2,:)=[1.5915    1.5372    1.8168    1.5945    2.0255    1.2382    1.2197    2.9079]*10^(-3);


fs5_array(1,:)=[10.6057    7.8691    9.1458    9.1001    9.1345   10.6533   10.5534    8.7336]*10^(-3);
fs5_array(2,:)=[20.4925   15.3622   13.3206   13.5449   12.8491   13.6289   17.5842   15.2465]*10^(-3);

ksw5_array(1,:)=[17.0711   17.5712   15.0185   15.6552   14.3134   11.9433   10.2896   12.5560];
ksw5_array(2,:)=[12.4536   18.4662   23.6659   16.1594   20.7696   25.0921   19.2268   23.6315];


T2s5_array(1,:)=[0.3943    0.5151    0.5085    0.4469    0.5101    0.5297    0.6021    0.6674]*10^(-3);
T2s5_array(2,:)=[0.4356    0.4806    0.4025    0.5114    0.5413    0.4144    0.4126    0.3881]*10^(-3);


value_R1obsfmap_array(1,:)=[0.4490    0.3932    0.4658    0.4378    0.4302    0.4427    0.4525    0.4773];
value_R1obsfmap_array(2,:)=[0.5621    0.5270    0.5577    0.5292    0.5206    0.5467    0.5578    0.5604];


value_R1fmap_array(1,:)=[0.4245    0.3619    0.4375    0.4137    0.4044    0.4163    0.4273    0.4516];
value_R1fmap_array(2,:)=[0.5221    0.4791    0.5089    0.4914    0.4822    0.5029    0.5084    0.5174];


value_pmfmap_array(1,:)=[ 0.0473    0.0540    0.0554    0.0441    0.0471    0.0487    0.0477    0.0509];
value_pmfmap_array(2,:)=[ 0.0946    0.1050    0.1150    0.0835    0.0829    0.1002    0.1166    0.1012];

value_kmfmap_array(1,:)=[  19.5676   21.8582   17.2545   30.1638   29.6267   31.1840   34.9704   32.2001];
value_kmfmap_array(2,:)=[ 14.1759   14.9835   14.2442   18.7022   16.1133   15.1452   14.1699   12.6700];

% value_R2_array(1)=28;
% value_R2_array(2)=24;
value_R2_array(1)=1/0.06;
value_R2_array(2)=1/0.06;







pulseduration=5;
i_SNR=5.0;

exciteduration=0.002;
exciteangle=90;

TR=0.02;


for iiii=1:1:256
gauss(iiii)=100;
end

rat_type=2;

max=2000;
step=50;
sep1=3.6*offppm;
sep2=3*offppm;
sep3=2*offppm
sep4=-1.6*offppm;
sep5=-3.3*offppm;

fs1=mean(fs1_array(rat_type, :),2);
% fs1=0.0015
fs2=mean(fs2_array(rat_type, :),2);
% fs2=0.003
fs3=0.0003;
%fs3=0.0005;k
fs4=0.003;
%fs4=0.0025
fs5=mean(fs5_array(rat_type, :),2);
%fs5=0.0065
fm=mean(value_pmfmap_array(rat_type, :),2);
%fm=0.1

ksw1=mean(ksw1_array(rat_type, :),2);
ksw2=mean(ksw2_array(rat_type, :),2);
ksw3=500;
ksw4=50;
ksw5=mean(ksw5_array(rat_type, :),2);
kmw=mean(value_kmfmap_array(rat_type, :),2);

R1S=1/1.5*betal;
R2S1=mean(1./T2s1_array(rat_type, :)*gamma,2);
R2S5=mean(1./T2s5_array(rat_type, :)*gamma,2);

R2S2=1/0.01;
R1W_obs=mean(value_R1obsfmap_array(rat_type, :)*betal_W,2);
R1W=mean(value_R1fmap_array(rat_type, :)*betal_W,2);
R2W=value_R2_array(rat_type)*gamma_W;
R1M=1/1.5*betal;

%k= -max:step:max;
offset= -max:step:max;
k=[-4000, -3500, -3000, -2500, offset, 2500, 3000,3500,4000]*delB0;






 cal_Lorentzian1=(fs1.*ksw1.*(tt.*42.6*2*pi).^2./((tt.*42.6*2*pi).^2+(R2S1+ksw1)*ksw1+ksw1./(R2S1+ksw1).*((k+sep1)*2*pi).^2));
 cal_Lorentzian2=mean((Zspectra_mor_amine_normal(:,1,rat_power,:))*delw1^2./((3*400*2*pi)^2./((1*42.6*2*pi).^2+(k'*2*pi).^2)).*((3*offppm*2*pi)^2./((delw1*42.6*2*pi).^2+(k'*2*pi).^2)),4)';
 cal_Lorentzian5=(fs5.*ksw5.*(tt.*42.6*2*pi).^2./((tt.*42.6*2*pi).^2+(R2S5+ksw5)*ksw5+ksw5./(R2S5+ksw5).*((k+sep5)*2*pi).^2));
 cal_Lorentzian6=mean((Zspectra_mor_MT_normal(:,1,3,:))*delw1^2,4)';
 cal_eff=R1W_obs.*((k)*2*pi).^2./((tt.*42.6*2*pi).^2+((k)*2*pi).^2)+R2W.*(tt.*42.6*2*pi).^2./((tt.*42.6*2*pi).^2+((k)*2*pi).^2);
 

 SS=R1W_obs./(cal_eff+cal_Lorentzian1./(1+fm)+cal_Lorentzian2./(1+fm)+cal_Lorentzian5./(1+fm)+cal_Lorentzian6).*(((k)*2*pi).^2./((tt.*42.6*2*pi).^2+((k)*2*pi).^2));

 
 for ii=1:89
    fZspectra_normal(ii,:,:,:)=Zspectra_normal(90-ii,:,:,:);
 end

 figure (1)
 hold on
plot(SS)
plot(squeeze(mean(fZspectra_normal(:,rat_power,:),3)))




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
 cal_Lorentzian1_cal=(fs1_cal.*ksw1_cal.*(tt.*42.6*2*pi).^2./((tt.*42.6*2*pi).^2+(R2S_cal+ksw1_cal)*ksw1_cal+ksw1_cal./(R2S_cal+ksw1_cal).*((k+sep1)*2*pi).^2));
 cal_Lorentzian2_normal_cal=fs2_cal*mean((Zspectra_mor_amine_normal(:,1,rat_power,:))*delw1^2./((3*400*2*pi)^2./((1*42.6*2*pi).^2+(k'*2*pi).^2)).*((3*offppm*2*pi)^2./((delw1*42.6*2*pi).^2+(k'*2*pi).^2)),4)';
 cal_Lorentzian5_cal=(fs5_cal.*ksw5.*(tt.*42.6*2*pi).^2./((tt.*42.6*2*pi).^2+(R2S5+ksw5)*ksw5+ksw5./(R2S5+ksw5).*((k+sep5)*2*pi).^2));
 cal_Lorentzian6_normal_cal=fm_cal*mean((Zspectra_mor_MT_normal(:,1,3,:))*delw1^2,4)';
 cal_eff_cal=R1W_cal_obs.*((k)*2*pi).^2./((tt.*42.6*2*pi).^2+((k)*2*pi).^2)+R2W_cal.*(tt.*42.6*2*pi).^2./((tt.*42.6*2*pi).^2+((k)*2*pi).^2);

 cal_Lorentzian2_tumor_cal=fs2_cal*mean((Zspectra_mor_amine_tumor(:,1,rat_power,:))*delw1^2./((3*400*2*pi)^2./((1*42.6*2*pi).^2+(k'*2*pi).^2)).*((3*offppm*2*pi)^2./((delw1*42.6*2*pi).^2+(k'*2*pi).^2)),4)';
 cal_Lorentzian6_tumor_cal=fm_cal*mean((Zspectra_mor_MT_tumor(:,1,3,:))*delw1^2,4)';

 sscal_normal = R1W_cal_obs./(cal_eff_cal+cal_Lorentzian1_cal./(1+fm_cal*fm)+cal_Lorentzian2_normal_cal./(1+fm_cal*fm)+cal_Lorentzian5_cal./(1+fm_cal*fm)+cal_Lorentzian6_normal_cal).*(((k)*2*pi).^2./((tt.*42.6*2*pi).^2+((k)*2*pi).^2));
 SS_cal_normal(:, ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1)=R1W_cal_obs./(cal_eff_cal+cal_Lorentzian1_cal./(1+fm_cal*fm)+cal_Lorentzian2_normal_cal./(1+fm_cal*fm)+cal_Lorentzian5_cal./(1+fm_cal*fm)+cal_Lorentzian6_normal_cal).*(((k)*2*pi).^2./((tt.*42.6*2*pi).^2+((k)*2*pi).^2));
 SS_cal_value_ref_normal=R1W_cal_obs./(cal_eff_cal+0./(1+fm_cal*fm)+cal_Lorentzian2_normal_cal./(1+fm_cal*fm)+cal_Lorentzian5_cal./(1+fm_cal*fm)+cal_Lorentzian6_normal_cal).*(((k)*2*pi).^2./((tt.*42.6*2*pi).^2+((k)*2*pi).^2));
 
 sscal_tumor = R1W_cal_obs./(cal_eff_cal+cal_Lorentzian1_cal./(1+fm_cal*fm)+cal_Lorentzian2_tumor_cal./(1+fm_cal*fm)+cal_Lorentzian5_cal./(1+fm_cal*fm)+cal_Lorentzian6_tumor_cal).*(((k)*2*pi).^2./((tt.*42.6*2*pi).^2+((k)*2*pi).^2));
 SS_cal_tumor(:, ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1)=R1W_cal_obs./(cal_eff_cal+cal_Lorentzian1_cal./(1+fm_cal*fm)+cal_Lorentzian2_tumor_cal./(1+fm_cal*fm)+cal_Lorentzian5_cal./(1+fm_cal*fm)+cal_Lorentzian6_tumor_cal).*(((k)*2*pi).^2./((tt.*42.6*2*pi).^2+((k)*2*pi).^2));
 SS_cal_value_ref_tumor=R1W_cal_obs./(cal_eff_cal+0./(1+fm_cal*fm)+cal_Lorentzian2_tumor_cal./(1+fm_cal*fm)+cal_Lorentzian5_cal./(1+fm_cal*fm)+cal_Lorentzian6_tumor_cal).*(((k)*2*pi).^2./((tt.*42.6*2*pi).^2+((k)*2*pi).^2));

 R1W_cal_matrix(ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1)=R1W_cal_obs;
 fm_cal_matrix(ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1)=fm_cal*fm;
 params_matrix(:, ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1) = [R1W_cal_obs; R2W_cal;R2S_cal; fs1_cal; fs2_cal; fs5_cal; fm_cal; ksw1_cal];
 mtr_normal =(1-sscal_normal) - (1-SS_cal_value_ref_normal);
 mtr_amp_normal(1,ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1) = mtr_normal(1,16);
 mtr_width_normal(1,ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1)= fwhm(mtr_normal,k);
 mtr_tumor = (1-sscal_tumor)-(1- SS_cal_value_ref_tumor);
 mtr_amp_tumor(1,ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1) = mtr_tumor(1,16);
 mtr_width_tumor(1,ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1)= fwhm(mtr_tumor,k);

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
matrix_input1(:,:)=reshape(SS_cal_normal(k_cut,:,:,:,:,:,:,:),  [length(k_cut') num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1]);   
matrix_input_all1(:,:)=reshape(SS_cal_normal,  [length(k') num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1]); 
matrix_input2(:,:)=reshape(SS_cal_tumor(k_cut,:,:,:,:,:,:,:),  [length(k_cut') num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1]);   
matrix_input_all2(:,:)=reshape(SS_cal_tumor,  [length(k') num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1]);
matrix_input = [matrix_input1, matrix_input2];
matrix_input_all = [matrix_input_all1, matrix_input_all2];

sz = size(matrix_input_all(:,1));
matrix_input_all_noise = matrix_input_all + (0.002+0.002*(i_SNR-1))*(rand(sz));
sz = size(matrix_input(:,1));
matrix_input_noise = matrix_input + (0.002+0.002*(i_SNR-1))*(rand(sz));


R1W_cal_matrix_output_all=reshape(R1W_cal_matrix,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1 1]);
fm_cal_matrix_output_all=reshape(fm_cal_matrix,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1 1]);
matrix_ip = [matrix_input; [R1W_cal_matrix_output_all.', R1W_cal_matrix_output_all.']];
matrix_ip_noise = [matrix_input_noise; [R1W_cal_matrix_output_all.', R1W_cal_matrix_output_all.']];

% Create Target Matrix
matrix_MTR_output1(:,1)=reshape(mtr_amp_normal,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1 1]);   
matrix_MTR_output1(:,2)=reshape(mtr_width_normal,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1 1]);
matrix_MTR_output2(:,1)=reshape(mtr_amp_tumor,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1 1]);   
matrix_MTR_output2(:,2)=reshape(mtr_width_tumor,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1 1]);
matrix_MTR_output(:,1)=[matrix_MTR_output1(:,1); matrix_MTR_output2(:,1)]; 
matrix_MTR_output(:,2)=[matrix_MTR_output1(:,2); matrix_MTR_output2(:,2)];

R1W_cal_matrix_output_all=[R1W_cal_matrix_output_all.', R1W_cal_matrix_output_all.'];
fm_cal_matrix_output_all=[fm_cal_matrix_output_all.', fm_cal_matrix_output_all.'];


R1W_cal_matrix_output_all2=[R1W_cal_matrix_output_all.', R1W_cal_matrix_output_all.'];
fm_cal_matrix_output_all2=[fm_cal_matrix_output_all.', fm_cal_matrix_output_all.'];

filename ='train_data';
save(filename,'matrix_input','matrix_input_noise','matrix_MTR_output');
