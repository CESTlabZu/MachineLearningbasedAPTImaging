function vals = pulsesolv1(w1, dw, ksw1,ksw2,ksw3,ksw4, ksw5,kmw, mnots1,mnots2,mnots3,mnots4,mnots5, mnotw, mnotm, R1S, R2S1,R2S2, R2S3,R2S4,  R2S5, R1W, R2W, R1M, R2M, sep1,sep2,sep3, sep4,sep5, init, duration,excitewait)
    fs1=mnots1/mnotw;
    fs2=mnots2/mnotw;
    fs3=mnots3/mnotw;
    fs4=mnots4/mnotw;
    fs5=mnots5/mnotw;
    fm=mnotm/mnotw;
    
    delter=2.3*400*2*pi;
    %ls=1./(pi*R2M)./(1+((dw-delter)./R2M).^2);
    %W=w1^2*sqrt(pi/2)/R2M*exp(-((dw-delter)/R2M)^2/2);


    % Lorentzian
    ls=1./(pi*R2M)./(1+((dw-delter)./R2M).^2);
    w=(w1)^2*pi*ls;  


   A=[-R2S1-ksw1,    -sep1-dw,      0,            fs1*ksw1,                     0,                           0,                                                             0,           0,           0,               0,            0,                  0,                   0,            0,          0,           0,           0,               0,           0,              0;                       
      sep1+dw,       -R2S1-ksw1,    w1,           0,                            fs1*ksw1,                    0,                                                             0,           0,           0,               0,            0,                  0,                   0,            0,          0,           0,           0,               0,           0,              0;
      0,             -w1,           -R1S-ksw1,    0,                            0,                           fs1*ksw1,                                                      0,           0,           0,               0,            0,                  0,                   0,            0,          0,           0,           0,               0,           0,              R1S*fs1;
      ksw1,           0,            0,            -R2W-fs1*ksw1-fs2*ksw2-fm*kmw-fs3*ksw3-fs4*ksw4-fs5*ksw5     ,  -dw,                         0,                           0,           ksw2,        0,               0,            ksw3,               0,                   0,            ksw4,       0,           0,           ksw5,            0,           0,              0;
      0,              ksw1,         0,            dw,                           -R2W-fs1*ksw1-fs2*ksw2-fm*kmw-fs3*ksw3-fs4*ksw4-fs5*ksw5,      w1,                          0,           0,          ksw2,             0,            0,                  ksw3,                0,            0,          ksw4,        0,           0,               ksw5,        0,              0;
      0,              0,            ksw1,         0,                            -w1,                          -R1W-fs1*ksw1-fs2*ksw2-fm*kmw-fs3*ksw3-fs4*ksw4-fs5*ksw5,     kmw,         0,           0,               ksw2,         0,                  0,                   ksw3,         0,          0,           ksw4,        0,               0,           ksw5,           R1W;
      0,              0,            0,            0,                            0,                           fm*kmw,                                                        -R1M-kmw-w   0,           0,               0,            0,                  0,                   0,            0,          0,           0,           0,               0,           0,              R1M*fm;
      0,              0,            0,            fs2*ksw2,                     0,                           0,                                                             0,           -R2S2-ksw2,  -sep2-dw,        0,            0,                  0,                   0,            0,          0,           0,           0,               0,           0,              0;           
      0,              0,            0,            0,                            fs2*ksw2,                    0,                                                             0,           sep2+dw,     -R2S2-ksw2,      w1,           0,                  0,                   0,            0,          0,           0,           0,               0,           0,              0;
      0,              0,            0,            0,                            0,                           fs2*ksw2,                                                      0,           0,           -w1,             -R1S-ksw2,    0,                  0,                   0,            0,          0,           0,           0,               0,           0,              R1S*fs2;
      0,              0,            0,            fs3*ksw3,                     0,                           0,                                                             0,           0,           0,               0,            -R2S3-ksw3,         -sep3-dw,            0,            0,          0,           0,           0,               0,           0,              0;
      0,              0,            0,            0,                            fs3*ksw3,                    0,                                                             0,           0,           0,               0,            sep3+dw,            -R2S3-ksw3,          w1,           0,          0,           0,           0,               0,           0,              0;
      0,              0,            0,            0,                            0,                           fs3*ksw3,                                                      0,           0,           0,               0,            0,                  -w1,                 -R1S-ksw3,    0,          0,           0,           0,               0,           0,              R1S*fs3;
      0,              0,            0,            fs4*ksw4,                     0,                           0,                                                             0,           0,           0,               0,            0,                  0,                   0,            -R2S4-ksw4, -sep4-dw,    0,           0,               0,           0,              0;
      0,              0,            0,            0,                            fs4*ksw4,                    0,                                                             0,           0,           0,               0,            0,                  0,                   0,            sep4+dw,    -R2S4-ksw4,  w1,          0,               0,           0,              0;
      0,              0,            0,            0,                            0,                           fs4*ksw4,                                                      0,           0,           0,               0,            0,                  0,                   0,            0,           -w1,        -R1S-ksw4,   0,               0,           0,              R1S*fs4;
      0,              0,            0,            fs5*ksw5,                     0,                           0,                                                             0,           0,           0,               0,            0,                  0,                   0,            0,            0,          0,          -R2S5-ksw5,      -sep5-dw,    0,              0;
      0,              0,            0,            0,                            fs5*ksw5,                    0,                                                             0,           0,           0,               0,            0,                  0,                   0,            0,            0,          0,          sep5+dw,         -R2S5-ksw5,  w1,             0;
      0,              0,            0,            0,                            0,                           fs5*ksw5,                                                      0,           0,           0,               0,            0,                  0,                   0,            0,            0,          0,          0,               -w1,         -R1S-ksw5,    R1S*fs5;
      0,              0,            0,            0,                            0,                           0,                                                             0,           0,           0,               0,            0,                  0,                   0,            0,            0,          0,          0,               0,           0,              0;];
      

    [V,D] = eig(A); % Calculate Eigenvectors and Eigenvalues
    ppp = real(V*diag(exp(diag(D*duration)))/V)*init; % Calculate Matrix Exponential
    vals=ppp;
end

