function a = Sw(w1, dw, ksw, R1W, R2W, M0s, M0w)

S11 = -((R1W + M0s/M0w*ksw) * sin(atan(w1/dw)) * sin(atan(w1/dw)) + (R2W + M0s/M0w*ksw) * cos(atan(w1/dw)) * cos(atan(w1/dw)));
S12 = -sqrt(w1*w1 + dw*dw);
S13 = -sin(atan(w1/dw)) * cos(atan(w1/dw)) * (R1W + M0s/M0w*ksw - R2W - M0s/M0w*ksw);
S21 = sqrt(w1 * w1 + dw * dw); 
S22 = -R2W - M0s/M0w*ksw; 
S23 = 0;
S31 = -sin(atan(w1/dw)) * cos(atan(w1/dw)) * (R1W + M0s/M0w*ksw - R2W - M0s/M0w*ksw); 
S32 = 0; 
S33 = -((R2W + M0s/M0w*ksw) * (sin(atan(w1/dw)) * sin(atan(w1/dw))) + (R1W + M0s/M0w*ksw) * (cos(atan(w1/dw)) * cos(atan(w1/dw))));

a = [ S11, S12, S13;
      S21, S22, S23;
      S31, S32, S33];