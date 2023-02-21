function a = Ss(w1, dw, ksw, sep, R1S, R2S)

S11 = -((R1S + ksw) * sin(atan(w1/(sep+dw))) * sin(atan(w1/(sep+dw))) + (R2S + ksw) * cos(atan(w1/(sep+dw))) * cos(atan(w1/(sep+dw))));
S12 = -sqrt(w1*w1 + (sep+dw)*(sep+dw));
S13 = -sin(atan(w1/(sep+dw))) * cos(atan(w1/(sep+dw))) * (R1S + ksw - R2S - ksw);
S21 = sqrt(w1 * w1 + (sep+dw)*(sep+dw)); 
S22 = -R2S - ksw; 
S23 = 0;
S31 = -sin(atan(w1/(sep+dw))) * cos(atan(w1/(sep+dw))) * (R1S + ksw - R2S - ksw); 
S32 = 0; 
S33 = -((R2S + ksw) * (sin(atan(w1/(sep+dw))) * sin(atan(w1/(sep+dw)))) + (R1S + ksw) * (cos(atan(w1/(sep+dw))) * cos(atan(w1/(sep+dw)))));

a = [ S11, S12, S13;
      S21, S22, S23;
      S31, S32, S33];