function a = Ts(w1, dw, ksw, sep)

T11 = cos( atan(w1/(sep+dw)) - atan(w1/dw) );
T12 = 0;
T13 = -sin( atan(w1/(sep+dw)) - atan(w1/dw) );
T21 = 0;
T22 = 1;
T23 = 0;
T31 = sin( atan(w1/(sep+dw)) - atan(w1/dw) );
T32 = 0;
T33 = cos( atan(w1/(sep+dw)) - atan(w1/dw) );

a = [ T11, T12, T13;
      T21, T22, T23;
      T31, T32, T33];
a = ksw * a;