function a = B(w1, dw, sep, R1S, R1W, M0s, M0w)

B11 = R1W * M0w * sin( atan(w1/dw) );
B21 = 0;
B31 = R1W * M0w * cos( atan(w1/dw) );
B41 = R1S * M0s * sin( atan(w1/(sep+dw)) );
B51 = 0;
B61 = R1S * M0s * cos( atan(w1/(sep+dw)) );

a = -[B11;
      B21;
      B31;
      B41;
      B51;
      B61];