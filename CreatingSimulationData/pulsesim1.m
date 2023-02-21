function  p = pulsesim1(dw, ksw1,ksw2,ksw3, ksw4,ksw5, kmw, mnots1,mnots2, mnots3,mnots4,mnots5, mnotw, mnotm, R1S, R2S,R2S5, R1W, R2W, R1M, R2M, sep1,sep2,sep3,sep4,sep5, duration, curve, angle, init, TR, exciteflag, excitewait, exciteduration, exciteangle)



for t=1:1:256
    w1=getsatpulse(curve, angle, t, duration, TR);
  %  w1=getsatpulse(curve, exciteangle, t, exciteduration, TR);
   [p,blah] = pulsesolv1(w1, dw, ksw1, ksw2,ksw3,ksw4,ksw5, kmw, mnots1,mnots2,mnots3,mnots4,mnots5, mnotw, mnotm, R1S, R2S,R2S5, R1W, R2W, R1M, R2M, sep1,sep2,sep3,sep4, sep5, init, duration/256);

   index = size(blah); 
   init = blah(index(1), :);   
end


p = blah(index(1), :);
end