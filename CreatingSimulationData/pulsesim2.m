function p = pulsesim2(dw, ksw, mnots, mnotw, R1S, R2S, R1W, R2W, sep, duration, curve, angle, init, TR, exciteflag, excitewait, exciteduration, exciteangle)
w1 = angle * sqrt(duration/TR);

[p,blah] = pulsesolv1(w1, dw, ksw, mnots, mnotw, R1S, R2S, R1W, R2W, sep, init, duration);  
index = size(blah); 
init = blah(index(1), :);

w1 = 0;

if (exciteflag == 1)
   [p,blah] = pulsesolv1(w1, dw, ksw, mnots, mnotw, R1S, R2S, R1W, R2W, sep, init, excitewait); 
   index = size(blah); 
   init = blah(index(1), :);
   
   w1 = exciteangle * sqrt(exciteduration/TR);
   
   [p,blah] = pulsesolv1(w1, dw, ksw, mnots, mnotw, R1S, R2S, R1W, R2W, sep, init, exciteduration);  
   index = size(blah); 
   init = blah(index(1), :);
   
   w1=0;
   [p,blah] = pulsesolv1(w1, dw, ksw, mnots, mnotw, R1S, R2S, R1W, R2W, sep, init, TR-(excitewait+duration+exciteduration)); 
else
   [p,blah] = pulsesolv1(w1, dw, ksw, mnots, mnotw, R1S, R2S, R1W, R2W, sep, init, TR-duration);
end
index = size(blah);
p = blah(index(1), :);
end