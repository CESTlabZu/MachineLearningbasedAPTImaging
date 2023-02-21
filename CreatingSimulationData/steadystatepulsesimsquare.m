function p = steadystatepulsesimsquare(dw, ksw, mnots, mnotw, R1S, R2S, R1W, R2W, sep, duration, curve, angle, TR, time, exciteflag, excitewait, exciteduration, exciteangle, spoil)
p = zeros(time,6);
init = [0;0;0;0;0;1];
t = 1;
steadycheck = 100000;
steadystate = 1;
while (steadystate == 1 && t < time)
   blah = pulsesim2(dw, ksw, mnots, mnotw, R1S, R2S, R1W, R2W, sep, duration, curve, angle, init, TR, exciteflag, excitewait, exciteduration, exciteangle); 
   %blah;
   p(t,:) = blah;
   if(spoil == 1)
     blah(4) = 0;
     blah(5) = 0;
     blah(1) = 0;
     blah(2) = 0;
   end
   init = blah;
   t = t+1;
   if ( mod(t,4) == 0 )
       if (floor(blah(6)*1000) == floor(steadycheck*1000 ))
           steadystate = 0;
           %t;
       end
       steadycheck = blah(6);
   end
   %p(t,:)= init;
end
p(t:time, :) = [];
%p = init;
end