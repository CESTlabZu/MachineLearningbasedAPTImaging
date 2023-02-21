function l = runsteadysimsquare(ksw, mnots, mnotw, R1S, R2S, R1W, R2W, sep, duration, curve, angle, TR, time, exciteflag, excitewait, exciteduration, exciteangle, k, spoil)
%l = zeros(11, 6);
%k= 35000:5000:85000;
%k= 100:100:1000;
l = zeros(size(k), 6);

for t=1:1:size(k)
    p = steadystatepulsesimsquare(k(t), ksw, mnots, mnotw, R1S, R2S, R1W, R2W, sep, duration, curve, angle, TR, time, exciteflag, excitewait, exciteduration, exciteangle, spoil);
    ind2 = size(p);
    l(t, :) = p(ind2(1), :)
end