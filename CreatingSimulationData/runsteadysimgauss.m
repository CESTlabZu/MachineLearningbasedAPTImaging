function l = runsteadysimgauss(ksw1,ksw2, ksw3,ksw4, ksw5,kmw, mnots1,mnots2, mnots3, mnots4, mnots5, mnotw, mnotm, R1S, R2S,R2S5, R1W, R2W, R1M, R2M, sep1,sep2,sep3,sep4, sep5, duration, curve, angle, TR, time, exciteflag, excitewait, exciteduration, exciteangle, k, spoil)
%l = zeros(11, 6);
%k= 35000:5000:85000;
%k= 100:100:1000;
l = zeros(89, 19);

for t=1:1:size(k)
    p = steadystatepulsesimgauss(k(t), ksw1,ksw2, ksw3,ksw4, ksw5, kmw, mnots1,mnots2, mnots3,mnots4, mnots5, mnotw, mnotm, R1S, R2S,R2S5, R1W, R2W, R1M, R2M, sep1,sep2, sep3,sep4,sep5, duration, curve, angle, TR, time, exciteflag, excitewait, exciteduration, exciteangle, spoil);
    ind2 = size(p);
    l(t, :) = p(ind2(1), :)
end