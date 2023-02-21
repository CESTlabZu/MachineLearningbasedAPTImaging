function [a,b] = pulsesolv1(w1, dw, ksw1,ksw2,ksw3,ksw4, ksw5,kmw, mnots1,mnots2,mnots3,mnots4,mnots5, mnotw, mnotm, R1S, R2S,R2S5, R1W, R2W, R1M, R2M, sep1,sep2,sep3, sep4,sep5, init, duration)
tspan = [0 duration];
y0 = init;
dw=dw;
delter=2.3*400*2*pi;
%ls=1./(pi*R2M)./(1+((dw-delter)./R2M).^2);
%W=w1^2*sqrt(pi/2)/R2M*exp(-((dw-delter)/R2M)^2/2);

ls=1./(pi*R2M)./(1+((dw-delter)./R2M).^2);
W=(w1)^2*pi*ls;  


%delter=0*400*5*2*pi;
%D=1000;
%lls=0;
%for ii=1:1:D+1
%    ls=sin(pi/2/D*(ii-1))*sqrt(2/pi)/abs(3*cos(pi/2/D*(ii-1))^2-1)/R2M*exp(-2*((dw-delter)/R2M/(3*cos(pi/2/D*(ii-1))^2-1))^2);
%   ls=ls+lls;
%       lls=ls;
%end
%ls=ls/(D+1);




%W=(w1)^2*pi*ls;  
[a,b] = ode45(@f, tspan, y0);
    function dydt = f(t,y)   
     dydt = [ ((-sep1 - dw)*y(2) - (R2S + ksw1)*y(1) + mnots1/mnotw*ksw1*y(4))
             ((sep1 + dw)*y(1) + w1*y(3) - (R2S + ksw1) * y(2) + mnots1/mnotw * ksw1 * y(5))
             (-w1*y(2) - R1S*(y(3) - mnots1) - ksw1*y(3) + mnots1/mnotw*ksw1*y(6))
             (-dw*y(5)-R2W*y(4)+ksw1*y(1)- mnots1/mnotw*ksw1*y(4)+ksw2*y(8)- mnots2/mnotw*ksw2*y(4))
             (dw*y(4) + w1*y(6)-R2W*y(5) + ksw1*y(2) - mnots1/mnotw*ksw1 *y(5)+ ksw2*y(9) - mnots2/mnotw*ksw2 *y(5))
             (-w1*y(5) - R1W*(y(6) - mnotw) + ksw1*y(3) - mnots1/mnotw*ksw1*y(6)+ ksw2*y(10) - mnots2/mnotw*ksw2*y(6)+ ksw3*y(13) - mnots3/mnotw*ksw3*y(6)+ ksw4*y(16) - mnots4/mnotw*ksw4*y(6)+ ksw5*y(19) - mnots5/mnotw*ksw5*y(6)+kmw*y(7)-mnotm/mnotw*kmw*y(6))
             (- R1M*(y(7) - mnotm) - kmw*y(7) + mnotm/mnotw*kmw*y(6)-W*y(7))
             ((-sep2 - dw)*y(9) - (1/0.01 + ksw2)*y(8) + mnots2/mnotw*ksw2*y(4))
             ((sep2 + dw)*y(8) + w1*y(10) - (1/0.01 + ksw2) * y(9) + mnots2/mnotw * ksw2 * y(5))
             (-w1*y(9) - R1S*(y(10) - mnots2) - ksw2*y(10) + mnots2/mnotw*ksw2*y(6))    
             ((-sep3 - dw)*y(12) - (1/0.01 + ksw3)*y(11) + mnots3/mnotw*ksw3*y(4))
             ((sep3 + dw)*y(11) + w1*y(13) - (1/0.01 + ksw3) * y(12) + mnots3/mnotw * ksw3 * y(5))
             (-w1*y(12) - R1S*(y(13) - mnots3) - ksw3*y(13) + mnots3/mnotw*ksw3*y(6))   
             ((-sep4 - dw)*y(15) - (1/0.001 + ksw4)*y(14) + mnots4/mnotw*ksw4*y(4))
             ((sep4 + dw)*y(14) + w1*y(16) - (1/0.001 + ksw4) * y(15) + mnots4/mnotw * ksw4 * y(5))
             (-w1*y(15) - R1S*(y(16) - mnots4) - ksw4*y(16) + mnots4/mnotw*ksw4*y(6))   
             ((-sep5 - dw)*y(18) - (R2S5 + ksw5)*y(17) + mnots5/mnotw*ksw5*y(4))
             ((sep5 + dw)*y(17) + w1*y(19) - (R2S5 + ksw5) * y(18) + mnots5/mnotw * ksw5 * y(5))
             (-w1*y(18) - R1S*(y(19) - mnots5) - ksw5*y(19) + mnots5/mnotw*ksw5*y(6))   
               ];
             
           
    end
end
