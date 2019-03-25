function f = sys(t,x,I)
global wI1 w11 w12 w13 w23 rx1 rx2 rx3 Kl1 Kl2 Kl3 tl1 tl2 tl3 Sxp Rlt Dx Dl mu

%x = [m1;m2;m3;p1;p2;p3];
n = 2; %General n
n23 = 10; %n protein 2 to 3
K = 1*10^(-12); %General K
Ki = 0.3*10^(-5); %Inducer K

a = round(t);
In = I(a+1,1);
m1 = x(1);
m2 = x(2);
m3 = x(3);
p1 = x(4);
p2 = x(5);
p3 = x(6);


fI1 = (In^n)/(Ki^n+In^n); %Inducer function
f12 = (p1^n)/(K^n+p1^n); %P1 to P2 function
f13 = (p1^n)/(K^n+p1^n); %P1 to P3 function
f23 = (p2^n23)/(K^n23+p2^n23); %P2 to P3 function

ui = (w11+wI1*fI1)/(1+w11+wI1*fI1);
up1 = (w11+w12*f12)/(1+w11+w12*f12);
up1p2 = (w11+w13*f13)/(1+w11+w13*f13+w23*f23);


rl1 = Kl1*Rlt*(m1/(Sxp*tl1+m1*tl1+m1))*10^3; %Translation rate for gene 1
rl2 = Kl2*Rlt*(m2/(Sxp*tl2+m2*tl2+m2))*10^3; %Translation rate for gene 2
rl3 = Kl3*Rlt*(m3/(Sxp*tl3+m3*tl3+m3))*10^3; %Translation rate for gene 3 


f = zeros(6,1);
f(1,1) = rx1*ui-m1*(Dx+mu);
f(2,1) = rx2*up1-m2*(Dx+mu);
f(3,1) = rx3*up1p2-m3*(Dx+mu);
f(4,1) = rl1-p1*(Dl+mu);
f(5,1) = rl2-p2*(Dl+mu);
f(6,1) = rl3-p3*(Dl+mu);

end