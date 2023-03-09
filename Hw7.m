%% Name: Sakthi Aadharsh Azhagar Gobinath Kavitha
%% Person ID: 50409103
%% Question 1
clear all
[t,h]=ode45(@(t,h) func(t,h),[0,10],[0,0]);
h
plot(t,h(:,1))
xlabel("Time (t)")
ylabel("Height (h)")
clear all
%% Question 2
% % clear all
syms x(z) 
u=8;
k=30;
Ka=5;
ca=1-x;
ca0=0.2;
l=6;
t=l/u;
tSpan=linspace(0,6,101);
%% a)
[z,x1]=ode45(@(z,x) nocat(z,x),tSpan,0);
figure(1)
plot(tSpan,x1);
%% b)
[z,x2]=ode45(@(z,x) cata1(z,x),tSpan,0);
i=1;
for z=linspace(0,6,101)
    a1(i)=a1_value(z,u);
    i=i+1;
end
figure(2)
plot(tSpan,x2)
hold on 
plot(tSpan,a1,"--")
hold off
%% c)
[z,x3,a2]=ode45(@(z,x) cata2(z,x),tSpan,0);
i=1;
for z=linspace(0,6,101)
    a2(i)=a2_value(z,u);
    i=i+1;
end
figure(3)
plot(tSpan,x3)
hold on 
plot(tSpan,a2,"--")
hold off
%% d)
[z,x4]=ode45(@(z,x) cata3(z,x),tSpan,0);
i=1;
for z=linspace(0,6,101)
    x=x4(i);
    a3(i)=a3_value(z,x,ca0,u);
    i=i+1;
end
figure(4)
plot(tSpan,x4)
hold on 
plot(tSpan,a3,"--")
hold off
%% Question 3
solveODEs
function solveODEs
ini=[100,0,0,423];
tSpan=[0:0.01:1];
[V,F]=ode45(@ODEs,tSpan,ini);
%% plot 1
T=F(:,4);
figure(1)
plot(V,T)
xlabel("V")
ylabel("T")
%% plot 2
figure(2)
plot(V,F(:,1))
hold on
plot(V,F(:,2))
hold on
plot(V,F(:,3))
hold off
xlabel("V")
ylabel("F")
legend("Fa","Fb","Fc")
end
function eq = ODEs(V,sys)
Fa=sys(1);
Fb=sys(2);
Fc=sys(3);
T=sys(4);
T0 = 373;
Cpa = 90;
Cpb = 90;
Cpc = 180;
k1a=10*exp(4000*((1/300)-(1/T)));
k2a=0.09*exp(9000*((1/300)-(1/T)));
ct0=0.1;
Ft=Fa+Fb+Fc;
ca=ct0*(Fa/Ft)*(T0/T);
r1=-k1a*ca;
r2=-k2a*ca^2;
dFadV=r1+r2;
dFbdV=-r1;
dFcdV=-r2/2;
dTdV=(4000*(373-T)+(-r1)*(20000)+(-r2)*(60000))/(Fa*Cpa+Fb*Cpb+Fc*Cpc);
eq(1,1)=dFadV;
eq(2,1)=dFbdV;
eq(3,1)=dFcdV;
eq(4,1) =dTdV;
end
%% Functions
function hdot=func(t,h)
L=150;%cm
R=0.25;%cm
mu=0.01;%poise (gcm-1s-1)
g=980.7;%cm2/s
rho=1;%g/cm3
fx=50;
hdot=zeros(2,1);
hdot(1)=h(2);
hdot(2)=g/L*(fx-(8*mu*L/(rho*g*R^2)*h(2)+h(1)));
end
function dxdz=nocat(z,x)
u=8;
k=30;
Ka=5;
ca=1-x;
ca0=0.2;
l=6;
t=z/u;
a=1;
dxdz=a*k*(1-x)/(u*(1+Ka*ca0*(1-x)));
end
function a1=a1_value(z,u)
t=z/u;
A_dash=12;
a1=1/(1+A_dash*t^1/2);
end
function [dxdz,a1]=cata1(z,x)
u=8;
k=30;
Ka=5;
ca=1-x;
ca0=0.2;
a1=a1_value(z,u);
dxdz=a1*k*(1-x)/(u*(1+Ka*ca0*(1-x)));
end
function a2=a2_value(z,u)
t=z/u;
kd2=17.5;
a2=1/(1+3*kd2*t)^1/3;
end
function [dxdz,a2]=cata2(z,x)
u=8;
k=30;
Ka=5;
ca=1-x;
ca0=0.2;
a2=a2_value(z,u);
dxdz=a2*k*(1-x)/(u*(1+Ka*ca0*(1-x)));
end
function a3=a3_value(z,x,ca0,u)
t=z/u;
kd3=140;
a3=1/(2*(0.5+kd3*ca0*x*t)^1/3);
end
function dxdz=cata3(z,x)
u=8;
k=30;
Ka=5;
ca=1-x;
ca0=0.2;
a3=a3_value(z,x,ca0,u);
dxdz=a3*k*(1-x)/(u*(1+Ka*ca0*(1-x)));
end
