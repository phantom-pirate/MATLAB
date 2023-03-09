%% Name:Sakthi Aadharsh Azhagar Gobinath Kavitha
%% PersonID :50409103
%% Question 2
syms u(t) v(t)
e1=diff(u)==2*u+v+4;
e2=diff(v)==-4*u-3*v+7;
sys=[e1;e2];
sol=dsolve(sys,u(0)==1,v(0)==2)
%% Question 3
syms Ca(t) Cb(t) k1 k2 Ca0
e1=diff(Ca(t))==-k1*Ca;
e2=diff(Cb(t))==k1*Ca-k2*Cb;
sys=[e1;e2];
sol=dsolve(sys,Ca(0)==Ca0,Cb(0)==0)
%% Question 6
syms y(x)
%% a)
c=x*diff(y,x)-y
sol=dsolve(c,y(2)==5)
%% b)
c=x*diff(diff(y,x),x)+diff(y)
sol=dsolve(c)
%% c)
Dy=diff(y,x);
c=x^2*diff(diff(y,x),x)+x*diff(y)-9*y
sol=dsolve(c,y(2)==1,Dy(2)==2)
%% d)
Dy=diff(y,x);
c=x^2*diff(diff(y,x),x)+x*diff(y)+y
sol=dsolve(c,y(1)==0,Dy(1)==2)
%% e)
c=x^3*diff(diff(diff(y,x),x),x)+x*diff(y)-y
sol=dsolve(c)
%% Question 8
clear all
syms y(x)
%% a)
c=diff(diff(y,x),x)-y-8*x;
sol=dsolve(c)
%% b)
c=diff(diff(y,x),x)-y-8*exp(x);
sol=dsolve(c)
%% c)
c=diff(diff(y,x),x)+4*diff(y,x)+4*y-2*exp(-2*x);
sol=dsolve(c)
%% Question 4
clear all
syms x1(t) x2(t) x3(t)
r=10;
V1=50;
V2=25;
V3=50;
e1=diff(x1)+r/V1*x1-r/V3*x3
e2=diff(x2)-r/V1*x1+r/V2*x2
e3=diff(x3)-r/V2*x2+r/V3*x3
sys=[e1;e2;e3];
sol=dsolve(sys)
%% Question 1
clear all
syms h(t)
L=150;%cm
R=0.25;%cm
mu=0.01;%poise (gcm-1s-1)
g=980.7%cm2/s
rho=1;%g/cm3
dh=diff(h,t)
c1=L/g*diff(diff(h,t),t)+8*mu*L/(rho*g*R^2)*diff(h,t)+h;%t<=0
c2=L/g*diff(diff(h,t),t)+8*mu*L/(rho*g*R^2)*diff(h,t)+h-50;%t>0
sol1=double(dsolve(c1,h(0)==0,dh(0)==0));
sol2=matlabFunction(dsolve(c2,h(0)==0,dh(0)==0));
i=1;
for t=linspace(-1,10)
    if t>0
        p(i)=sol2(t);
    else
        p(i)=sol1;
    end
    i=i+1;
end
t=linspace(-1,10);
p
plot(t,p)
xlabel("Time (t)")
ylabel("Height (h)")




