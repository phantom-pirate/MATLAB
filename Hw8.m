%% Name:Sakthi Aadharsh Azhagar Gobinath Kavitha
%% PersonID :50409103
%% Question 1
clear all
L=150;%cm
R=0.25;%cm
mu=0.01;%poise (gcm-1s-1)
g=980.7;%cm2/s
rho=1;%g/cm3
fx=50;
A=[0,1;-(8*mu*L/(rho*g*R^2))/(L/g),-(1/(L/g))];
h0 = [0,0];
tspan = [0 10];
[t,h] = ode45(@(t,h)func(t,h),tspan,h0);
figure(1)
plot(h(:,1),h(:,2))
xlabel('h1');
ylabel('h2')
figure(2)
hold on
for h10=0:2:10
    for h20=0:2:10
        h0=[h10,h20];
        tspan=[0 10];
        [t,h]=ode45(@(t,h)func(t,h),tspan,h0);
        plot(h(:,1),h(:,2))            
    end
end
hold off
xlabel('h1');
ylabel('h2');
figure(3)
plot(real(eig(A)),imag(eig(A)),'o')
grid on
axis([-3 3 -3 3])
ax = gca;
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xlabel('Re'), ylabel('Im')
%% b)
%% Question 2
clear all
%% a)
syms x(t) y(t)
a = [2 -1; 1 -3];
b=[0;0]
sol=a\(-b);
e=eig(a);
if e(1)/e(2)==1
    if(e(1)>0)
        disp("Unstable")
    elseif(e(1)>0)
        disp("Stable")
    end
elseif and(isreal(e(1)),isreal(e(2)))
    if and(e(1)>0,e(2)>0)
        disp("Unstable Node")
    elseif e(1)*e(2)<0
        disp("Saddle Point")
    elseif and(e(1)<0,e(2)<0)
        disp("Stable Node")
    end
else
    if real(e(1))==0
        disp("Center")
    elseif real(e(1))>0
        disp("Unstable spiral")
    elseif real(e(1))<0
        disp("stable spiral")
    end
end
%% b)
clear all
syms x(t) y(t)
a = [1 -1; 1 3];
b=[0;-4];
sol=a\(-b)
e=eig(a);
if e(1)/e(2)==1
    if(e(1)>0)
        disp("Unstable")
    elseif(e(1)<0)
        disp("Stable")
    end
elseif and(isreal(e(1)),isreal(e(2)))
    if and(e(1)>0,e(2)>0)
        disp("Unstable Node")
    elseif e(1)*e(2)<0
        disp("Saddle Point")
    elseif and(e(1)<0,e(2)<0)
        disp("Stable Node")
    end
else
    if real(e(1))==0
        disp("Center")
    elseif real(e(1))>0
        disp("Unstable spiral")
    elseif real(e(1))<0
        disp("stable spiral")
    end
end
%% c)
clear all
syms x(t) y(t)
a = [1 -2; 1 -1];
b=[3;2];
sol=a\(-b)
e=eig(a);
if e(1)/e(2)==1
    if(e(1)>0)
        disp("Unstable")
    elseif(e(1)<0)
        disp("Stable")
    end
elseif and(isreal(e(1)),isreal(e(2)))
    if and(e(1)>0,e(2)>0)
        disp("Unstable Node")
    elseif e(1)*e(2)<0
        disp("Saddle Point")
    elseif and(e(1)<0,e(2)<0)
        disp("Stable Node")
    end
else
    if real(e(1))==0
        disp("Center")
    elseif real(e(1))>0
        disp("Unstable spiral")
    elseif real(e(1))<0
        disp("stable spiral")
    end
end
%% d)
clear all
syms x(t) y(t)
a = [2 -2; 1 4];
b=[-4;3];
sol=a\(-b)
e=eig(a);
if e(1)/e(2)==1
    if(e(1)>0)
        disp("Unstable")
    elseif(e(1)<0)
        disp("Stable")
    end
elseif and(isreal(e(1)),isreal(e(2)))
    if and(e(1)>0,e(2)>0)
        disp("Unstable Node")
    elseif e(1)*e(2)<0
        disp("Saddle Point")
    elseif and(e(1)<0,e(2)<0)
        disp("Stable Node")
    end
else
    if real(e(1))==0
        disp("Center")
    elseif real(e(1))>0
        disp("Unstable spiral")
    elseif real(e(1))<0
        disp("stable spiral")
    end
end

%% Question 3
%% a)
clear all
syms f1(x,y) f2(x,y) 
f1=1-y^2;
f2=x+2*y;
s=solve([f1,f2],[x,y])
a=matlabFunction(jacobian([f1;f2],[x y]),"Vars",[x y]);
%% at the point 1
e=vpa(eig(a(s.x(1),s.y(1))))
disp("Metastable, Saddle point")
%% at the point 2
e=eig(a(s.x(2),s.y(2)))
disp("Unstable, Oscillatory")
%% b)
clear all
syms f1(x,y) f2(x,y) 
f1=2-4*x-15*y;
f2=4-x^2;
s=solve([f1,f2],[x,y])
a=matlabFunction(jacobian([f1;f2],[x y]),"Vars",[x y]);
%% at the point 1
e=vpa(eig(a(s.x(1),s.y(1))))
disp("Stable,Oscillatory")
%% at the point 2
e=eig(a(s.x(2),s.y(2)))
disp("Metastable, Saddle point")
%% c)
clear all
syms f1(x,y) f2(x,y) 
f1=x-2*y;
f2=4*x-x^3;
s=solve([f1,f2],[x,y])
a=matlabFunction(jacobian([f1;f2],[x y]),"Vars",[x y]);
%% at the point 1
e=vpa(eig(a(s.x(1),s.y(1))))
disp("Non stable,Oscillatory")
%% at the point 2
e=vpa(eig(a(s.x(2),s.y(2))))
disp("Metastable, Saddle point")
%% at the point 3
e=vpa(eig(a(s.x(3),s.y(3))))
disp("Metastable, Saddle point")
%% d)
clear all
syms f1(x,y) f2(x,y) 
f1=x-y-x^2+x*y;
f2=-y-x^2;
s=solve([f1,f2],[x,y])
disp("Critical point 1")
[s.x(1),s.y(1)]
a=matlabFunction(jacobian([f1;f2],[x y]),"Vars",[x y]);
%% at the point 1
e=vpa(eig(a(s.x(1),s.y(1))))
disp("Non stable,Oscillatory")
%% at the point 2
e=vpa(eig(a(s.x(2),s.y(2))))
disp("Metastable, Saddle point")
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