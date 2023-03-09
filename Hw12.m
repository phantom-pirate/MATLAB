%% Name:Sakthi Aadharsh Azhagar Gobinath Kavitha
%% PersonID :50409103
%% Question 1
clear all
h=80;
dx=1.2/100;
dy=dx;
k=15;
e=2*10^6;
Ts=25;
q=5000;
syms T1 T2 T3 T4 T5 T6 T7 T8 T9
%node 1
n1=h*dx/2*(Ts-T1)+k*dy/2*(T2-T1)/dx+k*dx/2*(T4-T1)/dx+e*dx*dy/4;
%node 2
n2=h*dx*(Ts-T2)+k*dy/2*(T3-T2)/dx+k*dx*(T5-T2)/dx+k*dy/2*(T1-T2)/dx+e*dx*dy/2;
%node 3
n3=h*(dx/2+dy/2)*(Ts-T3)+k*dx/2*(T6-T3)/dy+k*dy/dx*(T2-T3)/2+e*dx*dy/4;
%node 4
n4=k*dy*(T5-T4)/dx+k*dx/2*(T1-T4)/dy+k*dx/2*(90-T4)/dy+e*(dx*dy)/2
%node 5
n5=k*dy*(T4-T5)/dx+k*dx*(T6-T5)/dy+k*dx*(90-T5)/dy+k*dx*(T2-T5)/dy+e*(dx*dy);
%node 6
n6=h*(dx/2+dy/2)*(Ts-T6)+k*dy/2*(T7-T6)/dx+k*dx*(90-T6)/dy+k*dy*(T5-T6)/dx+k*dx/2*(T3-T6)/dy+e*3*dx*dy/4;
%node 7
n7=h*dx*(Ts-T7)+k*dy/2*(T8-T7)/dx+k*dx*(90-T7)/dy+k*dy/2*(T6-T7)/dx+e*dx*dy/2;
%node 8
n8=h*dx*(Ts-T8)+k*dy/2*(T9-T8)/dx+k*dx*(90-T8)/dy+k*dy/2*(T7-T8)/dx+e*dx*dy/2;
%node 9
n9=h*dx/2*(Ts-T9)+q*dy/2+k*dx/2*(90-T9)/dy+k*dy/2*(T8-T9)/dx+e*dx*dy/4;
%system of equation
[A,B]=equationsToMatrix([n1, n2, n3, n4, n5, n6, n7, n8, n9], [T1, T2, T3, T4, T5, T6, T7, T8, T9]);
X=vpa(linsolve(A,B));
%% Question 2
clear all
k=15;
h=80;
dx=1.2/100;
dy=dx;
l=dx;
q=5000;
alpha=3.2*10^-6;
e=2*10^6;
Ts=25;
dt=10;
t=alpha*dt/l^2;
i=1;
T1(i)=25;
T2(i)=25;
T3(i)=25;
T4(i)=25;
T5(i)=25;
T6(i)=25;
T7(i)=25;
T8(i)=25;
T9(i)=25;
hr=[60,180,300,600,3600];
for j=hr
    for i=1:j/dt
        %node 1
        T1(i+1)=(1-4*t-2*t*h*l/k)*T1(i)+2*t*(T2(i)+T4(i)+h*l/k*Ts+e*l^2/(2*k));
        %node 2
        T2(i+1)=(1-4*t-2*t*h*l/k)*T2(i)+t*(T1(i)+T3(i)+2*T5(i)+2*h*l/k*Ts+e*l^2/k);
        %node 3
        T3(i+1)=(1-4*t-4*t*h*l/k)*T3(i)+2*t*(T2(i)+T6(i)+2*h*l/k*Ts+e*l^2/(2*k));
        %node 4
        T4(i+1)=(1-4*t)*T4(i)+t*(T1(i)+2*T5(i)+90+e*l^2/k);
        %node 5
        T5(i+1)=(1-4*t)*T5(i)+t*(T2(i)+T4(i)+T6(i)+90+e*l^2/k);
        %node 6
        T6(i+1)=(1-4*t-4*t*h*l/(3*k))*T6(i)+t/3*(2*T3(i)+4*T5(i)+2*T7(i)+4*90+4*h*l/k*Ts+3*3*l^2/k);
        %node 7
        T7(i+1)=(1-4*t-2*t*h*l/k)*T7(i)+t*(T6(i)+T8(i)+2*90+2*h*l/k*Ts+e*l^2/k);
        %node 8
        T8(i+1)=(1-4*t-2*t*h*l/k)*T8(i)+t*(T7(i)+T9(i)+2*90+2*h*l/k*Ts+e*l^2/k);
        %node 9
        T9(i+1)=(1-4*t-2*t*h*l/k)*T9(i)+2*t*(T8(i)+90+q*l/k+h*l/k*Ts+e*l^2/(2*k)); 
    end
    T3(i)
end
%% Question 3
clear all
N=100 ;
T0(1,1)=100;
T0(2:N+1,1)=25;

tsp= [0 20];
x0=T0(2:N);
[t,T]=ode45(@(t,x) rodtemp(t,x),tsp,x0);
plot(t,T)
xlabel("Time")
ylabel("Temperature")


function fval=rodtemp(t,x)
%parameters
Ts=25;
lambda=424;
rho=10500;
cp=236;
k1=5.25*10^-7;
a = 0;
b = 1;
N = 101;
dx = (b-a)/(N);
alpha=lambda/(rho*cp*dx^2);
%Getting Temperature
N=length(x)+1;
T(1)=100;
T(2:N)=x;
T(N+1)=Ts;

%define dTdt
dTdt=zeros(N+1,1);
for i=2:N
    dTdt(i,1)=alpha*(T(i+1)-2*T(i)+T(i-1))-k1*(T(i)-Ts); 
end
fval=dTdt(2:N);
end
