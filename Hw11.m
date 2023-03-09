%% Name:Sakthi Aadharsh Azhagar Gobinath Kavitha
%% PersonID :50409103
%% Question 2
clc
clear all
[x,c] = themal_conduc;
plot(x,c)
hold on
%% Question 3
clear all
[x,c] = themal_co
plot(x,c)
hold off
%% Question 4
N=101;
t=linspace(0,2*180,N);
dt=0.5;%assumption
theta_gues=0.7*cos(t)+0.5*sin(t);
theta_final=fsolve(@(y)pendulum(y,dt,N),theta_gues);
plot(t,theta_gues)
hold on
plot(t,theta_final)
legend("linearized", "numerically solved")
hold off
ylabel("theta")
xlabel("time")
%% Question 5
clear all
syms c(x)
e=diff(diff(c,x),x)-c;
ana=matlabFunction(dsolve(e,c(0)==0,c(1)==0.5));
dx=0.01;
x=0:dx:1;
N=1/dx+1;
a=zeros(N);
f=zeros(N,1);
ca=0;
cb=0.5;
alpha=-2-dx^2;
j=1;
for i=3:N-1
    a(i,i-1)=1;
    a(i,i)=alpha;
    a(i,i+1)=1;
    f(i)=0;
end
i=2;
a(i,i)=alpha;
a(i,i+1)=1;
f(i)=-ca;
i=N-1;
a(i,i)=alpha;
a(i,i-1)=1;
f(i)=-cb;
a= a(2:N-1,2:N-1); % creating a matrix of only the unknowns for c2:cn-1
f= f(2:N-1);
y = a\f;
y=[ca;y;cb]
figure(1)
plot(x,y,"o")
hold on
plot(x,ana(x))
hold off
f

for dx=linspace(0.001,0.5)
    x=0:dx:1;
    N=length(x);
    a=zeros(N);
    f=zeros(N,1);
    ca=0;
    cb=0.5;
    alpha=-2-dx^2;
    for i=3:N-1
        a(i,i-1)=1;
        a(i,i)=alpha;
        a(i,i+1)=1;
        f(i)=0;
    end
    i=2;
    a(i,i)=alpha;
    a(i,i+1)=1;
    f(i)=-ca;
    i=N-1;
    a(i,i)=alpha;
    a(i,i-1)=1;
    f(i)=-cb;
    a= a(2:N-1,2:N-1); % creating a matrix of only the unknowns for c2:cn-1
    f= f(2:N-1);
    y = a\f;
    y=transpose([ca;y;cb]);
    n(j)=norm((y-ana(x)),2)/norm(y,2);
    j=j+1;
end
figure(2)
dx=linspace(0.001,0.5)
loglog(dx,n)
%% Question 6
clear all
clc

dz=0.01;
N=1/dz+1;

x_ini=ones(N,1);
y_ini=ones(N,1);
ini_gues=[x_ini,y_ini];
func=@(x,y) coupled(x,y,dz,N)
sol=fsolve(func,ini_gues)
function sys=coupled(x,y,dz,N)

e1=zeros(1,N);
e2=zeros(1,N);
e1(1)=x(3);%sample boundary to be changed into dirichlet BC
e2(1)=y(3)-2*x(i)*dz;%sample boundary to be changed into dirichlet BC
for i = 2: N-1
 e1(i)=x(i+1)-2*x(i)+x(i-1)+dz^2*x(i)*exp(y(i))
 e2(i)=y(i+1)-2*y(i)+y(i-1)+dz^2*(x(i)+y(i)^2)
end

e1(1)=x(N)-0.8;
e2(1)=y(N)-1;
sys=[e1 e2];
end
%% function
function [x,c] = themal_co
a=0; 
ca=1;
b = 1; 
cb = 5;
dx = 0.01
n=1/dx+1;
x =dx;
A = zeros(n);
f = zeros(n,1);
A(1,1)=1;
f(1)=ca;
for i = 2:n-1     
    kim1=(1+1/4*((x-dx)+x)^2);
    kip1=(1+1/4*((x+dx)+x)^2);
    ki=(1+1/4*(((x-dx)+x)^2)+((x+dx)+x)^2);
    A(i,i-1)=(ki/dx^2-(kip1-kim1)/(4*dx^2));
    A(i,i)=-2*ki/dx^2;
    A(i,i+1)=(ki/dx^2+(kip1-kim1)/(4*dx^2));
    f(i) = 0;   
    x=x+dx;
end
A(n,n) =1;
f(n) = cb;

A = sparse(A);
c = A\f;
x=a:dx:b;
x = x';
end

function [x,c] = themal_conduc
a=0; 
ca=1;
b = 1; 
cb = 5;
dx = 0.01
n=1/dx+1;
x =dx;
A = zeros(n);
f = zeros(n,1);
A(1,1)=1;
f(1)=ca;
for i = 2:n-1     
    kim1=(1+(x-dx)^2);
    kip1=(1+(x+dx)^2);
    ki=(1+x^2);
    A(i,i-1)=(ki/dx^2-(kip1-kim1)/(4*dx^2));
    A(i,i)=-2*ki/dx^2;
    A(i,i+1)=(ki/dx^2+(kip1-kim1)/(4*dx^2));
    f(i) = 0;   
    x=x+dx;
end
A(n,n) =1;
f(n) = cb;

A = sparse(A);
c = A\f;
x=a:dx:b;
x = x';
end

function R=pendulum(y,dt,N)
R=zeros(1,N);
alpha=0.7;
beta=alpha;
R(1)=y(1)-alpha;
for i = 2: N-1
R(i)=(y(i+1)-2*y(i)+y(i-1))+(dt^2)*sin(y(i));
end
R(N)=y(N)-beta;
end
