%% Homework 2
%% Name and Person ID : Sakthi Aadharsh Azhagar Gobinath Kavitha [50409103]
%% 1)

%% Picard's Method
clear all
funtol=0.0001;
numi=1000;
k=0;
for x=[-0.1, 0.6, 1.99,2.01]
    xi=x;
    k=0;
    gap=0;        
    while 1        
        f=x^2-3*x+2;        
        if abs(f)<funtol
            sprintf("For the initial value of %d, the system has converged at %d for the value of x = %d",xi,k,x)      
            break;
        elseif k==numi
                sprintf("For the initial value of %d, the system has not converged at %d",xi,k)      
                break; 
        else
            k=k+1;
            x=x+f; 
        end                       
    end
end
syms x
f=x^2-3*x+2;  
fplot(f,[-0.5,2.5])
%% 2)
%% a)
syms x
%% Given Data
rf=10;
k=0.4655;
%% Stopping Criteria and Initial guess
funtol=0.0001;
numi=1000;
a=0.67;% value must be between 0 and 1
%% Input the function
fa=(((rf^2)*(x^3)/(20*(1-x)^2))-k);
df=matlabFunction(diff(fa),"var",{x});
%% Newton raphson method
while 1
    f=(((rf^2)*(a^3)/(20*(1-a)^2))-k);   
    a=a-f/df(a);       
    if abs(f/df(a))<funtol     
        break;
    elseif k>numi
        break;
    end 
i=i+1;
end
sprintf("The value obtained via Newton-Raphson Method is %d",a)

%% b) fzero function
clear x
j=1;

 for i=linspace(0,1)
     if i<1
         x(j)=fzero(@(x) (((rf^2)*(x^3)/(20*(1-x)^2))-k),i); 
         j=j+1;
     end
 end
 sprintf("The value obtained via fzero Method is %d",(uniquetol(x,funtol)))
%% 3)
clear all
%% Initialise
R=0.032;
k=1.9*10^-4;
syms r
D=4*10^-7;
%% Relations
phi=R*sqrt(k/D);
gamma=r/R;
%% Input the function
f=matlabFunction((1/gamma)*(exp(-phi*gamma-1))-0.5);
%% plot
fplot(f)
% From the plot the initial value is taken to be 0.1
%% Solution
solution=fzero(f,0.1)

%% 4)
%system of linear and non linear system to model the working of microscale liver and lungs .
%% a) Linear system of equations
clear all
%lung
%llu=(2-1.5*R)*Clu-0.5*R*Cl-779.3+780*R;
%liver
%lli=Clu-Cl-76;
i=1;
for R=linspace(0.6,0.9)
    la=[(2-1.5*R),-0.5*R;1,-1];
    lb=[-(-779.3+780*R);76];
    lx=la\lb;
    Clung(i)=lx(1);
    Cliver(i)=lx(2);
    i=i+1;
end
R=linspace(0.6,0.9);
plot(R,Clung,"-",R,Cliver,"--")
legend("Concentration in lung","Concentration in liver")
title("Plot of Naphthalene Concentration in lungs and liver vs Recycle ratio")

%% b) Non Linear system of equation:
clear all
%% initialise
r=[0.6,0.8,0.9];
R=r(1);
syms x y 
%lung
fy=780*(1-R)+R*(0.5*x+1.5*y)-(8.75*y)/(2.1+y)*0.08-2*y;
%liver
fx=0.5*y-(118*x)/(7+x)*0.322-0.5*x;
f=[fx;fy];
a=[0;0];
k=1;
 while 1
    j=jacobian(f,[x y]);
    b=a;
    jf=matlabFunction(j*f);   
    a=a-jf(a(1),a(2));
    k=k+1;
    if and(abs(a(1)-b(1))<0.0001,abs(a(2)-b(2))<0.001)
        break;
    elseif k-1>1000
        break;
    end
 end
sol_nl(:,1)=a;
a=[0;0];
R=r(2);
 while 1
    j=jacobian(f,[x y]);
    b=a;
    jf=matlabFunction(j*f);   
    a=a-jf(a(1),a(2));
    k=k+1;
    if and(abs(a(1)-b(1))<0.0001,abs(a(2)-b(2))<0.001)
        break;
    elseif k-1>1000
        break;
    end
 end
sol_nl(:,2)=a;
a=[0;0];
R=r(3);
 while 1
    j=jacobian(f,[x y]);
    b=a;
    jf=matlabFunction(j*f);   
    a=a-jf(a(1),a(2));
    k=k+1;
    if and(abs(a(1)-b(1))<0.0001,abs(a(2)-b(2))<0.001)
        break;
    elseif k-1>1000
        break;
    end
 end
sol_nl(:,3)=a
i=1;
for R=[0.6,0.8,0.9]
    la=[(2-1.5*R),-0.5*R;1,-1];
    lb=[-(-779.3+780*R);76];
    lx=la\lb;
    Clung(i)=lx(1);
    Cliver(i)=lx(2);
    i=i+1;
end
sol_l=[Cliver;Clung];
rel_err=abs(sol_l-sol_nl);
R=[0.6,0.8,0.9];
plot(R,rel_err(1,:),"-",R,rel_err(2,:),"--")
legend("Relative Error in Concentration in lung","Relative Error in Concentration in liver")
title("Plot of Relative Error in Naphthalene Concentration in lungs and liver vs Recycle ratio")
%% c) fsolve Method
i=1;
for R=[0.6,0.8,0.9]
    f= @(x) [0.5*x(2)-(118*x(1))/(7+x(1))*0.322-0.5*x(1);780*(1-R)+R*(0.5*x(1)+1.5*x(2))-(8.75*x(2))/(2.1+x(2))*0.08-2*x(2)];
    sol_n(:,i)=fsolve(f,[0,0]);
    i=i+1;
end
rel_err=abs(sol_l-sol_n);
R=[0.6,0.8,0.9];
plot(R,rel_err(1,:),"-",R,rel_err(2,:),"--")
legend("Relative Error in Concentration in lung","Relative Error in Concentration in liver")
title("Plot of Relative Error in Naphthalene Concentration in lungs and liver vs Recycle ratio")


