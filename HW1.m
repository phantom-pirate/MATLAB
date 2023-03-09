%% Homework 1
%% Name and Person ID : Sakthi Aadharsh Azhagar Gobinath Kavitha [50409103]
%% Question 2
%% 2)
clear all
%% System of equations
a=[4,2,-1,1;0,2,2,1;0,0,2,4;0,0,0,8]
b=[0;1;-1;4]
%% Classification of the system
if det(a)~=0
    if rank(a)==rank([a b])
        if rank(a)==length(a)          
            x=a\b
            sprintf("The system has an unique solution (%f, %f, %f, %f)",x)
        end
    end
elseif rank(a)~=rank([a b])
    sprintf("The system has no solution")
else
    sprintf("The system has infinitely many solution")    
end
%% Question 4
%% (a)
clear all
mat_a=[1,1,0,-2,0,0;1,0,-1,0,0,1;-1,0,0,1,1,-1;1,0,0,0,-2,0;0,1,0,0,0,-2];
%% Rank of the matrix
r_a=rank(mat_a)
%% Length of the matrix
l_a=length(mat_a)
%% Reduced Row Echelon Rorm
a=rref(mat_a)
%% (b)
clear all
mat_b=[1,1,-2,0,0;1,0,-1,-1,1;1,0,0,-2,0;0,1,0,0,-2];
%% Rank of the matrix
r_b=rank(mat_b)
%% Length of the matrix
l_b=length(mat_b)
%% Reduced Row Echelon Rorm
b=rref(mat_b)
%% (c)
clear all
mat_c=[1,1,-1,0,0,0,0,0,0,0;0,0,1,1,-1,-1,0,0,0,0;0,0,0,0,0,1,-1,-2,-1,0;0,0,0,1,-1,0,0,0,1,-1;1,0,0,1,0,1,0,0,0,0]
%% Rank of the matrix
r_c=rank(mat_c)
%% Length of the matrix
l_c=length(mat_c)
%% Reduced Row Echelon Rorm
c=rref(mat_c)
%% Question 5
clear all
%% (a) Representation of the system of equation in matrix form 
a=[100,125,125,62.5;80,110,120,25;140,80,120,100;90,104.8,60,137.33]
b=[6625;5290;7300;6539]
%% (b) Determinant of the system of matrices
d=det(a)
%% (c) Roots of the system of equation
x=a\b
%% Question 6
clear all
%% (a) Representation of the system of equation in matrix form 
a=[0.07,0.18,0.15,0.24;0.04,0.24,0.1,0.65;0.54,0.42,0.54,0.1;0.35,0.16,0.21,0.01]
b=70*[0.15;0.25;0.4;0.2]
%% (b) Determinant of the system of matrices
d=det(a)
%% (c) Roots of the system of equation
x=a\b
%% Question 3
clear all
% Number of reactions can be calculated from the length and rank of matrix
r=[2,0,2,1,1,0,2,3,4;6,1,5,3,4,2,4,8,10];
% rank of the matrix
r_r=rank(r);
% length of matrix
l_r=length(r);
%number of reactions
n=l_r-r_r;
%% Question 1
%% (a)
clear all
syms nmax
% values of variable "x":
x=linspace(-1,1,101);
%Summation terms:
f=((cospi(nmax)-1)/(nmax*pi)^2*cospi(nmax*x))+(cospi(nmax))/(nmax*pi)*sinpi(nmax*x);
%% For nmax = 10
y=1/4+vpasum(f,nmax,1,10);
%% For nmax = 100
y2=1/4+vpasum(f,nmax,1,100);
%% For nmax = 1000
y3=1/4+vpasum(f,nmax,1,1000);
%% Plotting the Curves:
plot(x,y,"m",x,y2,"r",x,y3,"b") 
xlabel("x");
ylabel("y(x)");
legend ("nmax=10","nmax=100","nmax=1000"); 
title("Plot of nmax");
grid on;
%% (b)
%values of yn
sfn=0;
for i=1:1000
    f1=((cospi(nmax)-1)/(nmax*pi)^2*cospi(nmax*x))+(cospi(nmax))/(nmax*pi)*sinpi(nmax*x);
    yn(i)=1/4+vpasum(f1,nmax,1,i);
end
for l=1:length(linspace(-1,1,101))
%values of "y*" (named as ys)
    if -1<=x(l)<=0
        v=-x(l);
        f1=((cospi(nmax)-1)/(nmax*pi)^2*cospi(nmax*v))+(cospi(nmax))/(nmax*pi)*sinpi(nmax*v);
        ys=1/4+vpasum(f1,nmax,1,1000);
    else
        x=0;
        f2=((cospi(nmax)-1)/(nmax*pi)^2*cospi(nmax*x))+(cospi(nmax))/(nmax*pi)*sinpi(nmax*x);
        ys=1/4+vpasum(f2,nmax,1,1000);
    end
%     E(nmax)=norm((yn-ys))
end
%E=norm((yn-ys));

