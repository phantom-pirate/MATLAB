%% Homework 4
%% Name and Person ID : Sakthi Aadharsh Azhagar Gobinath Kavitha [50409103]
clc
%% Question 1
%% (a)
clear all
%% input the function
syms x 
f=(x-1)^2*log(x);
x0=1;
%% taking the second derivative
h=matlabFunction(hessian(f),"Var",{x});
%% Second derivative test  
D=eig(h(x0))
%% Classification
if D>0
    disp("Local Maximum");
elseif D<0
    disp("Local Minimum");
else
    %% Attempting first derivative test
    g=matlabFunction(gradient(f),"Var",{x});
    if (sign(g(x0+0000005))==sign(g(x0-0.0000005)))
         disp("Inflection Point");
    elseif sign(g(x0-0.5))==-1
        if sign(g(x0+0.5))==1
             disp("Local Minimum");
        else
             disp("Local Maximum");        
        end
    end
end
%% Plotting the graph
fplot(f,[0.5,1.5])
%% (b)
clear all
%% input the function
syms x 
f=3*(x-1)^4+5;
x0=1;
%% taking the second derivative
h=matlabFunction(hessian(f),"Var",{x});
%% Second derivative test  
D=eig(h(x0))
%% Classification
if D>0
    disp("Local Maximum");
elseif D<0
    disp("Local Minimum");
else
    %% Attempting first derivative test
    g=matlabFunction(gradient(f),"Var",{x});
    if (sign(g(x0+0.000005))==sign(g(x0-0.0000005)))
         disp("Inflection Point");
    elseif sign(g(x0-0.5))==-1
        if sign(g(x0+0.5))==1
             disp("Local Minimum");
        else
             disp("Local Maximum");        
        end
    end
end
%% Plotting the graph
fplot(f,[0.5,1.5])
%% (c)
clear all
%% input the function
syms x 
f=(x-1)^2*(x-2)*(x);
x0=1;
%% taking the second derivative
h=matlabFunction(hessian(f),"Var",{x});
%% Second derivative test  
D=eig(h(x0))
%% Classification
if D>0
    disp("Local Maximum");
elseif D<0
    disp("Local Minimum");
else
    %% Attempting first derivative test
    g=matlabFunction(gradient(f),"Var",{x});
    if (sign(g(x0+0000005))==sign(g(x0-0.0000005)))
         disp("Inflection Point");
    elseif sign(g(x0-0.5))==-1
        if sign(g(x0+0.5))==1
             disp("Local Minimum");
        else
             disp("Local Maximum");        
        end
    end
end
%% Plotting the graph
fplot(f,[0.5,1.5])
%% (d)
clear all
%% input the function
syms x 
f=(x+1)*(x-3)*(1-x);
x0=1;
%% taking the second derivative
h=matlabFunction(hessian(f),"Var",{x});
%% Second derivative test  
D=eig(h(x0))
%% Classification
if D>0
    disp("Local Maximum");
elseif D<0
    disp("Local Minimum");
else
    %% Attempting first derivative test
    g=matlabFunction(gradient(f),"Var",{x});
    if (sign(g(x0+0.000005))==sign(g(x0-0.0000005)))
         disp("Inflection Point");
    elseif sign(g(x0-0.5))==-1
        if sign(g(x0+0.5))==1
             disp("Local Minimum");
        else
             disp("Local Maximum");        
        end
    end
end
%% Plotting the graph
fplot(f,[0.5,1.5])
%% (e)
clear all
%% input the function
syms x 
f=sin(5*(x-1)^4);
x0=1;
%% taking the second derivative
h=matlabFunction(hessian(f),"Var",{x});
%% Second derivative test  
D=eig(h(x0))
%% Classification
if D>0
    disp("Local Maximum");
elseif D<0
    disp("Local Minimum");
else
    %% Attempting first derivative test
    g=matlabFunction(gradient(f),"Var",{x});
    if (sign(g(x0+0000005))==sign(g(x0-0.0000005)))
         disp("Inflection Point");
    elseif sign(g(x0-0.5))==-1
        if sign(g(x0+0.5))==1
             disp("Local Minimum");
        else
             disp("Local Maximum");        
        end
    end
end
%% Plotting the graph
fplot(f,[0.5,1.5])
%% (f)
clear all
%% input the function
syms x 
f=exp(8*(x-1)^5);
x0=1;
%% taking the second derivative
h=matlabFunction(hessian(f),"Var",{x});
%% Second derivative test  
D=eig(h(x0))
%% Classification
if D>0
    disp("Local Maximum");
elseif D<0
    disp("Local Minimum");
else
    %% Attempting first derivative test
    g=matlabFunction(gradient(f),"Var",{x});
    if (sign(g(x0+0.000005))==sign(g(x0-0.0000005)))
         disp("Inflection Point");
    elseif sign(g(x0-0.5))==-1
        if sign(g(x0+0.5))==1
             disp("Local Minimum");
        else
             disp("Local Maximum");        
        end
    end
end
%% Plotting the graph
fplot(f,[0.5,1.5])
%% (g)
clear all
%% input the function
syms x 
f=(1-x)*sin((x^2-1)^3);
x0=1;
%% taking the second derivative
h=matlabFunction(hessian(f),"Var",{x});
%% Second derivative test  
D=eig(h(x0))
%% Classification
if D>0
    disp("Local Maximum");
elseif D<0
    disp("Local Minimum");
else
    %% Attempting first derivative test
    g=matlabFunction(gradient(f),"Var",{x});
    if (sign(g(x0+0000005))==sign(g(x0-0.0000005)))
         disp("Inflection Point");
    elseif sign(g(x0-0.5))==-1
        if sign(g(x0+0.5))==1
             disp("Local Minimum");
        else
             disp("Local Maximum");        
        end
    end
end
%% Plotting the graph
fplot(f,[0.5,1.5])
%% (h)
clear all
%% input the function
syms x 
f=exp(-(log(x))^3);
x0=1;
%% taking the second derivative
h=matlabFunction(hessian(f),"Var",{x});
%% Second derivative test  
D=eig(h(x0))
%% Classification
if D>0
    disp("Local Maximum");
elseif D<0
    disp("Local Minimum");
else
    %% Attempting first derivative test
    g=matlabFunction(gradient(f),"Var",{x});
    if (sign(g(x0+0.000005))==sign(g(x0-0.0000005)))
         disp("Inflection Point");
    elseif sign(g(x0-0.5))==-1
        if sign(g(x0+0.5))==1
             disp("Local Minimum");
        else
             disp("Local Maximum");        
        end
    end
end
%% Plotting the graph
fplot(f,[0.5,1.5])

%% Question 2

%% (a)
clear all
%% Input the function
syms x
f=exp(-x^3);
fplot(f,[-3,-2])
%% finding the critical points
g=matlabFunction(gradient(f),"vars",{x})
r=fzero(g,0)
%% Second Derivative test
h=matlabFunction(hessian(f),"vars",{x});
D=eig(h(r))
if D>0
    disp("Local Maximum");
elseif D<0
    disp("Local Minimum");
else
    if (sign(g(r+0.000005))==sign(g(r-0.0000005)))
         disp("Inflection Point");
    elseif sign(g(r-0.0000005))==-1
        if sign(g(r+0.0000005))==1
             disp("Local Minimum");
        else
             disp("Local Maximum");        
        end
    end
end
%% plot
s=linspace(-1.5,0.5,201);
for i=1:length(s)      
    fs(i)=exp(-1*s(i)^3);
end
for i=1:length(s)      
    gs(i)=g(s(i)); 
end
plot(s,fs,"-",s,gs,"--")
legend("Function","First Derivative")
title("Plot of Functiona and its Derivtive")

%% (b)
clear all
%% Input the function
syms x
f=1/(x^2-4*x+5);
fplot(f,[-3,3])
%% finding the critical points
g=matlabFunction(gradient(f),"vars",{x})
r=fzero(g,0)
%% Second Derivative test
h=matlabFunction(hessian(f),"vars",{x});
D=eig(h(r))
if D>0
    disp("Local Maximum");
elseif D<0
    disp("Local Minimum");
else
    if (sign(g(r+0.000005))==sign(g(r-0.0000005)))
         disp("Inflection Point");
    elseif sign(g(r-0.0000005))==-1
        if sign(g(r+0.0000005))==1
             disp("Local Minimum");
        else
             disp("Local Maximum");        
        end
    end
end
%% plot
s=linspace(-3,3,201);
for i=1:length(s)      
    fs(i)=exp(-1*s(i)^3);
end
for i=1:length(s)      
    gs(i)=g(s(i)); 
end
plot(s,fs,"-",s,gs,"--")
legend("Function","First Derivative")
title("Plot of Functiona and its Derivtive")

%% (c)
clear all
%% Input the function
syms x
f=exp(12*x-x^3);
fplot(f,[0,4])
%% finding the critical points
g=matlabFunction(gradient(f),"vars",{x})
r=fzero(g,0)
%% Second Derivative test
h=matlabFunction(hessian(f),"vars",{x});
D=eig(h(r))
if D>0
    disp("Local Maximum");
elseif D<0
    disp("Local Minimum");
else
    if (sign(g(r+0.000005))==sign(g(r-0.0000005)))
         disp("Inflection Point");
    elseif sign(g(r-0.0000005))==-1
        if sign(g(r+0.0000005))==1
             disp("Local Minimum");
        else
             disp("Local Maximum");        
        end
    end
end
%% plot
s=linspace(0,4,201);
for i=1:length(s)      
    fs(i)=exp(-1*s(i)^3);
end
for i=1:length(s)      
    gs(i)=g(s(i)); 
end
plot(s,fs,"-",s,gs,"--")
legend("Function","First Derivative")
title("Plot of Functiona and its Derivtive")

%% (d)
clear all
%% Input the function
syms x
f=-3*(log(x)^3);
fplot(f,[0,3])
%% finding the critical points
g=matlabFunction(gradient(f),"vars",{x})
r=fsolve(g,0.25)
%% Second Derivative test
h=matlabFunction(hessian(f),"vars",{x});
D=eig(h(r))
if D>0
    disp("Local Maximum");
elseif D<0
    disp("Local Minimum");
else
    if (sign(g(r+0.000005))==sign(g(r-0.0000005)))
         disp("Inflection Point");
    elseif sign(g(r-0.0000005))==-1
        if sign(g(r+0.0000005))==1
             disp("Local Minimum");
        else
             disp("Local Maximum");        
        end
    end
end
%% plot
s=linspace(0,3,201);
for i=1:length(s)      
    fs(i)=exp(-1*s(i)^3);
end
for i=1:length(s)      
    gs(i)=g(s(i)); 
end
plot(s,fs,"-",s,gs,"--")
legend("Function","First Derivative")
title("Plot of Functiona and its Derivtive")

clc

%% (e)
clear all
%% Input the function
syms x
f=log(x)^4;
fplot(f,[0,4])
%% finding the critical points
g=matlabFunction(gradient(f),"vars",{x})
r=fsolve(g,0.25)
%% Second Derivative test
h=matlabFunction(hessian(f),"vars",{x});
D=eig(h(r))
if D>0
    disp("Local Maximum");
elseif D<0
    disp("Local Minimum");
else
    if (sign(g(r+0.000005))==sign(g(r-0.0000005)))
         disp("Inflection Point");
    elseif sign(g(r-0.0000005))==-1
        if sign(g(r+0.0000005))==1
             disp("Local Minimum");
        else
             disp("Local Maximum");        
        end
    end
end
%% plot
s=linspace(0,4,201);
for i=1:length(s)      
    fs(i)=exp(-1*s(i)^3);
end
for i=1:length(s)      
    gs(i)=g(s(i)); 
end
plot(s,fs,"-",s,gs,"--")
legend("Function","First Derivative")
title("Plot of Functiona and its Derivtive")


%% (f)
clear all
%% Input the function
syms x
f=exp(-sin(x));
fplot(f,[-1,5])
%% finding the critical points
g=matlabFunction(gradient(f),"vars",{x})
r=fzero(g,0)
%% Second Derivative test
h=matlabFunction(hessian(f),"vars",{x});
D=eig(h(r))
if D>0
    disp("Local Maximum");
elseif D<0
    disp("Local Minimum");
else
    if (sign(g(r+0.000005))==sign(g(r-0.0000005)))
         disp("Inflection Point");
    elseif sign(g(r-0.0000005))==-1
        if sign(g(r+0.0000005))==1
             disp("Local Minimum");
        else
             disp("Local Maximum");        
        end
    end
end
%% plot
s=linspace(-1,5,201);
for i=1:length(s)      
    fs(i)=exp(-1*s(i)^3);
end
for i=1:length(s)      
    gs(i)=g(s(i)); 
end
plot(s,fs,"-",s,gs,"--")
legend("Function","First Derivative")
title("Plot of Functiona and its Derivtive")

%% (g)
clear all
%% Input the function
syms x
f=x^2*exp(-x);
fplot(f,[-3,-2])
%% finding the critical points
g=matlabFunction(gradient(f),"vars",{x})
r=fzero(g,0)
%% Second Derivative test
h=matlabFunction(hessian(f),"vars",{x});
D=eig(h(r))
if D>0
    disp("Local Maximum");
elseif D<0
    disp("Local Minimum");
else
    if (sign(g(r+0.000005))==sign(g(r-0.0000005)))
         disp("Inflection Point");
    elseif sign(g(r-0.0000005))==-1
        if sign(g(r+0.0000005))==1
             disp("Local Minimum");
        else
             disp("Local Maximum");        
        end
    end
end
%% plot
s=linspace(-1.5,0.5,201);
for i=1:length(s)      
    fs(i)=exp(-1*s(i)^3);
end
for i=1:length(s)      
    gs(i)=g(s(i)); 
end
plot(s,fs,"-",s,gs,"--")
legend("Function","First Derivative")
title("Plot of Functiona and its Derivtive")
%% (h)
clear all
%% Input the function
syms x
f=exp(4*sqrt(x)-x^2);
fplot(f,[0,5])
%% finding the critical points
g=matlabFunction(gradient(f),"vars",{x})
r=fzero(g,1)
%% Second Derivative test
h=matlabFunction(hessian(f),"vars",{x});
D=eig(h(r))
if D>0
    disp("Local Maximum");
elseif D<0
    disp("Local Minimum");
else
    if (sign(g(r+0.000005))==sign(g(r-0.0000005)))
         disp("Inflection Point");
    elseif sign(g(r-0.0000005))==-1
        if sign(g(r+0.0000005))==1
             disp("Local Minimum");
        else
             disp("Local Maximum");        
        end
    end
end
%% plot
s=linspace(-1.5,0.5,201);
for i=1:length(s)      
    fs(i)=exp(-1*s(i)^3);
end
for i=1:length(s)      
    gs(i)=g(s(i)); 
end
plot(s,fs,"-",s,gs,"--")
legend("Function","First Derivative")
title("Plot of Functiona and its Derivtive")

%% Question 3
%% a
clear all
%% input the function
syms x y
f=2*x^2+x*y+y^2+7*y+8;
%% Find the critical points
l=matlabFunction(gradient(f),"Vars",{[x y]});
x0=[0 0];
r=fsolve(l,x0)
%% Second Derivative Test
h=matlabFunction(hessian(f),"Vars",{[x y]});
E=eig(h(r))
%% Classification
fxx=matlabFunction(diff(f,x,2),"Vars",{[x y]});
fxx=fxx(x0);
fyy=matlabFunction(diff(f,y,2),"Vars",{[x y]});
fyy=fyy(x0);
fxy=matlabFunction(diff(diff(f,x),y),"Vars",{[x y]});
fxy=fxy(x0);
D=fxx*fyy-fxy^2;
if D>0&fxx<0
    disp("Local Minimum");
else if D>0&fxx>0
         disp("Local Maximum");
    else D<0
        disp("saddle point")
    end
end

%% b
clear all
%% input the function
syms x y
f=x^2+x*y+y^2+6*y;
%% Find the critical points
l=matlabFunction(gradient(f),"Vars",{[x y]});
x0=[0 0];
r=fsolve(l,x0)
%% Second Derivative Test
h=matlabFunction(hessian(f),"Vars",{[x y]});
E=eig(h(r))
%% Classification
fxx=matlabFunction(diff(f,x,2),"Vars",{[x y]});
fxx=fxx(x0);
fyy=matlabFunction(diff(f,y,2),"Vars",{[x y]});
fyy=fyy(x0);
fxy=matlabFunction(diff(diff(f,x),y),"Vars",{[x y]});
fxy=fxy(x0);
D=fxx*fyy-fxy^2;
if D>0&fxx<0
    disp("Local Minimum");
else if D>0&fxx>0
         disp("Local Maximum");
    else D<0
        disp("saddle point")
    end
end

%% c
clear all
%% input the function
syms x y
f=x^2+4*x*y+y^2+6*x+8;
%% Find the critical points
l=matlabFunction(gradient(f),"Vars",{[x y]});
x0=[0 0];
r=fsolve(l,x0)
%% Second Derivative Test
h=matlabFunction(hessian(f),"Vars",{[x y]});
E=eig(h(r))
%% Classification
fxx=matlabFunction(diff(f,x,2),"Vars",{[x y]});
fxx=fxx(x0);
fyy=matlabFunction(diff(f,y,2),"Vars",{[x y]});
fyy=fyy(x0);
fxy=matlabFunction(diff(diff(f,x),y),"Vars",{[x y]});
fxy=fxy(x0);
D=fxx*fyy-fxy^2;
if D>0&fxx<0
    disp("Local Minimum");
else if D>0&fxx>0
         disp("Local Maximum");
    else D<0
        disp("saddle point")
    end
end
%% d 
clear all
%% input the function
syms x y
f=-6*x*y+x^5;
%% Find the critical points
l=matlabFunction(gradient(f),"Vars",{[x y]});
x0=[0 0];
r=fsolve(l,x0)
%% Second Derivative Test
h=matlabFunction(hessian(f),"Vars",{[x y]});
E=eig(h(r))
%% Classification
fxx=matlabFunction(diff(f,x,2),"Vars",{[x y]});
fxx=fxx(x0);
fyy=matlabFunction(diff(f,y,2),"Vars",{[x y]});
fyy=fyy(x0);
fxy=matlabFunction(diff(diff(f,x),y),"Vars",{[x y]});
fxy=fxy(x0);
D=fxx*fyy-fxy^2;
if D>0&fxx<0
    disp("Local Minimum");
else if D>0&fxx>0
         disp("Local Maximum");
    else D<0
        disp("saddle point")
    end
end

%% e
clear all
%% input the function
syms x y
f=(2*x-y)*x;
%% Find the critical points
l=matlabFunction(gradient(f),"Vars",{[x y]});
x0=[0 0];
r=fsolve(l,x0)
%% Second Derivative Test
h=matlabFunction(hessian(f),"Vars",{[x y]});
E=eig(h(r))
%% Classification
fxx=matlabFunction(diff(f,x,2),"Vars",{[x y]});
fxx=fxx(x0);
fyy=matlabFunction(diff(f,y,2),"Vars",{[x y]});
fyy=fyy(x0);
fxy=matlabFunction(diff(diff(f,x),y),"Vars",{[x y]});
fxy=fxy(x0);
D=fxx*fyy-fxy^2;
if D>0&fxx<0
    disp("Local Minimum");
else if D>0&fxx>0
         disp("Local Maximum");
    else D<0
        disp("saddle point")
    end
end
%% f
clear all
%% input the function
syms x y
f=exp(-(x^2+y^2+1));
%% Find the critical points
l=matlabFunction(gradient(f),"Vars",{[x y]});
x0=[0 0];
r=fsolve(l,x0)
%% Second Derivative Test
h=matlabFunction(hessian(f),"Vars",{[x y]});
E=eig(h(r))
%% Classification
fxx=matlabFunction(diff(f,x,2),"Vars",{[x y]});
fxx=fxx(x0);
fyy=matlabFunction(diff(f,y,2),"Vars",{[x y]});
fyy=fyy(x0);
fxy=matlabFunction(diff(diff(f,x),y),"Vars",{[x y]});
fxy=fxy(x0);
D=fxx*fyy-fxy^2;
if D>0&fxx<0
    disp("Local Minimum");
else if D>0&fxx>0
         disp("Local Maximum");
    else D<0
        disp("saddle point")
    end
end

%% g
clear all
%% input the function
syms x y
f=log(x^2+x*y+y^2+4);
%% Find the critical points
l=matlabFunction(gradient(f),"Vars",{[x y]});
x0=[0 0];
r=fsolve(l,x0)
%% Second Derivative Test
h=matlabFunction(hessian(f),"Vars",{[x y]});
E=eig(h(r))
%% Classification
fxx=matlabFunction(diff(f,x,2),"Vars",{[x y]});
fxx=fxx(x0);
fyy=matlabFunction(diff(f,y,2),"Vars",{[x y]});
fyy=fyy(x0);
fxy=matlabFunction(diff(diff(f,x),y),"Vars",{[x y]});
fxy=fxy(x0);
D=fxx*fyy-fxy^2;
if D>0&fxx<0
    disp("Local Minimum");
else if D>0&fxx>0
         disp("Local Maximum");
    else D<0
        disp("saddle point")
    end
end
%% h
clear all
%% input the function
syms x y
f=exp(x^2-y^2);
%% Find the critical points
l=matlabFunction(gradient(f),"Vars",{[x y]});
x0=[0 0];
r=fsolve(l,x0)
%% Second Derivative Test
h=matlabFunction(hessian(f),"Vars",{[x y]});
E=eig(h(r))
%% Classification
fxx=matlabFunction(diff(f,x,2),"Vars",{[x y]});
fxx=fxx(x0);
fyy=matlabFunction(diff(f,y,2),"Vars",{[x y]});
fyy=fyy(x0);
fxy=matlabFunction(diff(diff(f,x),y),"Vars",{[x y]});
fxy=fxy(x0);
D=fxx*fyy-fxy^2;
if D>0&fxx<0
    disp("Local Minimum");
else if D>0&fxx>0
         disp("Local Maximum");
    else D<0
        disp("saddle point")
    end
end
%% i
clear all
%% input the function
syms x y
f=2*x*y-sin(x*y);
%% Find the critical points
l=matlabFunction(gradient(f),"Vars",{[x y]});
x0=[0 0];
r=fsolve(l,x0)
%% Second Derivative Test
h=matlabFunction(hessian(f),"Vars",{[x y]});
E=eig(h(r))
%% Classification
fxx=matlabFunction(diff(f,x,2),"Vars",{[x y]});
fxx=fxx(x0);
fyy=matlabFunction(diff(f,y,2),"Vars",{[x y]});
fyy=fyy(x0);
fxy=matlabFunction(diff(diff(f,x),y),"Vars",{[x y]});
fxy=fxy(x0);
D=fxx*fyy-fxy^2;
if D>0&fxx<0
    disp("Local Minimum");
else if D>0&fxx>0
         disp("Local Maximum");
    else D<0
        disp("saddle point")
    end
end
%% j
clear all
%% input the function
syms x y
f=x^3-x*y+x^5;
%% Find the critical points
l=matlabFunction(gradient(f),"Vars",{[x y]});
x0=[0 0];
r=fsolve(l,x0)
%% Second Derivative Test
h=matlabFunction(hessian(f),"Vars",{[x y]});
E=eig(h(r))
%% Classification
fxx=matlabFunction(diff(f,x,2),"Vars",{[x y]});
fxx=fxx(x0);
fyy=matlabFunction(diff(f,y,2),"Vars",{[x y]});
fyy=fyy(x0);
fxy=matlabFunction(diff(diff(f,x),y),"Vars",{[x y]});
fxy=fxy(x0);
D=fxx*fyy-fxy^2;
if D>0&fxx<0
    disp("Local Minimum");
else if D>0&fxx>0
         disp("Local Maximum");
    else D<0
        disp("saddle point")
    end
end
%% k
clear all
%% input the function
syms x y
f=log(1+x^2+y^2);
%% Find the critical points
l=matlabFunction(gradient(f),"Vars",{[x y]});
x0=[0 0];
r=fsolve(l,x0)
%% Second Derivative Test
h=matlabFunction(hessian(f),"Vars",{[x y]});
E=eig(h(r))
%% Classification
fxx=matlabFunction(diff(f,x,2),"Vars",{[x y]});
fxx=fxx(x0);
fyy=matlabFunction(diff(f,y,2),"Vars",{[x y]});
fyy=fyy(x0);
fxy=matlabFunction(diff(diff(f,x),y),"Vars",{[x y]});
fxy=fxy(x0);
D=fxx*fyy-fxy^2;
if D>0&fxx<0
    disp("Local Minimum");
else if D>0&fxx>0
         disp("Local Maximum");
    else D<0
        disp("saddle point");
    end
end
%% l
clear all
%% input the function
syms x y
f=x^3-x*y+x^2*y-2;
%% Find the critical points
l=matlabFunction(gradient(f),"Vars",{[x y]});
x0=[0 0];
r=fsolve(l,x0)
%% Second Derivative Test
h=matlabFunction(hessian(f),"Vars",{[x y]});
E=eig(h(r))
%% Classification
fxx=matlabFunction(diff(f,x,2),"Vars",{[x y]});
fxx=fxx(x0);
fyy=matlabFunction(diff(f,y,2),"Vars",{[x y]});
fyy=fyy(x0);
fxy=matlabFunction(diff(diff(f,x),y),"Vars",{[x y]});
fxy=fxy(x0);
D=fxx*fyy-fxy^2;
if D>0&fxx<0
    disp("Local Minimum");
else if D>0&fxx>0
         disp("Local Maximum");
    else D<0
        disp("saddle point");
    end
end
%% m
clear all
%% input the function
syms x y
f=1/(x*y+y^2+y+1);
%% Find the critical points
l=matlabFunction(gradient(f),"Vars",{[x y]});
x0=[0 0];
r=fsolve(l,x0)
%% Second Derivative Test
h=matlabFunction(hessian(f),"Vars",{[x y]});
E=eig(h(r))
%% Classification
fxx=matlabFunction(diff(f,x,2),"Vars",{[x y]});
fxx=fxx(x0);
fyy=matlabFunction(diff(f,y,2),"Vars",{[x y]});
fyy=fyy(x0);
fxy=matlabFunction(diff(diff(f,x),y),"Vars",{[x y]});
fxy=fxy(x0);
D=fxx*fyy-fxy^2;
if D>0&fxx<0
    disp("Local Minimum");
else if D>0&fxx>0
         disp("Local Maximum");
    else D<0
        disp("saddle point");
    end
end
%% n
clear all
%% input the function
syms x y
f=2*x*y+x^8;
%% Find the critical points
l=matlabFunction(gradient(f),"Vars",{[x y]});
x0=[0 0];
r=fsolve(l,x0)
%% Second Derivative Test
h=matlabFunction(hessian(f),"Vars",{[x y]});
E=eig(h(r))
%% Classification
fxx=matlabFunction(diff(f,x,2),"Vars",{[x y]});
fxx=fxx(x0);
fyy=matlabFunction(diff(f,y,2),"Vars",{[x y]});
fyy=fyy(x0);
fxy=matlabFunction(diff(diff(f,x),y),"Vars",{[x y]});
fxy=fxy(x0);
D=fxx*fyy-fxy^2;
if D>0&fxx<0
    disp("Local Minimum");
else if D>0&fxx>0
         disp("Local Maximum");
    else D<0
        disp("saddle point");
    end
end

%% Question 4
%% a
clear all
%% Input the function
syms x y z
f=x^2+z^2+3*x*y;
%% finding critical points
l=matlabFunction(gradient(f),"Vars",{[x y z]});
x0=[0 0 0];
r=fsolve(l,x0)
h=matlabFunction(hessian(f),"Vars",{[x y z]});
E=eig(h(r))
%% Classification
fxx=matlabFunction(diff(f,x,2),"Vars",{[x y z]});
fxx=fxx(x0);
fyy=matlabFunction(diff(f,y,2),"Vars",{[x y z]});
fyy=fyy(x0);
fzz=matlabFunction(diff(f,y,2),"Vars",{[x y z]});
fzz=fzz(x0);
fxy=matlabFunction(diff(diff(f,x),y),"Vars",{[x y z]});
fxy=fxy(x0);
fyz=matlabFunction(diff(diff(f,y),z),"Vars",{[x y z]});
fyz=fyz(x0);
fxz=matlabFunction(diff(diff(f,x),z),"Vars",{[x y z]});
fxz=fxz(x0);
m=[fxx fxy fxz;fxy fyy fyz;fxz fyz fzz];
D=det(m);
if D>0&&fxx<0
    disp("Local Minimum");
else if D>0&&fxx>0
         disp("Local Maximum");
    else D<0
        disp("saddle point")
    end
end
%% b
clear all
syms x y z
%% Input the function
f=exp(2*x^2+z*x+3-5*z^2);
%% finding critical points
l=matlabFunction(gradient(f),"Vars",{[x y z]});
x0=[0 0 0];
r=fsolve(l,x0)
%% Second Derivative Test
h=matlabFunction(hessian(f),"Vars",{[x y z]});
E=eig(h(r))
%% Classification
fxx=matlabFunction(diff(f,x,2),"Vars",{[x y z]});
fxx=fxx(x0);
fyy=matlabFunction(diff(f,y,2),"Vars",{[x y z]});
fyy=fyy(x0);
fzz=matlabFunction(diff(f,y,2),"Vars",{[x y z]});
fzz=fzz(x0);
fxy=matlabFunction(diff(diff(f,x),y),"Vars",{[x y z]});
fxy=fxy(x0);
fyz=matlabFunction(diff(diff(f,y),z),"Vars",{[x y z]});
fyz=fyz(x0);
fxz=matlabFunction(diff(diff(f,x),z),"Vars",{[x y z]});
fxz=fxz(x0);
m=[fxx fxy fxz;fxy fyy fyz;fxz fyz fzz];
D=det(m);
if D>0&&fxx<0
    disp("Local Minimum");
else if D>0&&fxx>0
         disp("Local Maximum");
    else D<0
        disp("saddle point")
    end
end
%% c
clear all
%% Input the function
syms x y z
f=x^2+y^2+z^2+x*y+x*z+y*z;
%% finding critical points
l=matlabFunction(gradient(f),"Vars",{[x y z]});
x0=[0 0 0];
r=fsolve(l,x0)
%% Second Derivative Test
h=matlabFunction(hessian(f),"Vars",{[x y z]});
E=eig(h(r))
%% Classification
fxx=matlabFunction(diff(f,x,2),"Vars",{[x y z]});
fxx=fxx(x0);
fyy=matlabFunction(diff(f,y,2),"Vars",{[x y z]});
fyy=fyy(x0);
fzz=matlabFunction(diff(f,y,2),"Vars",{[x y z]});
fzz=fzz(x0);
fxy=matlabFunction(diff(diff(f,x),y),"Vars",{[x y z]});
fxy=fxy(x0);
fyz=matlabFunction(diff(diff(f,y),z),"Vars",{[x y z]});
fyz=fyz(x0);
fxz=matlabFunction(diff(diff(f,x),z),"Vars",{[x y z]});
fxz=fxz(x0);
m=[fxx fxy fxz;fxy fyy fyz;fxz fyz fzz];
D=det(m);
if D>0&&fxx<0
    disp("Local Minimum");
else if D>0&&fxx>0
         disp("Local Maximum");
    else D<0
        disp("saddle point")
    end
end
%% d
clear all
%% Input the function
syms x y z
f=log(x^2+y^2+z^2+1);
l=matlabFunction(gradient(f),"Vars",{[x y z]});
x0=[0 0 0];
r=fsolve(l,x0)
%% Second Derivative Test
h=matlabFunction(hessian(f),"Vars",{[x y z]});
E=eig(h(r))
%% Classification
fxx=matlabFunction(diff(f,x,2),"Vars",{[x y z]});
fxx=fxx(x0);
fyy=matlabFunction(diff(f,y,2),"Vars",{[x y z]});
fyy=fyy(x0);
fzz=matlabFunction(diff(f,y,2),"Vars",{[x y z]});
fzz=fzz(x0);
fxy=matlabFunction(diff(diff(f,x),y),"Vars",{[x y z]});
fxy=fxy(x0);
fyz=matlabFunction(diff(diff(f,y),z),"Vars",{[x y z]});
fyz=fyz(x0);
fxz=matlabFunction(diff(diff(f,x),z),"Vars",{[x y z]});
fxz=fxz(x0);
m=[fxx fxy fxz;fxy fyy fyz;fxz fyz fzz];
D=det(m);
if D>0&&fxx<0
    disp("Local Minimum");
else if D>0&&fxx>0
         disp("Local Maximum");
    else D<0
        disp("saddle point")
    end
end
%% f
clear all
%% Input the function
syms w x y z
f=exp(w^2+z^2-x*w-y*z);
%% finding critical points
l=matlabFunction(gradient(f),"Vars",{[w x y z]});
x0=[0 0 0 0];
r=fsolve(l,x0)
%% Second Derivative Test
h=matlabFunction(hessian(f),"Vars",{[w x y z]});
E=eig(h(r))
%% Classification
fww=matlabFunction(diff(f,w,2),"Vars",{[w x y z]});
fww=fww(x0);
fxx=matlabFunction(diff(f,x,2),"Vars",{[w x y z]});
fxx=fxx(x0);
fyy=matlabFunction(diff(f,y,2),"Vars",{[w x y z]});
fyy=fyy(x0);
fzz=matlabFunction(diff(f,y,2),"Vars",{[w x y z]});
fzz=fzz(x0);
fxy=matlabFunction(diff(diff(f,x),y),"Vars",{[w x y z]});
fxy=fxy(x0);
fyz=matlabFunction(diff(diff(f,y),z),"Vars",{[w x y z]});
fyz=fyz(x0);
fxz=matlabFunction(diff(diff(f,x),z),"Vars",{[w x y z]});
fxz=fxz(x0);
fwx=matlabFunction(diff(diff(f,x),w),"Vars",{[w x y z]});
fwx=fwx(x0);
fwy=matlabFunction(diff(diff(f,y),w),"Vars",{[w x y z]});
fwy=fwy(x0);
fwz=matlabFunction(diff(diff(f,z),w),"Vars",{[w x y z]});
fwz=fwz(x0);
m=[fww fwx fwy fwz;fwx fxx fxy fxz;fwy fxy fyy fyz;fwz fxz fyz fzz];
D=det(m);
if D>0&&fxx<0
    disp("Local Minimum");
else if D>0&&fxx>0
         disp("Local Maximum");
    else D<0
        disp("saddle point")
    end
end
%% g
clear all
%% Input the function
syms w x y z
f=7-4*(w^2+x^2+y^2+z^2)+4*x*w+2*w*z;
%% finding critical points
l=matlabFunction(gradient(f),"Vars",{[w x y z]});
x0=[0 0 0 0];
r=fsolve(l,x0)
%% Second Derivative Test
h=matlabFunction(hessian(f),"Vars",{[w x y z]});
E=eig(h(r))
%% Classification
fww=matlabFunction(diff(f,w,2),"Vars",{[w x y z]});
fww=fww(x0);
fxx=matlabFunction(diff(f,x,2),"Vars",{[w x y z]});
fxx=fxx(x0);
fyy=matlabFunction(diff(f,y,2),"Vars",{[w x y z]});
fyy=fyy(x0);
fzz=matlabFunction(diff(f,y,2),"Vars",{[w x y z]});
fzz=fzz(x0);
fxy=matlabFunction(diff(diff(f,x),y),"Vars",{[w x y z]});
fxy=fxy(x0);
fyz=matlabFunction(diff(diff(f,y),z),"Vars",{[w x y z]});
fyz=fyz(x0);
fxz=matlabFunction(diff(diff(f,x),z),"Vars",{[w x y z]});
fxz=fxz(x0);
fwx=matlabFunction(diff(diff(f,x),w),"Vars",{[w x y z]});
fwx=fwx(x0);
fwy=matlabFunction(diff(diff(f,y),w),"Vars",{[w x y z]});
fwy=fwy(x0);
fwz=matlabFunction(diff(diff(f,z),w),"Vars",{[w x y z]});
fwz=fwz(x0);
m=[fww fwx fwy fwz;fwx fxx fxy fxz;fwy fxy fyy fyz;fwz fxz fyz fzz];
D=det(m);
if D>0&&fxx<0
    disp("Local Minimum");
else if D>0&&fxx>0
         disp("Local Maximum");
    else D<0
        disp("saddle point")
    end
end