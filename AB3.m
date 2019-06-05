function [h,t,y]=AB3(a,b,y0,N)

%% Input 

% a           starting time
% b           final time
% y0          initial value
% N           # of time interval


%% Output

% h           step size
% t           time grid points
% y           numerical solution


%% 
h=(b-a)/N;               % step size        
t=linspace(a,b,N+1);     % time grid points   
y=zeros(1,N+1);          % numerical solution 
y(1)=y0;                 % initial value

%% Initiating AB3 we use an ERK(3,3) and other options (ERK1 and ERK2)
for j=1:2 %ERK3
    %ff(j)=f(t(j),y(j)); Comapct implementation
    k1=f(t(j),y(j));
    %k2=f(t(j)+0.5*h,y(j)+0.5*h*k1);
    %k3=f(t(j)+0.75*h,y(j)+0.75*h*k2);
    %y(j+1)=y(j)+h*(2/9*k1+3/9*k2+4/9*k3);
    k2=f(t(j)+2/3*h,y(j)+2/3*h*k1);%OptiRK2
    y(j+1)=y(j)+h*(1/4*k1+3/4*k2); %OptiRK2
    %y(j+1)=y(j)+h*k1;
end

%% The AB3 method
for j=3:N
    %ff(j)=f(t(j),y(j)); Comapct implementation
    %y(j+1)=y(j)+h*(23/12*ff(j)-16/12*ff(j-1)+5/12*ff(j-2)); Comapct implementation
    y(j+1)=y(j)+h*(23/12*f(t(j),y(j))-16/12*f(t(j-1),y(j-1))+5/12*f(t(j-2),y(j-2)));
end


%% The vector field f, i.e. the RHS of y'(t)=f(t,y(t))
%% You have to modify it appropriately for this case 
function rhs=f(t,y)
rhs=y+t;
