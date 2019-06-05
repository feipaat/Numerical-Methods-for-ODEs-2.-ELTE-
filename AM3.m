function [h,t,y]=AM3(a,b,y0,N)
%%
%
% Predictor-corrector based implementation of AM3 by using
% RK3 and AB3 methods
%

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

%% Initiating AM3 we use an ERK(3,3)
for j=1:2 %ERK3
    %ff(j)=f(t(j),y(j)); Comapct implementation
    k1=f(t(j),y(j));
    k2=f(t(j)+0.5*h,y(j)+0.5*h*k1);
    k3=f(t(j)+0.75*h,y(j)+0.75*h*k2);
    y(j+1)=y(j)+h*(2/9*k1+3/9*k2+4/9*k3);
    %y(j+1)=y(j)+h*k1;
end

%% The predictor-corrector implementation
for j=3:N
    ypred(j+1)=y(j)+h*(23/12*f(t(j),y(j))-16/12*f(t(j-1),y(j-1))+5/12*f(t(j-2),y(j-2)));
    %y(j+1)=y(j)+h*(5/12*f(t(j+1),ypred(j+1))+8/12*f(t(j),y(j))-1/12*f(t(j-1),y(j-1))); %PEC version
    ypp(j+1)=y(j)+h*(5/12*f(t(j+1),ypred(j+1))+8/12*f(t(j),y(j))-1/12*f(t(j-1),y(j-1))); %PECE version
    y(j+1)=y(j)+h*(5/12*f(t(j+1),ypp(j+1))+8/12*f(t(j),y(j))-1/12*f(t(j-1),y(j-1))); %PECE version
    %ypp(j+1)=y(j)+h*(5/12*f(t(j+1),y(j+1))+8/12*f(t(j),y(j))-1/12*f(t(j-1),y(j-1))); %PECEC version
    %y(j+1)=y(j)+h*(5/12*f(t(j+1),ypp(j+1))+8/12*f(t(j),y(j))-1/12*f(t(j-1),y(j-1))); %PECECE version
    error_esti(j+1)=1/10*(y(j+1)-ypred(j+1));
end
%error_esti %for demonstration

%% The vector field f, i.e. the RHS of y'(t)=f(t,y(t))
%% You have to modify it appropriately for this case 
function rhs=f(t,y)
rhs=t+y;
