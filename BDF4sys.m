function [h,t,y]=BDF4sys(a,b,y0,N)
%%
%
% Implementation of BDF4 by using Octave's built-in Newton's solver (fsolve)
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
m=length(y0);            % dimension of the system
y=zeros(m,N+1);          % numerical solution 
y(:,1)=y0;               % initial value

%% Initiating BDF4 we use RK4
for j=1:3
	k1=f(t(j),y(:,j));
	k2=f(t(j)+0.5*h,y(:,j)+0.5*h*k1);
        k3=f(t(j)+0.5*h,y(:,j)+0.5*h*k2);
	k4=f(t(j)+h,y(:,j)+h*k3);
        y(:,j+1)=y(:,j)+h*(1/6*k1+1/3*k2+1/3*k3+1/6*k4);
end

%% The BDF4 method by using fsolve
for j=4:N
	yy=fsolve(@(yy) yy-48/25*y(:,j)+36/25*y(:,j-1)-16/25*y(:,j-2)+3/25*y(:,j-3)-h*12/25*f(t(j+1),yy), y(:,j));
	y(:,j+1)=yy;
end


%% The vector field f, i.e. the RHS of y'(t)=f(t,y(t))
%% You have to modify it appropriately for this case 
function rhs=f(t,y)
rhs=zeros(2,1);
rhs(1)=y(2);
rhs(2)=-0.3*y(2)-sin(y(1));
