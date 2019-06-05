function [h,t,y]=BDF2(a,b,y0,N)
%%
%
% Implementation of BDF2 by using Octave's built-in Newton's solver (fsolve)
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

%% Initiation BDF2 we use EE
y(2)=y(1)+h*f(t(1),y(1));

%% The BDF2 method by using fsolve
for j=2:N
        yy=fsolve(@(yy) yy-4/3*y(j)+1/3*y(j-1)-h*2/3*f(t(j+1),yy), y(j));
	y(j+1)=yy;
end

%% The vector field f, i.e. the RHS of y'(t)=f(t,y(t))
%% You have to modify it appropriately for this case 
function rhs=f(t,y)
rhs=t+y;
