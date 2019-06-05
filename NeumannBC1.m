function [h,y]=NeumannBC1(N)
%% Input 

% N           # of time interval


%% Output

% h           step size
% y           numerical solution


%%
a=0; b=1; %interval
alpha=1; beta=exp(1); % boundary conditions

%Grid
h=(b-a)/(N+1);
% A_h
e=ones(N+2,1);
A_h=(1/h^2)*spdiags([e -2*e e],[-1,0,1],N+2,N+2);
A_h(1,1)=-1/h; A_h(1,2)=1/h;
A_h(N+2,N+1)=0; A_h(N+2,N+2)=1;

%the rhs function f
x=linspace(a,b,N+2)';
f=exp(x); f(1)=alpha; f(N+2)=beta;
% numerical solution
y=A_h\f;
