% Generates the discretization of the problem
%
%   u''+lambda*k1(x)*u+f(u)k2(x)*u=0
%
% on domain -1,1 with u(-1)=u(1)=0
% with n discretization points where f(u)=alpha(u)/beta(u) and
%   alpha(u)=int_-1^1 g(x)u(x);
%   beta(u)=u'(0)
%   g(x)=exp(-gamma x^2) (gaussian with scaling gamma)
%
% Returns:  A,B,a,b,xv in the discretization
%    Au+lambda*B*u+(a'*u)/(b'*u)*C*u=0
% where xv is the x-values of the elements in u to be used in plotting.
%
% Example: [A,B,C,a,b,xv]=gen_pde(100,@(x) 1.2+1.1*tanh(5*x), @(x) 0.5*cos(x),20)
function [A,B,C,a,b,xv]=gen_pde(n,k1,k2,gamma)

if (mod(n,2) ~= 0)
    error("n must be even");
end

h=2/(n+1);
e = ones(n,1);
A = spdiags([e -2*e e], -1:1, n, n)/h^2;
% only interior points
xv=linspace(-1,1,n+2)';
xv=xv(2:end-1);

g=@(x) exp(-gamma*(x).^2);
% Trapezoidal rule for alpha:
a=g(xv)*h;
a(1)=a(1)/2; % Half of border points
a(end)=a(end)/2; % Half of border points


% Finite difference for u'(0) to compute b
n1=n/2;
b=zeros(n,1);
b(n1)=-1/(2*h);
b(n1+1)=1/(2*h);

B=spdiags(k1(xv),0,n,n);
C=spdiags(k2(xv),0,n,n);
