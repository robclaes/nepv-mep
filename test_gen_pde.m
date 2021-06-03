clear;
close all;

rng(1);
n=100;
g = randn(n,1);
k1 = @(x) 1+0.5*tanh(5*x);
k2 = @(x) 1+0.5*cos(pi*x);
[A,B,C,r,s,xv]=gen_pde(n,k1,k2,10); %
B=full(B);
C=full(C);
W1 = {-A,B,C};
W2 = {-A-g*r',B,C-g*s'}; % With eigenvector nonlinearity
Wn = {W1{:}; W2{:}}; % With eigenvector nonlinearity
fn = @(x) (r'*x)/(s'*x); % With eigenvector nonlinearity
tol = 1e-14;
sigma= 1;
x0 = rand(n,1); x0 = x0/norm(x0);

[xn,ln,histn] = sylvinv(Wn,sigma,x0,tol,fn,5);
[xnn,lnn,mnn,histnn] = resinv_symm(Wn,ln,fn(xn),randn(n,1),xn,tol,fn,20); % Can RIS reduce residu?

figure(1);
plot(xv,-xn,"--"); hold on
plot(xv,xnn/norm(xnn));hold off;
xlabel("x");ylabel("u(x)");
legend("II","RIS");

figure(2);
semilogy(1:5,histn.resnorm); hold on;
semilogy(6:25,histnn.resnormnl); hold off;
xlabel("Iteration k");
ylabel("||\rho||_2");
legend("II", "RIS");