%% Problem specification

clear; close all;
rng(2);
n=5;
A=round(10*randn(n,n));
B=round(10*randn(n,n));
C=round(10*randn(n,n));
g=randn(n,1);
r=randn(n,1);
s=randn(n,1);

f=@(x) (r'*x)/(s'*x);
AA={-A,B,C};
BB={-A-g*r',B,C-g*s'};
W = {AA{:};BB{:}};


%% Solve problem using operator determinants

[V,D,symmind] = eigopdet(W);
figure(1);
nonsymmind = setdiff(1:n^2,symmind);
plot(diag(D(nonsymmind,nonsymmind)),"kx"); hold on;
plot(diag(D(symmind,symmind)),"ko");
axis equal;
xlabel("Re"); ylabel("Im");
xlim([-7 7]);
ylim([-5 5]);

%% Setup problem for resinv 

v1 = randn(n,1);
v2 = v1;
TOL = 1e-14;

%% Solve problem using resinv
%  All methods converge to the desired eigenvalue
rng(0);
pert = randn(n,1);
pert = 1e-1*pert/norm(pert);
ind = 8;
X = reshape(V(:,ind),n,n);
x10 = X(:,1)+pert;
x10 = x10/(v1'*x10);
x20 = x10;
m0 = f(x10);
pertl = randn() + 1i*randn();
pertl = pertl/abs(pertl)*1e-1;
l0 = D(ind,ind)+pertl; 
figure(1);plot(real(l0),imag(l0),"k.","MarkerSize",12);
confac1symm = convergence_factor(diag(D(symmind,symmind)),l0);
confac1 = convergence_factor(diag(D),l0);
theo_conv_symm1 = confac1symm.^(1:14);
[x1,x2,l,m,hist1] = resinv(W, l0, m0,v1,v2,x10,x20,TOL,f);
[x1s,ls,ms,hist1s] = resinv_symm(W, l0, m0,v1,x10,TOL,f);
[xx,ll,hist] = sylvinv(W,l0,x10,TOL,f);


figure;
semilogy(hist1.resnormnl,"k."); hold on;
semilogy(hist1s.resnormnl,"ko");
semilogy(hist.resnorm,"k+");
semilogy(theo_conv_symm1/theo_conv_symm1(1)*hist.resnorm(1),"k");
legend("RI", "RIS", "II","Theoretical rate II");
xlabel("Iteration k");
ylabel("||\rho||_2");


%% Solve problem using resinv
%  (Symmetric) resinv converges to symmetric (not closest) eigenvalue
%  Sylvester converges to closest symmetric eigenvalue
rng(3);
x10 = randn(n,1)+1i*randn(n,1);
x10 = x10/(v1'*x10);
x20 = x10;
m0 = f(x10);
l0 = -2i; 
figure(1);
plot(real(l0),imag(l0),"ks",'MarkerFaceColor',[0 0 0]);
legend("Nonsymmetric eigenvalues", "Symmetric eigenvalues","\sigma_1","\sigma_2");
confac2_symm = convergence_factor(diag(D(symmind,symmind)),l0);
theo_conv_symm2 = confac2_symm.^(1:60);
[x1,x2,l,m,hist1] = resinv(W, l0, m0,v1,v2,x10,x20,TOL,f);
[x1s,ls,ms,hist1s] = resinv_symm(W, l0, m0,v1,x10,TOL,f);
[xx,ll,hist] = sylvinv(W,l0,x10,TOL,f);

figure;
semilogy(hist1.resnormnl,"k."); hold on;
semilogy(hist1s.resnormnl,"ko");
semilogy(hist.resnorm,"k+");
semilogy(theo_conv_symm2/theo_conv_symm2(1)*hist.resnorm(1),"k");
legend("RI", "RIS", "II","Theoretical rate II");
xlabel("Iteration k");
ylabel("||\rho||_2");




function f=convergence_factor(values, sigma)
    dis = abs(values-sigma);
    sortdis = sort(dis);
    f= sortdis(1)/sortdis(2);
end