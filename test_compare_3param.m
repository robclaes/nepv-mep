%% Problem specification
clear; close all;
rng(0);

n=10;
A=round(10*randn(n,n));
B=round(10*randn(n,n));
C1=round(10*randn(n,n));
r1=randn(n,1); 
s1=randn(n,1); 
g1=randn(n,1);
C2=round(10*randn(n,n));
r2=randn(n,1);
s2=randn(n,1);
g2=randn(n,1);

f1=@(x) (r1'*x)/(s1'*x); 
f2=@(x) (r2'*x)/(s2'*x);

AA={-A,B,C1,C2};
BB={-A-g1*r1',B,C1-g1*s1',C2};
CC={-A-g2*r2',B,C1,C2-g2*s2'};
Wi = {AA{:};BB{:};CC{:}};

%% Use operator determinants to find symmetric solutions
[V,D,symmind] = eigopdet3(Wi);

%% Set g2 = g1 (not possible for operator determinant approach
g2 = g1;
AA={-A,B,C1,C2};
BB={-A-g1*r1',B,C1-g1*s1',C2};
CC={-A-g2*r2',B,C1,C2-g2*s2'};
W = {AA{:};BB{:};CC{:}};
% W = Wi;


%% Show local convergence
rng(0);
v1 = randn(n,1);
v2 = v1; 
v3 = v1;
pert = randn(n,1)+1i*randn(n,1);
pert = 1e-2*pert/norm(pert);
x10 = V(1:n,symmind(150)) + pert;
x10 = x10/(v1'*x10);
x20 = x10; 
x30 = x10;
TOL = 1e-14;

pertl = randn() + 1i*randn();
pertl = pertl/abs(pertl)*1e-2;
l0 = D(symmind(150),symmind(150))+pertl;
[~,closestind] = min(abs(diag(D)-l0));
closest = D(closestind,closestind);
[~,closestsymmind] = min(abs(diag(D(symmind,symmind)) - l0));
closestsymm = D(symmind(closestsymmind),symmind(closestsymmind));
m10 = f1(x10);
m20 = f2(x10);

[x1,x2,x3,l,m1,m2,hist1] = resinv3(W,l0,m10,m20,v1,v2,v3,x10,x20,x30,TOL,f1,f2);
[x1s,ls,m1s,m2s,hist1s] = resinv_symm3(W,l0,m10,m20,v1,x10,TOL,f1,f2);
[xi,li,histi] = inverse_iteration(Wi,l0,x10,TOL,{f1, f2});

confac = convergence_factor(diag(D(symmind,symmind)),l0);
confacc = convergence_factor(diag(D),l0);
theo_conv = confac.^(1:11);

figure
semilogy(hist1.resnormnl,"k.");hold on; 
semilogy(hist1s.resnormnl,"ko");
semilogy(histi.resnorm,"k+");
semilogy(theo_conv/theo_conv(1)*histi.resnorm(1),"k");
legend("RI", "RIS","II", "Theoretical rate II");
xlabel("Iteration k");
ylabel("||\rho||_2");


%% Symmetric finds a (not nearest) solution, nonsymmetric does not
rng(1); % 1
v1 = randn(n,1);
v2 = v1; 
v3 = v1;
x10 = randn(n,1) + 1i*randn(n,1);
x10 = x10/(v1'*x10);
x20 = x10; 
x30 = x10;
TOL = 1e-14;

% l0 = -1;
l0 = -4-3i;
m10 = f1(x10);
m20 = f2(x10);


[x1,x2,x3,l,m1,m2,hist1] = resinv3(W,l0,m10,m20,v1,v2,v3,x10,x20,x30,TOL,f1,f2);
[x1s,ls,m1s,m2s,hist1s] = resinv_symm3(W,l0,m10,m20,v1,x10,TOL,f1,f2);
[xi,li,histi] = inverse_iteration(Wi,l0,x10,TOL,{f1, f2});

confac = convergence_factor(diag(D(symmind,symmind)),l0);
confac2 = convergence_factor(diag(D),l0);
theo_conv = confac.^(1:82);


figure
semilogy(hist1.resnormnl,"k.");hold on; 
semilogy(hist1s.resnormnl,"ko");
semilogy(histi.resnorm,"k+");
semilogy(theo_conv/theo_conv(1)*histi.resnorm(1),"k");
legend("RI", "RIS","II", "Theoretical rate II");
xlabel("Iteration k");
ylabel("||\rho||_2");


function f=convergence_factor(values, sigma)
    dis = abs(values-sigma);
    sortdis = sort(dis);
    f= sortdis(1)/sortdis(2);
end