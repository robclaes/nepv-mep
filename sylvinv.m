function [x,l,hist] = sylvinv(W,sigma,x0,tol,f,max_iter)
    % Sylvester inverse iteration solver for nonlinear eigenvalue problems.
    % W is a 2x3 cell array where the first column corresponds to the left
    % side of the 2EP.
    % f is a function handle represting the function f(x) = r'*x / s'*x.
    if nargin<6
        max_iter = 100;
    end
    t0 = tic();
    dW23 = decomposition(W{2,3});
    dW13 = decomposition(W{1,3});
    A_lyap = dW23\(-W{2,1}+sigma*W{2,2});
    B_lyap = -(dW13\(-W{1,1}+sigma*W{1,2})).';
    M1 = (dW23\W{1,2});
    M2 = (dW13\W{2,2}).';
    C_lyap = @(X) M1*X - X*M2;
    
    hist.resnorm = NaN(max_iter,1);
    hist.resnorm2 = NaN(max_iter,1);
    hist.time_count = NaN(max_iter,1);
    hist.l = NaN(max_iter,1);
    k=1;
    X = x0*x0.';
    
    flag = 'real';
    if ~isreal(A_lyap) || ~isreal(B_lyap)
        flag = 'complex'; % Need complex Schur form
    end
    [QA,TA] = schur(A_lyap,flag);
    [QB,TB] = schur(B_lyap,flag);
    
    while (k==1 || (hist.resnorm(k-1)>tol && k<=max_iter))
        Xn = sylvester_iter(QA, TA, QB, TB, -C_lyap(X));
        X = (Xn+Xn.')/(2*norm(Xn,'fro'));
        X = X/ norm(X(:));
        [hist.l(k),hist.resnorm(k)] = rayleigh(W,X,f);
        hist.resnorm2(k) = sylvres(W,hist.l(k)-1,X);
        hist.time_count(k) = toc(t0);
        k=k+1;
    end
    x = X(:,1)/norm(X(:,1));
    l = hist.l(k-1);
end

function [lambda,res] = rayleigh(W,X,f)
    D1X = -W{2,1}*X*W{1,3}.'+W{2,3}*X*W{1,1}.';
    D0X = W{2,3}*X*W{1,2}.'-W{1,2}*X*W{1,3}.';
    lambda = (X(:)'*D1X(:)) / (X(:)'*D0X(:));
    x = X(:,1)/norm(X(:,1));
    res = norm((-W{1,1} + lambda*W{1,2} + f(x)*W{1,3})*x);
end

function res= sylvres(W,sigma,X)
    X = X/norm(X(:,1));
    f = (-W{2,1}+sigma*W{1,2})*X*W{1,3}.'  - W{2,3}*X*(-W{1,1}+sigma*W{1,2}).'  -  W{2,3}*X*W{1,2}.' + W{1,2}*X*W{1,3}.';
    res = norm(f(:));
end
