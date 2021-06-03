function [x,l,hist] = inverse_iteration(W, sigma, x0, tol, f)
    n = size(W{1,1},1);
    m = size(W,1)-1;
    Delta0 = opdet(W,0);
    Delta1 = opdet(W,1);
    t0 = tic();
    left = decomposition(Delta1-sigma*Delta0,'qr');
    hist.resnorm = NaN(100,1);
    hist.l = NaN(100,1);
    hist.time_count = NaN(100,1);
    k=1;
    x = x0;
    for i = 2:size(W,1)
        x = kron(x,x0);
    end
    while (k==1 || (hist.resnorm(k-1)>tol && k<150))
        x = left\(Delta0*x);
        x = symmetrize(x,n,m);
        [hist.resnorm(k),hist.l(k)] = nlresidual(W,Delta0,Delta1,x,f);
        hist.time_count(k) = toc(t0);
        k = k+1;
    end
    x = x(1:n);
    l = hist.l(k-1);
end

function [x] = symmetrize(x,n,m)
    X = reshape(x,n*ones(1,m+1));
    per = perms(1:(m+1));
    Y = permute(X,per(1,:));
    for i = 2:size(per,1)
        Y = Y + permute(X,per(i,:));
    end
    x = Y(:);
    x = x/norm(x);
end

function [res,lambda] = nlresidual(W,D0,D1,x,f)
    n = size(W{1,1},1);
    lambda = x'*D1*x/(x'*D0*x);
    T = -W{1,1} + lambda*W{1,2};
    xx = x(1:n)/norm(x(1:n));
    for i = 3:size(W,2)
        T = T + f{i-2}(xx)*W{1,i};
    end
    res = norm(T*xx);
end