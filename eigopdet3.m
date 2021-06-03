function [V,D, symmind] = eigopdet3(W)
    % V, D describe all the eigenvalues and eigenvectors of the operator
    % determinant eigenvalue problem of the given 3EP in W.
    % symmind contains the indices of the symmetric eigenvectors in V.
    Delta0 = opdet(W,0);
    Delta1 = opdet(W,1);
    [V,D] = eig(Delta1,Delta0);
    n = round(size(V,2)^(1/3));
    symmind = [];
    for ind = 1:size(V,2)
        X = reshape(V(:,ind),n,n,n);
        Xslice = X(:,:,1);
        if norm(Xslice-Xslice.')<sqrt(eps)
            symmind = [symmind ind];
        end
    end
end

