function [V,D, symmind] = eigopdet(W)
    % V, D describe all the eigenvalues and eigenvectors of the operator
    % determinant eigenvalue problem of the given 2EP in W.
    % symmind contains the indices of the symmetric eigenvectors in V.
    Delta0 = opdet(W,0);
    Delta1 = opdet(W,1);
    [V,D] = eig(Delta1,Delta0);
    n = sqrt(size(V,2));
    symmind = [];
    for ind = 1:size(V,2)
        X = reshape(V(:,ind),n,n);
        if norm(X-X.')<sqrt(eps)
            symmind = [symmind ind];
        end
    end
end

