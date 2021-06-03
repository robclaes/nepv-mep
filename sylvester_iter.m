function [X] = sylvester_iter(QA, TA, QB, TB, C)
    % Specialized sylvester function for use in iterative methods where the
    % matrices A and B remain constant.
    % This functions solves A*X + X*B = C, similarly as sylvester(A,B,C)
    % does. Note that this is different from lyap(A,B,C).
    % Assumes a Schur decomposition of A and B is already done and given in
    % the factors QA, TA, QB and TB.
    CC = QA'*C*QB;
    X = matlab.internal.math.sylvester_tri(TA, TB, CC, 'I', 'I', 'notransp');
    X = QA*X*QB';
end