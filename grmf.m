function [F] = grmf(W1,W2,Y,k,Iteration_max,lamda_L,lamda_1,lamda_2,p_nearest_neighbor)
%tju cs, bioinformatics. This program is coded by reference follow:
%ref:
%[1] Zheng X, Ding H, Mamitsuka H, et al. 
%	Collaborative matrix factorization with multiple similarities for predicting drug-target interactions[C]
%	ACM SIGKDD International Conference on Knowledge Discovery and Data Mining. ACM, 2013:1025-1033.
%and
%[2] Shen Z, Zhang Y H, Han K, et al. 
%    miRNA-Disease Association Prediction with Collaborative Matrix Factorization[J]. Complexity, 2017, 2017(9):1-9.
%
%[3] Ezzat A, Zhao P, Wu M, et al. Drug-Target Interaction Prediction with Graph Regularized Matrix Factorization[J]. 
%    IEEE/ACM Transactions on Computational Biology & Bioinformatics, 2016, PP(99):1-1.
%
% Different from the above articles, our program did not use multiple similarities matrix (ref [1]) and WKNKN (ref [2])


%Graph Regularized Matrix Factorization (GRMF)
%This program is used to Collaborative filtering. 
% W1 : the kernel of object 1, (m-by-m)
% W2 : the kernel of object 2, (n-by-n)
% Y  : binary adjacency matrix, (m-by-n)
% k  : the k is the dimension of the feature spaces
% Iteration_max  : the Iteration_maxis the max numbers of Iteration
%lamda_1 (It is recommended to set 0.001),
%lamda_2 (It is recommended to set 0.001), 
%lamda_L (It is recommended to set 0.1) : these are regularization coefficients of kernel W1, kernel W2, A, B.

fprintf('Graph Regularized Matrix Factorization\n'); 
F=[];

[m,n]=size(Y);
%initial value of A and B
fprintf('initial value of A and B\n'); 
[U1,S_k,V1] = svds(Y,k);
A = U1*(S_k^0.5);  %(m-by-k)
B = V1*(S_k^0.5);  %(n-by-k)

%objective function:
% min    ||Y - AB'||^2 + lamda_l*(||A||^2 + ||B||^2) + lamda_w1*tr(A'*Lw1*A) + lamda_w2*tr(B'*Lw2*B)
% 
% the ||x|| is F norm



%1.Sparsification of the similarity matrices
%1.1
fprintf('Sparsification of the similarity matrices 1\n');
N_1 = preprocess_PNN(W1,p_nearest_neighbor);
S_1 = N_1.*W1;
d_1 = sum(S_1);
D_1 = diag(d_1);
L_D_1 = D_1 - S_1;

%1.2
fprintf('Sparsification of the similarity matrices 2\n');
N_2 = preprocess_PNN(W2,p_nearest_neighbor);
S_2 = N_2.*W2;
d_2 = sum(S_2);
D_2 = diag(d_2);
L_D_2 = D_2 - S_2;

%sloving the problem by alternating least squares
fprintf('Sloving by Alternating least squares\n');
for i=1:Iteration_max
	A = (Y*B - lamda_1*L_D_1*A)/(B'*B + lamda_L*eye(k));
	B = (Y'*A - lamda_2*L_D_2*B)/(A'*A + lamda_L*eye(k));
	%A = (Y*B - lamda_1*W1*A)/(B'*B + lamda_L*eye(k) );
	%B = (Y'*A - lamda_2*W2*B)/(A'*A + lamda_L*eye(k) );
	%A = (Y*B + lamda_1*W1*A)/(B'*B + lamda_L*eye(k) + lamda_1*A'*A);
	%B = (Y'*A + lamda_2*W2*B)/(A'*A + lamda_L*eye(k) + lamda_2*B'*B);

end

%reconstruct Y*
fprintf('Reconstruct Y*\n');
F = A*B';
end

function S=preprocess_PNN(S,p)
%preprocess_PNN sparsifies the similarity matrix S by keeping, for each
%drug/target, the p nearest neighbors and discarding the rest.
%
% S = preprocess_PNN(S,p)

    NN_mat = zeros(size(S));

    % for each drug/target...
    for j=1:length(NN_mat)
        row = S(j,:);                           % get row corresponding to current drug/target
        row(j) = 0;                             % ignore self-similarity
        [~,indx] = sort(row,'descend');         % sort similarities descendingly
        indx = indx(1:p);                       % keep p NNs
        NN_mat(j,indx) = S(j,indx);             % keep similarities to p NNs
        NN_mat(j,j) = S(j,j);                   % also keep the self-similarity (typically 1)
    end

    % symmetrize the modified similarity matrix
    S = (NN_mat+NN_mat')/2;

end