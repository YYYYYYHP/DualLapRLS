function [F] = nrlmf(W1,W2,Y,k,Iteration_max,lamda_1,lamda_2,c,alpha_1,beta_2,learn_rate,p_nearest_neighbor)
%tju cs, bioinformatics. This program is coded by reference follow:
%ref:
%[1] Liu Y, Wu M, Miao C, et al. 
%    Neighborhood Regularized Logistic Matrix Factorization for Drug-Target Interaction Prediction.[J]. 
%    Plos Computational Biology, 2016, 12(2):e1004760.
%
%
%


%Neighborhood Regularized Logistic Matrix Factorization (LMF)
%This program is used to Collaborative filtering. 
% W1 : the kernel of object 1, (m-by-m)
% W2 : the kernel of object 2, (n-by-n)
% Y  : binary adjacency matrix, (m-by-n)
% k  : the k is the dimension of the feature spaces
% Iteration_max  : the Iteration_maxis the max numbers of Iteration
%lamda_1 (0.125), lamda_2 (0.125), alpha_1 (0.0001) beta_2(0.0001) : the are regularization coefficients of kernel W1, kernel W2, U, V.
%c (3) 
%learn_rate: the learn rate of Gradient decline (0.001)

F=[];
fprintf('Neighborhood Regularized Logistic Matrix Factorization\n'); 



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


%%2 Regularized Logistic Matrix Factorization
[m,n]=size(Y);
%initial value of U and V 
%U = normrnd(0,1/sqrt(k),m,k);  %(m-by-k)
%V = normrnd(0,1/sqrt(k),k,n);  %(k-by-n)
%V = V';

[U1,S_k,V1] = svds(Y,k);
U = U1*(S_k^0.5);  %(m-by-k)
V = V1*(S_k^0.5);  %(n-by-k)

[P] = sigm_v(U*V');
%objective function:

%sloving the problem by Gradient decline
for i=1:Iteration_max
	delta_U = P*V + (c - 1)*(Y.*P)*V - c*Y*V + (lamda_1*eye(m) + alpha_1*L_D_1)*U;
	delta_V = P'*U + (c - 1)*(Y'.*P')*U - c*Y'*U + (lamda_2*eye(n) + beta_2*L_D_2)*V;
	U = U - learn_rate*delta_U;
	V = V - learn_rate*delta_V;
	[P] = sigm_v(U*V'); 
end

%reconstruct Y*
F = U*V';

[F] = sigm_v(F);

end


function norm_a = sigm_v(a)

	norm_a = (exp(a))./(1 + exp(a));

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