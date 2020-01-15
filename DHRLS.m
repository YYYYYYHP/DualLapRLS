function [LapA] = DHRLS(W1,W2,y, beta,lambda1,lambda2,knn,iterMax)
%dual Hypergraph regularized least square model
%tju cs, bioinformatics. This program is recoded by reference follow:
%ref:
           
% W1 : the kernel of object 1, (m-by-m)
% W2 : the kernel of object 2, (n-by-n)
% y  : binary adjacency matrix, (m-by-n)
%lambda1: Regularized item for W1 (0.9)
%lambda2: Regularized item for W2 (0.9)

[num_1,num_2] = size(y);

%1.
%1.1
%Hypergraph Laplacian Regularized

L_D_11 = construct_Hypergraphs_knn(W1, knn);


%1.2
%Hypergraph Laplacian Regularized

L_D_22 = construct_Hypergraphs_knn(W2, knn);

%1.3 Solve Laplacian Regularized least square
fprintf('Dual Hypergraph regularized least square model\n');
LapA =[];

%Matrix initialization
Ta_d = randn(num_1,num_2);
Ta_t = randn(num_2,num_1);
for ii=1:iterMax

	Ta_d = pinv(W1*W1 + beta*eye(num_1) + lambda1*W1*L_D_11*W1 )*(2*W1*y - W1* Ta_t'*W2);
	Ta_t = pinv(W2*W2 + beta*eye(num_2) + lambda2*W2*L_D_22*W2 )*(2*W2*y' -W2* Ta_d'*W1);
    fprintf('---------------iteration :%d\n',ii);
end
%LapA =W1*Ta_d;
%LapA = (W2*Ta_t)';
%LapA = (W1*Ta_d + (W2*Ta_t)')/2;
LapA = W1*Ta_d+(W2*Ta_t)';
end


