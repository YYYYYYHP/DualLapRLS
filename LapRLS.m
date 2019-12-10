function [LapA] = LapRLS(W1,W2,y, lambda,p_nearest_neighbor,type)
%tju cs, bioinformatics. This program is recoded by reference follow:
%ref:
%[1] Xia Z, Wu L Y, Zhou X, et al. 
%      Semi-supervised drug-protein interaction prediction from heterogeneous biological spaces[J]. 
%           Bmc Systems Biology, 2010, 4(S2):1-16.
% W1 : the kernel of object 1, (m-by-m)
% W2 : the kernel of object 2, (n-by-n)
% y  : binary adjacency matrix, (m-by-n)
%lambda: Regularized item (4)
%p_nearest_neighbor: the p nearest neighbor samples (26)

%Network of Laplacian Regularized Least Square
[num_1,num_2] = size(y);
%1.Sparsification of the similarity matrices
%1.1
%fprintf('Caculating nearest neighbor graph 1\n');
%N_1 = neighborhood_Com(W1,p_nearest_neighbor);
%fprintf('Sparsification of the similarity matrix 1\n');
%S_1 = N_1.*W1;
S_1 = W1;
d_1 = sum(S_1);
D_1 = diag(d_1);
L_D_1 = D_1 - S_1;

d_tmep_1=eye(num_1)/(D_1^(1/2));
L_D_11 = d_tmep_1*L_D_1*d_tmep_1;

%1.2
%fprintf('Caculating nearest neighbor graph 2\n');
%N_2 = neighborhood_Com(W2,p_nearest_neighbor);
%fprintf('Sparsification of the similarity matrix 2\n');
%S_2 = N_2.*W2;
S_2 = W2;
d_2 = sum(S_2);
D_2 = diag(d_2);
L_D_2 = D_2 - S_2;

d_tmep_2=eye(num_2)/(D_2^(1/2));
L_D_22 = d_tmep_2*L_D_2*d_tmep_2;

fprintf('Laplacian Regularized Least Square\n');
%A_1 = (W1/(W1 + lambda*L_D_1*W1))*y;
%A_2 = (W2/(W2 + lambda*L_D_2*W2))*y';

A_1 = W1*pinv(W1 + lambda*L_D_11*W1)*y;
A_2 = W2*pinv(W2 + lambda*L_D_22*W2)*y';

LapA =[];
if type==1
	LapA = (A_1 + A_2')/2;
else
	A_2 = A_2';
	LapA=A_1.*(A_1>=A_2)+A_2.*(A_1<A_2);
end