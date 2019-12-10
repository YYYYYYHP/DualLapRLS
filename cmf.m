function [F] = cmf(W1,W2,Y,k,Iteration_max,lamda_1,lamda_2,lamda_L)
%tju cs, bioinformatics. This program is coded by reference follow:
%ref:
%[1] Zheng X, Ding H, Mamitsuka H, et al. 
%	Collaborative matrix factorization with multiple similarities for predicting drug-target interactions[C]
%	ACM SIGKDD International Conference on Knowledge Discovery and Data Mining. ACM, 2013:1025-1033.
%and
%[2] Shen Z, Zhang Y H, Han K, et al. 
%    miRNA-Disease Association Prediction with Collaborative Matrix Factorization[J]. Complexity, 2017, 2017(9):1-9.
%
% Different from the above articles, our program did not use multiple similarities matrix (ref [1]) and WKNKN (ref [2])


%Collaborative Matrix Factorization (CMF)
%This program is used to Collaborative filtering. 
% W1 : the kernel of object 1, (m-by-m)
% W2 : the kernel of object 2, (n-by-n)
% Y  : binary adjacency matrix, (m-by-n)
% k  : the k is the dimension of the feature spaces
% Iteration_max  : the Iteration_maxis the max numbers of Iteration
%lamda_1 (1), lamda_2 (1), lamda_L (1) : the are regularization coefficients of kernel W1, kernel W2, A, B.


F=[];

[m,n]=size(Y);
%initial value of A and B 
%A = abs(rand(m,k));  %(m-by-k)
%B = abs(rand(k,n));  %(k-by-n)
%B = B';

[U1,S_k,V1] = svds(Y,k);
A = U1*(S_k^0.5);  %(m-by-k)
B = V1*(S_k^0.5);  %(n-by-k)

%objective function:
% min    ||Y - AB'||^2 + lamda_L*(||A||^2 + ||B||^2) + lamda_1*||W1 - AA'||^2 + lamda_2*||W2 - BB'||^2
% A,B>0
% the ||x|| is F norm
%sloving the problem by alternating least squares
for i=1:Iteration_max
	A = (Y*B + lamda_1*W1*A)/(B'*B + lamda_L*eye(k) + lamda_1*A'*A);
	B = (Y'*A + lamda_2*W2*B)/(A'*A + lamda_L*eye(k) + lamda_2*B'*B);
	%A = (Y*B - lamda_1*W1*A)/(B'*B + lamda_L*eye(k) );
	%B = (Y'*A - lamda_2*W2*B)/(A'*A + lamda_L*eye(k) );
end

%reconstruct Y*
F = A*B';
