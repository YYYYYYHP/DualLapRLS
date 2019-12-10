function [LapA] = DLapRLS(W1,W2,y, lambda1,lambda2,iterMax,pre)
%Dual Laplacian regularized least square 
%tju cs, bioinformatics. This program is recoded by reference follow:
%ref:
%[1]  Zheng Xia, Ling-Yun Wu, Xiaobo Zhou, Stephen TC Wong 
%            Semi-supervised drug-protein interaction prediction from heterogeneous biological spaces[J]. 
%                BMC Systems Biology, 2010, 4(2):s6.
%
%[2] Ghosh A, Sekhar C C. 
%             Label Correlation Propagation for Semi-supervised Multi-label Learning[J]. 2017. 
% 
%[3] Chen G, Song Y, Wang F, et al. 
%          Semi-supervised Multi-label Learning by Solving a Sylvester Equation[C]// 
%            Siam International Conference on Data Mining, SDM 2008, April 24-26, 2008, Atlanta, Georgia, Usa. DBLP, 2008:410-419.
%                
% W1 : the kernel of object 1, (m-by-m)
% W2 : the kernel of object 2, (n-by-n)
% y  : binary adjacency matrix, (m-by-n)
%lambda1: Regularized item for W1 (0.9)
%lambda2: Regularized item for W2 (0.9)

%Network of Laplacian Regularized Least Square
[num_1,num_2] = size(y);
%0 preprocesses the interaction matrix Y 
if pre==1
	fprintf('Y preprocesses\n');
   y = preprocess_Y(y,W1,W2,3,0.5);
else
   fprintf('No preprocesses\n');
end

%1.
%1.1

S_1 = W1;
d_1 = sum(S_1);
D_1 = diag(d_1);
L_D_1 = D_1 - S_1;
%Laplacian Regularized

L_D_11 = (D_1^(-0.5))*L_D_1*(D_1^(-0.5));


%1.2

S_2 = W2;
d_2 = sum(S_2);
D_2 = diag(d_2);
L_D_2 = D_2 - S_2;
%Laplacian Regularized

L_D_22 = (D_2^(-0.5))*L_D_2*(D_2^(-0.5));

%1.3 Solve Laplacian Regularized least square
fprintf('Dual Laplacian Regularized Least Square\n');
LapA =[];

%Matrix initialization
Ta_d = randn(num_1,num_2);
Ta_t = randn(num_2,num_1);
for ii=1:iterMax

	Ta_d=(W1*W1 + lambda1*L_D_11 )\(W1*(2*y - Ta_t'*W2));
	temp_t = (2*y - W1*Ta_d)*W2/(lambda2*L_D_22 + W2*W2);
	Ta_t  = temp_t';

end
%LapA =W1*Ta_d;
%LapA = (W2*Ta_t)';
LapA = (W1*Ta_d + (W2*Ta_t)')/2;

end







function similarities_N = neighborhood_Com(similar_m,kk)

similarities_N=zeros(size(similar_m));

mm = size(similar_m,1);

for ii=1:mm
	
	for jj=ii:mm
		iu = similar_m(ii,:);
		iu_list = sort(iu,'descend');
		iu_nearest_list_end = iu_list(kk);
		
		ju = similar_m(:,jj);
		ju_list = sort(ju,'descend');
		ju_nearest_list_end = ju_list(kk);
		if similar_m(ii,jj)>=iu_nearest_list_end & similar_m(ii,jj)>=ju_nearest_list_end
			similarities_N(ii,jj) = 1;
			similarities_N(jj,ii) = 1;
		elseif similar_m(ii,jj)<iu_nearest_list_end & similar_m(ii,jj)<ju_nearest_list_end
			similarities_N(ii,jj) = 0;
			similarities_N(jj,ii) = 0;
		else
			similarities_N(ii,jj) = 0.5;
			similarities_N(jj,ii) = 0.5;
		end
	
	
	end


end

end




function newY=preprocess_Y(Y,Sd,St,K,eta)
%preprocesses the interaction matrix Y by replacing each
%of the 0's (i.e. presumed non-interactions) with a continuous value
%between 0 and 1. For each 0, the K nearest known drugs are used to infer
%a value, the K nearest known targets are used to infer another value, and
%then the average of the two values is used to replace that 0.
 % decay values to be used in weighting similarities later
    eta = eta .^ (0:K-1);
	newY =[];
    y2_new1 = zeros(size(Y));
    y2_new2 = zeros(size(Y));

    empty_rows = find(any(Y,2) == 0);   % get indices of empty rows
    empty_cols = find(any(Y)   == 0);   % get indices of empty columns

    % for each drug i...
    for i=1:length(Sd)
        drug_sim = Sd(i,:); % get similarities of drug i to other drugs
        drug_sim(i) = 0;    % set self-similiraty to ZERO

        indices  = 1:length(Sd);    % ignore similarities 
        drug_sim(empty_rows) = [];  % to drugs of 
        indices(empty_rows) = [];   % empty rows

        [~,indx] = sort(drug_sim,'descend');    % sort descendingly
        indx = indx(1:K);       % keep only similarities of K nearest neighbors
        indx = indices(indx);   % and their indices

        % computed profile of drug i by using its similarities to its K
        % nearest neighbors weighted by the decay values from eta
        drug_sim = Sd(i,:);
        y2_new1(i,:) = (eta .* drug_sim(indx)) * Y(indx,:) ./ sum(drug_sim(indx));
    end

    % for each target j...
    for j=1:length(St)
        target_sim = St(j,:); % get similarities of target j to other targets
        target_sim(j) = 0;    % set self-similiraty to ZERO

        indices  = 1:length(St);        % ignore similarities 
        target_sim(empty_cols) = [];    % to targets of
        indices(empty_cols) = [];       % empty columns

        [~,indx] = sort(target_sim,'descend');  % sort descendingly
        indx = indx(1:K);       % keep only similarities of K nearest neighbors
        indx = indices(indx);   % and their indices

        % computed profile of target j by using its similarities to its K
        % nearest neighbors weighted by the decay values from eta
        target_sim = St(j,:);
        y2_new2(:,j) = Y(:,indx) * (eta .* target_sim(indx))' ./ sum(target_sim(indx));
    end

    % average computed values of the modified 0's from the drug and target
    % sides while preserving the 1's that were already in Y 
    newY = max(Y,(y2_new1 + y2_new2)/2);
	

end