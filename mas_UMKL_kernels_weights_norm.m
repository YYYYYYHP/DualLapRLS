function [w] = mas_UMKL_kernels_weights_norm(Kernels_list,lamda_1,lamda_2)
%Unsupervised multiple kernel learning based on maximizes the average similarity 
%tju cs, bioinformatics. This program is recoded by reference follow:
%ref:
%[1] Mariette J, Villa-Vialaneix N. Unsupervised multiple kernel learning for heterogeneous data integration[J]. 
%              Bioinformatics, 2018.
%
%[2] Lavit, C., Escoufier, Y., Sabatier, R., and Traissac, P. (1994). The act (statis method). 
%               Computational Statistics & Data Analysis, 18(1), 97 – 119.
%

fprintf('Unsupervised multiple kernel learning based on maximizes the average similarity \n');
num_kernels = size(Kernels_list,3);

w = zeros(num_kernels,1);




M = zeros(num_kernels,num_kernels);

for i=1:num_kernels
	for j=1:num_kernels
		kk1 = Kernels_list(:,:,i);
		kk2 = Kernels_list(:,:,j);
		mm = trace(kk1'*kk2);
		m1 = trace(kk1*kk1');
		m2 = trace(kk2*kk2');
		M(i,j) = mm/(sqrt(m1*m2));
	end
end

N_U = size(Kernels_list,1);
l=ones(N_U,1);
H = eye(N_U) - (l*l')/N_U;

M = zeros(num_kernels,num_kernels);

for i=1:num_kernels
	for j=1:num_kernels
		kk1 = H*Kernels_list(:,:,i)*H;
		kk2 = H*Kernels_list(:,:,j)*H;
		%a1 = kk1(:);
		%a2 = kk2(:);
		
		%mm = dot( a1,a2 )/( sqrt( sum( a1.*a1 ) ) * sqrt( sum( a2.*a2 ) ) );
		mm = trace(kk1'*kk2);
		m1 = trace(kk1*kk1');
		m2 = trace(kk2*kk2');
		M(i,j) = mm/(sqrt(m1)*sqrt(m2));
		%M(i,j) = (mm);
	end
end
d_1 = sum(M);
D_1 = diag(d_1);
LapM = D_1 - M;

v = randn(num_kernels,1);
%v(1:num_kernels) = 1/num_kernels;

falpha = @(v)obj_function(v,M,LapM,lamda_1,lamda_2);
        
        % Optimal v
[x_alpha, fval_alpha] = optimize_weights(v, falpha);
% cvx_begin
%     variable v(num_kernels,1);
%     minimize( norm(v'*K_values,2)*TEMP_V1 - v'*TEMP_V2);
% 	%subject to
% 	v >= 0;
% 	sum(v)==1;
% cvx_end

w = x_alpha;


end

function [J] = obj_function(w,Ma,Lm,regcoef1,regcoef2)
    
    J = -1*(w'*Ma*w) + regcoef1*w'*Lm*w  + regcoef2*(norm(w,2))^2;
end

function [x, fval] = optimize_weights(x0, fun)
    n = length(x0);
    Aineq   = [];
    bineq   = [];
    Aeq     = ones(1,n);
    beq     = 1;
    LB      = zeros(1,n);
    UB      = ones(1,n);

    options = optimoptions('fmincon','Algorithm','interior-point', 'Display', 'notify');
    [x,fval] = fmincon(fun,x0,Aineq,bineq,Aeq,beq,LB,UB,[],options);
end




