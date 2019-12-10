function [ y,K_COM1,K_COM2 ] = loaddata( task )
%LOADDATA 此处显示有关此函数的摘要
%   此处显示详细说明
if strcmp(task, 'DTI')
    regcoef = 10000;
    dataname = 'ic';
    %dataname = 'e';
    %dataname = 'nr';
    %dataname = 'gpcr';
    [y,l1,l2] = loadtabfile(['data/interactions/' dataname '_admat_dgc.txt']);
    k1_paths = {['data/kernels/' dataname '_simmat_proteins_sw-n.txt'],...
        ['data/kernels/' dataname '_simmat_proteins_go.txt'],...
        ['data/kernels/' dataname '_simmat_proteins_ppi.txt'],...
        };
    K1 = [];
    for i=1:length(k1_paths)
        [mat, labels] = loadtabfile(k1_paths{i});
        mat = process_kernel(mat);
        K1(:,:,i) = Knormalized(mat);
    end
    k2_paths = {['data/kernels/' dataname '_simmat_drugs_simcomp.txt'],...
                ['data/kernels/' dataname '_simmat_drugs_sider.txt'],...
                };
    K2 = [];
    for i=1:length(k2_paths)
        [mat, labels] = loadtabfile(k2_paths{i});
        mat = process_kernel(mat);
        K2(:,:,i) = Knormalized(mat);
    end
    % 2. multiple kernel
    [weight_v1] = FKL_weights(K1,y,1,regcoef);
    K_COM1 = combine_kernels(weight_v1, K1);
    [weight_v2] = FKL_weights(K2,y,2,regcoef);
    K_COM2 = combine_kernels(weight_v2, K2);
elseif strcmp(task, 'GDI')
    y = load('./data2/GDI_matrix.txt');
    fprintf('---------------get Y \n')
    K_COM1 = load('./data2/Kdiseases_TKA_200.txt');
    fprintf('---------------get K_COM1 \n')
    K_COM2 = load('./data2/Kgenes_TKA_200.txt');
    fprintf('---------------get K_COM2 \n')
end

end

