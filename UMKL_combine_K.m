lamda_1 = [0.1,0.01,0.001];
lamda_2 = [0.1,0.01,0.001];
k1_paths = {['data2/kernels/disease_simmat_overlap.txt'],...
            ['data2/kernels/disease_simmat_semantic.txt'],...
            };
K1_list = [];
for i=1:length(k1_paths)
    K1_list(:,:,i) = load(k1_paths{i});
end
k2_paths = {['data2/kernels/gene_simmat_overlap.txt'],...
            ['data2/kernels/gene_simmat_gocc.txt'],...
            ['data2/kernels/gene_simmat_ppicos.txt'],...
            ['data2/kernels/gene_simmat_sw.txt'],...
            };
K2_list = [];
for i=1:length(k2_paths)
    K2_list(:,:,i) = load(k2_paths{i});
end
% 2. multiple kernel
% for l1 = lamda_1
%     for l2 = lamda_2
%         [weight_v1] = mas_UMKL_kernels_weights_norm(K1_list,lamda_1,lamda_2)
%         K_COM1 = combine_kernels(weight_v1, K1_list);
%         [weight_v2] = mas_UMKL_kernels_weights_norm(K2_list,lamda_1,lamda_2)
%         K_COM2 = combine_kernels(weight_v2, K2_list);
%         w = [lamda_1,lamda_2,weight_v1,weight_v2]
%         fid=fopen('./weights.txt','at+'); 
%         for i = 1:length(w)
%             fprintf(fid,'%f ',w(i));  
%             if i == length(w)
%                 fprintf(fid,'\n');  
%             else
%                 fprintf(fid,',');
%             end
%         end
%         fclose(fid); 
%     end
% end
y = load('./data2/interactions/GDI_matrix.txt');
[weight_v1] = cka_kernels_weights(K1_list,y,1)
K_COM1 = combine_kernels(weight_v1, K1_list);
[weight_v2] = cka_kernels_weights(K2_list,y,2)
K_COM2 = combine_kernels(weight_v2, K2_list);
save_results('./data2/Kdiseases_cka.txt',K_COM1)
save_results('./data2/Kgenes_cka.txt',K_COM2)