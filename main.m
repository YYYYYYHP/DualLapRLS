clear
seed = 12345678;
nfolds = 10; 
% task = 'DTI';
task = 'GDI';

method = 'DLapRLS';
% method = 'DLapRLS_C';
% method = 'DLapRLS_C1';
% method = 'CMF';
% method = 'LapRLS';

one_kernel = 1;

globa_true_y_lp=[];
globa_predict_y_lp=[];

results = [];

[ y,K_COM1,K_COM2,K1_list,K2_list ] = loaddata( task ,one_kernel);

if one_kernel
    K_COM1 = K1_list(:,:,1);
    K_COM2 = K2_list(:,:,1);
end

if strcmp(method , 'DLapRLS')
    lamda_1 = [4,2,1,2^-1,2^-2];
    lamda_2 = [4,2,1,2^-1,2^-2];
    iter_max = 15;
    WKNKN = 0;
    for l1 = lamda_1
        for l2 = lamda_2
            fprintf('-- lamda_1: %f - lamda_2: %f \n', l1,l2)
            [ mean_aupr,mean_auc,globa_true_y_lp,globa_predict_y_lp ] = runNFolds(method,y,K_COM1,K_COM2,0,nfolds,WKNKN,l1,l2,iter_max);
            results = cat(1,results,[l1,l2,mean_aupr,mean_auc]);
            save_results(['./results/DLapRLS/' task '/mean_NWKNKN_onek_result.txt'],results);
        end
    end
  
elseif strcmp(method , 'DLapRLS_C')
    alpha1 = 0.5;
    lamda_1 = [4,2,1,2^-1,2^-2];
    lamda_2 = [4,2,1,2^-1,2^-2];
    iter_max = 15;
    WKNKN = 0;
    for a1 = alpha1
        a2 = 1-a1;
        for l1 = lamda_1
            for l2 = lamda_2
                [ mean_aupr,mean_auc,globa_true_y_lp,globa_predict_y_lp ] = runNFolds(method,y,K_COM1,K_COM2,1,nfolds,WKNKN,a1,a2,l1,l2,iter_max);
                results = cat(1,results,[a1,a2,l1,l2,mean_aupr,mean_auc]);
                save_results(['./results/DLapRLS_C/' task '/firstFold_result.txt'],results);
            end
        end
    end
    
elseif strcmp(method , 'DLapRLS_C1')
%     alpha1 = [0.7,0.8,0.9,1.1];
%     beta = [0.01,0.001,0.0001];
    alpha1 = 1;
    beta = 0.01;
    lamda_1 = [4,2,1,2^-1,2^-2];
    lamda_2 = [4,2,1,2^-1,2^-2];
    iter_max = 10;
    WKNKN = 0;
    for a1 = alpha1
        a2 = 2-a1;
        for b = beta
            for l1 = lamda_1
                for l2 = lamda_2
                    fprintf('-METHOD:%s - alpha1: %f - alpha2:%f - beta:%f - lamda_1 :%f - lamda_2 :%f  \n',method,a1,a2,b,l1,l2)
                    [ mean_aupr,mean_auc,globa_true_y_lp,globa_predict_y_lp ] = runNFolds(method,y,K_COM1,K_COM2,1,nfolds,WKNKN,a1,a2,b,l1,l2,iter_max);
                    results = cat(1,results,[a1,a2,b,l1,l2,mean_aupr,mean_auc]);
                    save_results(['./results/DLapRLS_C1/' task '/firstFold_result.txt'],results);
                end
            end
        end
    end
    
elseif strcmp(method , 'LapRLS')
    lamda_1 = [2^-3,2^-2,2^-1,1,2,2^2,2^3];
    WKNKN = 0;
    for l1 = lamda_1
        fprintf('-METHOD:%s  - lamda_1: %f\n', method,l1)
        [ mean_aupr,mean_auc,globa_true_y_lp,globa_predict_y_lp ] = runNFolds(method,y,K_COM1,K_COM2,0,nfolds,WKNKN,l1);
        results = cat(1,results,[l1,mean_aupr,mean_auc]);
        save_results(['./results/LapRLS/' task '/GDI_mean_result.txt'],results); %2^-4,
    end
    
elseif strcmp(method , 'CMF')
    K = 500:100:900;
    if strcmp(task,'DTI')
        K = [50,60,70,80,90,100];
    end
    lamda_1 = [2^-2,2^-1,0,1,2,2^2,2^3];
    lamda_2 = [2^-2,2^-1,0,1,2,2^2,2^3];
    lamda_L = [2^-2,2^-1,0,1,2,2^2,2^3];
    iter_max = 40;
    WKNKN = 0;
    for k = K
        for l1 = lamda_1
            for l2 = lamda_2
                for lL = lamda_L
                    fprintf('-METHOD:%s - K:%d - lamda_1 :%f - lamda_2 :%f - lamda_L: %f\n',method,k,l1,l2,lL)
                    [ mean_aupr,mean_auc,globa_true_y_lp,globa_predict_y_lp ] = runNFolds(method,y,K_COM1,K_COM2,0,nfolds,WKNKN,k,l1,l2,lL,iter_max);
                    results = cat(1,results,[l1,l2,mean_aupr,mean_auc]);
                    save_results(['./results/CMF/' task '/firstFold_result.txt'],results);
                end
            end
        end
    end
end
