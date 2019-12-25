clear
seed = 12345678;
nfolds = 10; 
% task = 'DTI';
task = 'GDI';

% method = 'DLapRLS';
% method = 'DLapRLS_C';
method = 'DLapRLS_C1';
% method = 'CMF';
% method = 'LapRLS';
% method = 'grmf'
% method = 'nrlmf'

one_kernel = 0;


globa_true_y_lp=[];
globa_predict_y_lp=[];

results = [];

if ~one_kernel
    [ y,K_COM1,K_COM2] = loaddata( task );
else
    fprintf('---------------one kernel -- disease semantic -- gene sw-- \n')
    y = load('./data2/interactions/GDI_matrix.txt');
    fprintf('---------------get Y \n')
    K_COM1 = load('data2/kernels/disease_simmat_semantic.txt');
    fprintf('---------------get K_COM1 \n')
    K_COM2 = load('data2/kernels/gene_simmat_sw.txt');
    fprintf('---------------get K_COM2 \n')
end

if strcmp(method , 'DLapRLS')
%     lamda_1 = [4,2,1,2^-1,2^-2];
%     lamda_2 = [4,2,1,2^-1,2^-2];
    lamda_1 = 1;
    lamda_2 = 0.25;
    iter_max = 15;
    WKNKN = 1;
    for l1 = lamda_1
        for l2 = lamda_2
            fprintf('-- lamda_1: %f - lamda_2: %f \n', l1,l2)
            [ mean_aupr,mean_auc,globa_true_y_lp,globa_predict_y_lp ] = runNFolds(method,y,K_COM1,K_COM2,0,nfolds,CVS,WKNKN,l1,l2,iter_max);
            results = cat(1,results,[l1,l2,mean_aupr,mean_auc]);
            %save_results(['./results/DLapRLS/' task '/mean_NWKNKN_onek_result.txt'],results);
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
                [ mean_aupr,mean_auc,globa_true_y_lp,globa_predict_y_lp ] = runNFolds(method,y,K_COM1,K_COM2,1,nfolds,CVS,WKNKN,a1,a2,l1,l2,iter_max);
                results = cat(1,results,[a1,a2,l1,l2,mean_aupr,mean_auc]);
                save_results(['./results/DLapRLS_C/' task '/firstFold_result.txt'],results);
            end
        end
    end
    
elseif strcmp(method , 'DLapRLS_C1')
%     alpha1 = 0:0.2:2;
%     beta = [0.01,0.001,0.0001];
    alpha1 = 1;
    beta = 0.01;
%     lamda_1 = [2,1,2^-1,2^-2];
%     lamda_2 = [2,1,2^-1,2^-2];
    lamda_1=1;
    lamda_2=0.25;
    iter_max = 10;
    WKNKN = 0;
    first = 0;
    CVS = 3;
    for a1 = alpha1
        a2 = 2-a1;
        for b = beta
            for l1 = lamda_1
                for l2 = lamda_2
                    fprintf('-METHOD:%s - alpha1: %f - alpha2:%f - beta:%f - lamda_1 :%f - lamda_2 :%f  \n',method,a1,a2,b,l1,l2)
                    [ mean_aupr,mean_auc,globa_true_y_lp,globa_predict_y_lp ] = runNFolds(method,y,K_COM1,K_COM2,first,nfolds,CVS,WKNKN,a1,a2,b,l1,l2,iter_max);
                    results = cat(1,results,[a1,a2,b,l1,l2,mean_aupr,mean_auc]);
                    save_results(['./results/DLapRLS_C1/' task '/1_025_GDI_mean_NWKNKN_CVS3_result.txt'],results);
                end
            end
        end
    end
    
elseif strcmp(method , 'LapRLS')
    lamda_1 = [2^-2,2^-1,1,2,2^2];
    WKNKN = 0;
    CVS = 1;
    first = 0;
    for l1 = lamda_1
        fprintf('-METHOD:%s  - lamda_1: %f\n', method,l1)
        [ mean_aupr,mean_auc,globa_true_y_lp,globa_predict_y_lp ] = runNFolds(method,y,K_COM1,K_COM2,first,nfolds,CVS,WKNKN,l1);
        results = cat(1,results,[l1,mean_aupr,mean_auc]);
        save_results(['./results/LapRLS/' task '/GDI_mean_CVS1_result.txt'],results); %2^-4,
    end
    
elseif strcmp(method , 'grmf')
    K = 700:100:900;
    lamda_1 = [0,10^-4,10^-3,10^-2,10^-1];
    lamda_2 = [0,10^-4,10^-3,10^-2,10^-1];
    lamda_L = [2^-2,2^-1,1,2];
    p_nearest = 5;
    iter_max = 10;
    WKNKN = 0;
    CVS = 1;
    first = 0;
    for k = K
        for l1 = lamda_1
            for l2 = lamda_2
                for lL = lamda_L
                    fprintf('-METHOD:%s - K:%d - lamda_1 :%f - lamda_2 :%f - lamda_L: %f\n',method,k,l1,l2,lL)
                    [ mean_aupr,mean_auc,globa_true_y_lp,globa_predict_y_lp ] = runNFolds(method,y,K_COM1,K_COM2,first,nfolds,CVS,WKNKN,k,l1,l2,lL,iter_max,p_nearest);
                    results = cat(1,results,[k,l1,l2,lL,mean_aupr,mean_auc]);
                    save_results(['./results/grmf/' task '/GDI_mean_CVS1_result.txt'],results); 
                end
            end
        end
    end
    
elseif strcmp(method , 'nrlmf')
    K = 700:100:900;
    lamda_1 = [2^-5,2^-4,2^-3,2^-2,2^-1,1];
    alpha_1 = [2^-5,2^-4,2^-3,2^-2,2^-1,1,2,4];
    beta_2 = [2^-5,2^-4,2^-3,2^-2,2^-1,1];
    c = 3;
    learn_rate = 0.001;
    p_nearest = 5;
    iter_max = 10;
    WKNKN = 0;
    CVS = 1;
    first = 0;
    for k = K
        for l1 = lamda_1
            l2 = l1;
            for a1 = alpha1
                for b2 = beta_2
                    fprintf('-METHOD:%s - K:%d - lamda_1 = lamda_2 :%f - alpha_1 :%f - beta_2: %f\n',method,k,l1,a1,b2)
                    [ mean_aupr,mean_auc,globa_true_y_lp,globa_predict_y_lp ] = runNFolds(method,y,K_COM1,K_COM2,first,nfolds,CVS,WKNKN,k,l1,l2,a1,b2,learn_rate,p_nearest,iter_max);
                    results = cat(1,results,[k,l1,l2,lL,mean_aupr,mean_auc]);
                    save_results(['./results/nrlmf/' task '/GDI_mean_CVS1_result.txt'],results); 
                end
            end
        end
    end
    
elseif strcmp(method , 'CMF')
    K = 700:100:900;
    if strcmp(task,'DTI')
        K = [50,60,70,80,90,100];
    end
    lamda_1 = [2^-2,2^-1,1,2];
    lamda_2 = [2^-2,2^-1,1,2];
    lamda_L = [2^-2,2^-1,1,2];
    iter_max = 40;
    WKNKN = 0;
    CVS = 3;
    first = 0;
    for k = K
        for l1 = lamda_1
            for l2 = lamda_2
                for lL = lamda_L
                    fprintf('-METHOD:%s - K:%d - lamda_1 :%f - lamda_2 :%f - lamda_L: %f\n',method,k,l1,l2,lL)
                    [ mean_aupr,mean_auc,globa_true_y_lp,globa_predict_y_lp ] = runNFolds(method,y,K_COM1,K_COM2,first,nfolds,CVS,WKNKN,k,l1,l2,lL,iter_max);
                    results = cat(1,results,[k,l1,l2,lL,mean_aupr,mean_auc]);
                    save_results(['./results/CMF/' task '/GDI_mean_CVS3_result.txt'],results);
                end
            end
        end
    end
end
