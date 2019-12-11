function [ mean_aupr,mean_auc,globa_true_y_lp,globa_predict_y_lp] = runNFolds(method,y,K_COM1,K_COM2,firstFold,nfolds,CVS,WKNKN,varargin)
    % split folds
    %     crossval_idx = crossvalind('Kfold', length(y(:)), nfolds);
    [num_D,num_G] = size(y);
    if CVS == 1
        crossval_idx = crossvalind('Kfold',y(:),nfolds);
    elseif CVS == 2
        KP = 1:1:num_D;
        crossval_idx = crossvalind('Kfold',KP,nfolds);
    elseif CVS == 3
        KP = 1:1:num_G;
        crossval_idx = crossvalind('Kfold',KP,nfolds);
    end
    
    fold_auc_fblm_ka = [];
    fold_aupr_fblm_ka = [];
    
    globa_true_y_lp=[];
    globa_predict_y_lp=[];

    for fold=1:nfolds
        tic;
        train_idx = find(crossval_idx~=fold);
        test_idx  = find(crossval_idx==fold);

        y_train = y;
        if CVS == 1
            y_train(test_idx) = 0;
        elseif CVS == 2
            y_train(test_idx,:) = 0;
        elseif CVS == 3
            y_train(:,test_idx) = 0;
        end

        if WKNKN
            y_train=preprocess_WKNKN(y_train,K_COM1,K_COM2,3,0.5);
        end
        
        if strcmp(method,'DLapRLS')
            %runNFolds(method,y,K_COM1,K_COM2,1,nfolds,l1,l2,iter_max)
            lamda_1 = varargin{1};
            lamda_2 = varargin{2};
            iter_max = varargin{3};
            [A_cos_com]  = DLapRLS(K_COM1,K_COM2,y_train, lamda_1,lamda_2,iter_max,0);
        elseif strcmp(method,'DLapRLS_C')
            % runNFolds(method,y,K_COM1,K_COM2,1,nfolds,a1,a2,l1,l2,iter_max);
            alpha_1 = varargin{1};
            alpha_2 = varargin{2};
            lamda_1 = varargin{3};
            lamda_2 = varargin{4};
            iter_max = varargin{5};
            [A_cos_com]  = DLapRLS_C(K_COM1,K_COM2,y_train,alpha_1,alpha_2,lamda_1,lamda_2,iter_max,0);
        elseif strcmp(method,'DLapRLS_C1')
            % runNFolds(method,y,K_COM1,K_COM2,1,nfolds,a1,a2,b,l1,l2,iter_max);
            alpha_1 = varargin{1};
            alpha_2 = varargin{2};
            beta = varargin{3};
            lamda_1 = varargin{4};
            lamda_2 = varargin{5};
            iter_max = varargin{6};
            [A_cos_com]  = DLapRLS_C1(K_COM1,K_COM2,y_train,alpha_1,alpha_2,beta,lamda_1,lamda_2,iter_max,0);
        elseif strcmp(method,'LapRLS')
            %runNFolds(method,y,K_COM1,K_COM2,1,nfolds,l1)
            lamda = varargin{1};
            [A_cos_com]  = LapRLS(K_COM1,K_COM2,y_train,lamda,0,1);
        elseif strcmp(method,'CMF')
            %runNFolds(method,y,K_COM1,K_COM2,1,nfolds,k,l1,l2,lL,iter_max)
            k =  varargin{1};
            lamda_1 =  varargin{2};
            lamda_2 =  varargin{3};
            lamda_L =  varargin{4};
            iter_max =  varargin{5};
            [A_cos_com] = cmf(K_COM1,K_COM2,y_train,k,iter_max,lamda_1,lamda_2,lamda_L);
        else
            %[A_cos_com]  = grtmf(K_COM1,K_COM2,y_train,l,l1,l2,item_max,k);
            %[A_cos_com] = test_TMF(K_COM1,K_COM2,y_train,k1,k2,k,2^l1,2^l1,2^l2,2^l2,item_max);
            %[A_cos_com] = C_TMF(K_COM1,K_COM2,y_train,k,iter_max,2^l,2^l1,2^l2);
            %A_cos_com] = TMF(K_COM1,K_COM2,y_train,k1,k2,2^l1);
            %[A_cos_com] = cmf_1(K_COM1,K_COM2,y_train,k,iter_max,2^l,2^l1,2^l2,0.01);
        end
        
    
        toc;
        
        %% 4. evaluate predictions
        yy=y;
        %yy(yy==0)=-1;
        %stats = evaluate_performance(y2(test_idx),yy(test_idx),'classification');
        if CVS == 1
            test_labels = yy(test_idx);
            predict_scores = A_cos_com(test_idx);
        elseif CVS == 2 
            test_labels1 = yy(test_idx,:); 
            test_labels = test_labels1(:);
            predict_scores = A_cos_com(test_idx,:);
            predict_scores = predict_scores(:);
        elseif CVS == 3 
            test_labels1 = yy(:,test_idx); 
            test_labels = test_labels1(:);
            predict_scores = A_cos_com(:,test_idx);
            predict_scores = predict_scores(:);
        end
        
        [X,Y,tpr,aupr_LGC_A_KA] = perfcurve(test_labels,predict_scores,1, 'xCrit', 'reca', 'yCrit', 'prec');

        [X,Y,THRE,AUC_LGC_KA,OPTROCPT,SUBY,SUBYNAMES] = perfcurve(test_labels,predict_scores,1);
        

        fprintf('-METHOD:%s -- FOLD %d - AUPR: %f - AUC: %f \n', method, fold, aupr_LGC_A_KA,AUC_LGC_KA)


        fold_aupr_fblm_ka=[fold_aupr_fblm_ka;aupr_LGC_A_KA];
        fold_auc_fblm_ka=[fold_auc_fblm_ka;AUC_LGC_KA];

        globa_true_y_lp=[globa_true_y_lp;test_labels];
        globa_predict_y_lp=[globa_predict_y_lp;predict_scores];
        if firstFold
            break;
        end
    end
    mean_aupr = mean(fold_aupr_fblm_ka)
    mean_auc = mean(fold_auc_fblm_ka)                   
end