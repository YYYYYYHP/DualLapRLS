clear
seed = 12345678;
nfolds = 10; 
% task = 'DTI';
task = 'GDI';

method = 'DLapRLS_C1';
[y,K_COM1,K_COM2] = loaddata( task );
alpha1 = 1;
beta = 0.01;
alpha2=1;

lamda_1=0.5;
lamda_2=0.5;
iter_max = 10;
WKNKN = 0;
first = 0;
CVS = 1;

[A_cos_com]  = DLapRLS_C1(K_COM1,K_COM2,y,alpha1,alpha2,beta,lamda_1,lamda_2,iter_max,0);
A_cos_com = A_cos_com - y;
[a,b]=sort(A_cos_com(:),'descend');
[row,col]=ind2sub(size(A),b);