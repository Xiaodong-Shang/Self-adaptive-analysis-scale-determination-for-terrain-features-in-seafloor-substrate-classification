clc
clear
%% data input
load xydep
Pnum=size(xydep,1);

%% Neighborhood determination and depth feature extraction
k_min=25;
k_max=45;
delta_k=2;
n1=200000;
n2=400000;
n3=600000;
n4=800000;
n5=1000000;
n6=1200000;
n7=1400000;
n8=1600000;
n9=Pnum;
%optNESS_EntropyBased
[DepFea1 opt_nn_size1] = optNESS_EntropyBased(xydep,k_min,k_max,delta_k,1,n1); % Neighborhood determination
[DepFea2 opt_nn_size2] = optNESS_EntropyBased(xydep,k_min,k_max,delta_k,n1+1,n2); % Neighborhood determination
[DepFea3 opt_nn_size3] = optNESS_EntropyBased(xydep,k_min,k_max,delta_k,n2+1,n3); % Neighborhood determination
[DepFea4 opt_nn_size4] = optNESS_EntropyBased(xydep,k_min,k_max,delta_k,n3+1,n4); % Neighborhood determination
[DepFea5 opt_nn_size5] = optNESS_EntropyBased(xydep,k_min,k_max,delta_k,n4+1,n5); % Neighborhood determination
[DepFea6 opt_nn_size6] = optNESS_EntropyBased(xydep,k_min,k_max,delta_k,n5+1,n6); % Neighborhood determination
[DepFea7 opt_nn_size7] = optNESS_EntropyBased(xydep,k_min,k_max,delta_k,n6+1,n7); % Neighborhood determination
[DepFea8 opt_nn_size8] = optNESS_EntropyBased(xydep,k_min,k_max,delta_k,n7+1,n8); % Neighborhood determination
[DepFea9 opt_nn_size9] = optNESS_EntropyBased(xydep,k_min,k_max,delta_k,n8+1,n9); % Neighborhood determination


