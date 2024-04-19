clc; clear; close all;

addpath(genpath('./lib'))

for tensor_num= 1
    switch tensor_num
    case 1
        load('foreman.mat');
        data_name = 'foreman';
    case 2
        load('hall.mat');
        data_name = 'hall';
    case 3 
        load('pavia.mat');
        data_name = 'pavia';
    case 4
        load('MRI.mat');
        data_name = 'MRI';
    case 5
        load('HSI.mat');
        data_name = 'HSI';
    case 6
        load('video_suzie.mat');
        data_name = 'video_suzie'; 
    end

%% Set sampling ratio
sr = 0.1;

%% Initial data
% normalized data
if max(X(:))>1
X=X/max(X(:));
end

Nway=[size(X,1), size(X,2), size(X,3)];
n1 = size(X,1); n2 = size(X,2); 
frames=size(X,3);

ratio = 0.01;%%fine tune
R=AdapN_Rank(X,ratio);
% R=Nway;
Y_tensorT = X;

p = round(sr*prod(Nway));
known = randsample(prod(Nway),p); data = Y_tensorT(known);
[known, id]= sort(known); data= data(id);
Y_tensor0= zeros(Nway); Y_tensor0(known)= data;

%% Initialization of the factor matrices X and A
for n = 1:3
    coNway(n) = prod(Nway)/Nway(n);
end
for i = 1:3
    Y0{i} = Unfold(Y_tensor0,Nway,i);
    Y0{i} = Y0{i}';
    X0{i}= rand(coNway(i), R(i));
    A0{i}= rand(R(i),Nway(i));
end
    
%% Initialization of the parameters
opts=[];
opts.maxit=2000;
opts.Ytr= Y_tensorT;
opts.tol=1e-4;
alpha=[1,1,1];
opts.alpha = alpha / sum(alpha);
rho=0.001;%%%
opts.rho1=rho;
opts.rho2=rho;
opts.rho3=rho;

opts.mu=1; %%%%%%%fine tune
opts.beta=1;%%%%%%%fine tune
opts.miter=300;
case0='LRTCDM';%
% case0='LRTCFM';%
%% Begin the comlpetion with LRTC-DM
        fprintf('\n');
        disp('Begin the comlpetion with LRTC-DM')
        t0= tic;
        [Y_LRTCDM, A, X, Out]= LRTC_DM(Y0, data, A0, X0,Y_tensor0, Nway, known, opts, case0);
        time= toc(t0);
        psnr0=PSNR(Y_LRTCDM, Y_tensorT);
end