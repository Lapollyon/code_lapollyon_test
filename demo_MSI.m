%% =================================================================
% This script runs the NL-FCTN decomposition-based TC method
%
% More details can be found in [1]
% [1] Wen-Jie Zheng, Xi-Le Zhao*, Yu-Bang Zheng*, and Zhi-Feng Pang, 
%     Nonlocal Patch-Based Fully-Connected Tensor Network Decomposition
%     for Multispectral Image Inpainting, IEEE Geoscience and 
%     Remote Sensing Letters, vol. 19, pp. 1-5, 2022.


%% =================================================================
%clc;
clear;
close all;
addpath(genpath('data'));
addpath(genpath('lib'));

EN_FCTN_TC        = 1;
EN_NL_FCTN_TC     = 1;
methodname  = {'Observed', 'FCTN-TC','NL-FCTN-TC'};
Mnum = length(methodname);

%% Load initial data
load('Omsi_4.mat')
X = Omsi_4;
if max(X(:))>1
    X = X/max(X(:));
end

%% evaluation indexes
Re_tensor  =  cell(Mnum,1);
psnr       =  zeros(Mnum,1);
ssim       =  zeros(Mnum,1);
sam        =  zeros(Mnum,1);
time       =  zeros(Mnum,1);

%% Sampling with random position
sample_ratio = 0.05;
fprintf('=== The sample ratio is %4.2f ===\n', sample_ratio);
Y_tensorT    = X;
Ndim         = ndims(Y_tensorT);
Nway         = size(Y_tensorT);
rand('seed',2);
index        = find(rand(prod(Nway),1)<sample_ratio);
Omega        = zeros(Nway);
Omega(index) = 1;
Omega        = Omega>0;
F            = zeros(Nway);
F(Omega)     = Y_tensorT(Omega);

%%
i  = 1;
Re_tensor{i} = F;
[psnr(i), ssim(i), sam(i)] = HSIQA(Y_tensorT*255, Re_tensor{i}*255);
enList = 1;

%% Perform  algorithms

%% Use FCTN_TC
i = i+1;
if EN_FCTN_TC
    % initialization of the parameters
    % Please refer to our paper to set the parameters
    opts=[];
    opts.max_R = [0,    7,   7;
                  0,    0,   7;
                  0,    0,   0];
    opts.R     = [0,    2,   2;
                  0,    0,   2;
                  0,    0,   0];
    opts.tol   = 1e-5;
    opts.maxit = 1000;
    opts.rho   = 0.1;
    opts.Xtrue = Y_tensorT;
    %%%%%
    fprintf('\n');
    disp(['performing ',methodname{i}, ' ... ']);
    t0= tic;
    [Re_tensor{i}, ~] = inc_FCTN_TC(F, Omega, opts);
    time(i) = toc(t0);
    [psnr(i), ssim(i), sam(i)] = HSIQA(Y_tensorT*255, Re_tensor{i}*255);
    enList = [enList,i];
end

%% Use NL_FCTN_TC
i = i+1;
if EN_NL_FCTN_TC
    opts               =   [];
    opts.patsize       =   6;
    opts.step          =   opts.patsize-1;
    opts.deNoisingIter =   1;
    opts.patnum        =   30;
    
    optsP              =   [];
    optsP.rho          =   0.01;
    optsP.maxit        =   100;
    optsP.R            =   [0,  6,  2,  4;
                            0,  0,  2,  4;
                            0,  0,  0,  2;
                            0,  0,  0,  0];
    fprintf('\n');
    disp(['performing ',methodname{i}, ' ... ']);
    t0= tic;
    Re_tensor{i} = NL_FCTN_TC(Re_tensor{i-1}, Y_tensorT, Omega, opts, optsP);
    time(i) = toc(t0);
    [psnr(i), ssim(i), sam(i)] = HSIQA(Y_tensorT*255, Re_tensor{i}*255);
    enList = [enList,i];
end

%% Show result
fprintf('\n');
fprintf('================== Result ==================\n');
fprintf(' %10.10s    %5.4s      %5.4s    %5.4s    \n', 'method', 'PSNR', 'SSIM', 'SAM');
for i = 1:length(enList)
    fprintf(' %10.10s    %5.4f    %5.4f    %5.4f    \n',...
        methodname{enList(i)},psnr(enList(i)), ssim(enList(i)), sam(enList(i)));
end
fprintf('================== Result ==================\n');
figure,
showMSIResult(Re_tensor,Y_tensorT,min(Y_tensorT(:)),max(Y_tensorT(:)),methodname,enList,1,Nway(3))