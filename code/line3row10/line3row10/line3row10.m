clear all; close all; clc; 
addpath('./Utils');
addpath('./Liblinear');
make; 
tic;
%Fixed random initialization
% s = RandStream('mt19937ar','Seed',1);
% RandStream.setGlobalStream(s);

TrnSize = 10000; 
ImgSize1 = 27; 
ImgSize2 = 20;
ImgFormat = 'gray'; %'color' or 'gray'

%% Loading data 
load('AR_27x20');  
load('labels');
addpath('./functions');
%% 444444 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A=[];
% for n = 1:length(X)
%     W=X{n};
%     a=ones(size(W,2),1)*n;
%     A=[A;a];
% end
% C=[];
% for m = 1:length(X)
%     W=X{m};
%     c = mat2imgcell(W,ImgSize1,ImgSize2,ImgFormat);
%     C=[C;c];
% end
% Labels=A;
% dataCell=C;
% for i = 1:length(dataCell)
%    A=dataCell{i};
%    B=A(1:4:ImgSize1,1:4:ImgSize2);
%    dataCell{i}=B;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%111111111111111111111111111111111111111111111111111111111111111111111&&
CTL.afa = [0.03 0];
time=0;
result=[];
for m = 1:1
    CTL.afa(1)=CTL.afa(1)+0.01;
    CTL.afa(2)=5e-3;
    for mm = 1:1
    CTL.afa(2)=CTL.afa(2)+1e-3;
    e=zeros(1,3);
    time=time+1;
%% %%%%%%%%%%%%%%训练集，测试集%%%%%%%%%%%%%%%%%%%
num=randperm(2600);

TrnData_ImgCell = dataCell(num(1:2000));
TrnLabels=labels(num(1:2000));

TestData_ImgCell=dataCell(num(2001:end));
TestLabels=labels(num(2001:end));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% ===== Reshuffle the training data =====
% Randnidx = randperm(size(mnist_train,1)); 
% mnist_train = mnist_train(Randnidx,:); 
% =======================================

nTestImg = length(TestLabels);
nTrnImg = length(TrnLabels);
%%  parameters
CTL.layers = 1;
CTL.delta=5e-2;
CTL.sizeD = [5,5,5; 5,5,5];
CTL.afa = [0.08 0.005];
CTL.numiter1 = [14 14];
CTL.numiter2 = [11 11];
CTL.stepD = [5e-14 4e-12];
CTL.stepZ1 = [8e-2 5e-2];% step_train
CTL.stepZ2 = [8e-2 5e-2];% step_exc
CTL.HistBlockSize = [3 3]; 
CTL.BlkOverLapRatio = 0;
CTL.Pyramid = [];
CTL.ImgSize1 = ImgSize1;
CTL.ImgSize2 = ImgSize2;
CTL.NumFilters = [CTL.sizeD(1,3) CTL.sizeD(2,3)];


%Fixed random initialization
% s = RandStream('mt19937ar','Seed',1);
% RandStream.setGlobalStream(s);


%% %%%%%%%%%%%%%%%%%%%%% 训练 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n ====== CTL Training ======= \n')
% TrnData_ImgCell = mat2imgcell(TrnData,ImgSize,ImgSize,ImgFormat); % convert columns in TrnData to cells 
% clear TrnData; 


[ftrain, D, Ola, OlD, OlZ]=CTLNet_train(TrnData_ImgCell,CTL);

for p = 1:CTL.layers
    x=(10:1:CTL.numiter1(p));
    subplot(CTL.layers,3,3*(p-1)+1)
    plot(x,Ola(p,x),'b');
    title('obj');
    hold on;
    subplot(CTL.layers,3,3*(p-1)+2)
    plot(x,OlD(p,x),'b');
    title('val_loss');
    hold on;
    subplot(CTL.layers,3,3*(p-1)+3)
    plot(x,OlZ(p,x),'b');
    title('val_z');
    hold on;
end

% [ftrain D BlkIdx] = PCANet_train(TrnData_ImgCell,CTL); % BlkIdx serves the purpose of learning block-wise DR projection matrix; e.g., WPCA
CTLNet_TrnTime = toc;
clear TrnData_ImgCell; 


fprintf('\n ====== Training Linear SVM Classifier ======= \n')

models = train(TrnLabels, ftrain', '-s 1 -q'); % we use linear SVM classifier (C = 1), calling libsvm library
LinearSVM_TrnTime = toc;
% clear ftrain; 

%% %%%%%%%%%%%%%%%% 特征提取，测试 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TestData_ImgCell = mat2imgcell(TestData,ImgSize,ImgSize,ImgFormat); % convert columns in TestData to cells 
% clear TestData; 

fprintf('\n ====== CTLNet Testing ======= \n')

nCorrRecog = 0;
RecHistory = zeros(nTestImg,1);


[ftest, Ola, OlD, OlZ]=CTLNet_ext(TestData_ImgCell, CTL, D);


for p = 1:CTL.layers
    x=(10:1:CTL.numiter2(p));
    subplot(CTL.layers,3,3*(p-1)+1)
    plot(x,Ola(p,x),'r');
    title('obj');
    hold on;
    subplot(CTL.layers,3,3*(p-1)+2)
    plot(x,OlD(p,x),'r');
    title('val_loss');
    hold on;
    subplot(CTL.layers,3,3*(p-1)+3)
    plot(x,OlZ(p,x),'r');
    title('val_z');
    hold on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for idx = 1:1:nTestImg
    ftest0=ftest(:,idx);
    [xLabel_est, accuracy, decision_values] = predict(TestLabels(idx),...
        sparse(ftest0'), models, '-q'); % label predictoin by libsvm
   
    if xLabel_est == TestLabels(idx)
        RecHistory(idx) = 1;
        nCorrRecog = nCorrRecog + 1;
    end
    
%     if 0==mod(idx,nTestImg/10)
%         fprintf('已测试样本数： %d ，精度： %.2f%%. \n',...
%             [idx 100*nCorrRecog/idx]); 
%     end 
    
    TestData_ImgCell{idx} = [];
    
end
Averaged_TimeperTest = toc/nTestImg;
Accuracy = nCorrRecog/nTestImg; 
ErRate = Accuracy;

%% Results display
fprintf('\n Results of CTLNet, followed by a linear SVM classifier');
% fprintf('\n     PCANet training sample: );
fprintf('\n      训练集样本总数: %.2f',nTrnImg);
fprintf('\n      测试集样本总数: %.2f',nTestImg);
fprintf('\n      CTL层数: %.2f',CTL.layers);
% fprintf('\n     Linear SVM training time: %.2f secs.', LinearSVM_TrnTime);
fprintf('\n     Testing Accuracy: %.2f%%', 100*ErRate);
% fprintf('\n     Average testing time %.2f secs per test sample. \n\n',Averaged_TimeperTest);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
e(1)=CTL.afa(1);
e(2)=CTL.afa(2);
e(3)=Accuracy;
result=[result;e];


    end
end
fprintf(datestr(now));
toc;

    