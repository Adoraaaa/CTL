function [f, Ola, OlD, OlZ] = CTLNet_ext(TrnData_ImgCell, CTL, D)


NumImg = length(TrnData_ImgCell);
ImgIdx = (1:NumImg)';

% D = cell(CTL.layers,1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [~, xpad] = PadFunc(InImg, psf_radius);

%D = cell(CTL.layers,K)
%K = CTL.sizeD(2,3);     %number of filters;

%% 训练%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    OutImgCell = TrnData_ImgCell; 
    f = cell(NumImg,1); % compute the PCANet training feature one by one 
    Ola=[]; OlD=[]; OlZ=[];
%%%%%%%%%%%%% cell转换为mat，作为每层的输入 %%%%%%%%%%%%%%%%%    
    for i = 1:CTL.layers
         L=length(OutImgCell);
         InImg=zeros(CTL.ImgSize1,CTL.ImgSize2,length(TrnData_ImgCell));
         for l = 1:L
             data=cell2mat(OutImgCell(l));
             data=data-repmat(mean(mean(data)), CTL.ImgSize1,CTL.ImgSize2);
             InImg(:,:,l)=data;
         end
%%%%%%%%%%%%%%%%% 从cell D取出每层的d %%%%%%%%%%%%%%%%%%%%%%%%
        size_kernel=CTL.sizeD(i,:);
        d = zeros(size_kernel); 
        for k = 1:CTL.sizeD(i,3)   %filter normalization
            d(:,:,k) = D{i,k};
        end
        
%%%%%%%%%%%%%%%%%% 输出特征为cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [OutImgCell, ImgIdx, Oa, OD, OZ]=CTL_ext_log(InImg, size_kernel, ImgIdx, d, CTL.afa(i), CTL.numiter2(i),...
            CTL.stepZ1(i), i, CTL.delta);
%         for n = 1:size_kernel(3)
%             D{i,n}=d(:,:,n); % 第n层训练出的d,排列存在cell D的第n行
%         end
        
        if i == CTL.layers
            for idx = 1:NumImg
                OutImgIndex = ImgIdx==idx; % select feature maps corresponding to image "idx"
                ImgIdx_i = ones(sum(OutImgIndex),1);
                OutImg_i = OutImgCell(OutImgIndex);
                
                [f{idx}] = HashingHist(CTL,ImgIdx_i,OutImg_i); % compute the feature of image "idx"
%                 OutImg(OutImgIndex) = cell(sum(OutImgIndex),1);
            end
        end
        f = sparse([f{:}]);
        
        Ola=[Ola;Oa];
        OlD=[OlD;OD];
        OlZ=[OlZ;OZ];      
    end
  
    
    
    
end