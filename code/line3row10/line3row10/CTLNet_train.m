function [f, D, Ola, OlD, OlZ] = CTLNet_train(TrnData_ImgCell, CTL)


NumImg = length(TrnData_ImgCell);
ImgIdx = (1:NumImg)';

D = cell(CTL.layers,1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [~, xpad] = PadFunc(InImg, psf_radius);

K = CTL.sizeD(2,3);     %number of filters
D = cell(CTL.layers,K);

%% ѵ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    OutImgCell = TrnData_ImgCell; 
    f = cell(NumImg,1); % compute the PCANet training feature one by one 
    Ola=[]; OlD=[]; OlZ=[];
%%%%%%%%%%%%% cellת��Ϊmat����Ϊÿ������� %%%%%%%%%%%%%%%%%    
    for i = 1:CTL.layers
         L=length(OutImgCell);
         InImg=zeros(CTL.ImgSize1,CTL.ImgSize2,length(TrnData_ImgCell));
         for l = 1:L
             data=cell2mat(OutImgCell(l));
             InImg(:,:,l)=data-repmat(mean(mean(data)), CTL.ImgSize1,CTL.ImgSize2);
         end
%%%%%%%%%%%%%%%%%%%%%%% ��ʼ��d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        size_kernel=CTL.sizeD(i,:);
        d = randn(size_kernel);    %no need to normalize
        d(:,:,1) = 1;  %set the first filter as a DC filter          
%         for k = 1:CTL.sizeD(i,3)   %filter normalization
%             d(:,:,k) = d(:,:,k) ./ (sqrt(size_kernel(1)*size_kernel(2))*norm(d(:,:,k),'fro'));
%         end
%%%%%%%%%%%%%%%%%%%% ÿ��ѵ�����������Ϊcell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [OutImgCell, d, ImgIdx, Oa,OD,OZ]=CTL_train_log(InImg, size_kernel, ImgIdx, d, CTL.afa(i), CTL.numiter1(i), CTL.stepD(i),...
            CTL.stepZ1(i), i, CTL.delta);
        for n = 1:size_kernel(3)
            D{i,n}=d(:,:,n); % ��n��ѵ������d,���д���cell D�ĵ�n��
        end
        
        if i == CTL.layers
            for idx = 1:NumImg
                OutImgIndex = ImgIdx==idx; % select feature maps corresponding to image "idx"
                ImgIdx_i = ones(sum(OutImgIndex),1);
                OutImg_i = OutImgCell(OutImgIndex);
                
                [f{idx}] = HashingHist(CTL,ImgIdx_i,OutImg_i); % compute the feature of image "idx"
%                 OutImgCell(OutImgIndex) = cell(sum(OutImgIndex),1);
            end
                    f = sparse([f{:}]);
        end
%         f = sparse([f{:}]);
        
        Ola=[Ola;Oa];
        OlD=[OlD;OD];
        OlZ=[OlZ;OZ];      
    end
  
    
    
    
end