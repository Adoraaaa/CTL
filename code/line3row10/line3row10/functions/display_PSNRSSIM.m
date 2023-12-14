function [PSNR,SSIM,info]=display_PSNRSSIM(X,Y,i,info)
for k=1:size(X,3)
    Xtemp=X(:,:,k);
    Ytemp=Y(:,:,k);
    info.PSNRoutput(k,i)=10*log10(1^2/mean((Xtemp(:)-Ytemp(:)).^2));
    info.SSIMoutput(k,i)=ssim(Xtemp,Ytemp,'DynamicRange',255);
end
PSNR=mean(info.PSNRoutput(:,i));
SSIM=mean(info.SSIMoutput(:,i));

fprintf('PSNR = %.2f, SSIM = %.2f\n', ...
    PSNR, SSIM);
return;