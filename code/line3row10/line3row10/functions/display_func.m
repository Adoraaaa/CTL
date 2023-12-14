function [] = display_func(iterate_fig, filter_fig, d, xfilt_d, z, x, psf_radius, iter)
sumxfilt_d=zeros(size(xfilt_d));
sumz=zeros(size(z));
for i=1:size(z,4)
    for j=1:size(z,3)
    sumxfilt_d(:,:,i)=sumxfilt_d(:,:,i)+xfilt_d(:,:,j,i)/size(z,3);
    sumz(:,:,i)=sumz(:,:,i)+z(:,:,j,i)/size(z,3);
    end
end
    figure(iterate_fig);clf
    subplot(4,6,1), imshow(x(:,:,2),[]); axis image; colormap gray; title(sprintf('Local iterate %d',iter));
    
    subplot(4,6,2), imshow(xfilt_d(:,:,1,2),[]); axis image; colormap gray; title('Filt. img');
    subplot(4,6,3), imshow(xfilt_d(:,:,6,2),[]); axis image; colormap gray;
    subplot(4,6,4), imshow(xfilt_d(:,:,11,2),[]); axis image; colormap gray;
    subplot(4,6,5), imshow(xfilt_d(:,:,16,2),[]); axis image; colormap gray;
    subplot(4,6,6), imshow(xfilt_d(:,:,21,2),[]); axis image; colormap gray;

    subplot(4,6,7), imshow(sumz(:,:,2),[]); axis image; colormap gray;
        
    subplot(4,6,8), imshow(z(:,:,1,2),[]); axis image; colormap gray; title('Spar. code'); 
    subplot(4,6,9), imshow(z(:,:,6,2),[]); axis image; colormap gray;
    subplot(4,6,10), imshow(z(:,:,11,2),[]); axis image; colormap gray;
    subplot(4,6,11), imshow(z(:,:,16,2),[]); axis image; colormap gray;
    subplot(4,6,12), imshow(z(:,:,21,2),[]); axis image; colormap gray;
    
    subplot(4,6,13), imshow(x(:,:,4),[]); axis image; colormap gray; title(sprintf('Local iterate %d',iter));
    
    subplot(4,6,14), imshow(xfilt_d(:,:,1,4),[]); axis image; colormap gray; title('Filt. img');
    subplot(4,6,15), imshow(xfilt_d(:,:,6,4),[]); axis image; colormap gray;
    subplot(4,6,16), imshow(xfilt_d(:,:,11,4),[]); axis image; colormap gray;
    subplot(4,6,17), imshow(xfilt_d(:,:,16,4),[]); axis image; colormap gray;
    subplot(4,6,18), imshow(xfilt_d(:,:,21,4),[]); axis image; colormap gray;

    subplot(4,6,19), imshow(sumz(:,:,4),[]); axis image; colormap gray;
     
    subplot(4,6,20), imshow(z(:,:,1,4),[]); axis image; colormap gray; title('Spar. code'); 
    subplot(4,6,21), imshow(z(:,:,6,4),[]); axis image; colormap gray;
    subplot(4,6,22), imshow(z(:,:,11,4),[]); axis image; colormap gray;
    subplot(4,6,23), imshow(z(:,:,16,4),[]); axis image; colormap gray;
    subplot(4,6,24), imshow(z(:,:,21,4),[]); axis image; colormap gray;
    drawnow;

    figure(filter_fig);clf
    sqr_k = ceil(sqrt(size(d,3)));
    pd = 1;
    d_disp = zeros( sqr_k * [size(d,1)+ pd, size(d,2) + pd] + [pd, pd]);
    for j = 0:size(d,3) - 1
        d_curr = d(:,:,j+1);
        d_disp( floor(j/sqr_k) * (size(d_curr,1) + pd) + pd + (1:size(d_curr,1)) , mod(j,sqr_k) * (size(d_curr,2) + pd) + pd + (1:size(d_curr,2)) ) = d_curr;
    end
    imagesc(d_disp); colormap gray; axis image; colorbar; 
    title(sprintf('Filter iterate %d',iter));
    drawnow;
        
return;