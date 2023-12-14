
function [] = displayx_func(iterate_fig, x, xre, iter)
    figure(iterate_fig);
    subplot(2,6,1), imshow(x(:,:,1),[]); axis image; colormap gray; title(sprintf('Local iterate %d',iter));
    subplot(2,6,2), imshow(x(:,:,3),[]); axis image; colormap gray; title('Filt. img');
    subplot(2,6,3), imshow(x(:,:,5),[]); axis image; colormap gray;
    subplot(2,6,4), imshow(x(:,:,7),[]); axis image; colormap gray;
    subplot(2,6,5), imshow(x(:,:,9),[]); axis image; colormap gray;
    subplot(2,6,6), imshow(x(:,:,10),[]); axis image; colormap gray;

    subplot(2,6,7), imshow(xre{1},[]); axis image; colormap gray;
    subplot(2,6,8), imshow(xre{3},[]); axis image; colormap gray; title('Spar. code'); 
    subplot(2,6,9), imshow(xre{5},[]); axis image; colormap gray;
    subplot(2,6,10), imshow(xre{7},[]); axis image; colormap gray;
    subplot(2,6,11), imshow(xre{9},[]); axis image; colormap gray;
    subplot(2,6,12), imshow(xre{10},[]); axis image; colormap gray;
    drawnow;
    
     figure(iterate_fig+1);
    subplot(2,2,1), imshow(x(:,:,1),[]); axis image; colormap gray; title(sprintf('Local iterate %d',iter));
    subplot(2,2,2), imshow(x(:,:,3),[]); axis image; colormap gray; title('Orig. img');
    
    subplot(2,2,3), imshow(xre{1},[]); axis image; colormap gray;
    subplot(2,2,4), imshow(xre{3},[]); axis image; colormap gray; title('Reco. img'); 
    drawnow;
return;