function [OutImg, OutImgIdx, Oa, OD, OZ]= CTL_ext_log(X,size_kernel,InImgIdx,D,afa,numiter,stepZ,layer,delta)
%% 参数%%%%%%%%%%%%%%%%%%%%%%
% if size_kernel(1)*size_kernel(2) > size_kernel(3)
%     error('The number of filters must be equal or larger than size of the filters.');
% end
K = size_kernel(3);     %number of filters
L = size(X,3);          %number of training images

%variables for filters
psf_radius = floor( [size_kernel(1)/2, size_kernel(2)/2] );

%dimensions of training images and sparse codes
% size_x = [size(X,1), size(X,2), L];
size_z = [size(X,1), size(X,2), K, L];

paddingchoice='zeros'; %options: 'circular' 'zeros'
[~, xpad] = PadFunc(X, psf_radius, paddingchoice);

%Initialization: sparse codes
Z = A_for_dk_xl(xpad, D, size_z);        
% Z = ProxSparse_log( Z, Z_old, afa, delta, stepZ);

ProxSparseL1 = @(u, a) sign(u) .* max( abs(u)-a, 0 );
%% 梯度下降更新D，Z
Oa=[]; OD=[]; OZ=[];
for i = 1:numiter
  Z_old=Z; D_old=D;
 %%%%%%%%更新Z%%%%%%%%%%%%%%%%%%%%%%%%%
    d_Z =Z - A_for_dk_xl(xpad,D,size_z);
    Z = Z_old - stepZ * d_Z;
    Z = ProxSparse_log( Z, Z_old, afa, delta, stepZ );
%     Z = ProxSparseL1(Z,afa);
 %%%%%%%%更新D%%%%%%%%%%%%%%%%%%%%%%%%
    Ah_kA_kd_p = zeros(size_kernel);
    for k = 1 : K
        A_kd_p=A_for_dk( xpad, D(:,:,k), size_x );
        Ah_kA_kd_p(:,:,k)=Ah_for_dk( xpad, A_kd_p - reshape(Z(:,:,k,:), [size_z(1),size_z(2),size_z(4)]), size_x, size_kernel ); 
    end
    Adnu = D - stepD * Ah_kA_kd_p;
    
    D = ProxFilterTightFrame( Adnu, size_kernel );

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [obj,f_d,f_z] = objectiveFunction( Z, -d_Z, afa );
    d_relErr = norm( D(:)-D_old(:) ) / norm( D(:) );
    z_relErr = norm( Z(:)-Z_old(:) ) / norm( Z(:) );
    s_z=nnz(Z);
    sparsity = 100*s_z/numel(Z);
    
    fprintf('layer: %d ,Iter %d,Sparsity %.2f %%, Obj %d, Obj_loss %d, Obj_Z %d, D Diff %5.5g, Z Diff %5.5g\n',...
        layer,i,sparsity, obj,f_d,f_z, d_relErr, z_relErr);
    
    Oa=[Oa obj];
    OD=[OD f_d];
    OZ=[OZ f_z];
    
end
OutImg = cell(K*L,1);
cnt=0;
for i = 1:L
    for j = 1:K
        cnt=cnt+1;
        OutImg{cnt}= Z(:,:,j,i);
    end
end

OutImgIdx = kron(InImgIdx,ones(K,1)); 






%% %%%%%%%%%%%%%%%%% D正交约束 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = ProxFilterTightFrame( Adnu, size_kernel )
%Solve orthogonal Procrustes problem

kernVec_size = size_kernel(1)*size_kernel(2);

AdNu = reshape(Adnu, [kernVec_size, size_kernel(3)]);
[U, ~, V] = svd(AdNu, 'econ');

D = sqrt(kernVec_size)^(-1) * U * V';

d = reshape( D, size_kernel );

return;

%%%%%%%%%%%%%%%%%%% Log-Reg %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Z = ProxSparse_log( z, z_old, afa, delta, stepZ )
ProxSparseL1 = @(u, a) sign(u) .* max( abs(u)-a, 0 );
signX=sign(z_old);
Z_DC=afa * signX.* (1-1./(abs(z_old)+delta)) ;
temp=z + stepZ * Z_DC;
Z = ProxSparseL1( temp, afa); 
return;

function [f_val,f_d,z2] = objectiveFunction( Z, A_dkxl, afa )

    %Dataterm
    f_d = 1/2 * norm(A_dkxl(:) ,2)^2;

    z2 = afa * norm(Z(:),2)^2;

    %Function val
    f_val = f_d + z2;
    
%     %Sparsity
%     sparsity = 100*f_z/numel(z);
    
% function [M, bpad] = PadFunc(b, psf_radius)
%     
% M = padarray(ones(size(b)), [psf_radius(1), psf_radius(2), 0], 0, 'both');    %mask
% %%%circular padding
% bpad = padarray(b, [psf_radius, psf_radius, 0], 'circular', 'both');
% %%%reflective padding
% % bpad = padarray(b, [psf_radius(1), psf_radius(2), 0], 'symmetric', 'both');     
% %%%%zero padding
% % bpad = padarray(b, [psf_radius, psf_radius, 0], 0, 'both');              
%     
% return;

%%%%%%%%%%%%%%%%%%%% Def: System Operators %%%%%%%%%%%%%%%%%%%%

% function Au = A_for_dk( xpad, u, size_x )
% 
% Au = zeros(size_x);
% for l = 1 : size_x(3)
%     %!!!NOTE: compensate rotating filter in conv2()
%     Au(:,:,l) = conv2( xpad(:,:,l), rot90(u,2), 'valid' );
% end
% 
% return;
% 
% 
% function Ahu = Ah_for_dk( xpad, u, size_x, size_kernel )
% 
% Ahu = zeros(size_kernel(1), size_kernel(2)); 
% for l = 1: size_x(3)
%     Ahu = Ahu + conv2( xpad(:,:,l), rot90(u(:,:,l),2), 'valid');    
% end
% 
% return;
% 
% 
% function x_filt = A_for_dk_xl( xpad, d, size_z )
% 
% x_filt = zeros(size_z);
% for l = 1 : size_z(4)
%     for k = 1 : size_z(3)
%         %!!!NOTE: compensate rotating filter in conv2()
%         x_filt(:,:,k,l) = conv2( xpad(:,:,l), rot90(d(:,:,k),2), 'valid' );
%     end
% end
% 
% return;





