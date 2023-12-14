function [OutImg,D,OutImgIdx,Oa,OD,OZ]= CTL_train_log(X,size_kernel,InImgIdx,D,afa,numiter,stepD,stepZ,layer,delta)
%% 参数%%%%%%%%%%%%%%%%%%%%%%
% if size_kernel(1)*size_kernel(2) > size_kernel(3)
%     error('The number of filters must be equal or larger than size of the filters.');
% end
K = size_kernel(3);     %number of filters
L = size(X,3);          %number of training images

%variables for filters
psf_radius = floor( [size_kernel(1)/2, size_kernel(2)/2] );

%dimensions of training images and sparse codes
size_x = [size(X,1), size(X,2), L];
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
    D_old=D; Z_old=Z;
 %%%%%%%%更新Z%%%%%%%%%%%%%%%%%%%%%%%%%
    d_Z =Z - A_for_dk_xl(xpad,D,size_z);
    Z = Z_old - stepZ * d_Z;
    Z = ProxSparse_log( Z, Z_old, afa, delta, stepZ );
%     Z = ProxSparseL1(Z,afa);
 %%%%%%%%更新D%%%%%%%%%%%%%%%%%%%%%%%%
%     d_D = zeros(size_kernel);
%     for k = 1 : K
%         A_kd_p=A_for_dk( xpad, D(:,:,k), size_x );
%         d_D(:,:,k)=Ah_for_dk( xpad, A_kd_p - reshape(Z(:,:,k,:), [size_z(1),size_z(2),size_z(4)]), size_x, size_kernel ); 
%     end
% %     for l = 1:L
% %         for k = 1:K
% %             d_D(:,:,k)=d_D(:,:,k) + conv2(xpad(:,:,l),rot90((conv2(xpad(:,:,l),D(:,:,k),'valid')-Z(:,:,k,l)),2),'valid');
% %         end
% %     end
%     Adnu = D - stepD * d_D;
%     
%     D = ProxFilterTightFrame( Adnu, size_kernel ); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %% Manifold update d_p
      for k = 1 : K
            param.size_kernel = [5, 5, 5];
            mx= param.size_kernel(1);
            my= param.size_kernel(2);
            n = param.size_kernel(3);
            d = randn(size_kernel);
            A_kl = @(u) A_for_dk_xl( xpad, u, size_z );
            d_old=d; 
            d_p=d;
            Omega_2D=(reshape(d_p,mx*my,n)');
            OOT=Omega_2D*Omega_2D';OTO=Omega_2D'*Omega_2D;
            %% update with d_p and Omega
                opts.record = 0;
                opts.mxitr  = 100;
                opts.xtol = 1e-5;
                opts.gtol = 1e-5;
                opts.ftol = 1e-5;
                opts.tau = 1e-5; opts.eta = 0.2;
                opts.rho = 1e-6; opts.nt = 5; opts.gamma = 0.85;
                opt.iscomplex=0; opts.record=1; opts.diter=2;
                psf_radius = floor( [mx-1, my-1] );
                cdsize=psf_radius-floor( [mx/2, my/2] );
                size_z = [size(X,1), size(X,2), K, L];
                size_zcd=size_z;
                tau_mo = 1;
                tau_sc = 1;
                weight = 1-eps;
                param1=10^(3/15)*1e-5;
                param.lambda = 1+eps;
           %% Initial function value and gradient
                % prepare for iterations
                G=zeros(param.size_kernel);
                for k = 1 : n
                    OmegakY=A_for_dk( xpad, d_p(:,:,k), size_x );
                    G(:,:,k)=Ah_for_dk( xpad, OmegakY- reshape(Z(:,:,k,:), [size_zcd(1),size_zcd(2),size_zcd(4)]), size_x, param.size_kernel );
                end
                G=reshape(G,[mx*my,n])';
                A_dkxl = A_kl(d); F = 1/2 * norm(A_dkxl(:) - Z(:), 2)^2;
                GX = G'*Omega_2D; out.nfe = 1;
                GXT = G*Omega_2D';  H = 0.5*(GXT - GXT');  RX = H*Omega_2D;
                dtX = G - Omega_2D*GX; nrmG  = norm(dtX, 'fro');
                Q = 1; Cval = F; tau = opts.tau;
                %% start searching
                for itr = 1 : opts.mxitr
                    XP = Omega_2D; FP = F; GP = G; dtXP = dtX;
                    % scale step size
                    nls = 1; deriv = opts.rho*nrmG^2; %deriv
                    while 1
                        % calculate G, F,
                        [Omega_2D, infX] = linsolve(eye(n) + tau*H, XP - tau*RX);
                        d_p=reshape(Omega_2D',param.size_kernel);
                        G=zeros(param.size_kernel);
                        for k = 1 : n
                            OmegakY=A_for_dk( xpad, d_p(:,:,k), size_x );
                            G(:,:,k)=Ah_for_dk( xpad, OmegakY- reshape(Z(:,:,k,:), [size_zcd(1),size_zcd(2),size_zcd(4)]), size_x, param.size_kernel );
                        end
                        G=reshape(G,[mx*my,n])';
                        A_dkxl = A_kl(d_p); F = 1/2 * norm(A_dkxl(:) - Z(:), 2)^2;
                        out.nfe = out.nfe + 1;
                        if F <= Cval - tau*deriv || nls >= 5
                            break;
                        end
                        tau = opts.eta*tau;          nls = nls+1;
                    end
                    GX = G'*Omega_2D;
                    GXT = G*Omega_2D';  H = 0.5*(GXT - GXT');  RX = H*Omega_2D;
                    dtX = G - Omega_2D*GX;    nrmG  = norm(dtX, 'fro');
                    
                    S = Omega_2D - XP;         XDiff = norm(S,'fro')/sqrt(n);
                    tau = opts.tau; FDiff = abs(FP-F)/(abs(FP)+1);
                    
                    YY = dtX - dtXP;     SY = abs(sum(sum(S.*YY)));
                    if mod(itr,2)==0; tau = sum(sum(S.*S))/SY;
                    else tau  = SY/sum(sum(YY.*YY)); end
                    
                    tau = max(min(tau, 1e20), 1e-20);
                    
                    if (opts.record >= 1)
                        if mod(itr,opts.diter)==1
                            fprintf('itr %4d, tao %3.2e F %4.3e nrmG %3.2e Xdiff %3.2e Fdiff %3.2e nls %2d\n', ...
                                itr, tau, F, nrmG, XDiff, FDiff, nls);
                        end
                    end
                    
                    crit(itr,:) = [nrmG, XDiff, FDiff];
                    mcrit = mean(crit(itr-min(opts.nt,itr)+1:itr, :),1);
                    if  ( XDiff < opts.xtol && FDiff < opts.ftol ) || nrmG < opts.gtol || all(mcrit(2:3) < 10*[opts.xtol, opts.ftol])
                        if itr <= 2
                            opts.xtol = 0.1*opts.xtol;
                            opts.gtol = 0.1*opts.gtol;
                        else
                            out.msg = 'converge';
                            break;
                        end
                    end
                    
                    Qp = Q; Q = opts.gamma*Qp + 1;
                end
                if (opts.record >= 1)
                    fprintf('itr %4d, tao %3.2e F %4.3e nrmG %3.2e Xdiff %3.2e Fdiff %3.2e nls %2d\n', ...
                        itr, tau, F, nrmG, XDiff, FDiff, nls);
                end
                out.feasi = norm(Omega_2D'*Omega_2D-eye(mx*my),'fro');
                if  out.feasi > 1e-5
                    Omega_2D = MGramSchmidt(Omega_2D);
                    out.nfe = out.nfe + 1;
                    out.feasi = norm(Omega_2D'*Omega_2D-eye(mx*my),'fro');
                end
                
                out.nrmG = nrmG;
                out.itr = itr;
                
                Omega_2D=normr(Omega_2D);
                d=reshape(Omega_2D',param.size_kernel);
                
                % UPDATE: Momentum coeff.
                tau_old = tau_mo;
                tau_mo = ( 1 + sqrt(1 + 4*tau_mo^2) ) / 2;
                %Extrapolation with momentum!
                w_d = min( (tau_old - 1)/tau_mo, weight*0.5*(param.lambda-1)/(param.lambda+1) );
                d_p = d + w_d .* (d - d_old);
                D=d_p;
         end
            
            
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
    
% function [bpad] = PadFunc(b, psf_radius)
%     
% M = padarray(ones(size(b)), [psf_radius(1), psf_radius(2), 0], 0, 'both');    %mask
% %%%circular padding
% % bpad = padarray(b, [psf_radius, psf_radius, 0], 'circular', 'both');
% %%%reflective padding
% % bpad = padarray(b, [psf_radius(1), psf_radius(2), 0], 'symmetric', 'both');     
% %%%%zero padding
% bpad = padarray(b, [psf_radius, psf_radius, 0], 0, 'both');              
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





