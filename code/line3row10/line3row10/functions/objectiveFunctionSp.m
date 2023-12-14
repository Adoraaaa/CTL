%%%%%%%%%%%%%%%%%%%% MISC %%%%%%%%%%%%%%%%%%%%

function [f_val, f_d, sparsity, spAY] = objectiveFunctionSp( z, A_dkxl,sptol, alpha )

    %Dataterm
    f_d = 1/2 * norm(A_dkxl(:) - z(:), 2)^2;
    %Regularizer
    f_z = nnz(z);

    %Function val
    f_val = f_d;% + alpha*f_z;
    
    %Sparsity
    sparsity = 100*f_z/numel(z);
    nzAY=sum(A_dkxl(:)>sptol);
    spAY= 100*nzAY/numel(z);
    
return;