%%%%%%%%%%%%%%%%%%%% MISC %%%%%%%%%%%%%%%%%%%%

function [f_val, f_d, sparsity, spZ, spZ1] = objectiveFunctionSptol( z, A_dkxl,sptol, alpha )

    %Dataterm
    f_d = 1/2 * norm(A_dkxl(:) - z(:), 2)^2;
    %Regularizer
    f_z = nnz(z);

    %Function val
    f_val = f_d;% + alpha*f_z;
    
    %Sparsity
    sparsity = 100*f_z/numel(z);
    nzZ=sum(abs(z(:))>sptol);
    nzZ1=sum(abs(z(:))>(sptol*10));
    N=numel(z);
    spZ= 100*nzZ/N;
    spZ1= 100*nzZ1/N;
    
return;