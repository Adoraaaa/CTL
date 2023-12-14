%%%%%%%%%%%%%%%%%%%% MISC %%%%%%%%%%%%%%%%%%%%

function [f_val, f_d, sparsity] = objectiveFunction( z, A_dkxl, afa )

    %Dataterm
    f_d = 1/2 * norm(A_dkxl(:) - z(:), 2)^2;
    %Regularizer
    f_z = afa * norm(A_dkxl(:) - z(:), 2)^2;

    %Function val
    f_val = f_d;% + alpha*f_z;
    
    %Sparsity
    sparsity = 100*nnz(z)/numel(z);
    
return;