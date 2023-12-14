function M =  majorMat_Ak_diag( xpad, size_x, size_kernel )
%Scaled identity majorization matrix in 
%Lem. 5.2 of DOI: 10.1109/TIP.2019.2937734

AtA_symb = zeros(size_kernel(1), size_kernel(2));
for l = 1 : size_x(3)
    P1x = xpad( 1 : 1+size_x(1)-1, 1 : 1+size_x(2)-1, l );
    for r2 = 1 : size_kernel(2)
        for r1 = 1 : size_kernel(1)
            Prx = xpad( r1 : r1+size_x(1)-1, r2 : r2+size_x(2)-1, l );
            AtA_symb(r1, r2) = AtA_symb(r1, r2) + Prx(:)' * P1x(:);
        end
    end
end

M = ( abs(AtA_symb(:))' * ones(size_kernel(1)*size_kernel(2),1) ) .* ...
    ones(size_kernel(1), size_kernel(1));

return;