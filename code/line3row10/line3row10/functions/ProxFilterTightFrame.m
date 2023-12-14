%%%%%%%%%%%%%%%%%%%% Def: Proximal Operators %%%%%%%%%%%%%%%%%%%%

function d = ProxFilterTightFrame( Adnu, size_kernel )
%Solve orthogonal Procrustes problem

kernVec_size = size_kernel(1)*size_kernel(2);

AdNu = reshape(Adnu, [kernVec_size, size_kernel(3)]);
[U, ~, V] = svd(AdNu, 'econ');

D = U * V';
D = normc(D)*sqrt(kernVec_size)^(-1) ;
d = reshape( D, size_kernel );

return;