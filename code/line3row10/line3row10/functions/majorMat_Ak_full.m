function M = majorMat_Ak_full( xpad, size_x, size_kernel, L, psf_radius, kern_xgrid, kern_ygrid )  
%Exact Hessian matrix in
%Prop. 5.1 of DOI: 10.1109/TIP.2019.2937734

%!!!NOTE: x-grid and y-grid are horizontal and vertical direction, resp.
%E.g. in matlab, X( y-grid indices, x-grid indices )
kern_xgrid_vec = kern_xgrid(:);
kern_ygrid_vec = kern_ygrid(:);

M = zeros( size_kernel(1)*size_kernel(2), size_kernel(1)*size_kernel(2) );

for k1 = 1 : size_kernel(1)*size_kernel(2)
    for k2 = 1 : k1
        
        k1x_coord = kern_xgrid_vec(k1);
        k1y_coord = kern_ygrid_vec(k1);
        
        k2x_coord = kern_xgrid_vec(k2);
        k2y_coord = kern_ygrid_vec(k2);
        
        for l = 1 : L               
            xpad_k1 = xpad( 1 + psf_radius(1) + k1y_coord : size_x(1) + psf_radius(1) + k1y_coord, ...
                1 + psf_radius(2) + k1x_coord : size_x(2) + psf_radius(2) + k1x_coord, l );

            xpad_k2 = xpad( 1 + psf_radius(1) + k2y_coord : size_x(1) + psf_radius(1) + k2y_coord, ...
                1 + psf_radius(2) + k2x_coord : size_x(2) + psf_radius(2) + k2x_coord, l );

            if k1 == k2
                M(k1, k2) = M(k1, k2) + (xpad_k1(:)'*xpad_k2(:))/2;
            else
                M(k1, k2) = M(k1, k2) + xpad_k1(:)'*xpad_k2(:);
            end            
        end
    end
end

M = M + M';

return;