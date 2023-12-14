function x_filt = A_for_dk_xl_shape( x, d, size_z, shape )

%x_filt = zeros(size_z);
for l = 1 : size_z(4)
    for k = 1 : size_z(3)
        %!!!NOTE: compensate rotating filter in conv2()
        x_filt(:,:,k,l) = conv2( x(:,:,l), rot90(d(:,:,k),2), shape );
    end
end

return;