function x_filt = A_for_dk_xl_2D( xpad, d, size_z )

x_filt = zeros(size_z);
for k = 1 : size_z(3)
    %!!!NOTE: compensate rotating filter in conv2()
    x_filt(:,:,k) = conv2( xpad(:,:), rot90(d(:,:,k),2), 'valid' );
end

return;