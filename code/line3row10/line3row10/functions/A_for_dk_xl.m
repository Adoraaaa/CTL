function x_filt = A_for_dk_xl( xpad, d, size_z )

%x_filt = zeros(size_z);
for l = 1 : size_z(4)
    for k = 1 : size_z(3)
        %!!!NOTE: compensate rotating filter in conv2()
        x_filt(:,:,k,l) = conv2( xpad(:,:,l), rot90(d(:,:,k),2), 'valid' );
    end
end

return;