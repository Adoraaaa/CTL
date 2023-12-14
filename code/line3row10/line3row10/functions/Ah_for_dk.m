function Ahu = Ah_for_dk( xpad, u, size_x, size_kernel )

Ahu = zeros(size_kernel(1), size_kernel(2)); 
for l = 1: size_x(3)
    Ahu = Ahu + conv2( xpad(:,:,l), rot90(u(:,:,l),2), 'valid');    
end

return;