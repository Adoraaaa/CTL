
%%%%%%%%%%%%%%%%%%%% Def: System Operators %%%%%%%%%%%%%%%%%%%%

function Au = A_for_dk( xpad, u, size_x )

%Au = zeros(size_x);
for l = 1 : size_x(3)
    %!!!NOTE: compensate rotating filter in conv2()
    Au(:,:,l) = conv2( xpad(:,:,l), rot90(u,2), 'valid' );
end

return;
