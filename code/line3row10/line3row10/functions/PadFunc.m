function [M, bpad] = PadFunc(b, psf_radius,paddingchoice)
    
M = padarray(ones(size(b)), [psf_radius(1), psf_radius(2), 0], 0, 'both');    %mask
switch paddingchoice
    case 'circular'
        %%%circular padding
        bpad = padarray(b, [psf_radius(1), psf_radius(2), 0], 'circular', 'both');
    case 'reflective'
        %%%reflective padding
        bpad = padarray(b, [psf_radius(1), psf_radius(2), 0], 'symmetric', 'both');
    case 'zeros'
        %%%%zero padding
        bpad = padarray(b, [psf_radius(1), psf_radius(2), 0], 0, 'both');
    case 'ones'
        %%ones padding
        bpad = padarray(b, [psf_radius(1), psf_radius(2), 0], 1, 'both');
end
return;