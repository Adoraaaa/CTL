function imN = mycol2im_set_nonoverlap_mfilter(patches, sz, n)

for i=1:size(patches,1)
    cnt = 1;
    im = zeros(sz + 2*(n-1));
    for dx = 1:n
        for dy = 1:n
            cur_sz = sz + 2*(n-1) - [dy-1, dx-1];
            rec = col2imstep(patches{i,cnt}, cur_sz, [n,n], [n,n]);
            im(dy:end, dx:end) = im(dy:end, dx:end) + rec;
            cnt = cnt + 1;
        end
    end
    
    imN(:,:,i) = im(n:end-n+1, n:end-n+1);
end
return

