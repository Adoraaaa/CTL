function [yout]=recovercol2im(X,szX,K)
% ========================================================
% ========================================================
M=szX(1); N=szX(2); 
n=sqrt(size(X,1)); 
yout=zeros(M,N); 
Weight=zeros(M,N); 
i=1; j=1;
if mod((M-n),K)==0
    for k=1:1:size(X,2)
        patch=reshape(X(:,k),[n,n]);
        yout(i:i+n-1,j:j+n-1)=yout(i:i+n-1,j:j+n-1)+patch;
        Weight(i:i+n-1,j:j+n-1)=Weight(i:i+n-1,j:j+n-1)+1;
        if i<M-n+1
            i=i+K;
        else
            i=1; j=j+K;
        end
    end
else
    for k=1:1:size(X,2)
        patch=reshape(X(:,k),[n,n]);
        yout(i:i+n-1,j:j+n-1)=yout(i:i+n-1,j:j+n-1)+patch;
        Weight(i:i+n-1,j:j+n-1)=Weight(i:i+n-1,j:j+n-1)+1;
        if i<M-2*n+1
            i=i+K;
        elseif i<M-n+1
            i=M-n+1;
        else
            i=1; 
            if j<M-2*n+1
                j=j+K;
            else
                j=M-n+1;
            end
        end
    end    
end
%yout=yout./Weight;
return
