%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  findDistanseBetweenDictionaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [num,ratio,totalDistances] = I_findDistanseBetweenDictionaries(original,new)
% first, all the column in oiginal starts with positive values.
catchCounter = 0;
totalDistances = 0;
%for i = 1:size(new,2)
%    new(:,i) = sign(new(1,i))*new(:,i);
%end
% Normalize
for i = 1:size(new,2)
    new(:,i)=new(:,i)/norm(new(:,i));
    original(:,i)=original(:,i)/norm(original(:,i));
end

for i = 1:size(original,2)
    d = original(:,i);%sign(original(1,i))*original(:,i);
    on=ones(1,size(original,2));
    distances =on-abs(sum ( (new.*repmat(d,1,size(new,2)))));%sum ( (new-repmat(d,1,size(new,2))).^2);
    [minValue,index] = min(distances);
    %errorOfElement = norm(new(:,index)-d,2);
    errorOfElement = 1-abs(new(:,index)'*d);%ColumnNormalised
    %errorOfElement = 1/size(original,2)-abs(new(:,index)'*d);%ColumnNormalisedtosqrt(M
    totalDistances = totalDistances+errorOfElement;
    catchCounter = catchCounter+(errorOfElement<0.01);
end
%------------????dicdist of Alain---------
% for i = 1:size(original,2)
%     original(:,i) = sign(original(1,i))*original(:,i);
% end
% for i = 1:size(original,2)
%     errorOfElement = 1-abs(new(:,i)'*original);
%     [errorOfElement,ind(i)] = min( errorOfElement);
%     totalDistances = totalDistances+errorOfElement;
%     catchCounter = catchCounter+(errorOfElement<0.01);
% end

ratio = 100*catchCounter/size(original,2);
num=catchCounter;