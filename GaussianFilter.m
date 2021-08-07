function [UG] = GaussianFilter(yin,Uin,std)

[nuy,nux,nuz] = size(Uin);
UG = zeros(nuy,nux,nuz);

stdarg = num2cell(std);
[stdx,stdy,stdz] = deal(stdarg{:});
kernel = 2*ceil(3*[stdx stdz])+1;

for idx = 1:nuy
    Uin(idx,:,:) = imgaussfilt(Uin(idx,:,:),[stdx stdz],'FilterSize',kernel);
end

for idx = 1:nuy
    y = reshape(yin,[nuy,1]);
    gaussian = 1/(stdy*sqrt(2*pi))*exp(-(y-y(idx)).^2/2/stdy^2);
    UG(idx,:,:) = trapz(yin,Uin.*gaussian,1);
end
