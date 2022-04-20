TS = [112,160,3,32];
Data = zeros(TS);
for ii = 1 : TS(end)
    imageName = ['c',num2str(ii-1),'.png'];
    IMin0 = double(imread(imageName));
    Data(:,:,:,ii) = IMin0;
end
save VideoData.mat Data
