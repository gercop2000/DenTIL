function M=getDensityMatrixCore(imgDim,partitions,centroids)

tileDim=imgDim/partitions;
M=[];
for i=1:tileDim:imgDim
    for j=1:tileDim:imgDim
        coords=centroids(centroids(:,1)>=i & centroids(:,1)<i+tileDim & ...
            centroids(:,2)>=j & centroids(:,2)<j+tileDim,:);
        M=[M;size(coords,1)];
    end
end

end