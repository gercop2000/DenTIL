function [ features,featureNames ] = getDenTILFeatures(imageOrTileArea,lympCentroids,nonLympCentroids,lympAreas)
%GETDENTILFEATURES Extract features associated to density of
%   tumor-infiltrating lymphotyes in a tile.
% Inputs:
%	imageOrTileArea		:	Tile to be analyzed. If the tile is not 
%							available or a faster result is desired, 
%							the tile area can be used instead.
% 	lympCentroids		:	Centroids of the lymphocytes 	
% 	nonLympCentroids	:	Centroids of the non lymphocytes
% 	lympAreas			: 	Areas of the lymphocytes

if numel(image)==1
	tileArea=imageOrTileArea;
else
	tileArea=getTissueArea(imageOrTileArea);
end

[features,featureNames ] = computeDenTILFeatures(tileArea,lympCentroids,nonLympCentroids,lympAreas)

end

function [features,featureNames ] = computeDenTILFeatures(tileArea,lympCentroids,nonLympCentroids,lympAreas)

numLymp=length(lympCentroids);
totNuclei=numLymp+length(nonLympCentroids);

%% Regular-density-based measures
%A=tileSize^2;
%totLympArea=sum(nucleiFeatures(prediction==1,1));
totLympArea=sum(lympAreas);

densLymp=numLymp/tileArea;
densAreaLymp=totLympArea/tileArea;
ratioLymp=numLymp/totNuclei;

%% Grouping-based measures
%groupingFactor=sum(1./pdist(lympCent))/length(lympCent);
groupingVector=getSumNodeWeightsThreshold(lympCentroids,'euclidean',.005);
normVect = normalizeVector(groupingVector,0);

maxGr=max(groupingVector);
minGr=min(groupingVector);
avgGr=mean(groupingVector);
stdGr=std(groupingVector);
medGr=median(groupingVector);
modeGr=mode(groupingVector);

numHighlyGroupedLymp=length(normVect(normVect>.5));

%% Convex-hull-based measures
if length(lympCentroids)>2
    [~,areaConvHull]=convhull([lympCentroids;nonLympCentroids]);
    [convHullLymp,areaConvHullLymp]=convhull(lympCentroids);
    convHullLymp=lympCentroids(convHullLymp,:);
    if length(nonLympCentroids)>2
        [convHullNonLymp,~]=convhull(nonLympCentroids);
        convHullNonLymp=nonLympCentroids(convHullNonLymp,:);
        intersArea=getIntersectedArea(convHullLymp,convHullNonLymp);
    else
        intersArea=0;
    end
    densLympConvHull=numLymp/areaConvHull;
    ratioConvHulls=areaConvHullLymp/areaConvHull;
else
    densLympConvHull=0;
    ratioConvHulls=0;
    intersArea=0;
end

%% Density-Matrix-based measures assuming the tile is a square
M=getDensityMatrixCore(sqrt(tileArea),5,lympCentroids);
M(M==0)=[];

maxM=max(M);
minM=min(M);
avgM=mean(M);
stdM=std(M);
medM=median(M);
modeM=mode(M);

%% compiling features

features=[densLymp,densAreaLymp,ratioLymp,maxGr,minGr,avgGr,...
    stdGr,medGr,modeGr,numHighlyGroupedLymp,densLympConvHull,...
    ratioConvHulls,intersArea,maxM,minM,avgM,stdM,medM,modeM,...
    ];

featureNames={'#Lymp/TissueArea','LympTotalArea/TissueArea',...
    '#Lymp/#TotalNuclei','MaxLympGroupingFactor','MinLympGroupingFactor',...
    'AvgLympGroupingFactor','StdLympGroupingFactor',...
    'MedianLympGroupingFactor','ModeLympGroupingFactor','NumHighlyGroupedLymp',...
    '#Lymp/TotalConvHullArea','LympConvHullArea/TotalConvHullArea',...
    'IntersectedAreaConvHullLymp&NonLymp','MaxDensityMatrixVal',...
    'MinDensityMatrixVal','AvgDensityMatrixVal','StdDensityMatrixVal',...
    'MedianDensityMatrixVal','ModeDensityMatrixVal',...
    };


end