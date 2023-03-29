function a = getTissueArea( I )
%GETTISSUEAREA Computes the image area corresponding to tissue

bw = rgb2gray(I)>220;
bw = bwareaopen(~bw, 1000);
bw = bwareaopen(~bw, 100);
a=sum(sum(~bw));

end