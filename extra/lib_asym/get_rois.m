function [outImg, oroi, iroi] = get_rois(inImg, roi, mar)
% ror is minX, minY, maxX, maxY
% oroi is boundary outer region of interest, reference is image
% iroi is boundary iner region of interest, reference is image
[height, width] = size(inImg);
minx = roi(1);
miny = roi(2);
maxx = roi(1)+roi(3)-1;
maxy = roi(2)+roi(4)-1;
iroi = [minx, miny, maxx, maxy];

nMiny= miny-mar;
if nMiny <1
    nMiny = 1;
end
nMaxy = maxy+mar;
if nMaxy > height
    nMaxy = height;
end
nMinx = minx-mar;
if nMinx <1
    nMinx = 1;
end                
nMaxx = maxx+mar;
if nMaxx > width
    nMaxx = width;
end
oroi = [nMinx, nMiny, nMaxx, nMaxy];

outImg = inImg(nMiny:nMaxy, nMinx:nMaxx);