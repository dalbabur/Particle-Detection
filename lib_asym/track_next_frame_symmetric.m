function [findFlag, data, trackingFlag] = track_next_frame_symmetric(currentPt, rst, rstFrameNum, trackingFlag, marx, mary)
currentFrame = currentPt(6);
nextFrame = currentFrame+1;
idx = find(rstFrameNum == nextFrame);
if isempty(idx)
    findFlag = 0;
    data = [];
    return
end

nextj = 0;
minDist = inf;
for i = 1:size(rst{idx,1},1)
    if trackingFlag{idx}(i) == 0
        if (rst{idx,1}(i,1) - currentPt(1)>0) && (rst{idx,1}(i,1) - currentPt(1)< marx)...
                && (abs(rst{idx,1}(i,2) - currentPt(2)) < mary) && (minDist > rst{idx,1}(i,1)-currentPt(1))
            minDist = rst{idx,1}(i,1)-currentPt(1);
            nextj = i;
        end
    end
end
if nextj == 0
    findFlag = 0;
    data = [];
    return
end

findFlag = 1;
data = [rst{idx}(nextj,:) nextFrame];
trackingFlag{idx}(nextj) = 1;