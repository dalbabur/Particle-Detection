clear all
close all
clc

%Open 11 cores
% s=matlabpool ('size');
% if s==0
%     matlabpool local 11
% else
%     matlabpool close
%     matlabpool local 11
% end
% s=gcp('nocreate');
% if isempty(s)
%     poolsize=0;
% else
%     poolsize=s.NumWorkers;
% end
% if poolsize==0
%     parpool('local',11)
% else
%     delete(s);
%     parpool('local',11)
% end

%% read data and set library
addpath lib_asym
importfile('mat/LinLUT.mat');


[fileName, pathName, ~] = uigetfile( ...
    {'*.cine','CINE-files (*.cine)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Select a movie file', ...
    'MultiSelect', 'on');

info=cineInfo(pathName,fileName); % changed a little bit
height = info.Height;
width = info.Width;
numFrames = info.NumFrames;

%Make a data directory
ImgDirName=[pathName '\' substr(fileName,0,-5),'\images'];
DataDirName=[pathName '\' substr(fileName,0,-5),'\data'];
mkdir(ImgDirName);
mkdir(DataDirName);


%% set variables
% for display
minVal =  1000;
maxVal = 3000;
marginLen = 10;
% for constructing background at least 20 frames
% ask user to define the range of backgroud images
prompt  = {'Starting frame # for BG (at least 20 frame clearance'};
def     = {'-80000'};
dlgtitle   = 'Frame range for the BG determination';
lines   = 1;
% options.Resize='on';
% options.WindowStyle='modal';
% options.Interpreter='tex';
answer  = inputdlg(prompt,dlgtitle,lines,def);

bgStartFrame = eval(answer{1})-info.startFrame+1;
bgEndFrame= bgStartFrame+20;

% bgStartFrame = 95; %80th frame (Current frame-initial frame)
% bgEndFrame = 115; %100th frame (Current frame-initial frame)
numTrainFrames = bgEndFrame-bgStartFrame+1;
% for finding frames containing objects
isObjThre = 500;
filterSize = 3;
% for removing small objects
dilationSize = 7;
smallObjSize = 50;
% for detection

% ask user to define the range of backgroud images
prompt  = {' numSamples','numIteration'};
def     = {'2000', '20'};
dlgtitle   = 'Please define interation parameters (increment of 1000 &10)';
lines   = 1;
% options.Resize='on';
% options.WindowStyle='modal';
% options.Interpreter='tex';
answer  = inputdlg(prompt,dlgtitle,lines,def);
numSamples = eval(answer{1}); % important parameter: increase up to 2000 -> it will consume 2 times of computational time (default:1000)
numIteration = eval(answer{2}); % important parameter: increase up to 20   -> it will consume 2 times of computational time (default:10)
% for detection
marx = 80;
mary = 30;

%% annotate region of interest
img=cineRead(pathName,fileName,1,info,LinLUT);
display_image(img, minVal,maxVal, 1); 
title('annotate region of interest')
roi = floor(getrect);                               % roi -> [xmin ymin width height]
[nImg,oroi, iroi] = get_rois(img, roi, marginLen);  % oroi, iroi ->[xmin ymin xmax ymax];

%% find background image (median image or you can use mean image)
tmp1 = uint16(zeros(height*width, numTrainFrames));
for i = 1 : numTrainFrames
    tmp2 = cineRead(pathName,fileName,bgStartFrame+i-1,info,LinLUT);
    tmp1(:,i) = tmp2(:);
end
bgImg = median(tmp1,2);
bgImg = reshape(bgImg, height, width);
display_image(bgImg, minVal, maxVal,1);
clear tmp1;clear tmp2;

%% annotate template
% find frames containing objects
flagObjFrame = false(numFrames,1);
bgImgRoi = bgImg(iroi(2): iroi(4), iroi(1): iroi(3));
bgImgORoi = bgImg(oroi(2): oroi(4), oroi(1): oroi(3));

%slice the variable for parfor
iroi2=iroi(2); iroi4=iroi(4);iroi1= iroi(1);iroi3= iroi(3);
oroi2=oroi(2); oroi4=oroi(4);oroi1=oroi(1);oroi3=oroi(3);
tic
parfor i = 1:numFrames
% for i = 1:numFrames
    i
    img = cineRead(pathName,fileName,i,info,LinLUT);
    img = img(iroi2:iroi4, iroi1:iroi3);
    forImg = img-bgImgRoi;
    if max(forImg(:)) > isObjThre
        flagObjFrame(i) = 1;
    end
%     display_image(forImg, 0, 500,2);
%     pause
end
toc
% % get rid of outliers (sporadic false detection) include inliers (sporadic missing detection)
% filter1 = ones(filterSize);
% for i = 1:numFrames
%     
% end

% annotate object
for i = 1:numFrames
    if flagObjFrame(i)
        img = cineRead(pathName,fileName,i,info,LinLUT);
        display_image(img, minVal, maxVal,2)
        button = questdlg('Can you see the object');
        if strcmp(button, 'Yes')
            roo = floor(getrect);
            [~, oroo, iroo] = get_rois(img, roo, marginLen);
            objImg = img(iroo(2):iroo(4), iroo(1):iroo(3));
            objbgImg = bgImg(iroo(2):iroo(4), iroo(1):iroo(3));
            break;
        end
    end
end
[objPlusEdge_t, objMinusEdge_t, objEdge_s] = find_edge(objImg, objbgImg);% objEdge_s-> saves edge
% change edge map to index
tmp1 = objPlusEdge_t(:);
tmp2 = 1:length(tmp1);
tmp3 = tmp2(tmp1);
[pey, pex] = ind2sub(size(objPlusEdge_t), tmp3');

tmp1 = objMinusEdge_t(:);
tmp2 = 1:length(tmp1);
tmp3 = tmp2(tmp1);
[ney, nex] = ind2sub(size(objMinusEdge_t), tmp3');

tmp1 = objEdge_s(:);
tmp2 = 1:length(tmp1);
tmp3 = tmp2(tmp1);
[ey, ex] = ind2sub(size(objEdge_s), tmp3');

% % annotate center of the object
% display_image(objImg, minVal, maxVal,3);
% title('point center of mass')
% [cx, cy] = getpts;
% cx = floor(cx); cy = floor(cy);

%% find center of mass and center of circumference
while(1)
    clearvars -global modeFlag;
    clearvars -global startPoints;
    clearvars -global endPoints;
    clearvars -global rectanglePoints;

    global modeFlag;
    global startPoints;
    global endPoints;
    global rectanglePoints;

    close all
    in = double(objImg);
    in = in-double(minVal);
    in = in/double(maxVal-minVal);
    figure('WindowButtonDownFcn',@buttonDown, 'WindowButtonUpFcn',@buttonUp);
    imshow(in)

    modeFlag = 0;
    disp('press any button when you are done')
    pause

    % calculate center of mass
    area = zeros(size(rectanglePoints,1),1);
    areaCenter = zeros(size(rectanglePoints,1),2);
    for i = 1:size(rectanglePoints,1)
        area(i) = sqrt((rectanglePoints(i,1)-rectanglePoints(i,3))^2+ (rectanglePoints(i,2)-rectanglePoints(i,4))^2)*...
            sqrt((rectanglePoints(i,3)-rectanglePoints(i,5))^2+ (rectanglePoints(i,4)-rectanglePoints(i,6))^2);
        areaCenter(i,:) = [(rectanglePoints(i,1)+rectanglePoints(i,5))/2, (rectanglePoints(i,2)+rectanglePoints(i,6))/2];
    end
    centerMass = floor([sum(area.*areaCenter(:,1))/sum(area) sum(area.*areaCenter(:,2))/sum(area)]);


    % calulate center of circumference
    rectPt = zeros(size(rectanglePoints,1)*4,2);
    for i = 1:size(rectanglePoints,1)
        rectPt(4*(i-1)+1,:) = [rectanglePoints(i,1) rectanglePoints(i,2)];
        rectPt(4*(i-1)+2,:) = [rectanglePoints(i,3) rectanglePoints(i,4)];
        rectPt(4*(i-1)+3,:) = [rectanglePoints(i,5) rectanglePoints(i,6)];
        rectPt(4*(i-1)+4,:) = [rectanglePoints(i,7) rectanglePoints(i,8)];
    end

    centerCur = zeros(1,2);
    minDist = 100000;
    for i = 1:size(objImg,1)
        for j = 1:size(objImg,2)
            maxDist = max(sqrt((rectPt(:,1)-j).^2 + (rectPt(:,2)-i).^2));
            if maxDist < minDist
                centerCur = [j,i];
                minDist = maxDist;
            end
        end
    end

    % display result
    close all
    display_image(objImg, minVal, maxVal,3);
    hold on
    plot(centerMass(1), centerMass(2), 'g*');
    plot(centerCur(1), centerCur(2), 'r.');
    for i = 1:size(rectanglePoints,1)
        plot([rectanglePoints(i,1), rectanglePoints(i,3), rectanglePoints(i,5), rectanglePoints(i,7), rectanglePoints(i,1)],...
            [rectanglePoints(i,2), rectanglePoints(i,4), rectanglePoints(i,6), rectanglePoints(i,8), rectanglePoints(i,2)]);
    end
    pause


    button = questdlg('Do you like the annotation');
    if strcmp(button, 'Yes')
        break;
    end
    
    clearvars -global modeFlag;
    clearvars -global startPoints;
    clearvars -global endPoints;
    clearvars -global rectanglePoints;
end



% annotate reference angle (direction where angle is 0)
cx = centerMass(1); cy = centerMass(2);
rectPt = floor(rectPt);

close all
display_image(objImg, minVal, maxVal,3);
hold on
plot(cx, cy, 'g*');
hold off
[directionPTx, directionPTy] = getpts; %User should click a point and hit enter to select the point
directionPTx = floor(directionPTx);
directionPTy = floor(directionPTy);
    
flipCenterCur = centerCur;
flipCenterCur(2) = size(objImg,1) - flipCenterCur(2);

flipCenterMass = centerMass;
flipCenterMass(2) = size(objImg,1) - flipCenterMass(2);

flipDirectionPTx = directionPTx;
flipDirectionPTy = size(objImg,1) - directionPTy;


%% detection
% coarse detection: gabor filter -> not good for manmade objects
% fine pose estimation: distance transformation, randomized method -> later
% if we have time more.
% local features -> may be proper for this application

tic
% for frames where there exist objects
rst = cell(numFrames,1);
roo4=roo(4);roo3=roo(3);roo2=roo(2);roo1=roo(1);
parfor i = 1:numFrames
% for i =1:numFrames
    i
    if flagObjFrame(i)
        img = cineRead(pathName,fileName,i,info,LinLUT);
        %img = img(oroi(2):oroi(4),oroi(1):oroi(3)); %command line for
        %regular for-loop
        img = img(oroi2:oroi4,oroi1:oroi3);
        %display_image(img, minVal, maxVal,2)
        
        % find region where there exist object
        forImg = (img-bgImgORoi);
        binImg = reshape((forImg(:)>isObjThre), size(forImg));
        
        % connect close regions
        binImg = bwmorph(binImg, 'dilate',dilationSize);
        
        % get rid of small object region
        CC = bwconncomp(binImg);
        flagBig = true(CC.NumObjects,1);
        for j = 1:CC.NumObjects
            if length(CC.PixelIdxList{j}) < smallObjSize
               flagBig(j) = 0;
            end
        end
        
        if sum(flagBig) == 0
            flagObjFrame(i) = 0;
        else
            rst{i,1} = zeros(sum(flagBig),5); % pose, flip, cost
        end
        % for all big object regions 
        cnt2 = 1;
        for j = 1:CC.NumObjects
            if flagBig(j)
                %find feature and describe it
                [h,w] = ind2sub(size(img),CC.PixelIdxList{j});
                minH = min(h)-floor(roo4/2);
                if minH <1
                    minH = 1;
                end
                maxH = max(h)+floor(roo4/2);
                if maxH > size(img,1);
                    maxH = size(img,1);
                end
                minW = min(w)-floor(roo3/2);
                if minW <1
                    minW = 1;
                end                
                maxW = max(w)+floor(roo3/2);
                if maxW > size(img,2);
                    maxW = size(img,2);
                end
                sceneImg = img(minH:maxH,minW:maxW);
                scenebgImg = bgImgORoi(minH:maxH,minW:maxW);
                [edge1, edge2, edge3] = find_edge(sceneImg, scenebgImg);
                [dist1, dist2, dist3] = distance_transform(edge1, edge2, edge3);
                
                % detect object with particle filter frame
                [pose, cost] = particle_filter_asymmetric(dist1, dist2, dist3, [pex, pey], [nex, ney], [ex, ey], [cx, cy], numSamples, numIteration, objImg);
%                 display_image(sceneImg, minVal, maxVal,2);
%                 hold on
%                 plot(pose(1), pose(2), '.')
                % change reference and save it

                pose(1) = pose(1)+minW-1;
                pose(2) = pose(2)+minH-1;

                rst{i,1}(cnt2,:) = [pose cost];
                cnt2 = cnt2+1;                
            end
        end
    end
end
toc
numObjFrames = sum(flagObjFrame);
tmpRst = cell(numObjFrames,1);
cnt = 1;
for i = 1:numFrames
    if flagObjFrame(i)
        tmpRst{cnt,1} = rst{i,1};
        cnt = cnt+1;
    end
end
rst = tmpRst;
tmp = 1:numFrames;
rstFrameNum = tmp(flagObjFrame);

%% get rid of objects closely detected
for i = 1:numObjFrames
    if size(rst{i,1},1) > 1
        tempLen = size(rst{i,1},1);
        flag = true(tempLen,1);
        
        % for each point
        for j = 1:tempLen-1
            if flag(j)
                ref = rst{i,1}(j,1:2);
                target = rst{i,1}(j+1:end,1:2);
                % find close points
                idx = sqrt((target(:,1)-ref(1)).^2+(target(:,2)-ref(2)).^2) < size(objImg,1);
                % if there are close points
                if sum(idx) ~=0
                    % find a point whose cost is minimum
                    minCost = rst{i,1}(j,5);
                    minIdx = j;
                    for k=j+1:tempLen
                        if minCost > rst{i,1}(k,5) && idx(k-j)
                            minCost = rst{i,1}(k,5);
                            minIdx = k;
                        end
                    end
                    
                    % get rid of the other points
                    idx = [1;idx];
                    for k = j:tempLen
                        if k ~=minIdx && idx(k-j+1)
                            flag(j) = 0;
                        end
                    end
                    
                end
                
            end
        end
        
        rst{i,1} = rst{i,1}(flag,:);
    end
end
% re order the detections in the same frame (small x goes first)
for i = 1:numObjFrames
    if size(rst{i,1},1) > 1
        tmp = rst{i,1};
        [~, idx] = sort(tmp(:,1));
        rst{i,1} = tmp(idx,:);
    end
end
    
% %% display results
% cnt = 1;
% for i = 1:numFrames
%     if flagObjFrame(i)
%         if ~isempty(rst{cnt,1})
%             img = cineRead(pathName,fileName,i,info,LinLUT);
%             img = img(oroi(2):oroi(4),oroi(1):oroi(3));
%             display_image(img, minVal, maxVal,2)
%             title(['frame number is ' num2str(i)])
%             hold on
%             for j = 1:size(rst{cnt,1},1)
%                 r = rst{cnt,1}(j,3);
%                 plot(rst{cnt,1}(j,1),rst{cnt,1}(j,2),'r*');
%                 plot([rst{cnt,1}(j,1); 30*cosd(r)+rst{cnt,1}(j,1)],[rst{cnt,1}(j,2); 30*sind(r)+rst{cnt,1}(j,2)],'r');
%                 text(rst{cnt,1}(j,1),rst{cnt,1}(j,2),['(' num2str(rst{cnt,1}(j,1)) ',' num2str(rst{cnt,1}(j,2)) ',' num2str(rst{cnt,1}(j,3)) ')' ' cost :' num2str(rst{cnt,1}(j,4))]);
%             end
%             hold off
%             pause
%         end
%         cnt = cnt+1;
%     end
% end
%% tracking
trackingFlag = cell(numObjFrames,1);
for i = 1:numObjFrames
    trackingFlag{i} = false(size(rst{i,1},1),1);
end

trackingRst = cell(numObjFrames,1);
cnt = 1;
for i = 1:numObjFrames
    for j = 1:size(rst{i,1},1)
        
        if trackingFlag{i}(j) == 0
            tmpRst = [rst{i,1}(j,:) rstFrameNum(i)]; %x, y, rot, flip, cost, frameNum
            trackingFlag{i}(j) = 1;
            while (1)
                [findFlag, data, trackingFlag] = track_next_frame_symmetric(tmpRst(end,:),rst, rstFrameNum, trackingFlag, marx, mary);
                if findFlag == 0
                    break;
                end
                tmpRst = [tmpRst; data];
            end
            trackingRst{cnt} = tmpRst;
            cnt = cnt+1;
        end
        
    end
end
%% calculate absolute position of the points
finalTrackingRst = cell(cnt-1,1);
for i = 1:cnt-1
    finalTrackingRst{i} = zeros(size(trackingRst{i},1), 6+2+2*size(rectPt,1)+1); %cMassX, cMassY, rot, rotRef, cost, frameNum, cCircumX, cCricumY, otherPoints, flip
    for j = 1:size(trackingRst{i},1)
        if trackingRst{i}(j,4) == 0
            trackingRst{i}(j,6) = trackingRst{i}(j,6)+info.startFrame-1;
            finalTrackingRst{i}(j,1:3) = trackingRst{i}(j,1:3); % save cMassX, cMassY, rot
            finalTrackingRst{i}(j,5:6) = trackingRst{i}(j,5:6); % save cost, frameNum
            finalTrackingRst{i}(j,end) = trackingRst{i}(j,4);   % flip or not
            
            % calculate absolute position of cCircumX and Y
            rot = trackingRst{i}(j,3);
            x = trackingRst{i}(j,1); y = trackingRst{i}(j,2);
            rMat = [cosd(rot) -sind(rot);
            sind(rot) cosd(rot)];
            tmpCenterCur = rMat*(centerCur-centerMass)';
            finalTrackingRst{i}(j,7:8) = tmpCenterCur'+[x y];

            % calculate rotRef
            dirPT = rMat*([directionPTx directionPTy]-centerMass)';
            finalTrackingRst{i}(j,4) = atan2(dirPT(2), dirPT(1))/pi*180;

            % calcuate absolute position of other points
            for k = 1:size(rectPt,1)
                tmpRectPt = rMat*(rectPt(k,:)-centerMass)';
                finalTrackingRst{i}(j,(k-1)*2+9:(k-1)*2+10) = tmpRectPt'+[x y];
            end
        else
            trackingRst{i}(j,6) = trackingRst{i}(j,6)+info.startFrame-1;
            finalTrackingRst{i}(j,1:3) = trackingRst{i}(j,1:3); % save cMassX, cMassY, rot
            finalTrackingRst{i}(j,5:6) = trackingRst{i}(j,5:6); % save cost, frameNum
            finalTrackingRst{i}(j,end) = trackingRst{i}(j,4);   % flip or not
            
            % calculate absolute position of cCircumX and Y
            rot = trackingRst{i}(j,3);
            x = trackingRst{i}(j,1); y = trackingRst{i}(j,2);
            rMat = [cosd(rot) -sind(rot);
            sind(rot) cosd(rot)];
            tmpCenterCur = rMat*(flipCenterCur-flipCenterMass)';
            finalTrackingRst{i}(j,7:8) = tmpCenterCur'+[x y];

            % calculate rotRef
            dirPT = rMat*([flipDirectionPTx flipDirectionPTy]-flipCenterMass)';
            finalTrackingRst{i}(j,4) = atan2(dirPT(2), dirPT(1))/pi*180;

            % calcuate absolute position of other points
            for k = 1:size(rectPt,1)
                tmpRectPt = rMat*([rectPt(k,1) size(objImg,1)-rectPt(k,2)]-flipCenterMass)';
                finalTrackingRst{i}(j,(k-1)*2+9:(k-1)*2+10) = tmpRectPt'+[x y];
            end
        end
        
    end
end

% update the positions
offsetX = oroi(1)-1;
offsetY = oroi(2)-1;
for i = 1:length(trackingRst)
    for j = 1:size(trackingRst{i},1)
        trackingRst{i}(j,1) = trackingRst{i}(j,1) + offsetX;
        trackingRst{i}(j,2) = trackingRst{i}(j,2) + offsetY;
        finalTrackingRst{i}(j,1) = finalTrackingRst{i}(j,1) + offsetX;
        finalTrackingRst{i}(j,2) = finalTrackingRst{i}(j,2) + offsetY;
        finalTrackingRst{i}(j,7) = finalTrackingRst{i}(j,7) + offsetX;
        finalTrackingRst{i}(j,8) = finalTrackingRst{i}(j,8) + offsetY;
        
        for k = 1:size(rectPt,1)
            finalTrackingRst{i}(j,(k-1)*2+9) = finalTrackingRst{i}(j,(k-1)*2+9)+offsetX;
            finalTrackingRst{i}(j,(k-1)*2+10) = finalTrackingRst{i}(j,(k-1)*2+10)+offsetY;
        end
    end
end


%% Kilho's display
%for i = 1:cnt-1
    %display_trcking_rst(trackingRst{i},bgImgORoi,oroi,pathName,fileName,info,LinLUT, ceil(size(objImg,1)/2), ceil(size(objImg,2)/2));
    %***display_trcking_rst_whole_frame_asymmetric(finalTrackingRst{i},bgImg,pathName,fileName,info,LinLUT, ceil(size(objImg,1)/2), ceil(size(objImg,2)/2));
    %display_trcking_rst_whole_frame_1(finalTrackingRst{i},bgImg,pathName,fileName,info,LinLUT, ceil(size(objImg,1)/2), ceil(size(objImg,2)/2));
    %i
    %pause
%end

%% Sj's display & save the overlayed image
% ask user to define the range of backgroud images
prompt  = {'# of consecutive images'};
def     = {'20'};
dlgtitle   = 'Threshold of consecutiveness';
lines   = 1;
% options.Resize='on';
% options.WindowStyle='modal';
% options.Interpreter='tex';
answer  = inputdlg(prompt,dlgtitle,lines,def);
CountThreshold=eval(answer{1});

CnT=1;
for i = 1:cnt-1
    if size(trackingRst{i},1)>CountThreshold
        %***display_trcking_rst_whole_frame(trackingRst{i},bgImg,pathName,fileName,info,LinLUT, ceil(size(objImg,1)/2), ceil(size(objImg,2)/2));
        %display_trcking_rst_whole_frame_asymmetric(trackingRst{i},bgImg,pathName,fileName,info,LinLUT, ceil(size(objImg,1)/2), ceil(size(objImg,2)/2));
        display_trcking_rst_whole_frame_asymmetric(finalTrackingRst{i},bgImg,pathName,fileName,info,LinLUT, ceil(size(objImg,1)/2), ceil(size(objImg,2)/2));
        if trackingRst{i}(1,6)<0
            CurrentFrameNumber=['_n',num2str(abs(trackingRst{i}(1,6)))];
        else
            CurrentFrameNumber=['_p',num2str(abs(trackingRst{i}(1,6)))];
        end
        %***saveas(gcf,[ImgDirName '\' substr(fileName,0,-5) CurrentFrameNumber],'jpg');
        saveas(gcf,[ImgDirName '\' substr(fileName,0,-5) CurrentFrameNumber '_FE2'],'jpg');
        subTrackingRst{CnT}=finalTrackingRst{i};
        CnT=CnT+1;
    end
end
%% save the result
%***save([DataDirName '\' fileName '.mat']);
%***save([DataDirName '\sub' substr(fileName,0,-5) '.mat']); %subpopulation that has consecutive detetion more than the threshold
% -----
% UPDATE BELOW SAVED FILE NAMES
particleType='AL';
save([DataDirName '\' fileName '_REV2.mat']);
save([DataDirName '\' particleType 'sub' substr(fileName,0,-5) '_REV2.mat']); %subpopulation that has consecutive detetion more than the threshold
%matlabpool close
poolobj = gcp('nocreate');
delete(poolobj);