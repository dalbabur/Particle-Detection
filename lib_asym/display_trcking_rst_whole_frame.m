function display_trcking_rst_whole_frame(data,bgImg,pathName,fileName,info,LinLUT, halfHeight, halfWidth)

rstImg = bgImg;
for i = 1:size(data,1)
    img = cineRead(pathName,fileName,(data(i,5)-info.startFrame)+1,info,LinLUT);
    forImg = img - bgImg;
    x = data(i,1); y = data(i,2);
    minH = y-halfHeight;
    if minH <1
        minH = 1;
    end
    maxH = y+halfHeight;
    if maxH > size(img,1);
        maxH = size(img,1);
    end
    minW = x-halfWidth;
    if minW <1
        minW = 1;
    end                
    maxW = x+halfWidth;
    if maxW > size(img,2);
        maxW = size(img,2);
    end
    rstImg(minH:maxH, minW:maxW) = rstImg(minH:maxH, minW:maxW)+forImg(minH:maxH, minW:maxW);
end
display_image(rstImg, 1000, 3000, 1000)
hold on
for i = 1:size(data,1)
    plot(data(i,1),data(i,2),'r*');
    plot([data(i,1); 30*cosd(data(i,3))+data(i,1)],[data(i,2); 30*sind(data(i,3))+data(i,2)],'r');
%     text(data(i,1),data(i,2),['(' num2str(round(data(i,1))) ',' num2str(round(data(i,2))) ',' num2str(round(data(i,3))) ',' num2str(data(i,4)) ',' num2str(data(i,5)) ')']);
end
plot(data(:,1),data(:,2));
hold off
%pause