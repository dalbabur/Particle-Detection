function buttonUp(src,eventdata)
global modeFlag;
global startPoints;
global endPoints;
global rectanglePoints;

modeFlag

if modeFlag == 0
    C = get (gca, 'CurrentPoint');
    hold on
    plot(C(1,1), C(1,2), '*');
    hold off
    endPoints = [endPoints; [C(1,1), C(1,2)]]
    hold on
    plot([startPoints(end,1);endPoints(end,1)], [startPoints(end,2);endPoints(end,2)]);
    hold off
    
elseif modeFlag == 1;
    % calculate 4 points of the rectangles
    C = get (gca, 'CurrentPoint');
    x = C(1,1); y = C(1,2);
    x1 = startPoints(end,1);
    y1 = startPoints(end,2);
    x2 = endPoints(end,1);
    y2 = endPoints(end,2);
    if (y2-y1) == 0
        points1 = [x ,y1];
        points2 = [x ,y2];
    elseif (x2-x1) ==0
        points1 = [x1, y];
        points2 = [x2, y];
    else
        m = (y2-y1)/(x2-x1);
        mNor = -1/m;
        A = [mNor, -1;
            m, -1];
        B = [mNor*x1-y1;
            m*x-y];
        points1 = inv(A)*B;
        
        A = [mNor, -1;
            m, -1];
        B = [mNor*x2-y2;
            m*x-y];
        points2 = inv(A)*B;
    end
    rectanglePoints = [rectanglePoints;
    x1, y1, x2, y2, points2', points1'];
    % draw the rectangles
    hold on
    plot([rectanglePoints(end,1), rectanglePoints(end,3), rectanglePoints(end,5), rectanglePoints(end,7), rectanglePoints(end,1)],...
        [rectanglePoints(end,2), rectanglePoints(end,4), rectanglePoints(end,6), rectanglePoints(end,8), rectanglePoints(end,2)]);
    hold off
    
end

if modeFlag == 0
    modeFlag = 1;
else
    modeFlag = 0;
end
