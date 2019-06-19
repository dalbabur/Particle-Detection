function buttonDown(src,eventdata)
global modeFlag;
global startPoints;

if modeFlag == 0;
    C = get (gca, 'CurrentPoint');
    hold on
    plot(C(1,1), C(1,2), '*');
    hold off
    startPoints = [startPoints; [C(1,1), C(1,2)]]
elseif modeFlag == 1;
    
end