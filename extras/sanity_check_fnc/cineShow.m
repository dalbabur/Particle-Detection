%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Droplet/Cell/Bead Detection
%
% Diego Alba 3/12/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cineShow(cine_folder,cine_file,f,info,LinLUT,...
        window_length,window_origin,window_height,prom)
    
    data = cineRead2(cine_folder,cine_file,f,info,LinLUT,...
        1:info.Height,window_length,window_origin);
    
    figure
    subplot(2,1,1)
    colormap('gray')
    imagesc(data);
    title(f)
    hold on 
    c={'b','r'};
    for i = 1:length(window_height)
        [~,l] = findpeaks(-(data(window_height(i),:)),'MinPeakProminence',prom,'MinPeakDistance',4,'MaxPeakWidth',30);
        if ~isempty(l), plot(l,window_height(i),[c{i},'*']), end
    end
    subplot(2,1,2)
    hold on
    for i = 1:length(window_height)
        findpeaks(-(data(window_height(i),:)),'MinPeakProminence',prom,'MinPeakDistance',4,'MaxPeakWidth',30);
    end
    
end
