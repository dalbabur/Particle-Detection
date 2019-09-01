%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Droplet/Cell/Bead Detection
%
% Diego Alba 3/12/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h,window_length,skips] = consistency_test(cine_folder,cine_file,frames,info,LinLUT,...
    window_height,window_length,window_origin,skip,pxpf,prom)

h = zeros(1200,length(window_length),length(skip));
skip_frames = floor(window_length/pxpf);
skips = floor(skip_frames'*skip);

for s = 1:length(skip)
    for w = 1:length(window_length)
        
        data = cineRead2(cine_folder,cine_file,frames,info,LinLUT,...
            window_height,window_length(w),window_origin,floor(skip_frames(w)*skip(s)));
        
        l = [];
        for p = 1:length(data(1,1,:))
            locss = [];
            for k = 1:length(data(:,1,1))
                [~,locs,~] = findpeaks(-(data(k,:,p)),'MinPeakProminence',prom,'MinPeakDistance',4);
                locss = [locss locs];
            end
            locss = sort(locss)+length(data(1,:,1))*(p-1);
            l = [l locss];
        end
        
        clocs = circshift(l,1);
        d = l(2:end)-clocs(2:end);
        
        dummy = (histcounts(d,'BinMethod','integers','Normalization','probability'));
        h(1:length(dummy),w,s) = dummy;
    end
end
h = squeeze(h);
% figure; 
% for i = 1:length(skip)
%     subplot(1,length(skip),i); 
%     imagesc(h(:,:,i))
%     xticks(1:length(window_length))
%     xticklabels(skips(:,i))
%     xlabel('# Skipped Frames')
%     
%     colorbar
%     a1Pos = get(gca,'Position');
%     axes('Position',[a1Pos(1) a1Pos(2)-.05 a1Pos(3) a1Pos(4)],'Color','none',...
%         'XTick',1:length(window_length),'XTickLabel',window_length','YTick',[],'YTickLabel',[])
%     xlabel('Window Length')
%     xlim([0.5 length(window_length)+0.5])
% end
end

