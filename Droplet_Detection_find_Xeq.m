%%
% Droplet Detection
%
% Diego Alba 3/12/2019
%%%%

cine_folder = 'A:\DACS Cell Sorting with Prof. Sangwon Kim\Images From Cell Sorting KO\Alignment videos\2019.07.31';
cine_file = 'Re=8 vid 1 part1.cine';


window_height = [1:64]; % vector of pixels at which to read data
window_length = 1280; % number of pixles 
window_origin = 0; % offset 
frames = [1  183]; % range of frames
cells = 1; % cells or beads
radius = 4;

n_locs = 15;
offset = 100;
x_locs = round(linspace(window_origin+offset, window_length-offset, n_locs));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
importfile('LinLUT.mat'); %a conversion between packed 10bit data to real 10bit data
info = cineInfo(cine_folder,cine_file);



%%
% Get sample frames 
tic
sample = cineRead2(cine_folder,cine_file,frames(1):frames(2),info, LinLUT, window_height, window_length, window_origin); 

if cells, sample = (sample + min(-sample(:))); else sample = -sample; end
toc
%% find wall
tic
all_pr = [];
minpwidth = 3;
for p = frames(1):frames(2)
    for k = 1:n_locs
        [~,~,~,pr] = findpeaks(-mean(sample(:,x_locs(k)-radius*2:x_locs(k)+radius*2,p),2),'MinpeakWidth',minpwidth); % only care about prominance first
        all_pr = [all_pr; pr];
    end
end
% plot prominance
figure
histogram(all_pr)
toc


%% set peak 
tic
prom = 1000;
pdist = 4;
toc
%%
all_locs = NaN(frames(2),n_locs,2);
for p = frames(1):frames(2)
    for k = 1:n_locs
        [~,locs] = findpeaks(-mean(sample(:,x_locs(k)-radius:x_locs(k)+radius,p),2),'MinPeakProminence',prom,'MinPeakDistance',pdist,'MinPeakWidth',minpwidth,'Annotate','extent');
        all_locs(p,k,1:length(locs)) = locs;
    end
end

wall = squeeze(round((nanmean(all_locs))));
toc
%%

% use findpeaks on data
all_pr = [];
all_w = [];
maxpwidth = radius*2;
for p = frames(1):frames(2)
    for k = 1:n_locs
        [~,~,w,pr] = findpeaks(mean(sample(wall(k,1):wall(k,2),x_locs(k)-radius:x_locs(k)+radius,p),2)); % only care about prominance first
        all_pr = [all_pr; pr];
        all_w = [all_w; w];
    end
end

% plot prominance
figure
histogram(all_pr)
figure
histogram(all_w)
toc


%% set peak features for rest of analysis
tic
prom = 900;
pdist = 4;
maxpwidth = 15;
toc

%%
% use findpeaks on data
tic
all_locs = zeros(frames(2),n_locs,5);
flag = zeros(frames(2),n_locs);
for p = frames(1):frames(2)
    for k = 1:n_locs
        [~,locs,w,pr] = findpeaks(mean(sample(wall(k,1):wall(k,2),x_locs(k)-radius:x_locs(k)+radius,p),2),'MinPeakProminence',prom,'MinPeakDistance',pdist,'MinPeakWidth',minpwidth,'MaxPeakWidth',maxpwidth,'Annotate','extent');
%         remove = zeros(1,length(w));
%         for i = 1:length(w)
%             if w(i) > 0.7*mean(wall(:,2)-wall(:,1))
%                 remove(i) = 1;
%             end
%         end
%         locs(logical(remove)) = [];
        all_locs(p,k,1:length(locs)) = locs;
        if ~isempty(locs)
            flag(p,k) = 1;
        end
    end
end
toc
[idx_p,idx_k] = find(flag);
for i = 1:length(idx_k)
    p = idx_p(i);
    k = idx_k(i);
    figure
    findpeaks(mean(sample(wall(k,1):wall(k,2),x_locs(k)-radius:x_locs(k)+radius,p),2),'MinPeakProminence',prom,'MinPeakDistance',pdist,'MinPeakWidth',minpwidth,'MaxPeakWidth',maxpwidth,'Annotate','extent');
end

if sum(sum((all_locs(:,:,1)))) == sum((all_locs(:)))
    disp('detected a maxiimum of 1 particle per region, all good')
else
    disp('detected more than 1 particle per region, not supported')
end
for k = 1:n_locs
    all_locs(:,k,:) = all_locs(:,k,:)+wall(k,1);
end
figure
subplot(1,2,1)
imagesc(all_locs(:,:,1))
subplot(1,2,2)
histogram(all_locs(all_locs>0))

%%
% calculate distance from the wall
dist = all_locs;
all_locs(all_locs==0) = NaN;
for k = 1:n_locs
    dist(:,k,:) = min(all_locs(:,k,:)-wall(k,1) ,wall(k,2) - all_locs(:,k,:));
end

figure
subplot(1,2,1)
imagesc(dist(:,:,1))
subplot(1,2,2)
histogram(dist)

