%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Droplet/Cell/Bead Detection
%
% Diego Alba 3/12/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cine_folder = 'C:\Users\Srivathsan Kalyan\OneDrive - Johns Hopkins University\WT Focusing Vids';
cine_file = 'Re=10_T1 part 2.cine';


window_height = 1:64; % vector of pixels at which to read data
window_length = 1280; % number of pixles 
window_origin = 0; % offset 
cells = 1; % cells or beads
radius = 3;
box = 100;
speed_part = 20;
n_locs = 15;
offset = 100;
x_locs = round(linspace(window_origin+offset, window_length-offset, n_locs));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
importfile('LinLUT.mat'); %a conversion between packed 10bit data to real 10bit data
info = cineInfo(cine_folder,cine_file);
frames = info.NumFrames;
frames=[1 frames];

%%
% Get sample frames 
tic
sample = cineRead2(cine_folder,cine_file,frames(1):frames(2),info, LinLUT, window_height, window_length, window_origin); 

if cells, sample = (sample + min(-sample(:))); else sample = -sample; end
toc
%% find wall
tic
all_pr = [];
minpwidth = 0;
for p = frames(1):frames(2)
    for k = 1:n_locs
        final_vals=sample(:,x_locs(k)-radius*2:x_locs(k)+radius*2,p);
        final_vals=final_vals';
        final=max(-final_vals)';
        [~,~,~,pr] = findpeaks(final,'MinpeakWidth',minpwidth); % only care about prominance first
        all_pr = [all_pr; pr];
    end
end
% plot prominance
figure
histogram(all_pr)
xlabel('Prominence')
ylabel('# of Occurances')
toc


%% set peak 
tic
prom = 1000; %Only Change this in this Section
pdist = 4; 
toc
%%
all_locs = NaN(frames(2),n_locs,2);
for p = 1:frames(2)
    for k = 1:n_locs
        final_vals=sample(:,x_locs(k)-radius*2:x_locs(k)+radius*2,p);
        final_vals=final_vals';
        final=min(-final_vals)';
        [~,locs] = findpeaks(final,'MinPeakProminence',prom,'MinPeakDistance',pdist,'MinPeakWidth',minpwidth,'Annotate','extent');
        all_locs(p,k,1:length(locs)) = locs;
    end
end

wall(:,:,1) = squeeze(round((nanmin(all_locs(:,:,1)))));
wall(:,:,2) = squeeze(round((nanmax(all_locs(:,:,2)))));
wall=squeeze(wall);
toc
%%

% use findpeaks on data
all_pr = [];
all_w = [];
maxpwidth = radius*2;
for p = frames(1):frames(2)
    for k = 1:n_locs
        final_vals=sample(wall(k,1):wall(k,2),x_locs(k)-radius*2:x_locs(k)+radius*2,p);
        final_vals=final_vals';
        final=max(final_vals)';
        [~,~,w,pr] = findpeaks(final); % only care about prominance first
        all_pr = [all_pr; pr];
        all_w = [all_w; w];
    end
end

% plot prominance
figure
histogram(all_pr)
hold on
title('Prominence of Each Peak')
figure
histogram(all_w)
hold on
title('Width of Each Peak')
toc


%% set peak features for rest of analysis
tic
prom = 1800;%Change this based on previous graph output
pdist = 4;
maxpwidth = 12;
toc

%%
% use findpeaks on data
tic
all_locs = zeros(frames(2),n_locs,5);
flag = zeros(frames(2),n_locs);
for p = frames(1):frames(2)
    for k = 1:n_locs
        final_vals=sample(wall(k,1):wall(k,2),x_locs(k)-radius*2:x_locs(k)+radius*2,p);
        final_vals=final_vals';
        final=max(final_vals)';
        [~,locs,w,pr] = findpeaks(final,'MinPeakProminence',prom,'MinPeakDistance',pdist,'MinPeakWidth',minpwidth,'MaxPeakWidth',maxpwidth,'Annotate','extent');
        all_locs(p,k,1:length(locs)) = locs;
        if ~isempty(locs)
            flag(p,k) = 1;
        end
    end
end
toc

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
xlabel('X-Locations')

subplot(1,2,2)
histogram(all_locs(all_locs>0))

[idx_p,idx_k] = find(flag);
for i = 1:length(idx_k)
    p = idx_p(i);
    k = idx_k(i);
    figure
    final_vals=sample(wall(k,1):wall(k,2),x_locs(k)-radius*2:x_locs(k)+radius*2,p);
    final_vals=final_vals';
    final=max(final_vals)';
    findpeaks(final,'MinPeakProminence',prom,'MinPeakDistance',pdist,'MinPeakWidth',minpwidth,'MaxPeakWidth',maxpwidth,'Annotate','extent');
    title(i)
end

delete = input('enter and array of indices to delete (in the form: [1, 17,25]): ');
idx_p(delete) = [];
idx_k(delete) = [];

particles = sub2ind(size(all_locs),idx_p,idx_k);

%%
% calculate distance from the wall
dist = all_locs(particles);
dist= min(dist-wall(idx_k,1),wall(idx_k,2)-dist);

dummy = zeros(frames(2),n_locs,5);
dummy(particles) = dist;

figure
subplot(1,2,1)
imagesc(dummy(:,:,1))
subplot(1,2,2)
histogram(dist)
all_locs(particles)
hold on
top_wall=mean(wall(:,1));
bottom_wall=mean(wall(:,2));
xlim([0;64])
top=[top_wall,0;top_wall,15];
bottom=[bottom_wall,0;bottom_wall,15];
plot(top(:,1),top(:,2))
plot(bottom(:,1),bottom(:,2))