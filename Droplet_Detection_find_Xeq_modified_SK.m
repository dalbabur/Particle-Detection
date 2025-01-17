%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Droplet/Cell/Bead Detection
%
% Diego Alba 3/12/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cine_folder = 'E:\Yuelin Particle Focusing Videos and Analysis Code\KO Alignment Vids';
cine_file = 'Re=8 vid 1 part3.cine';
window_height = 1:64; % vector of pixels at which to read data
window_length = 1280; % number of pixles 
window_origin = 0; % offset 
cells = 1; % cells or beads
radius = 3; %distance around locations for finding particles
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
sample = cineRead2(cine_folder,cine_file,frames(1):frames(2),info,...
    LinLUT, window_height, window_length, window_origin);%Read the whole cine file
if cells, sample = (sample + min(-sample(:))); else sample = -sample; end
toc
%% find wall
tic
all_pr = [];
minpwidth = 0;
for p = frames(1):2:frames(2) %Read/detect the wall at every frame, could 
    % change only to 1 frame
    for k = 1:n_locs %Check to see how the wall varies with x location
        final_vals=sample(:,x_locs(k)-radius*2:x_locs(k)+radius*2,p);  %Segment the file
        final_vals=final_vals';
        final=mean(final_vals)';
        [~,~,~,pr] = findpeaks(final,'MinpeakWidth',minpwidth); % only care about prominance first
        all_pr = [all_pr; pr];
    end
end
% plot prominance
figure
histogram(all_pr)
toc
%This will determine the wall prominance. We expect this to vary from video
%to video, but have seen that using ~1000-1300 is a good starting point
%% set peak 
tic
prom = 1200; %Only Change this in this Section
pdist = 4; 
toc
%%
all_locs = NaN(frames(2),n_locs,2); %initialize the wall vector, should all
%get replaced
for p = 2 %Only determine wall location using 1 frame, can change the frame
    %number if not optimal
    for k = 1:n_locs %Go through each location, and determine where the 
        %wall is at that location, for distance from wall
        final_vals=sample(:,x_locs(k)-radius*2:x_locs(k)+radius*2,p);
        final_vals=final_vals';
        final=max(-final_vals)';
        [~,locs] = findpeaks(final,'MinPeakProminence',prom,'MinPeakDistance',pdist,'MinPeakWidth',minpwidth,'Annotate','extent');
        all_locs(p,k,1:length(locs)) = locs;
        %store locations of peaks/walls as all_locs
    end
end
%get rid of all the nans by using the squeeze functions, should be left
%with a 2x15 array
wall(:,:,1) = squeeze(round((nanmin(all_locs(:,:,1)))));
wall(:,:,2) = squeeze(round((nanmax(all_locs(:,:,2)))));
wall=squeeze(wall);
toc
%if you rerun this, you MUST clear wall
%% Find possible locations for the particle, by using find peaks, and Thresholding peaks
% use findpeaks on data
all_pr = [];
all_w = [];
maxpwidth = radius*2;
for p = frames(1):2:frames(2)
    for k = 1:n_locs
        final_vals=sample(wall(k,1):wall(k,2),x_locs(k)-radius*2:x_locs(k)+radius*2,p);
        %only searching for peaks between the pixel coordinates of the wall
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
%use this section as the peaks for the particles. For a bright particle, we
%usually see the peak prominance >1800, sometimes can go a little lower.
%However, there is a shadow from the wall which could show up if you go too
%low.
%% set peak features for rest of analysis
tic
prom = 1800;%Change this based on previous graph output
pdist = 4;
maxpwidth = 20;
toc
%%
% use findpeaks on data
tic
all_locs = zeros(frames(2),n_locs,5);
flag = zeros(frames(2),n_locs);
for p = frames(1):2:frames(2)
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
%This section should actually find the particle from your video/data. If
%you set the prominance value too low, it will find more than 1
%particle/location, which means you can rerun it after using a higher
%prominance value. 
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
%comment this part out when you're done trouble shooting, this displays the
%graph at every instance of detection
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

delete = input('enter and array of indices to delete (in the form: [1, 17,25]): where the numbers correspond to the detection from the graph on the left');
idx_p(delete) = []; %use this to find locations of particle for point plotting
idx_k(delete) = []; %use this to find locations of particle for point plotting

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

