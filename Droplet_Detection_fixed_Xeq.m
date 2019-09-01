%%
% Droplet Detection
%
% Diego Alba 3/12/2019
%%%%

cine_folder = 'T:\Encapsulation\Cell-Encapsulation-Videos-DA\cells\wide channel';
cine_file = '24mil_cellspml-30ulpmin-20x.cine';


window_height = [23 38]; % vector of pixels at which to read data
window_length = 640; % number of pixles 
window_origin = 0; % offset 
frames = [1  576538]; % range of frames
cells = 1; % cells or beads
avg = {};
% from previous analysis
mean_vel = []; %18.57135909; % if known, else leave as [];
% std_vel = 9.003709145;
% n = 22618;
% prom = 150;
% pdist = 4;
% minpwidth = 0;
% maxpwidth = 7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
importfile('LinLUT.mat'); %a conversion between packed 10bit data to real 10bit data
info = cineInfo(cine_folder,cine_file);



%%
% Get sample frames to calculate particle velocity
% Also check peak prominance (important for particle detection)
if isempty(mean_vel)
tic
fr = floor(info.NumFrames);
skip_frames = floor(fr/1000);

sample = zeros(length(window_height),window_length,length(1:skip_frames:fr),2);

% load data
for i = 1:2 % if FPS too high, mean_vel is going to be wrong. Compare i*15 rather that i+1
    sample(:,:,:,i) = cineRead2(cine_folder,cine_file,[i:skip_frames:fr],info,LinLUT,...
        window_height,window_length,window_origin);
end  
if cells, sample = (sample + min(-sample(:))); else sample = -sample; end

if ~isempty(avg)
    sample2 = zeros(2,window_length,length(1:skip_frames:fr),2);
    for i = 1:2
        sample2(i,:,:,:) = mean(sample(avg{i},:,:,:));
    end
    clear sample
else
    sample2 = sample;
    clear sample
end


% use findpeaks on data
allpr = [];
for p = 1:length(sample2(1,1,:,1))
    for k = 1:length(sample2(:,1,1,1))
        [~,~,~,pr] = findpeaks((sample2(k,:,p,1))); % only care about prominance first
        allpr = [allpr pr];
    end
end

% plot prominance
figure
histogram(allpr)
toc
end
%% set peak features for rest of analysis
tic
prom = 450;
pdist = 4;
minpwidth = 0;
maxpwidth = 10;
toc

%% check peak features for rest of analysis

cine_frame = 400;
check = cineRead2(cine_folder,cine_file,cine_frame,info,LinLUT,...
        window_height,window_length,window_origin);
if cells, check = (check + min(-check(:))); else, check = -check; end  
figure
hold on
findpeaks((check(1,:)),'MinPeakProminence',prom,'MinPeakDistance',pdist,'MinPeakWidth',minpwidth,'MaxpeakWidth',maxpwidth','Annotate','extents')
findpeaks((check(2,:)),'MinPeakProminence',prom,'MinPeakDistance',pdist,'MinPeakWidth',minpwidth,'MaxpeakWidth',maxpwidth','Annotate','extents')

figure
check2 = cineRead2(cine_folder,cine_file,cine_frame,info,LinLUT,...
        1:info.Height,window_length,window_origin);
imagesc(check2)
hold on 
plot([0,640],[window_height(1),window_height(1)])
plot([0,640],[window_height(2),window_height(2)])

     
%%
% Calculate particle velocity 
if isempty(mean_vel)
tic
allt = NaN(80,length(sample2(1,1,:,1)),2);
e = 0;
for p = 1:length(sample2(1,1,:,1))
    xlocs = [];
    for k = 1:length(sample2(:,1,1,1))
        [~,t1] = findpeaks((sample2(k,:,p,1)),'MinPeakProminence',prom,'MinPeakDistance',pdist,'MinPeakWidth',minpwidth,'MaxpeakWidth',maxpwidth);
        [~,t2] = findpeaks((sample2(k,:,p,2)),'MinPeakProminence',prom,'MinPeakDistance',pdist,'MinPeakWidth',minpwidth,'MaxpeakWidth',maxpwidth);
        for i = 1:length(t2)
            dummy = t2(i) - t1;
            dummy = min(dummy(dummy>0));
            if ~isempty(dummy)
                if dummy > 50, allt(i,p,k) = NaN; e=e+1;
                else, allt(i,p,k) = dummy;  
                end
            else, e=e+1; 
            end
        end
    end
end

n = sum(~isnan(allt(:)))
error = e/(e+n)
avg = (e+n)/length(sample2(1,1,:,1))
mean_vel = nanmean(allt(:))
std_vel = nanstd(allt(:))

% plot velocity distribution
figure
histogram((allt(~isnan(allt(:)))),'BinMethod','integer','Normalization','probability')
xlabel('Particle Velocity (px/frame)')
ylabel(['Probability, n = ',num2str(n)])
title('Particle Velocity Distribution')
% xticks(linspace(0,50,10))
% xticklabels(round(linspace(0,50,10)/1.66*info.frameRate/10^6,3))

toc

end
%%
% now that velocity is known, calculate how many frames it would take for a
% particle to move an entire window length. skip that many frames when
% loading the data

skip_frames = round(window_length/mean_vel,0);
xfr = frames(1):skip_frames:frames(2);
tic
data = cineRead2(cine_folder,cine_file,xfr,info,LinLUT,...
    window_height,window_length,window_origin);
toc
if cells, data = (data + min(-data(:))); else data = -data; end


%%
% now find location of each particle, and keep track of overall position

tic
npeaks = zeros(length(data(1,1,:)),length(data(:,1,1)));
allws = NaN(length(data(1,1,:)),length(data(:,1,1)),80);
all_x = [];
all_y = [];
test = zeros(window_length,length(data(1,1,:)));

u = 0;
for p = 1:length(data(1,1,:))
    xlocs = [];
    ylocs = [];
    for k = 1:length(data(:,1,1))
        [~,locs,ws] = findpeaks((data(k,:,p)),'MinPeakProminence',prom,'MinPeakDistance',pdist,'MinPeakWidth',minpwidth,'MaxpeakWidth',maxpwidth');
        npeaks(p,k) = length(locs);
        allws(p,k,1:length(ws)) = ws;
        xlocs = [xlocs locs];
        ylocs = [ylocs window_height(k)*ones(1,length(xlocs))];
    end
    test(xlocs,p) = 1;
    [xlocs, idx] = sort(xlocs);
    ylocs = ylocs(idx);
    xlocs = xlocs+length(data(1,:,1))*(p-1);
    if length(xlocs) ~= length(unique(xlocs)), u = u+1; end
    all_x = [all_x xlocs];
    all_y = [all_y ylocs];
end

% to find distance between particles, subtract overall positions
clocs = circshift(all_x,1);
old_d = all_x(2:end)-clocs(2:end);

d = zeros(1,length(all_x)-1);
d_type = zeros(1,length(all_x)-1);
for i = 1:(length(all_x)-1)
    d(i) = norm([all_x(i+1) all_y(i+1)] - [all_x(i) all_y(i)]);
    if all_y(i+1) == all_y(i), d_type(i) = 1; end
end
d_type = logical(d_type);

% bin the data
[old_h, old_e] = (histcounts(old_d,'BinMethod','integers','Normalization','probability'));
[h,e1] = (histcounts(d,'BinMethod','integers','Normalization','count'));
[h1,e2] = (histcounts(d(d_type),'BinMethod','integers','Normalization','count'));
[h0,e3] = (histcounts(d(~d_type),'BinMethod','integers','Normalization','count'));
h1 = [zeros(1,e1(1)-.5) h1];
h0 = [zeros(1,e1(1)-.5) zeros(1,e3(1)-e2(1)) h0];
h = h./length(d); h1 = h1./length(d); h0 = h0./length(d);
toc

% plot location density for each frame
figure
imagesc(test)

% plot total particle count per pixel
figure
plot(sum(test'))

% plot histogram for peak width ~ particle size
figure
histogram(allws(~isnan(allws(:))),50,'Normalization','probability')
xlabel('Particle Size (peak width, px)')
ylabel(['Probability, n = ',num2str(length(d)+1)])
title('Particle Size Distribution')

% plot distance histogram
figure
hold on
plot(old_h(2:end),'k--'),plot(h1),plot(h0)
xlim([0 length(h)])
xticks(linspace(0,length(h),25))
% xticklabels(floor(linspace(0,length(h),25)/1.66))
xlabel('Particle Distance (px)')
ylabel(['Probability, n = ',num2str(length(d))])
title('Particle Distance Distribution')


% calculate and plot train size
dmax = 10*5;
train = find(d>dmax);
train = train - circshift(train,1);
train(1)=d(1);
figure
histogram(train,'BinMethod','integer','Normalization','probability')
xlabel('Particle Train Size (# beads)')
ylabel(['Probability, n = ',num2str(length(train))])
title(['Particle Train Size Distribution, dmax = ',num2str(dmax/1.66),' um'])



%%
% save data and process to excel file
tic
[~,dummy] = xlsread('analyzed.xlsx');
[row, ~] = size(dummy);
row = row+1;

cinesplit = strsplit(cine_file,'-');

xlswrite('analyzed.xlsx',{date},'Sheet1',['A',num2str(row)])
xlswrite('analyzed.xlsx',{cine_folder},'Sheet1',['B',num2str(row)])
xlswrite('analyzed.xlsx',{cine_file},'Sheet1',['C',num2str(row)])
xlswrite('analyzed.xlsx',cinesplit(1),'Sheet1',['D',num2str(row)])
xlswrite('analyzed.xlsx',cinesplit(2),'Sheet1',['E',num2str(row)])
xlswrite('analyzed.xlsx',cinesplit(4),'Sheet1',['F',num2str(row)])
xlswrite('analyzed.xlsx',info.frameRate,'Sheet1',['G',num2str(row)])
xlswrite('analyzed.xlsx',{num2str(window_height)},'Sheet1',['H',num2str(row)])
xlswrite('analyzed.xlsx',window_origin,'Sheet1',['I',num2str(row)])
xlswrite('analyzed.xlsx',window_length,'Sheet1',['J',num2str(row)])
xlswrite('analyzed.xlsx',{num2str(frames)},'Sheet1',['K',num2str(row)])
xlswrite('analyzed.xlsx',prom,'Sheet1',['L',num2str(row)])
xlswrite('analyzed.xlsx',pdist,'Sheet1',['M',num2str(row)])
xlswrite('analyzed.xlsx',maxpwidth,'Sheet1',['N',num2str(row)])
xlswrite('analyzed.xlsx',mean_vel,'Sheet1',['O',num2str(row)])
xlswrite('analyzed.xlsx',std_vel,'Sheet1',['P',num2str(row)])
xlswrite('analyzed.xlsx',n,'Sheet1',['Q',num2str(row)])

xlswrite('analyzed.xlsx',sum(test,2),'test',[char(row+63),'1'])
xlswrite('analyzed.xlsx',all_x(:),'all_x',[char(row+63),'1'])
xlswrite('analyzed.xlsx',all_y(:),'all_y',[char(row+63),'1'])
xlswrite('analyzed.xlsx',d(:),'d',[char(row+63),'1'])
xlswrite('analyzed.xlsx',d_type(:),'d_type',[char(row+63),'1'])
xlswrite('analyzed.xlsx',h1(:),'h1',[char(row+63),'1'])
xlswrite('analyzed.xlsx',h0(:),'h0',[char(row+63),'1'])
toc