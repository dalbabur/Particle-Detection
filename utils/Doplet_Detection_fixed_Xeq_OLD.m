%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Droplet/Cell/Bead Detection
%
% Diego Alba 3/12/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cine_folder = 'H:\Encapsulation\Cell-Encapsulation-Videos-DA\cells\wide channel';
cine_file = '33.7mil_cellpml-22.5ulpmin-20x-2.cine';


window_height = [26 38]; % vector of pixels at which to read data
window_length = 1280; % number of pixles 
window_origin = 0; % offset 
frames = [9  576528]; % range of frames

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
importfile('LinLUT.mat'); %a conversion between packed 10bit data to real 10bit data
info = cineInfo(cine_folder,cine_file);



%%
% Get sample frames to calculate particle velocity
% Also check peak prominance (important for particle detection)

tic
fr = floor(info.NumFrames);
skip_frames = floor(fr/1000);

sample = zeros(length(window_height),window_length,length(1:skip_frames:fr),2);

% load data
for i = 1:2
    sample(:,:,:,i) = cineRead2(cine_folder,cine_file,[i:skip_frames:fr],info,LinLUT,...
        window_height,window_length,window_origin);
end
sample = (sample + min(-sample(:)));    

% use findpeaks on data
allpr = [];
for p = 1:length(sample(1,1,:,1))
    for k = 1:length(sample(:,1,1,1))
        [~,~,~,pr] = findpeaks((sample(k,:,p,1))); % only care about prominance first
        allpr = [allpr pr];
    end
end

% plot prominance
figure
histogram(allpr)
toc

%% set peak features for rest of analysis
tic
prom = 750;
pdist = 4;
minpwidth = 0;
maxpwidth = 7;
toc

%% check peak features for rest of analysis

cine_frame = 288268;
check = cineRead2(cine_folder,cine_file,cine_frame,info,LinLUT,...
        window_height,window_length,window_origin);
check = (check + min(-check(:)));    
findpeaks((check(1,:)),'MinPeakProminence',prom,'MinPeakDistance',pdist,'MinPeakWidth',minpwidth,'MaxpeakWidth',maxpwidth','Annotate','extents')


% check2 = cineRead2(cine_folder,cine_file,cine_frame,info,LinLUT,...
%         1:64,window_length,window_origin);
% 
% imagesc(check2)
%     
%%
% Calculate particle velocity 

tic
allt = NaN(80,length(sample(1,1,:,1)),2);
e = 0;
for p = 1:length(sample(1,1,:,1))
    locss = [];
    for k = 1:length(sample(:,1,1,1))
        [~,t1] = findpeaks((sample(k,:,p,1)),'MinPeakProminence',prom,'MinPeakDistance',pdist,'MinPeakWidth',minpwidth,'MaxpeakWidth',maxpwidth);
        [~,t2] = findpeaks((sample(k,:,p,2)),'MinPeakProminence',prom,'MinPeakDistance',pdist,'MinPeakWidth',minpwidth,'MaxpeakWidth',maxpwidth);
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
avg = (e+n)/length(sample(1,1,:,1))
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


%%
% now that velocity is known, calculate how many frames it would take for a
% particle to move an entire window length. skip that many frames when
% loading the data

skip_frames = round(window_length/mean_vel,0);
xfr = frames(1):skip_frames:frames(2);
tic
data = cineRead(cine_folder,cine_file,xfr(1),info,LinLUT);
toc
% data = (data + min(-data(:)));    


%%
tic
allc = [];
allr = [];
for i = 1:length(xfr)
A = cineRead(cine_folder,cine_file,xfr(i),info,LinLUT);
SE3 = strel('disk', 5);
B=imtophat(A,SE3);
% imshow(B)
contrastAdjusted = imadjust(B);
% imshow(contrastAdjusted)
% imshow(B)
% imshow(C)
[centers,radii] = imfindcircles(B,[1 5],'ObjectPolarity','bright', ...
          'Sensitivity',0.90);
% imagesc(A)
% viscircles(centers,radii)
allc = [allc; centers+1280*(i-1)];
allr = [allr; radii];
end
toc

%%
% now find location of each particle, and keep track of overall position

tic
npeaks = zeros(length(data(1,1,:)),length(data(:,1,1)));
allws = NaN(length(data(1,1,:)),length(data(:,1,1)),80);
l = [];
test = zeros(window_length,length(data(1,1,:)));

u = 0;
for p = 1:length(data(1,1,:))
    locss = [];
    for k = 1:length(data(:,1,1))
        [~,locs,ws] = findpeaks((data(k,:,p)),'MinPeakProminence',prom,'MinPeakDistance',pdist,'MinPeakWidth',minpwidth,'MaxpeakWidth',maxpwidth');
        npeaks(p,k) = length(locs);
        allws(p,k,1:length(ws)) = ws;
        locss = [locss locs];
    end
    test(locss,p) = 1;
    locss = sort(locss)+length(data(1,:,1))*(p-1);
    if length(locss) ~= length(unique(locss)), u = u+1; end
    l = [l locss];
end

% to find distance between particles, subtract overall positions
clocs = circshift(l,1);
d = l(2:end)-clocs(2:end);

% bin the data
h = (histcounts(d,'BinMethod','integers','Normalization','probability'));
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
plot(h)
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
xlswrite('analyzed.xlsx',l(:),'Sheet2',[char(row+63),'1'])
xlswrite('analyzed.xlsx',d(:),'Sheet3',[char(row+63),'1'])
xlswrite('analyzed.xlsx',sum(test')','Sheet4',[char(row+63),'1'])
