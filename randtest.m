dist = xlsread('analyzed.xlsx','d');
d_type = xlsread('analyzed.xlsx','d_type');
[~,~,info] = xlsread('analyzed.xlsx','Sheet1');
%%
W = [10 10 10 10 10 10 10 40 10 10 10 40 10 10];  % um
pxconv = [10/16 10/16 10/16 10/16 10/16 10/16 10/16 40/28 10/16 10/16 10/16 40/28 10/16 10/16]; % um/px

fnan = ~isnan(dist);
con = [];
for j = [9,10,11,13,14,12]
    dummy = [];
    n_frames = str2num(info{j,11});
    n_frames = n_frames(2)-n_frames(1);
    con = [con sum(~isnan(dist(:,j-1)))/(info{j,10}*pxconv(j)*W(j)*35*10^(-9)*n_frames/(info{j,10}/info{j,15}))];
    for i = 1:500
        dummy(:,i) = randperm(info{j,10}*round(n_frames/(info{j,10}/info{j,15})),sum(~isnan(dist(:,j-1))));
    end
    randlocs = sort(dummy);
    dummy = circshift(randlocs,[1,0]);
    randdist = randlocs(2:end,:)-dummy(2:end,:);
    dummy = (histcounts(randdist,'BinMethod','integer','Normalization','probability'));
    randhist(j,1:length(dummy)) = dummy;
end

d1 = d_type;
d1(isnan(d1)) = 0;
d0 = -(d_type-1);
d0(isnan(d0)) = 0;
d1 = logical(d1); d0 = logical(d0);

%%
dj = dist(:,8:10);
[dummy,e] = (histcounts(dj(fnan(:,8:10)),'BinMethod','integers','Normalization','count'));
[dummy1,e1] = (histcounts(dj(d1(:,1:3)),'BinMethod','integers','Normalization','count'));
[dummy0,e0] = (histcounts(dj(d0(:,1:3)),'BinMethod','integers','Normalization','count'));
dummy = [zeros(1,e(1)-.5) dummy];
dummy1 = [zeros(1,e(1)-.5) dummy1];
dummy0 = [zeros(1,e(1)-.5) zeros(1,e0(1)-e1(1)) dummy0];
dummy = dummy./(sum(sum(fnan(:,8:10))));
dummy1 = dummy1./(sum(sum(fnan(:,8:10))));
dummy0 = dummy0./(sum(sum(fnan(:,8:10))));

normed(1,1:500) = dummy(1:500)./mean(randhist(8:10,1:500));
normed(2,1:500) = dummy1(1:500)./mean(randhist(8:10,1:500));
normed(3,1:500) = dummy0(1:500)./mean(randhist(8:10,1:500));

dj = dist(:,12:13);
[dummy,e] = (histcounts(dj(fnan(:,12:13)),'BinMethod','integers','Normalization','count'));
[dummy1,e1] = (histcounts(dj(d1(:,5:6)),'BinMethod','integers','Normalization','count'));
[dummy0,e0] = (histcounts(dj(d0(:,5:6)),'BinMethod','integers','Normalization','count'));
dummy = [zeros(1,e(1)-.5) dummy];
dummy1 = [zeros(1,e(1)-.5) dummy1];
dummy0 = [zeros(1,e(1)-.5) zeros(1,e0(1)-e1(1)) dummy0];
dummy = dummy./(sum(sum(fnan(:,12:13))));
dummy1 = dummy1./(sum(sum(fnan(:,12:13))));
dummy0 = dummy0./(sum(sum(fnan(:,12:13))));

normed(4,1:500) = dummy(1:500)./mean(randhist(12:13,1:500));
normed(5,1:500) = dummy1(1:500)./mean(randhist(12:13,1:500));
normed(6,1:500) = dummy0(1:500)./mean(randhist(12:13,1:500));


dj = dist(:,11);
[dummy,e] = (histcounts(dj(fnan(:,11)),'BinMethod','integers','Normalization','count'));
[dummy1,e1] = (histcounts(dj(d1(:,4)),'BinMethod','integers','Normalization','count'));
[dummy0,e0] = (histcounts(dj(d0(:,4)),'BinMethod','integers','Normalization','count'));
dummy = [zeros(1,e(1)-.5) dummy];
dummy1 = [zeros(1,e(1)-.5) dummy1];
dummy0 = [zeros(1,e(1)-.5) zeros(1,e0(1)-e1(1)) dummy0];
dummy = dummy./(sum(sum(fnan(:,11))));
dummy1 = dummy1./(sum(sum(fnan(:,11))));
dummy0 = dummy0./(sum(sum(fnan(:,11))));

normed(7,1:500) = dummy(1:500)./(randhist(11,1:500));
normed(8,1:500) = dummy1(1:500)./(randhist(11,1:500));
normed(9,1:500) = dummy0(1:500)./(randhist(11,1:500));
%% 
% plot w/ normalized x (check pixels/particle in PCC) (3um are like 4px,
% cells are like 7px)

figure
hold on
plot(linspace(0,200/(3/pxconv(3)),200),normed(4,1:200))
plot(linspace(0,200/(3/pxconv(4)),200),normed(1,1:200))
plot(linspace(0,200/(10/pxconv(8)),200),normed(7,1:200))
plot([0 200],[1 1],'k--')
xlim([0 12])
xticks(1:15)
legend('0.055 %v/v - 3um beads, 10um channel','0.086 %v/v - 3um beads, 10um channel','0.45 %v/v - 10um cells, 40um channel',' random distribution')
xlabel('Particle Distance (# diameters)')
ylabel('Normalization Factor')

%%
% calculate and plot train size
dmax = 5*5;
train = find(randdist(9,:)>dmax);
train = train - circshift(train,1);
train(1)=d(1);
figure
histogram(train,'BinMethod','integer','Normalization','probability')
xlabel('Particle Train Size (# beads)')
ylabel(['Probability, n = ',num2str(length(train))])
title(['Particle Train Size Distribution, dmax = ',num2str(dmax/1.66),' um'])
