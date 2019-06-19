
l = data(~isnan(data(:,1)),1);
for j = 1:length(l)
    sort(abs(l-l(j)));  
end

%%
l1 = l(~isnan(l(:,1)),1);
n = 30;
ntot = length(l1);
dist = NaN(ntot-1, n);
hs = NaN(ntot-1,n);

for i = 1:n
    clocs = circshift(l1,i);
    dist(1:(end-i+1),i) = l1(i+1:end)-clocs(i+1:end);
    dummy = histcounts(dist(:,i),'BinMethod','integer');
    hs(1:length(dummy),i) = dummy;
end
histogram(dist(:),'BinMethod','integer')


%%
randh = NaN(6,3100);
normed = randh;
normed2 = normed;
n = round(14205*1.76*[1:6],0);
for j = 1:6
dummy = [];
for i = 1:500
%     dummy(:,i) = randperm(3198000,sum(~isnan(data(:,j))));
    dummy(:,i) = randperm(3198000,n(j));
end
rando = sort(dummy);
dummy = circshift(rando,[1,0]);
rand = rando(2:end,:)-dummy(2:end,:);
dummy = (histcounts(rand,'BinMethod','integer','Normalization','probability'));
randh(j,1:length(dummy)) = dummy;
% dummy = (histcounts(data2(:,j),'BinMethod','integers','Normalization','probability'));
% normed(j,1:length(dummy)) = dummy./randh(j,1:length(dummy));
end

%%
figure
hold on
for i = 1:6
plot(1:100,i*s*exp(-i*s*(1:100)))
end