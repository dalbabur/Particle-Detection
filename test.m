figure
k=1;
for i= 1:2
for j = 1:2
subplot(2,2,k)
imagesc(corr(squeeze(all_h(:,i,:)),squeeze(all_h(:,j,:))));
title(skips(j))
caxis([0.9 1]); k=k+1;
end
end

all_h2 = reshape(all_h,1200,25);
figure
first = corr(all_h2);
second = corr(first);
imagesc(second)

f = 26;
figure
findpeaks(-sample(1,:,f,1),'MinPeakProminence',prom,'MinPeakDistance',4)
hold on
findpeaks(-sample(1,:,f,2),'MinPeakProminence',prom,'MinPeakDistance',4)

for i = 1:length(d)
   ds(i,:) = d(i)-d;
end