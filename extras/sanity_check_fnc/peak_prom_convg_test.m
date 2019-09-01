prom = floor(linspace(120,190,20));
h = zeros(1200,length(prom));

for i = 1:length(prom)
    dummy = dtest(data,prom(i));
    h(1:length(dummy),i) = dummy;
end

figure
imagesc(corr(h))
yticks(1:length(prom))
yticklabels(prom)
xticks(1:length(prom))
xticklabels(prom)
colorbar