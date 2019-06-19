function h = dtest(data,prom)

npeaks = zeros(length(data(1,1,:)),length(data(:,1,1)));
allws = zeros(size(npeaks));
l = [];

u = 0;
for p = 1:length(data(1,1,:))
    locss = [];
    for k = 1:length(data(:,1,1))
        [~,locs,ws] = findpeaks(-(data(k,:,p)),'MinPeakProminence',prom);
        npeaks(p,k) = length(locs);
        allws(p,k) = mean(ws);
        locss = [locss locs];
    end
    locss = sort(locss)+length(data(1,:,1))*(p-1);
    if length(locss) ~= length(unique(locss)), u = u+1; end
    l = [l locss];
end

clocs = circshift(l,1);
d = l(2:end)-clocs(2:end);
        
h = (histcounts(d,'BinMethod','integers','Normalization','probability'));

end