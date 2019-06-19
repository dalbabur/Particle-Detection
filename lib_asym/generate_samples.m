function samples = generate_samples(tmpCost, samples)
numSamples = length(tmpCost);
xstd = std(samples(:,1));
ystd = std(samples(:,2));
rotstd = std(samples(:,3));

tmpCost = 1./(tmpCost+0.01);
nCost = ceil(tmpCost/sum(tmpCost)*numSamples);
[num, I] = sort(nCost, 'descend');
nSamples = zeros(numSamples,3);
cnt1 = 1;
cnt2 = 1;
flag = 0;
while (1)
    nS = num(cnt2);
    tmp1 = cnt1+nS-1;
    if  tmp1 == numSamples
        flag = 1;
    elseif tmp1 > numSamples        
        tmp1 = numSamples;
        flag = 1;
        nS = tmp1-cnt1+1;
    end
    nSamples(cnt1:tmp1,1) = xstd/5*randn(nS,1)+samples(I(cnt2),1);
    nSamples(cnt1:tmp1,2) = ystd/5*randn(nS,1)+samples(I(cnt2),2);
    nSamples(cnt1:tmp1,3) = rem(rotstd/5*randn(nS,1)+samples(I(cnt2),3), 360);
    cnt1 = cnt1+nS;
    cnt2 = cnt2+1;
    if flag
        break;
    end
end
samples = nSamples;