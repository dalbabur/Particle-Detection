function p = greedy_matching(d1, d2)
len1 = size(d1,1);
len2 = size(d2,1);
p = zeros(len1,2);

for i = 1:len1
    p(i,1) = i;
    tmp = repmat(d1(i,:), len2,1);
    [~, p(i,2)] = min(sum((d2-tmp).^2,2));
end
