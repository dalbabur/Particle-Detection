function [pose, cost] = particle_filter(dist1, dist2, dist3, e1, e2, e3, cMass, numSamples, numIteration)
% pose is row, col tranlation and rotation
[height1, width1] = size(dist1);


% initial random samples
samples = [rand(numSamples, 1)*width1, rand(numSamples, 1)*height1, 360*rand(numSamples,1)];
minCost = inf;
pose    = samples(:,1);
cost    = inf;
for i = 1:numIteration
%     figure(1)
%     imshow(dist1)
%     hold on
%     plot(samples(:,1), samples(:,2), '.')
%     hold off
%     pause
    % calculate cost
    tmpCost = cal_cost(dist1, dist2, dist3, e1, e2, e3, samples, cMass);
    % find min and save it
    [tempMin, idx] = min(tmpCost);
    if minCost > tempMin
        pose = samples(idx,:);
        cost = tempMin;
    end
    
    % generate sample based on the cost (gaussian)
    if i ~= numIteration
        if min(tmpCost) == inf
            samples = [rand(numSamples, 1)*width1, rand(numSamples, 1)*height1, 360*rand(numSamples,1)];
        else
            samples = generate_samples(tmpCost, samples);
        end
    end
end