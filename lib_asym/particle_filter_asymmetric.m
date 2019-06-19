function [pose, cost] = particle_filter_asymmetric(dist1, dist2, dist3, e1, e2, e3, cMass, numSamples, numIteration, objImg)
% pose is row, col tranlation and rotation
[height1, width1] = size(dist1);
[h, ~] = size(objImg);
flipE1 = e1;
flipE2 = e2;
flipE3 = e3;
flipE1(:,2) = h - flipE1(:,2);
flipE2(:,2) = h - flipE2(:,2);
flipE3(:,2) = h - flipE3(:,2);
fCMass      = cMass;
fCMass(2)   = h - fCMass(2);

% initial random samples
samples = [rand(2*numSamples, 1)*width1, rand(2*numSamples, 1)*height1, 360*rand(2*numSamples,1) randi(2,[2*numSamples,1])-1];
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
    tmpCost = cal_cost_asymmetric(dist1, dist2, dist3, e1, e2, e3, samples, cMass, fCMass, flipE1, flipE2, flipE3);
    % find min and save it
    [tempMin, idx] = min(tmpCost);
    if minCost > tempMin
        pose = samples(idx,:);
        cost = tempMin;
    end
    
    % generate sample based on the cost (gaussian)
    if i ~= numIteration
        if min(tmpCost) == inf
            samples = [rand(numSamples, 1)*width1, rand(numSamples, 1)*height1, 360*rand(numSamples,1) randi(2,[2*numSamples,1])-1];
        else
            samples = generate_samples_asymmetric(tmpCost, samples);
        end
    end
end