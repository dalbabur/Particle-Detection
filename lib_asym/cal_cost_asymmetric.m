function cost = cal_cost_asymmetric(dist1, dist2, dist3, e1, e2, e3, samples, cMass, fCMass, ef1, ef2, ef3)
numSamples = size(samples,1);
cost  = zeros(numSamples,1);

[height, width] = size(dist1);
e1(:,1) = e1(:,1)-cMass(1);
e1(:,2) = e1(:,2)-cMass(2);

e3(:,1) = e3(:,1)-cMass(1);
e3(:,2) = e3(:,2)-cMass(2);

lene1 = size(e1,1);
lene3 = size(e3,1);


ef1(:,1) = ef1(:,1)-fCMass(1);
ef1(:,2) = ef1(:,2)-fCMass(2);

ef3(:,1) = ef3(:,1)-fCMass(1);
ef3(:,2) = ef3(:,2)-fCMass(2);


for i = 1:numSamples
    
    if samples(i,4) == 0
        count = 0;
        x = samples(i,1); y=samples(i,2); rot=samples(i,3);
        rMat = [cosd(rot) -sind(rot);
        sind(rot) cosd(rot)];
        % for positive edge (time domain)
        tmp1 = rMat*e1';
        nx = round(tmp1(1,:)+x);
        ny = round(tmp1(2,:)+y);
        for j = 1:lene1
            if nx(j) >= 1 && ny(j) >= 1 && nx(j) <= width && ny(j) <= height
                cost(i) = cost(i) + dist1(ny(j), nx(j));
                count = count+1;
            end
        end


        % for spatial edge
        tmp1 = rMat*e3';
        nx = round(tmp1(1,:)+x);
        ny = round(tmp1(2,:)+y);
        for j = 1:lene3
            if nx(j) >= 1 && ny(j) >= 1 && nx(j) <= width && ny(j) <= height
                cost(i) = cost(i) + dist3(ny(j), nx(j));
                count = count+1;
            end
        end
        if count >= (lene3+lene1)/2
            cost(i) = cost(i)/count;
        else
            cost(i) = inf;
        end
    else
        count = 0;
        x = samples(i,1); y=samples(i,2); rot=samples(i,3);
        rMat = [cosd(rot) -sind(rot);
        sind(rot) cosd(rot)];
        % for positive edge (time domain)
        tmp1 = rMat*ef1';
        nx = round(tmp1(1,:)+x);
        ny = round(tmp1(2,:)+y);
        for j = 1:lene1
            if nx(j) >= 1 && ny(j) >= 1 && nx(j) <= width && ny(j) <= height
                cost(i) = cost(i) + dist1(ny(j), nx(j));
                count = count+1;
            end
        end


        % for spatial edge
        tmp1 = rMat*ef3';
        nx = round(tmp1(1,:)+x);
        ny = round(tmp1(2,:)+y);
        for j = 1:lene3
            if nx(j) >= 1 && ny(j) >= 1 && nx(j) <= width && ny(j) <= height
                cost(i) = cost(i) + dist3(ny(j), nx(j));
                count = count+1;
            end
        end
        if count >= (lene3+lene1)/2
            cost(i) = cost(i)/count;
        else
            cost(i) = inf;
        end        

    end
end