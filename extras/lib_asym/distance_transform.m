function [dist1, dist2, dist3] = distance_transform(edge1, edge2, edge3)
dist1 = bwdist(edge1, 'euclidean');
dist2 = bwdist(edge2, 'euclidean');
dist3 = bwdist(edge3, 'euclidean');