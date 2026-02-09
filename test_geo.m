% Example usage with pdist2:
X = [40.71,-74.01; 48.86,2.35 ; 52.30, 13.25]; % Matrix of points 
Y = [40.71,-74.01; 48.86,2.35]; % Matrix of points
D = pdist2(X, Y, @geodist)
DA = pdist2(X,Y)

% % A custom function to calculate geodesic distance between coordinates (conceptual example)
% function D2 = geodesic_pdist(ZI, ZJ)
%     % ZI is 1-by-2, ZJ is m2-by-2 (lat, lon pairs)
%     m2 = size(ZJ, 1);
%     D2 = zeros(m2, 1);
%     for k = 1:m2
%         % Use the Mapping Toolbox distance function inside the loop
%         % This iterates for each pair, which is less efficient than vectorized calls
%         D2(k) = distance(ZI(1), ZI(2), ZJ(k, 1), ZJ(k, 2), referenceEllipsoid('WGS84'));
%     end
% end
