% A custom function to calculate geodesic distance between coordinates
function D = geodist(X1, X2)
    wgs84 = wgs84Ellipsoid("km");
    % ZI is 1-by-2, ZJ is m2-by-2 (lat, lon pairs)

    % [lat1, lat2] = meshgrid(ZI, ZI);
    % [lon1, lon2] = meshgrid(ZJ, ZJ);

    [lat1, lat2] = meshgrid(X1(:,1), X2(:,1));
    [lon1, lon2] = meshgrid(X1(:,2), X2(:,2));

    % Calculate all distances
    [D, ~] = distance(lat1, lon1, lat2, lon2, wgs84);
end