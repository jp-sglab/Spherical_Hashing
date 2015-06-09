

% xData : matrix ( number of data * dimension )

bit = 64;                      % binary code length
nTrain = 20000;                 % number of training samples

% compute training set with random sampling
[nData, dim] = size(xData);
R = randperm( nData );  xTrain = xData(R(1:nTrain),:);

% compute centers and radii of hyper-spheres with training set
[centers, radii] = SphericalHashing( xTrain , bit );

% compute distances from centers
dData = distMat( xData , centers );

% compute binary codes for data points
th = repmat( radii' , size( dData , 1 ) , 1 );
bData = zeros( size(dData) );
bData( dData <= th ) = 1;
bData = compactbit(bData);
