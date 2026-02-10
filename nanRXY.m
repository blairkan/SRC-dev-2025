function [Rxy,Rxx,Ryy,Ryx] = nanRXY(X,Y)
% BK: Compute all covariance matrices in one big calculation. Concat input
% X and Y matrices and then return specific inputs as different row/column
% regions of the output covariance matrix. 

% BK: This is not doing regularization, rather specifying row/column
% indices for cutting RXY into final function outputs
D=size(X,1);

% BK: Size of concat input matrix is [time x (nXFeatures+nYFeatures)]. Size
% of output matrix is square nXFeatures+nYFeatures (think X'X) temporal
% covariance. 
RXY=nancov([X.' Y.'],'pairwise');

% Upper left
Rxx=RXY(1:D,1:D);

% Lower right
Ryy=RXY(D+1:end,D+1:end);

% Upper right
Rxy=RXY(1:D,D+1:end);

% Lower left
if nargout==4, Ryx=RXY(D+1:end,1:D); end

% D=size(X,1);
% RXY=nancov([X.' Y.'],'pairwise');
% Rxx1=RXY(1:D,1:D);
% Rxx2=RXY(2*D+1:3*D,2*D+1:3*D);
% 
% Ryy=RXY(D+1:2*D,D+1:);
% Rxy=RXY(1:D,D+1:end);