function [H, W, U, V, A, dC, R] = computeCCA(X, Y, Kx, Ky)
% [H, W, U, V, A, dC, R] = computeCCA(X, Y, Kx, Ky)
% -----------------------------------------------------------
% Blair dev version 2026
% Regularized canonical correlation
% 
% Inputs (required)
% - X: Stimulus feature convolution matrix: [time x lags] with intercept
%   column of 1s at the end
% - Y: EEG data: [time x electrode]
%
% Inputs (optional)
% - Kx: Regularization parameter for the stimulus. If empty or not
%   specified, will default to the rank of covariance matrix Rxx.
% - Ky: Regularization parameter for the EEG. If etmpty or not specified,
%   will default to the rank of the covariance matrix Ryy. 
%
% Outputs
% - H: Set of temporal filters [lag x component] 
%   --> Used to be "A"
% - W: Set of spatial filters [electrode x component] 
%   --> Used to be "B"
% - U: Temporally filtered stimulus feature [time x component] 
%   --> transpose of previous output
% - V: Spatially filtered EEG [time x component] 
%   --> transpose of previous output
% - A: Forward model of W for topoplots 
%   --> New
% - dC: Eigenvalues of the eigenvalue decomposition
%   --> New
% - R: Struct of covariance matrices Rxx, Ryy, Rxy, Ryx 
%   --> struct instead of returning only two of them as separate variables

% History
% - 2/3/2026: Blair - adapted from myCanonCorr

%% Check inputs
if nargin<2, error('JD: at least two argument required'); end
if ~ismatrix(X), error('JD: X must have two dimensions'); end
if ~ismatrix(Y), error('JD: Y must have two dimensions'); end
if size(X,1)>size(X,2), X=X.'; warning('JD: transposing X'); end
if size(Y,1)>size(Y,2), Y=Y.'; warning('JD: transposing Y'); end

%% Compute covariance matrices
[Rxy,Rxx,Ryy,Ryx] = nanRXY(X,Y);

% Store covariance matrices in output struct
R.Rxy = Rxy; R.Rxx = Rxx; R.Ryy = Ryy; R.Ryx = Ryx;

% Provide default Kx, Ky values if needed
if nargin<3 || isempty(Kx), Kx=rank(Rxx); end
if nargin<4 || isempty(Ky), Ky=rank(Ryy); end

%% Compute temporal filters H 
Rxxnsq=regSqrtInv(Rxx,Kx); % regularized Rxx^(-1/2)
M = Rxxnsq*Rxy*regInv(Ryy,Ky)*Ryx*Rxxnsq; % textbook
M = (M+M')/2; % fix nummerical precision asymmetric
[C,Dc]=eig(M);  % textbook
[dC,indx]=sort(diag(Dc),'descend'); % NEW - save also first output as dC
C=C(:,indx(1:min(Kx,Ky))); % dump zero eigenvalue dimensions 
H=Rxxnsq*C; % invert coordinate transformation

%% Compute spatial filters W
Ryynsq=regSqrtInv(Ryy,Ky); % regularized Ryy^(-1/2)
D=Ryynsq*Ryx*Rxxnsq*C;
W=Ryynsq*D;


%[A,Da]=eig(regInv(Rxx,Kx)*Rxy*regInv(Ryy,Ky)*Ryx);
%[B,Db]=eig(regInv(Ryy,Ky)*Ryx*regInv(Rxx,Kx)*Rxy);

if sum(~isreal(H(:))) | sum(~isreal(W(:)))
    warning('Imaginary components. Something is wrong!'); 
    H=real(H); W=real(W);
end

%% Compute spatial filter forward models A 
A = Ryy * W * inv(W' * Ryy * W);

%% Project input data through respective filters

%%% Blair note: Figure out later why the dimensions are like this:
% - H: Lag x component
% - X: lag x time
% - W: Channel x component
% - Y: Channel x time
U = X.' * H;
V = Y.' * W;

return

% nVars=min(size(U,1),size(V,1));
% rhos=zeros(nVars,1);
% pvals=zeros(nVars,1);
% for n=1:nVars
%     n
%     [r,p]=corrcoef(U(n,:),V(n,:));
%     rhos(n)=r(1,2);
%     pvals(n)=p(1,2);
% end
% 
% return










