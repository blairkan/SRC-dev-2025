function [H,W,rhos,pvals,U,V,Rxx,Ryy] = computeCCA(X,Y,Kx,Ky)
%[H,W,rhos,pvals,U,V] = myCanonCorr(X,Y,Kx,Ky) 
% -----------------------------------------------------------
% Blair dev version 2026
% Regularized canonical correlation
% 
% Inputs
% - X: Stimulus feature convolution matrix: [time x lags] with intercept
%   column of 1s at the end
% - Y: EEG data: [time x electrode]
% - Kx: Regularization parameter for the stimulus
% - Ky: Regularization parameter for the EEG
%
% Outputs
% - H: Set of temporal filters [lag x component] 
%   --> Used to be "A"
% - W: Set of spatial filters [electrode x component] 
%   --> Used to be "B"
% - A: Forward model of W for topoplots 
%   --> New
% - U: Temporally filtered stimulus feature [time x component] 
%   --> transpose of previous output
% - V: Spatially filtered EEG [time x component] 
%   --> transpose of previous output
% - R: Struct of covariance matrices Rxx, Ryy, Rxy, Ryx 
%   --> struct instead of returning only two of them as separate variables

% History
% - 2/3/2026: Blair - adapted from myCanonCorr

if nargin<2, error('JD: at least two argument required'); end
if ~ismatrix(X), error('JD: X must have two dimensions'); end
if ~ismatrix(Y), error('JD: Y must have two dimensions'); end
if size(X,1)>size(X,2), X=X.'; warning('JD: transposing X'); end
if size(Y,1)>size(Y,2), Y=Y.'; warning('JD: transposing Y'); end


[Rxy,Rxx,Ryy,Ryx] = nanRXY(X,Y);
if nargin<4 || isempty(Ky), Ky=rank(Ryy); end
if nargin<3 || isempty(Kx), Kx=rank(Rxx); end

% compute A
Rxxnsq=regSqrtInv(Rxx,Kx); % regularized Rxx^(-1/2)
M = Rxxnsq*Rxy*regInv(Ryy,Ky)*Ryx*Rxxnsq; % textbook
M = (M+M')/2; % fix nummerical precision asymmetric
[C,Dc]=eig(M);  % textbook
[~,indx]=sort(diag(Dc),'descend'); 
C=C(:,indx(1:min(Kx,Ky))); % dump zero eigenvalue dimensions 
H=Rxxnsq*C; % invert coordinate transformation

% compute B
Ryynsq=regSqrtInv(Ryy,Ky); % regularized Ryy^(-1/2)
D=Ryynsq*Ryx*Rxxnsq*C;
W=Ryynsq*D;


%[A,Da]=eig(regInv(Rxx,Kx)*Rxy*regInv(Ryy,Ky)*Ryx);
%[B,Db]=eig(regInv(Ryy,Ky)*Ryx*regInv(Rxx,Kx)*Rxy);

if sum(~isreal(H(:))) | sum(~isreal(W(:)))
    warning('Imaginary components. Something is wrong!'); 
    H=real(H); W=real(W);
end

U=H.'*X;
V=W.'*Y;

nVars=min(size(U,1),size(V,1));
rhos=zeros(nVars,1);
pvals=zeros(nVars,1);
for n=1:nVars
    n
    [r,p]=corrcoef(U(n,:),V(n,:));
    rhos(n)=r(1,2);
    pvals(n)=p(1,2);
end

return

X=randn(5,1000);
Y=randn(5,1000);
[H,W,rhos,pvals,U,V] = myCanonCorr(X,Y,5,5);
rhos
pvals










