function [Xscore, Yscore, Xweights, Yweights, betas, R] = pls_ycgosu(X, Y, ncomp)

 originalY = Y;
 originalX = X;
 
 Y0 = Y - mean(Y);X0 = X - mean(X);
 
 [Xrow,Xcol] = size(X0);
 [Yrow,Ycol] = size(Y0);
 
 if Xrow ~= Yrow
     error('Size mismatch!')
 end
 
 if nargin < 2
     error('Not enough inputs!')
 elseif nargin == 2
     numcomp = min(Xrow-1, Xcol);
 elseif nargin == 3
     if ~(isnumeric(ncomp))
         error('Last argument input should be a positive integer!')
     elseif ncomp == 0
         numcomp = min(Xrow-1, Xcol);
     else
         numcomp = ncomp;
     end
 elseif nargin > 4
     error('Too much inputs!')
 end
 
X = X0;
Y = Y0; 

% preallocation
Xscore= zeros(Xrow, numcomp);
Yscore = zeros(Xrow, numcomp);
Xweights = zeros(Xcol, numcomp);
Yweights = zeros(Ycol, numcomp);
b = zeros(numcomp, 1);

for i = 1:numcomp
    
    [W, ~, Q] = svd(X' * Y, 'econ');
    t = X * W(:, 1); % w
    u = Y * Q(:, 1); % q, if Y is univariate, first u is equal to Y(mean centered).
    Xweights(:, i) = W(:, 1);
    Yweights(:, i) = Q(:, 1);
    
    Xscore(:, i) = t;
    Yscore(:, i) = u;
    b(i, :) = u'*t / (t'*t); % isn't corr(u, t) more proper? and what would increasing b mean?
    
    X = X - t * (t' * X) / (t' * t);
    Y = Y - t * (t' * Y) / (t' * t);
    
    % X is "deflated" in each for loop. It means information on the space
    % which has maximum covariance with Y and X are removed. Consequently,
    % data after deflation are moved to orhtogonal space of the previous space.
    % So, cov(Xscore) yields diagonal matrix. 
    
end

% but we can also derive combination coefficients without deflating the data.
% say every t is some combination of X0(mean centered X).
% so we can write it as T = XR, where R is the combination matrix of X
% that has the combination coefficients of X.
% and also we have T matrix(Xscore) above.

T = Xscore;
R = pinv(X0) * T;
% obtained R is similar to Xweights but slightly different in that R' * R
% does not yield identity matrix. 

betas = R * diag(b) * Yweights';
% now we can predict the data with this beta.
% beta has p x m set of multivariate regression coefficients.

betas = [mean(originalY, 1) - mean(originalX, 1) * betas;betas];
% adding intercepts, only needed for reconstructing originalY. 
% originalY = [ones(n,1) originalX]*beta + Yresiduals
% Y0 = X0*beta(2:end,:) + Yresiduals
              
end

% main difference between simpls(already implemented in matlab. refer to function
% 'plsregress'.) and NIPALS algorithm is in process of deflation.
% In the NIPALS algorithm, X and Y are both regressed out respect to obtained t,
% which is a weighted value of X.( t = X * Xweights. and Xweights here contains unit vector which are pointing the direction
% that maximize the covariance between X*Xweights and Y*Yweights)

% However, in the simpls algorithm, X' * Y are regressed out respect to P which
% is a regression coeffients of X respect to t. (P = (X'* t)/(t'*t))
