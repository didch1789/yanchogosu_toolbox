function [Xscore, Yscore, Xweights, Yweights, betas, R, numFeatures] = pls_ycgosu_v2(X, Y, varargin)
% X, Y
% 
% numcomp 0 will start autocheck

originalY = Y;
originalX = X;
 
Y0 = Y - mean(Y);X0 = X - mean(X);
 
[Xrow,Xcol] = size(X0);
[Yrow,Ycol] = size(Y0);

if Xrow ~= Yrow
 error('Size mismatch!')
end

autocheck = 0;
numcomp = min(Xrow-1, Xcol);
 
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'numcomponents', 'numcomps', 'numcomp'}
                numcomp = varargin{i+1};
                numFeatures = numcomp;
                autocheck = 0;
                if numcomp == 0
                    numcomp = min(Xrow-1, Xcol);
                    autocheck = 1;
                end
        end
    end
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
    u = Y * Q(:, 1); % q, if Y is univariate, first u is equal to Y(mean centered) and Q is always just 1.
    Xweights(:, i) = W(:, 1);
    Yweights(:, i) = Q(:, 1);
    
    Xscore(:, i) = t;
    Yscore(:, i) = u;
    b(i, :) = u'*t / (t'*t); % b is the regression coefficents. t as the independent, u as the dependent variable.
    
    Ylast = Y;
    X = X - t * (t' * X) / (t' * t); % regress out t
    Y = Y - t * (t' * Y) / (t' * t);
    
    if autocheck
        if mean(var(Ylast, 0, 1) - var(Y, 0, 1)) < 0.01
            Xscore = Xscore(:, 1:i);
            Yscore = Yscore(:, 1:i);
            Xweights = Xweights(:, 1:i);
            Yweights = Yweights(:, 1:i);
            b = b(1:i, :);
            numFeatures = i;
            break
        end
    end
    
    % X is "deflated" in each for loop. It means information on the space
    % which has maximum covariance with Y and X are removed. Consequently,
    % data after deflation are moved to orhtogonal space of the t.
    % So, cov(Xscore) yields diagonal matrix. 
    
    % since Y is also deflated, var(Y, 0, 1) gets smaller in the iteration process which means, once Y is
    % projected on to t, projected Y gets smaller variance and I just wanted to stop before it has very small variance.
    % if Y has no variance, it would be useless to derive svd(X'*Y). So I
    % added autocheck to stop the loop if the variance difference between Ylast and Y is
    % small.
    
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
% In the NIPALS algorithm, both X and Y are regressed out respect to obtained t,
% which is a weighted value of X.( t = X * Xweights. and Xweights here contains unit vector which are pointing the direction
% that maximize the covariance between X and Y)
% However, in the simpls algorithm, X' * Y are regressed out respect to P which
% is a regression coeffients of X respect to t. (P = (X'* t)/(t'*t))

% one special note to add...!
% getting the full components possible (like extracting 9 comps from 10 x 100 matrix)
% will yield same results with full components possbile from PCA since and
% orthogonal 9 components would retain all the variance in the original 10 x 100 matrix.
