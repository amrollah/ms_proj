% =========================================================================
% Coded by Youngmin Ha (y.ha.1@research.gla.ac.uk) on 21 Jan 2015
%   by referring to Thayananthan, Arasanathan. "Template-based pose 
%   estimation and tracking of 3D hand motion." Cambridge, UK: Department 
%   of Engineering, University of Cambridge (2005).
% =========================================================================
% FUNCTION mrvr performs original multivariate relevance vector regression
%   developed by Thayananthan.
% ******************************** Input **********************************
% Phi:      matrix, design matrix
% T:        column vector or matrix, target (i.e. true function + noise)
% maxIters:	scalar, maximum iteration number of EM algorithm
% tolerance:scalar, tolerance value to check convergence of EM algorithm
% ******************************** Output *********************************
% used:     column vector, index of relevance vectors
% alpha: 	column vector, inverse variance of weight
% Mu:       column vector or matrix, mean of weight
% invSigma: matrix, inverse of covariance matrix of weight
% Omega:    scalar or matrix, estimated covariance matrix of noise
% nIters:   scalar, number of iterations of EM algorithm
% =========================================================================

function [used, alpha, Mu, invSigma, Omega, nIters] = mrvr(Phi, T, ...
    maxIters, tolerance)
%% constants
N = size(T,1); % # of training samples
V = size(T,2); % # of output dimensions

assert(size(Phi,1) == N, 'unexpected matrix size of Phi')
assert(size(Phi,2) == N + 1, 'unexpected matrix size of Phi')

%% initialisation
alpha = inf(N+1,1);
beta = 1./(.1.*var(T)); % beta = sigma.^-2
mask = logical(alpha < inf); % false means alpha=inf, true means alpha<inf

%% EM algorithm
isConverged = false;
for iterNum = 1:maxIters
    % update s', q', s, and q
    sPrime = nan(N+1,1);
    qPrime = nan(N+1,V);
    s = nan(N+1,V);
    q = nan(N+1,V);
    allAlphaInf = all(~mask);
    for i = 1:N+1
        phiSq = Phi(:,i)'*Phi(:,i);
        if allAlphaInf % if all alpha are inf
            for j = 1:V
                sPrime(i,j) = beta(j).*phiSq;
                qPrime(i,j) = beta(j).*(Phi(:,i)'*T(:,j));
            end
        else
            for j = 1:V
                temp = (beta(j).^2).*(Phi(:,i)'*PhiSigmaPhi(:,:,j));
                sPrime(i,j) = beta(j).*phiSq - temp*Phi(:,i);
                qPrime(i,j) = beta(j).*(Phi(:,i)'*T(:,j)) - temp*T(:,j);
            end
        end
        
        if ~mask(i) % isinf(alpha(i))
            s(i,:) = sPrime(i,:);
            q(i,:) = qPrime(i,:);
        else
            temp = alpha(i)./(alpha(i) - sPrime(i,:));
            s(i,:) = temp.*sPrime(i,:);
            q(i,:) = temp.*qPrime(i,:);
        end
    end
    
    % select i which maximizes deltaL(i)
    deltaL = nan(N+1,1);
    task = repmat('non', N+1, 1);
    alphaNew = nan(N+1,1);
    for i = 1:N+1
        % find roots of polynomial equation
        polyCoeff = zeros(V+1,2.*V+1);
        temp = V;
        for j = 1:V     
            temp = conv(temp, [1, 2.*s(i,j), s(i,j).^2]); % [alpha^2, alpha^1, alpha^0]
        end
        polyCoeff(1,:) = temp;
        for j = 1:V
            temp = [-1, -s(i,j) - q(i,j).^2, 0]; % [alpha^2, alpha^1, alpha^0]
            for k = 1:V
                if k ~= j
                    temp = conv(temp, [1, 2.*s(i,k), s(i,k).^2]); % [alpha^2, alpha^1, alpha^0]
                end
            end
            polyCoeff(j+1,:) = temp;
        end
        sumPolyCoeff = sum(polyCoeff);
        if sumPolyCoeff(1) ~= 0
            error('The order of polynomial should be 2V-1')
        end
        polyRoots = roots(sumPolyCoeff);
        
        % find positive real roots
        posRealRoots = polyRoots(real(polyRoots) > 0 & imag(polyRoots) == 0); 
        nPosRealRoots = length(posRealRoots);
        
        % update task and calculate deltaL
        switch nPosRealRoots
            case 0 % if # of positive real roots is 0
                alphaNew(i) = inf;
                [task(i,:), deltaL(i)] = calcDeltaL(mask(i), alpha(i), ...
                    alphaNew(i), sPrime(i,:), qPrime(i,:), s(i,:), q(i,:));
            case 1 % if # of positive real roots is 1
                alphaNew(i) = posRealRoots;
                [task(i,:), deltaL(i)] = calcDeltaL(mask(i), alpha(i), ...
                    alphaNew(i), sPrime(i,:), qPrime(i,:), s(i,:), q(i,:));
            otherwise % if # of positive real roots is more than 1
                taskRoots = repmat('non', nPosRealRoots, 1);
                deltaLRoots = nan(nPosRealRoots, 1);
                for k = 1:nPosRealRoots
                    [taskRoots(k,:), deltaLRoots(k)] = calcDeltaL(mask(i), alpha(i), ...
                        posRealRoots(k), sPrime(i,:), qPrime(i,:), s(i,:), q(i,:));
                end
                [deltaL(i), ind] = max(deltaLRoots);
                alphaNew(i) = posRealRoots(ind);
                task(i,:) = taskRoots(ind,:);
        end
    end
    if all(isnan(deltaL)) % if all tasks are non
        break
    elseif all(deltaL(~isnan(deltaL)) == -inf)
        warning('all values of deltaL are -inf')
        i = randi(N+1);
    else
        [~, i] = max(deltaL);
    end

    % estimate beta, where beta = sigma.^-2
    if iterNum ~= 1
        for j = 1:V
            beta(j) = (N - nMask + alpha(mask)'*diag(Sigma(:,:,j))) ...
            ./sum((T(:,j) - PhiMask*Mu(:,j)).^2);
        end
    end
    
    % update alpha(i)
    switch task(i,:)
        case 'est'
            changeLogAlpha = log(alpha(i)./alphaNew(i));
            alpha(i) = alphaNew(i);
            
            % check convergence criteria
            if abs(changeLogAlpha) < tolerance
                if all(isinf(alphaNew(~mask))) % ~mask means "out of the model"
                    isConverged = true;
                end
            end
        case 'add'
            alpha(i) = alphaNew(i);
            mask(i) = true;
        case 'del'
            alpha(i) = inf;
            mask(i) = false;
        otherwise
            error('task should be one of est, add, and del')
    end
    
    % update PhiMask, Sigma, and Mu
    PhiMask = Phi(:,mask);
    nMask = sum(mask);
    invSigma = zeros(nMask,nMask,V);
    Sigma = zeros(nMask,nMask,V);
    Mu = zeros(nMask,V);
    PhiSigmaPhi = zeros(N,N,V);
    PhiMaskPhiMask = PhiMask'*PhiMask;
    for j = 1:V
        invSigma(:,:,j) = beta(j).*PhiMaskPhiMask + diag(alpha(mask));
        assert(all(all(invSigma(:,:,j) == invSigma(:,:,j)')), ...
            'invSigma should be symmetric')

        % the below code is used instead of Sigma = inv(invSigma)
        % by referring to SB2_FullStatistics.m developed by Tipping
        % http://www.vectoranomaly.com/downloads/downloads.htm
        Upper = chol(invSigma(:,:,j));
        invU = inv(Upper);
        Sigma(:,:,j) = invU*invU';

        Mu(:,j) = beta(j).*(Sigma(:,:,j)*PhiMask'*T(:,j));
        PhiSigmaPhi(:,:,j) = PhiMask*Sigma(:,:,j)*PhiMask';
    end
        
    if isConverged
        break
    end    
end
nIters = iterNum;
disp(['1. # of iterations = ' num2str(iterNum)])

% index of relevance vectors
used = find(mask);

%% estimated covariance matrix of noise
% OmegaTilde is sample covariance matrix
yPredict = Phi(:,used)*Mu;
diff = T - yPredict;
OmegaTilde = diff'*diff./(N - 1);

% convert OmegaTilde to orrelation matrix Rhat
Rhat = corrcov(OmegaTilde);

% Omega is estimated covariance matrix of noise
Dhat = zeros(V);
Dhat(logical(eye(V))) = 1./sqrt(beta);
Omega = Dhat*Rhat*Dhat;

    %% ====================================================================
    % FUNCTION calcDeltaL calculates change in marginal likelihood
    % ****************************** Input ********************************
    % mask:     scalar, ~isinf(alpha)
    % alpha:    scalar, old value of alpha
    % alphaNew: scalar, new value of alpha
    % sPrime:   scalar
    % qPrime:   scalar
    % s:        scalar
    % q:        scalar
    % ****************************** Output *******************************
    % task:     string, one of 'est', 'add', 'del', and 'non'
    % deltaL: 	deltaL, change in marginal likelihood
    % =====================================================================
    function [task, deltaL] = calcDeltaL(mask, alpha, alphaNew, ...
            sPrime, qPrime, s, q)
        
        if ~isinf(alphaNew) % if i should be in model
            if mask % ~isinf(alpha(i)) % if i was already in model
                % reestimate alpha_i
                task = 'est';
                tempCalc = 1./alphaNew - 1./alpha;
                qPrimeSq = qPrime.^2;
                deltaL = sum(qPrimeSq./(sPrime + 1./tempCalc) ...
                    - log(1 + sPrime.*tempCalc));
            else
                % add i
                task = 'add';
                tempCalc = alphaNew + s;
                deltaL = sum((q.^2)./tempCalc + log(alphaNew./tempCalc));
            end
        elseif mask % ~isinf(alpha(i)) % if i was already in model
            % delete i
            task = 'del';
            qPrimeSq = qPrime.^2;
            deltaL = sum(qPrimeSq./(sPrime - alpha) ...
                - log(1 - sPrime./alpha));
        else
            task = 'non';
            deltaL = nan;
        end
        
        if imag(deltaL)
            % deltaL has a imaginary number if inside of log is negative,
            %   but it is numerical error. Therefore, replace it with -inf
            deltaL = -inf;
            warning('deltaL has an imaginary number')
        end
    end
end