function [ws_loo,ys_loo,Az_loo,Az_perm_dist] = logist_fast(X,y,l2_lambda,mode,varargin)
% FUNCTION [ws_loo,ys_loo,Az_loo,Az_perm_dist] = logist_fast(X,y,l2_lambda,'loo'[, cv])
%	*** OR ***
% FUNCTION [ws_loo,ys_loo,Az_loo,Az_perm_dist] = logist_fast(X,y,l2_lambda,'loo_and_bootstrap',cv,yperms)
%
%	INPUTS:
% 		X is the feature data (matrix size DxN -- D = # features, N = # trials)
% 		y are the binary (0/1) truth labels (vector size Nx1)
% 		l2_lambda is the l2 regularization value (e.g., 1e-6)
%		mode can be either 'loo' or 'loo_and_bootstrap'
%		If mode == 'loo', then run cross-validation only (no permutation testing)
%			the folds can be specified by an optional cv structure (otherwise assumed to be LOO)
%			the cv structure has the following fields:
%				cv.numFolds: the number of cross-validation folds
%				cv.incTrials: a (cv.numFolds x 1) cell array.
%					cv.incTrials{i} contains the indices of the trials to be included in the ith fold
%				cv.outTrials: a (cv.numFolds x 1) cell array.
%					cv.outTrials{i} contains the indices of the trials that are held out in the ith fold
%		If mode == 'loo_and_bootstrap', then permutation testing will also be run
%			the cv structure specifies the cross-validation folds (see above)
%			yperms is an N x numperms matrix that specifies the label permutations for the (numperms) permutations
%
%	OUTPUTS:
%		ws_loo: (D+1) x cv.numFolds matrix, giving the LR weights for each cv fold
%		ys_loo: N x 1 vector, giving the 
%		Az_loo: scalar:  the area-under-ROC curve
%		If mode=='loo_and_bootstrap' then the Az permutation distribution is also returned in Az_perm_dist:
%			Az_perm_dist: numperms x 1 vector giving the Az values for each of numperms permutations
%
%	EXAMPLE USAGE:
%		D = 50; % Number of variables
%		N = 100; % Number of trials
%		X = rand(D,N);
%		y = round(rand(N,1));
%		l2_lambda = 1e-1;
%		mode='loo_and_bootstrap';
%		% Setup the cross-validation structure (LOO)
%			cv = [];
%			cv.numFolds = N; % LOO cross-validation
%			cv.incTrials = cell(cv.numFolds,1);
%			cv.outTrials = cell(cv.numFolds,1);
%			for j=1:cv.numFolds; cv.outTrials{j} = j; cv.incTrials{j} = setdiff(1:N,j); end;
%		% Setup the label permutations
%			numperms = 100; % The number of permutations to run
%			yperms = zeros(N,numperms);
%			for j=1:numperms; yperms(:,j) = y(randperm(N)); end;
%		% Now call the function
%		[ws_loo,ys_loo,Az_loo,Az_perm_dist] = logist_fast(X,y,l2_lambda,mode,cv,yperms);
%
% Reference:
% @article{conroy2012,
% author = {Bryan R Conroy and Paul Sajda},
% title = {Fast, Exact Model Selection and Permutation Testing
%			for l2-Regularized Logistic Regression},
% conference = {AISTATS},
% year = {2012}
% }
%
% Author: Bryan Conroy (bc2468@columbia.edu, 2012)
%
%123456789012345678901234567890123456789012345678901234567890123456789012
%
% Copyright (C) 2012 Bryan Conroy
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

X = double(X);
y = double(y);
l2_lambda = double(l2_lambda);
conv_tol = 1e-11;
chunkSize = 5000;

looMode = 'average';

% Add all one's feature -- to estimate the bias
[D,N] = size(X);
X = [X;ones(1,N)];
D = D+1;

numBatches = 1;
if strcmp(mode,'loo')
    if isempty(varargin)
        cv = [];
        cv.numFolds = N;
        cv.outTrials = cell(1,cv.numFolds);
        for p=1:cv.numFolds
            cv.outTrials{p} = p;
        end
    else
        cv = varargin{1};
    end
    
    yperms = assertVec(y,'col');    
elseif strcmp(mode,'loo_and_bootstrap')
    cv = varargin{1};
    yperms = [assertVec(y,'col'),varargin{2}];
else    
    error('Unknown mode');
end
    
numperms = size(yperms,2);
% Make chunkSize the next closest multiple of cv.numFolds
chunkSize = cv.numFolds*max(1,round(chunkSize/cv.numFolds));
P = cv.numFolds*numperms;
if P > chunkSize
    numBatches = ceil(P/chunkSize);
    P = chunkSize;
    disp(['There are too many permutations to run in one batch.  This analysis will be run in ',num2str(numBatches), ' batches.']);pause(1e-9);
end
    
[y,zeroOutInds,perms] = getCurrentYAndZeroInds(1,yperms,cv,chunkSize);

ws_loo = zeros(D,cv.numFolds);
if strcmp(looMode,'average')
    ys_loo = zeros(cv.numFolds,1);
else
    ys_loo = zeros(N,1);
end
Az_perm_dist = zeros(numperms-1,1);

% Set the regularization matrix L
% Set the last diagonal element to zero to not penalize the bias
L = diag([ones(D-1,1);0]);

[Q,Z] = qr(X,0);
if size(Q,2) < D
    % Then we're better off using a low-rank representation for X
      
    % Now we need to find a basis for components y for which L*y is in the
    % range of U.  Also the component y must be orthogonal to U (U'*y=0)
    [Qperp,~] = qr(L'*Q,0);
    Qperp = (eye(size(Q,1))-Q*Q')*Qperp;
    [Qtmp,~] = qr([Q,Qperp],0);
    Q = Qtmp(:,1:size(Q,2));
    Qperp = Qtmp(:,(size(Q,2)+1):end);
    
    B = Q'*L*Qperp;
    D = Qperp'*L*Qperp;
    DinvBT = D\(B'); % DinvBT = D^(-1)*B';
    
    additiveMatrixFactor = 2*l2_lambda*(Q'*L*Q - B*DinvBT);
    
    D = size(Q,2);
    finalInvertMat = Q - Qperp*DinvBT;
    
    error('Didn''t reset Z here!');
else
    additiveMatrixFactor = 2*l2_lambda*L;
    finalInvertMat = eye(D);
    Z = X;
end

XT = X';

ws = zeros(D,P);
% Pre-allocate some of the big matrices
ws_last = zeros(D,P);
wsprev = ws;
onesMat = ones(N,P);
RTREP = zeros(N,P);
RTMAT = sparse(N,N);
RT_minus_R = zeros(N,P);
yvals = zeros(N,P);
mus = zeros(N,P);
R = zeros(N,P);
Minv_ls = zeros(D,P);

% OK now let's get started with the iterations
iter = 0;
locsNotFullyConverged = 1:P;
batchNum = 1;
while 1
    iter = iter + 1;
    disp(['Working on batch #',num2str(batchNum),' out of ',num2str(numBatches),'; Iteration #',num2str(iter),'...']);pause(1e-9);

    ws_last(:,:) = ws;

    yvals(:,locsNotFullyConverged) = XT*ws(:,locsNotFullyConverged);
	mus(:,locsNotFullyConverged) = 1./(1+exp(-yvals(:,locsNotFullyConverged)));

    R(:,locsNotFullyConverged) = mus(:,locsNotFullyConverged).*(1-mus(:,locsNotFullyConverged));
    R(zeroOutInds) = 0;
                
	RT = max(R(:,locsNotFullyConverged),[],2);%%median(R(:,locsNotFullyConverged),2);
	if max(RT) == 0
	    break;
    end
        
	neg_ls_without_X(:,locsNotFullyConverged) = (R(:,locsNotFullyConverged).*yvals(:,locsNotFullyConverged) + y(:,locsNotFullyConverged) - mus(:,locsNotFullyConverged));
    neg_ls_without_X(zeroOutInds) = 0;
%    neg_ls = X*neg_ls_without_X;

    RTMAT = spdiags(RT,0,RTMAT);
    M = X*(RTMAT*XT) + additiveMatrixFactor;
    MinvX = M\X;
	if strncmp(lastwarn,'Matrix is close to singular or badly scaled.',44)
	    warning('Bad conditioning. Suggest reducing subspace.')
	    singularwarning=1;
	    break;
    end
                
	% Iteration is the following:
	% ws = Minv*N*wsprev + Minv*neg_ls
	Minv_ls(:,locsNotFullyConverged) = MinvX*neg_ls_without_X(:,locsNotFullyConverged);
    RTREP(:,locsNotFullyConverged) = RTMAT*onesMat(:,locsNotFullyConverged);
%    RTREP = repmat(RT,1,length(locsNotFullyConverged));
	RT_minus_R(:,locsNotFullyConverged) = RTREP(:,locsNotFullyConverged) - R(:,locsNotFullyConverged);        
        
    
    performLowRankUpdates = 0;
    ratio_thresh = 0.5;
    maxRank = 1;%round(0.05*D);

    hasLowRankCorrections = zeros(1,P);
    lowRankLocs = cell(1,P);
    lowRankCorrectionMatsLeft = cell(1,P);
    lowRankCorrectionMatsRight = cell(1,P);
    lowRankCorrectionSizes = zeros(1,P);
    if performLowRankUpdates == 1
        RT_ratio = zeros(N,P);    
        RT_ratio(:,locsNotFullyConverged) = abs(RT_minus_R(:,locsNotFullyConverged))./RTREP(:,locsNotFullyConverged);
        RT_ratio(zeroOutInds) = 0;
        %disp('Zeroing out RT_ratio(zeroOutInds) for now.  May not want to do this in the future to speed up the algorithm.');pause(1e-9);
        colinds = max(RT_ratio(:,locsNotFullyConverged)) > ratio_thresh;
        lowRankColIndsUnique = locsNotFullyConverged(colinds);    
        avgLowRankNum = 0;
        if ~isempty(lowRankColIndsUnique)
            XTMinvX = XT*MinvX;
            for p=lowRankColIndsUnique
                lowRankLocs{p} = find(RT_ratio(:,p) > ratio_thresh);
                if isempty(lowRankLocs{p}); 
                    continue; 
                end;
    
                if length(lowRankLocs{p}) > maxRank
                    % Then make sure that we get everywhere that RT_ratio == 1
                    numToKeep = max(maxRank,length(find(RT_ratio(lowRankLocs{p},p) == 1)));
                    [~,inds] = sort(RT_ratio(lowRankLocs{p},p),'descend');
                    lowRankLocs{p} = lowRankLocs{p}(inds(1:numToKeep));
                end

                lowRankCorrectionSizes(p) = length(lowRankLocs{p});
                avgLowRankNum = avgLowRankNum + length(lowRankCorrectionSizes(p));

                hasLowRankCorrections(p) = 1;
                lowRankCorrectionMatsLeft{p} = MinvX(:,lowRankLocs{p});
                lowRankCorrectionMatsRight{p} = (diag(1./RT_minus_R(lowRankLocs{p},p)) - XTMinvX(lowRankLocs{p},lowRankLocs{p}))\XT(lowRankLocs{p},:);%XT(lowRankLocs{p},:)*MinvX(:,lowRankLocs{p})

                RT_minus_R(lowRankLocs{p},p) = 0;
                RT_ratio(lowRankLocs{p},p) = 0;        
            end
        end
    else
        %disp('Not performing low rank updates.  May want to do this in the future to speed up the algorithm.');pause(1e-9);
    end
    hasLowRankCorrectionsCell = num2cell(hasLowRankCorrections);
%    avgLowRankNum/length(locsNotFullyConverged)
    
	locsNotConverged = locsNotFullyConverged;
	ninnerloops = 0;
	while 1
        ninnerloops = ninnerloops + 1;
        %fprintf(1,['\tPerforming inner loop #',num2str(ninnerloops),'...\n']);pause(1e-9);
        wsprev(:,:) = ws;

        %XTwsprev = XT*wsprev(:,locsNotConverged);
        ws(:,locsNotConverged) = MinvX*(RT_minus_R(:,locsNotConverged).*(XT*wsprev(:,locsNotConverged))) + Minv_ls(:,locsNotConverged);

        % Now do the low-rank corrections
        lowRankLocsCurr = locsNotConverged(find(hasLowRankCorrections(locsNotConverged)));
        if ~isempty(lowRankLocsCurr)
            wsc = mat2cell(ws(:,lowRankLocsCurr),D,ones(1,length(lowRankLocsCurr)));
            ws(:,lowRankLocsCurr) = cell2mat(cellfun(@lowRankCorrection,wsc,lowRankCorrectionMatsLeft(lowRankLocsCurr),lowRankCorrectionMatsRight(lowRankLocsCurr),hasLowRankCorrectionsCell(lowRankLocsCurr),'UniformOutput',false));
        end

        tol_measure = (1/D)*sum((ws(:,locsNotConverged)-wsprev(:,locsNotConverged)).^2)./sum(ws(:,locsNotConverged).^2);
        locsNew = find(tol_measure > conv_tol);
        if isempty(locsNew)
            break;
        end
        locsNotConverged = locsNotConverged(locsNew);
    end
    fprintf(1,'\tNumber of inner loops: %d\n',ninnerloops);pause(1e-9);
        
    locsNotFullyConverged = find((1/D)*sum((ws_last-ws).^2)./sum(ws_last.^2) > conv_tol);
    
    if isempty(locsNotFullyConverged)
        % Save out the converged results
        eInd = 0;
        while 1
            sInd = eInd + 1;
            if sInd > length(perms)
            	disp(['For debugging purposes:  sInd=',num2str(sInd),'; length(perms)=',num2str(length(perms))]);pause(1e-9);
                break
            end
            eInd = sInd + cv.numFolds - 1;
            currInds = sInd:eInd;
            if max(abs(perms(currInds)-mean(perms(currInds))))
                error('Something''s very wrong here!');
            end
            permInd = perms(sInd);

            if strcmp(looMode,'average')
                ys_loo_perm = zeros(cv.numFolds,1);
                truth = zeros(cv.numFolds,1);
            else
                ys_loo_perm = zeros(N,1);
                truth = zeros(N,1);
            end
                
            for foldNum = 1:cv.numFolds
                if strcmp(looMode,'average')
                    ys_loo_perm(foldNum) = mean(XT(cv.outTrials{foldNum},:)*ws(:,currInds(foldNum)));
                    truth(foldNum) = yperms(cv.outTrials{foldNum}(1),permInd);
                else
                    ys_loo_perm(cv.outTrials{foldNum}) = XT(cv.outTrials{foldNum},:)*ws(:,currInds(foldNum));
                    truth(cv.outTrials{foldNum}) = yperms(cv.outTrials{foldNum},permInd);
                end
            end
            % Now compute the roc value and store it
            Az_val = rocarea(ys_loo_perm,truth);
            
            if permInd == 1
                % Then this is the LOO run -- save out the ws
                ws_loo(:,:) = ws(:,currInds);
                ys_loo = ys_loo_perm;
                Az_loo = Az_val;
            else
                Az_perm_dist(permInd-1) = Az_val;
            end
        end
        
        if batchNum < numBatches        
            % We need to re-fill
            iter = 0;
            batchNum = batchNum + 1;
            ws = zeros(D,P);
            [y,zeroOutInds,perms,locsNotFullyConverged] = getCurrentYAndZeroInds(batchNum,yperms,cv,chunkSize);
        else
            break;
        end
    end
end
ws = finalInvertMat*ws;

end


function [y,zeroOutInds,perms,locsNotFullyConverged] = getCurrentYAndZeroInds(batchNum,yperms,cv,chunkSize)
    [N,numperms] = size(yperms);    
    P = min(chunkSize,cv.numFolds*numperms);
    
    numpermsPerBatch = P/cv.numFolds;
    if abs(numpermsPerBatch - round(numpermsPerBatch)) > 0
        error('Something''s not right');
    end
    % Now determine the set of permutations to run based on batchNum
    permsStart = (batchNum-1)*numpermsPerBatch + 1;
    permsEnd = min(numperms, permsStart + numpermsPerBatch - 1);
    permsCurr = permsStart:permsEnd;
    numpermsCurr = length(permsCurr);
    
    perms = reshape(repmat(assertVec(permsCurr,'row'),cv.numFolds,1),numpermsCurr*cv.numFolds,1);
    
    y = zeros(N,P);
    eInd = 0;
    for n=permsCurr
        sInd = eInd + 1;
        eInd = sInd + cv.numFolds - 1;
        y(:,sInd:eInd) = repmat(yperms(:,n),1,cv.numFolds);
    end
    
    locsNotFullyConverged = 1:eInd;
    
    totalOutTrials = 0;
    for p=1:cv.numFolds
        totalOutTrials = totalOutTrials + length(cv.outTrials{p});
        cv.outTrials{p} = assertVec(cv.outTrials{p},'row');
    end
    
    zeroOutInds = zeros(1,totalOutTrials*numpermsCurr,'uint32');
    eInd = 0;
    % Compute the indices to be zeroed out in the NxP gradient matrix
    for foldNum=1:cv.numFolds
        for n=1:numpermsCurr
            sInd = eInd + 1;
            eInd = sInd + length(cv.outTrials{foldNum}) - 1;
            zeroOutInds(sInd:eInd) = cv.outTrials{foldNum} + N*(foldNum+(n-1)*cv.numFolds - 1);%sub2ind([N,P],assertVec(cv.outTrials{foldNum},'row'),(foldNum+(n-1)*cv.numFolds)*ones(1,length(cv.outTrials{foldNum})));
        end
    end
end

