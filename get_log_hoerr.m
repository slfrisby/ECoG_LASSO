function o = get_log_hoerr(varargin)

    % by Saskia (based on ger_hoerr.r, by Tim).
    % addpath('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/glmnet-matlab-master');
    addpath(genpath('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/WISC_MVPA'));
    
    % Inputs:
    % - X: X matrix of data (items x features).
    % - Y: labels vector (items x 1). E.g. Y = repmat([zeros(50,1);ones(50,1)],4,1);
    % - cvind: cross-validation index (used to fix which stimuli go into
    % the training set and which go into the test set - this ensures that
    % results replicate exactly). 
    % - ho: number of held-out items from each category (note that there are 2
    % categories so the actual holdout set is twice this size). Default: 5
    % a - alpha (1 for LASSO), Default: 1
    % - returnModels: whether to output fitted model. Default: True
    % - acrossRun: draw the test set from one fMRI run and the training set
    % from the other runs (stimuli that appear in the training set do not
    % appear in the test set). Default: True
    % - permstruct: If set, conducts a permutation test. Matrix dictating 
    % the order in which to scramble rows of data for permutation testing.
    % CURRENTLY IMPLEMENTED ONLY FOR ACROSSRUN = 0.
    p = inputParser;
    addRequired(p,'X');
    addRequired(p,'Y');
    addParameter(p,'cvind',[])
    addParameter(p,'ho',5);
    addParameter(p,'a',1);
    addParameter(p,'returnModels',1);
    addParameter(p,'acrossRun',1);
    addParameter(p,'permstruct',[])
    parse(p,varargin{:});
    
    X = p.Results.X;
    Y = p.Results.Y;
    cvind = p.Results.cvind;
    ho = p.Results.ho;
    a = p.Results.a;
    returnModels = p.Results.returnModels;
    acrossRun = p.Results.acrossRun;
    permstruct = p.Results.permstruct;
    
    % set parameters
    nitems = size(X,1);
    nreps = nitems/100;
    % find number of cross-validation folds
    if acrossRun
        nfolds = nitems/(ho*2);
        % keep track of which run the holdout set should be drawn from on
        % each fold
        runidx = [];
        for i = 1:nreps
            runidx = [runidx; i*ones((ho*2),1)];
        end
    elseif ~acrossRun
        nfolds = nitems/(ho*2*nreps);
    end
    % keep track of which folds should have the same holdout sets (this
    % does not affect within-run decoding; for across-run decoding, the
    % same holdout sets are repeated for each run)
    foldidx = repmat((1:(ho*2))',nreps,1);
    
    Y = Y(1:nitems);
    
    % configure outputs
    output = zeros(nfolds,1);
    if returnModels
        mo = {};
    end
    predictedcoords = zeros(nitems,1);
    if ~isempty(permstruct)
        permoutput = zeros(nfolds,100);
        permpredictedcoords = zeros(nitems,100);
    end
    
    % if cvind is not set, scramble order of living and nonliving. This step makes a vector of the
    % numbers 1:100, but 1-50 stay in the first half and 51-100 stay in the
    % second half
    if isempty(cvind)
        ford = [randperm(50),50+randperm(50)];
    end
    
    % loop through folds
    for i = 1:nfolds
        
        % if cvind is set manually 
        if ~isempty(cvind)
            s = cvind(1:nitems/nreps) ~= i;
        elseif isempty(cvind)
            % make a vector of 'true's
            s = logical(ones(nitems/nreps,1));
            % find the indices of the current holdout set.             
            % For the first fold this
            % is the 1st to 5th and 51st to 55th stimuli in ford, etc. 
            hwin = ((foldidx(i)-1)*ho + 1):(foldidx(i)*ho);
            hwin = [hwin,hwin+50];
            % find the indices of these stimuli within s and replace them with
            % 'false'
            hoind = ford(hwin);
            s(hoind) = 0;
        end
    
        if acrossRun
            % separate training and test data
            currentRun = logical(ones(size(X,1),1));
            currentRun(100*(runidx(i)-1)+1:100*runidx(i)) = 0;
            xtrn = X(currentRun,:);
            ytrn = Y(currentRun);
            xtst = X(~currentRun,:);
            ytst = Y(~currentRun);
            % make train and test versions of s
            strn = repmat(s,nreps-1,1);
            stst = s;
            % keep only the stimuli marked true in the training set
            xtrn = xtrn(strn,:);
            ytrn = ytrn(strn);
            % keep only the stimuli marked false in the test set
            xtst = xtst(~stst,:);
            ytst = ytst(~stst);   
        elseif ~acrossRun
            % put the stimuli marked as true into the training set
            xtrn = X(s,:);
            ytrn = Y(s);
            % put the stimuli marked as false into the test set
            xtst = X(~s,:);
            ytst = Y(~s);
        end
    
        % set options structure for training. This includes setting the alpha
        % parameter to 1 (LASSO; default) and lambda to the default range
        % used by WISC_MVPA for SOSLASSO.
        options = glmnetSet;
        options.lambda = linspace(3,0.2);

        % train the model on the training set
        m = cvglmnet(xtrn,ytrn,'binomial',options,'deviance',9);
    
        % test the model on the test set
        pred = cvglmnetPredict(m,xtst,'lambda_min');
        % get predicted coordinates
        if acrossRun
            coordidx = currentRun;
            coordidx(~coordidx) = s;
            predictedcoords(~coordidx) = pred;
        elseif ~acrossRun
            predictedcoords(~s) = pred;
        end

        % calculate the error under the ROC curve
        [~,~,~,AUC] = perfcurve(ytst,pred,1);
        % store model
        mo{i} = m;
        % save and print the AUC
        output(i) = AUC;
        disp(['Fold ',num2str(i),' hold-out AUC: ',num2str(AUC)]);

        % if permutation test
        if ~isempty(permstruct)
            % for perms 1:100
            for perm = 1:100
                % scramble the data
                Xp = X(permstruct(:,perm),:);
                % put the stimuli marked as true into the training set
                xptrn = Xp(s,:);
                ytrn = Y(s);
                % put the stimuli marked as false into the test set
                xptst = Xp(~s,:);
                ytst = Y(~s);

                % train perm model 
                mp = cvglmnet(xptrn,ytrn,'binomial',options,'deviance',9);
                % assess model using the value of lambda selected by the
                % REAL model
                permpred = cvglmnetPredict(mp,xptst,m.lambda_min);
                permpredictedcoords(~s,perm) = permpred;
                % calculate the error under the ROC curve
                [~,~,~,AUC] = perfcurve(ytst,permpred,1);
                permoutput(i,perm) = AUC;
            end
        end
    end

    % prepare output
    o.output = output;
    o.predictedcoords = predictedcoords;
    if returnModels
        o.models = mo;
    end
    if ~isempty(permstruct)
        o.permoutput = permoutput;
        o.permpredictedcoords = permpredictedcoords;
    end
end


