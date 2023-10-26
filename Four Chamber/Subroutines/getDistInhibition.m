function [f_orco,time2State,exitTrajectory,timeOfState,kinOfState,biasOfState,dFSmooth,FSmooth] = ...
    getDistInhibition(f_orco,crossingtype,plotFig)
%close all;%clear;
if isempty(f_orco)
    load('DataModel/Orco Retinal_10Dec202060sAvg_32_25_3.mat', 'f_orco');
end
spd = f_orco.calcSpd;

% calculate smoothed change in firing rate
dF = f_orco.calcDeltaFR;%gradient(f_orco.spk);
dF(:,end-1:end) = repmat(dF(:,end-2),1,2);
hist = ceil(0.2.*f_orco.fs)-1;% 200 ms period
dFSmooth = dF(:,1:end-hist);
FSmooth = f_orco.spk(:,1:end-hist);
for i = 2:hist+1
    dFSmooth = dF(:,i:end-(hist+1)+i)+dFSmooth;
    FSmooth = f_orco.spk(:,i:end-(hist+1)+i)+FSmooth;
end
dFSmooth = [zeros(size(dFSmooth,1),hist),dFSmooth]./(hist+1);% average
FSmooth = [zeros(size(FSmooth,1),hist),FSmooth]./(hist+1);% average

% get tracks where the firing rate doesn't change and is either during
% inhibition or at inside steady state
InsideSS = 18.55;
thresh = 50;
exitTrajectory = cell(f_orco.nFly,1);
spkHist = [];%zeros(1,121);
dFHist = [];%zeros(1,121);
for fly = 1:f_orco.nFly
    if strcmpi(crossingtype,'enter')
        % inside steady state
        stableFdF = (abs(FSmooth(fly,:) - InsideSS)<1) & (abs(dFSmooth(fly,:)) <1);
    elseif strcmpi(crossingtype,'exit')
        % inhibition
        stableFdF = (FSmooth(fly,:) == 0) & (abs(dFSmooth(fly,:)) <1);
    end
    stableFdF(1:f_orco.fs) = false;
    [startNdx,endNdx,type] = startEndSeq(stableFdF);
    exitTrajectory{fly,1} = [startNdx(type==true)';endNdx(type==true)']';
    
    % get the firing rate and change in firing rate for 2 s before entering
    % the track to 2 s after entering the track (for plotting later)
    tmpSpk = FSmooth(fly,:);
    tmpdSpk = dFSmooth(fly,:);
    if sum(type==true)>0
        tmp = startNdx(type==true)'+[-60:60];
        badNdx = tmp>10800;
        
        tmpSpk = tmpSpk(min(tmp,10800));
        tmpdSpk = tmpdSpk(min(tmp,10800));
        tmpSpk(badNdx) = nan;
        tmpdSpk(badNdx) = nan;
        
        spkHist = [spkHist;tmpSpk];
        dFHist = [dFHist;tmpdSpk];
    end
end
% get the average firing rate and change in firing rate for +/-2 seconds
% relative to start of the tracks
spkHist = nanmean(spkHist);
dFHist = nanmean(dFHist);

%%

distType = {'logN','logN','beta','beta','logN','exponential'};
zGrid = f_orco.model.params{1, 1}.KNN.Time;
dz = f_orco.model.params{1, 1}.KNN.ratio(3);

if strcmpi(crossingtype,'enter')
    sce = 'insideSSKin';
elseif strcmpi(crossingtype,'exit')
    sce = 'inhibitionKin';
end
fprintf('Generating %s distributions\n',sce)

time2State = cell(1,4);
timeOfState = cell(1,4);
for params = 1:6
    m = f_orco.model.params{1,params};
    if params<3
        MaxVal = 25;% max average speed is 25 mm/s
    elseif params<5
        MaxVal = 360;% maximum average curvature is 360 degrees
    else
        MaxVal = 5;% maximum average duration is 5 seconds
    end
    
    % get tracks where the firing rate doesn't change and is either during
    % inhibition or at inside steady state
    allTime = m.allTime;
    allFly = m.allFly;
    if strcmpi(crossingtype,'enter')
        % inside steady state
        stableFdF = (abs(m.allSpk - InsideSS)<1) & (abs(m.alldSpk) <1);
    elseif strcmpi(crossingtype,'exit')
        % inhibition
        stableFdF = (m.allSpk == 0) & (abs(m.alldSpk) <1);
    end
    
    % get variables
    allTimeFE = m.allTimeFE(stableFdF);
    allTime = allTime(stableFdF);
    allFly = allFly(stableFdF);
    inhibitionCurv = m.allDat(stableFdF);
    turnBias = m.turnBias(stableFdF);
    timeOfState{params} = [m.allFly(stableFdF), m.allTime(stableFdF)-1];
    kinOfState{params} = inhibitionCurv;
    biasOfState{params} = sign(turnBias-0.5);
    
    % get the time since start of inside steady state or inhibition for
    % each locomotor track
    t_sinceStableFdF = [];
    for fly = 1:f_orco.nFly
        if strcmpi(crossingtype,'enter')
            flyStableFdF = (abs(FSmooth(fly,:) - InsideSS)<1) & (abs(dFSmooth(fly,:)) <1);
        elseif strcmpi(crossingtype,'exit')
            flyStableFdF = (FSmooth(fly,:) == 0) & (abs(dFSmooth(fly,:)) <1);
        end
        flyStableFdF(1:f_orco.fs) = false;
        [startNdx,endNdx,type] = startEndSeq(flyStableFdF);
        startNdx(type==0) = [];
        endNdx(type==0) = [];
        if ~isempty(startNdx)
            [closestPtBef,~,idxB] = findBeforeAfter(allTime(allFly==fly),startNdx-1,'before');
            t_sinceStableFdF = [t_sinceStableFdF;allTime(allFly==fly)-closestPtBef-1];
        end
    end
    t_sinceStableFdF = t_sinceStableFdF./f_orco.fs;
    
    % set up parameter and turn bias fields to save in the fly object
    f_orco.model.params{1,params}.(sce).a = cell(1,numel(zGrid));
    f_orco.model.params{1,params}.(sce).b = cell(1,numel(zGrid));
    f_orco.model.params{1,params}.(sce).TurnBias = cell(1,numel(zGrid));
    if params == 1 || params == 3
        f_orco.model.TurnBias.(sce){1,1} = cell(1,numel(zGrid));
    elseif params == 2 || params == 4
        f_orco.model.TurnBias.(sce){1,2} = cell(1,numel(zGrid));
    end
    
    % loop through each time slice since first entry
    for slice = 1:numel(zGrid)
        % get the locomotor tracks within the time window of the time slice
        currInhCurv = inhibitionCurv((allTimeFE>zGrid(slice)-dz) & (allTimeFE<zGrid(slice)+dz));
        currTurnBias = turnBias((allTimeFE>zGrid(slice)-dz) & (allTimeFE<zGrid(slice)+dz));
        currT_sinceStableFdF = t_sinceStableFdF((allTimeFE>zGrid(slice)-dz) & (allTimeFE<zGrid(slice)+dz));
        
        % get the kinematics and fit to model
        if strcmpi(distType(params),'beta')
            xx = [0:MaxVal/72:MaxVal]./MaxVal;% beta is bounded by 0 and 1
        else
            xx = [0:MaxVal/72:MaxVal];
        end
        
        % set up sliding window for time since start of steady
        % state/inhibition
        window = 0.5;%seconds
        nOverlap = 0.3;%seconds
        tStart = [-0.5:(window-nOverlap):5-window];
        tEnd = tStart+window;
        Y = zeros(numel(xx),numel(tStart));
        if strcmpi(distType(params),'exponential')
            phatAll = nan(numel(tStart),1);% exponential has one parameter
        else
            phatAll = nan(numel(tStart),2);% logN and beta has 2 parameters
        end
        mu = nan(numel(tStart),1);
        avgBias = nan(numel(tStart),1);
        % loop through each window and fit to distribution
        for t = 1:numel(tStart)
            currCurv = currInhCurv(currT_sinceStableFdF>=tStart(t) & currT_sinceStableFdF<tEnd(t));
            currBias = currTurnBias(currT_sinceStableFdF>=tStart(t) & currT_sinceStableFdF<tEnd(t));
            minCurv = min(currInhCurv).*30;
            % only fit to distribution if there are over 10 tracks
            if numel(currCurv)>10
                if strcmpi(distType(params),'logN')
                    phat = lognfit(currCurv);
                    Y(:,t) = lognpdf(xx,phat(1),phat(2));
                    phatAll(t,:) = phat;
                end
                if strcmpi(distType(params),'beta')
                    phat = betafit(currCurv./MaxVal);
                    Y(:,t) = betapdf(xx,phat(1),phat(2));
                    phatAll(t,:) = phat;
                end
                if strcmpi(distType(params),'exponential')
                    phat = expfit(currCurv-minCurv./f_orco.fs);
                    Y(:,t) = exppdf(xx-minCurv./f_orco.fs,phat(1));
                    phatAll(t,:) = phat;
                end
                mu(t) = mean(currCurv);
                avgBias(t) = sum(currBias)./numel(currBias);
            end
        end
        
        % if over 3 time windows had enough data to fit distributions to
        if sum(~isnan(phatAll))>3
            
            % Use a spline fit on the fitted distributions to model how the
            % distribution of the kinematic evolves across time
            tt = (tStart+tEnd)./2;tt(isnan(phatAll(:,1))) = [];
            order = 1;
            s = cell(1,size(phatAll,2));
            for nParam = 1:size(phatAll,2)
                s{nParam} = spline(tt,phatAll(~isnan(phatAll(:,nParam)),nParam));
                s{nParam} = fnxtr(s{nParam},order);
            end
            % use spline fit on the turn optimality to model how the
            % optimality evolves across time
            s3 = spline(tt,avgBias(~isnan(phatAll(:,1)),1));
            s3 = fnxtr(s3,order);
            
            % add model parameters to the fly object (time slice parameters)
            f_orco.model.params{1,params}.(sce).tt = zGrid./f_orco.fs;
            f_orco.model.params{1,params}.(sce).dt = dz;
            
            % add model parameters to the fly object (kinematics spline fits)
            f_orco.model.params{1,params}.(sce).a{1,slice} = @(t) ppval(s{1},t);
            if ~strcmpi(distType(params),'exponential')
                f_orco.model.params{1,params}.(sce).b{1,slice} = @(t) ppval(s{2},t);
            else
                f_orco.model.params{1,params}.(sce).b{1,slice} = @(t) ppval(s{1},t);
            end
            if ~strcmpi(distType(params),'beta')
                f_orco.model.params{1,params}.(sce).maxVal = 1;
            else
                f_orco.model.params{1,params}.(sce).maxVal = MaxVal;
            end
            f_orco.model.params{1,params}.(sce).shiftVal = minCurv;
            f_orco.model.params{1,params}.(sce).fitFun = distType(params);
            f_orco.model.params{1,params}.(sce).TurnBias{1,slice} = @(t) ppval(s3,t);
            % add model parameters to the fly object (optimality spline fits)
            if strcmpi(m.state.state,'sharp turns')
                f_orco.model.TurnBias.(sce){1,1}{1,slice} = @(t) ppval(s3,t);
            elseif strcmpi(m.state.state,'curved walks')
                f_orco.model.TurnBias.(sce){1,2}{1,slice} = @(t) ppval(s3,t);
            elseif strcmpi(m.state.state,'stops')
                f_orco.model.TurnBias.(sce){1,3}{1,slice} = @(t) ppval(s3,t);
            else
                disp('Incorrect movement state, must be "sharp turn", "curved walks", or "stops"')
            end
            % plot analysis (would not recommend when considering a lot of
            % time slices because you will get hundreds of figures)
            if plotFig
                % generate pdf heatmaps and sampe from the pdf (used for
                % plotting only)
                x = [0:0.01:3];
                Y2 = zeros(numel(xx),numel(x));
                sampTime = rand(10000,1).*3;
                if strcmpi(distType(params),'logN')
                    for t = 1:numel(x)
                        Y2(:,t) = lognpdf(xx,ppval(s{1},x(t)),ppval(s{2},x(t)));
                    end
                    sampCurv = lognrnd(ppval(s{1},sampTime),ppval(s{2},sampTime));
                end
                if strcmpi(distType(params),'beta')
                    for t = 1:numel(x)
                        Y2(:,t) = betapdf(xx,ppval(s{1},x(t)),ppval(s{2},x(t)));
                    end
                    sampCurv = betarnd(ppval(s{1},sampTime),ppval(s{2},sampTime)).*MaxVal;
                end
                if strcmpi(distType(params),'exponential')
                    for t = 1:numel(x)
                        Y2(:,t) = exppdf(xx-minCurv./f_orco.fs,ppval(s{1},x(t)));
                    end
                    sampCurv = exprnd(ppval(s{1},sampTime))+minCurv./f_orco.fs;
                end
                muSamp = nan(numel(tStart),1);
                for t = 1:numel(tStart)
                    if strcmpi(distType(params),'beta')
                        currCurv = sampCurv(sampTime>=tStart(t) & sampTime<tEnd(t))./MaxVal;
                    else
                        currCurv = sampCurv(sampTime>=tStart(t) & sampTime<tEnd(t));
                    end
                    if numel(currCurv)>10
                        muSamp(t) = mean(currCurv);
                    end
                end
                % plot analysis
                plottingFunc(f_orco,s,tt,phatAll,xx,stableFdF,MaxVal,Y,Y2,sampCurv,mu,sampTime,x,m,params,muSamp,tStart,tEnd,currT_sinceStableFdF,currInhCurv,currTurnBias,avgBias,spkHist,dFHist,thresh,slice)
            end
        end
    end
    
    % set any time slices (since first entry) that do not have enough data
    % points to fit a spline fit to as equal to the closest time slice that
    % does have enough data
    emptyZGrid = cellfun(@(x) isempty(x),f_orco.model.params{1,params}.(sce).a);
    notEmptyZGrid = find(~emptyZGrid);
    emptyZGrid = find(emptyZGrid);
    if ~isempty(emptyZGrid) && ~isempty(notEmptyZGrid)
        [closestPtBef,closestPtAft,idxB] = findBeforeAfter(emptyZGrid,notEmptyZGrid,'both');
        dToNonEmptyGrid = abs([closestPtBef-emptyZGrid; closestPtAft-emptyZGrid]);
        [~,ndx] = min(dToNonEmptyGrid);
        closestZGrid = zeros(1,size(idxB,1));
        for k = 1:size(idxB,1)
            closestZGrid(k) = idxB(k,ndx(k));
        end
        closestZGrid = notEmptyZGrid(closestZGrid);
        
        f_orco.model.params{1,params}.(sce).TurnBias(1,emptyZGrid) = ...
            f_orco.model.params{1,params}.(sce).TurnBias(1,closestZGrid);
        f_orco.model.params{1,params}.(sce).a(1,emptyZGrid) = ...
            f_orco.model.params{1,params}.(sce).a(1,closestZGrid);
        f_orco.model.params{1,params}.(sce).b(1,emptyZGrid) = ...
            f_orco.model.params{1,params}.(sce).b(1,closestZGrid);
        if params == 1 || params == 3
            f_orco.model.TurnBias.(sce){1,1}(1,emptyZGrid) = ...
                f_orco.model.TurnBias.(sce){1,1}(1,closestZGrid);
        elseif params == 2 || params == 4
            f_orco.model.TurnBias.(sce){1,2}(1,emptyZGrid) = ...
                f_orco.model.TurnBias.(sce){1,2}(1,closestZGrid);
        end
    end
end


end

function [] = plottingFunc(f_orco,s,tt,phatAll,xx,firstState,MaxVal,Y,Y2,sampCurv,mu,sampTime,x,m,params,muSamp,tStart,tEnd,firstStateTime,inhibitionCurv,turnBias,avgBias,spkHist,dFHist,thresh,slice)
m2 = f_orco.model.params{1,params}.inhibitionKin;
% plotting functions
figure;set(gcf,'Position',[2 42 838 924])
for nParam = 1:numel(s)
    subplot(3,2,nParam);
    plot(tt,phatAll(~isnan(phatAll(:,nParam)),nParam));hold on;
    plot([0:0.01:3],ppval(s{nParam},[0:0.01:3]));legend({'emp','spline'})
    xlim([-0.5 5]);
    xlabel('seconds');ylabel(['param ' num2str(nParam)])
end
subplot(3,2,3);
imagesc((tStart+tEnd)./2,xx.*m2.maxVal,Y);set(gca,'YDir','normal');hold on;
colormap(flipud(hot));
title('Empirical fit')
subplot(3,2,4);
imagesc(x,xx.*m2.maxVal,Y2);set(gca,'YDir','normal');hold on;
colormap(flipud(hot));
title('Spline fit')
subplot(3,2,5);
scatter(sampTime(1:1000),sampCurv(1:1000));title('1000 samples')
subplot(3,2,6);
plot((tStart+tEnd)./2,mu);hold on;plot((tStart+tEnd)./2,muSamp.*m2.maxVal)
legend({'emp mean','sampled mean'})
for sbplt = 3:6
    subplot(3,2,sbplt);xlim([-0.5 5]);ylim([0 MaxVal]);
    xlabel('seconds');ylabel(m.state.kin);
end
suptitle(['After peak in df ' m.state.state]);


tt = (-60:60);
figure;set(gcf,'Position',[2 42 838 924])
subplot(3,1,1);
yyaxis left;imagesc((tStart+tEnd)./2,xx.*m2.maxVal,Y);set(gca,'YDir','normal');hold on;
colormap(flipud(hot))
scatter(firstStateTime,inhibitionCurv,'MarkerEdgeAlpha',.5);
plot((tStart+tEnd)./2,mu,'-c','LineWidth',2);
xlabel('seconds');ylabel(m.state.kin);
ylim([0 MaxVal]);xlim([-2 5])
yyaxis right;plot(tt./f_orco.fs,dFHist,'k','Linewidth',2);%hold on;
plot([-2 5],-[thresh thresh],'k')
xlabel('seconds');ylabel('average delta firing rate');
ylim([-150 150])
subplot(3,1,2);
yyaxis left;hold on;
scatter(firstStateTime,turnBias.*0.9+rand(numel(turnBias),1).*0.1,'MarkerEdgeAlpha',.5);
plot((tStart+tEnd)./2,avgBias,'-c','LineWidth',2);
plot((0:0.1:3),f_orco.model.params{1,params}.inhibitionKin.TurnBias{slice}(0:0.1:3),'--k','LineWidth',2);
xlabel('seconds');ylabel('turn optimality');
ylim([0 1]);xlim([-2 5])
yyaxis right;plot(tt./f_orco.fs,spkHist,'k','Linewidth',2);hold on;
xlabel('seconds');ylabel('average firing rate');
subplot(3,1,3);
scatter(m.alldSpk,m.allSpk,[],[0.5 0.5 0.5]);hold on;
scatter(m.alldSpk(firstState),m.allSpk(firstState),'r*');
xlim([-180 180]);
xlabel('df spks/s/s');ylabel('f (spks/s)');
end