function [f_orco] = getTransitionRate2(f_orco,plotFig)

n = -5;%-2
deltaf = 1;
deltaDf = 0.1;
dSpk = f_orco.calcDeltaFR;

% xGrid = [0:deltaf:50];
% yGrid = [-5:deltaDf:5];
xGrid = (-5:deltaDf:5).*f_orco.fs;
yGrid = 0:deltaf:50;
[XX,YY] = meshgrid(xGrid,yGrid);

f = [-deltaf./2:deltaf:50+deltaf./2];
df = [-5-deltaDf./2:deltaDf:5+deltaDf./2].*f_orco.fs;
nState = max(f_orco.states.ndx(:));
allF = cell(nState);alldF = cell(nState);
fAtDec = cell(nState);dfAtDec = cell(nState);
baseLineStateTot = zeros(nState,1);baseLineStateTrans = zeros(nState,1);
for fly = 1:f_orco.nFly
    tmpFlySpk = f_orco.spk(fly,:);
    tmpFlyDSpk = dSpk(fly,:);
    %fe = 1;%find(tmpFlySpk~=tmpFlySpk(1) & tmpFlyDSpk~=0,1);
    fe = find(tmpFlySpk~=tmpFlySpk(1) & tmpFlyDSpk~=0,1);
    
    [startNdx,endNdx,type] = startEndSeq(f_orco.states.ndx(fly,:));
    
    baseLineTrans = endNdx(endNdx<fe);
    baseLineTransType = type(endNdx<fe);
    startNdx(endNdx<fe) = [];
    type(endNdx<fe) = [];
    endNdx(endNdx<fe) = [];
    prevState = type(1:end-1);
    nextState = type(2:end);
    
    stateTot = f_orco.states.ndx(fly,fe:numel(tmpFlySpk));
    
    fAtDecision = mean(tmpFlySpk(max(endNdx'+[n:0],1)),2);
    fTotal = mean(tmpFlySpk(max([fe:numel(tmpFlySpk)]'+[n:0],1)),2);
    
    dfAtDecision = mean(tmpFlyDSpk(max(endNdx'+[n:0],1)),2);
    dfTotal = mean(tmpFlyDSpk(max([fe:numel(tmpFlyDSpk)]'+[n:0],1)),2);
    
    for stateBef = 1:nState
        for stateAft = 1:nState
            currNdx = [(prevState==stateBef & nextState==stateAft) false];
            fAtDec{stateBef,stateAft} = [fAtDec{stateBef,stateAft}; fAtDecision(currNdx)];
            dfAtDec{stateBef,stateAft} = [dfAtDec{stateBef,stateAft}; dfAtDecision(currNdx)];
        end
        allF{stateBef} = [allF{stateBef}; fTotal(stateTot==stateBef)];
        alldF{stateBef} = [alldF{stateBef}; dfTotal(stateTot==stateBef)];
    end
    baseLineStateTrans = baseLineStateTrans+histcounts(baseLineTransType,[1:5])';
    baseLineStateTot = baseLineStateTot+histcounts(f_orco.states.ndx(fly,1:fe-1),[1:5])';
end

for stateBef = 1:nState
    DecisionPointsFDf3{stateBef,1} = histcounts2(cell2mat(fAtDec(stateBef,:)'),cell2mat(dfAtDec(stateBef,:)'),f,df);
    TotalPointsFDf3{stateBef,1} = histcounts2(cell2mat(allF(stateBef,:)'),cell2mat(alldF(stateBef,:)'),f,df);
    for stateAft = 1:nState
        DecisionPointsFDf2{stateBef,stateAft} = histcounts2(fAtDec{stateBef,stateAft},dfAtDec{stateBef,stateAft},f,df);
    end
end

for stateBef = 1:nState
    tmpF_trans = cell2mat(fAtDec(stateBef,:)');
    tmpdF_trans = cell2mat(dfAtDec(stateBef,:)');
    tmpF_all = cell2mat(allF(stateBef,:)');
    tmpdF_all = cell2mat(alldF(stateBef,:)');
    
    a = histcounts(tmpF_trans(abs(tmpdF_trans)<1),[0:2:50]);
    b = histcounts(tmpF_all(abs(tmpdF_all)<1),[0:2:50]);
    a(b<15) = nan;
    
    a2 = histcounts(tmpdF_trans(abs(tmpF_trans-4.701)<2),[-150:10:150]);
    b2 = histcounts(tmpdF_all(abs(tmpF_all-4.701)<2),[-150:10:150]);
    a2(b2<15) = nan;
    
    a2 = histcounts(tmpdF_trans(abs(tmpF_trans-4.701)<2),[flip(-2.^(0:0.5:7)), 2.^(0:0.5:7)]);
    b2 = histcounts(tmpdF_all(abs(tmpF_all-4.701)<2),[flip(-2.^(0:0.5:7)), 2.^(0:0.5:7)]);
end

x = [flip(-2.^(0:0.5:7)), 2.^(0:0.5:7)];
a2 = histcounts(tmpdF_trans(abs(tmpF_trans-4.701)<2),[flip(-2.^(0:0.5:7)), 2.^(0:0.5:7)]);
b2 = histcounts(tmpdF_all(abs(tmpF_all-4.701)<2),[flip(-2.^(0:0.5:7)), 2.^(0:0.5:7)]);

[flip(-2.^(0:2:9)), 2.^(0:2:9)]

tmpF_trans(tmpdF_trans==0);
tmpF_all(tmpdF_all==0);

figure;
a = histcounts(tmpF_trans(abs(tmpdF_trans)<1),[0:2:50]);
b = histcounts(tmpF_all(abs(tmpdF_all)<1),[0:2:50]);

tmp = DecisionPointsFDf3{2,1}./TotalPointsFDf3{2,1};
baselineRate = baseLineStateTrans./baseLineStateTot;

%tmp = squeeze(sum(DecisionPointsFDf,1)./sum(TotalPointsFDf,1));
% figure;surf(XX',YY',tmp)

i = 1;
a = DecisionPointsFDf3{i,1};
b = TotalPointsFDf3{i,1};
if plotFig
    plotFigs(f_orco,DecisionPointsFDf3,TotalPointsFDf3,XX,YY,n,baselineRate,plotFig);
end

% tmp = a./b;
% tmp(b<60) = nan;
end

function [] = plotFigs(f_orco,DecisionPointsFDf3,TotalPointsFDf3,XX,YY,n,baselineRate,plotFig)
xlims = [min(XX(:)) max(XX(:))];
ylims = [min(YY(:)) max(YY(:))];

currFig = get(gcf,'Number');
k = 1;
figure(currFig+1);set(gcf,'Position',[2 42 824 924])
figure(currFig+2);set(gcf,'Position',[2 42 824 924])
figure(currFig+3);set(gcf,'Position',[2 42 824 924])
for i = 1:3
    a = DecisionPointsFDf3{i,1};
    b = TotalPointsFDf3{i,1};
    
    nK = 9;
    K = (1/nK)*ones(nK);
    %tmp2 = conv2(tmp,K,'same');
    tmp2 = conv2(a,K,'same')./conv2(b,K,'same');
    tmp2(conv2(b,K,'same')<5) = nan;

    marginaldF = sum(conv2(a,K,'same'),1)./sum(conv2(b,K,'same'),1);
    marginalF = sum(conv2(a,K,'same'),2)./sum(conv2(b,K,'same'),2);

    % marginaldF = sum(a,1)./sum(b,1);
    % marginalF = sum(a,2)./sum(b,2);

    if plotFig
        %%%
        figure(currFig+1);
        subplot(3,2,k)
        surf(XX,YY,tmp2)
        xlabel('df (spk/s^2)');ylabel('f (hz)');zlabel('# of transitions/ # of data points')
        xlim(xlims);ylim(ylims);
        title(['Baseline rate: ' num2str(round(baselineRate(i),2))])
        subplot(3,2,k+1)
        surf(XX,YY,1./tmp2./f_orco.fs)
        xlabel('df (spk/s^2)');ylabel('f (hz)');zlabel('duration (s) before transition')
        xlim(xlims);ylim(ylims);
        title(['Baseline dur: ' num2str(round(1./baselineRate(i))./f_orco.fs,2)])
        suptitle([num2str(nK) ' Grid Conv, ' num2str(abs(n-1)./f_orco.fs.*1000) ' ms average'])

        %%%
        figure(currFig+2);
        subplot(3,2,k)
        plot(YY(:,1),marginalF)
        xlabel('f (hz)');ylabel('# of transitions/ # of data points');
        xlim(ylims);
        title(['Baseline rate: ' num2str(round(baselineRate(i),2))])
        subplot(3,2,k+1)
        plot(YY(:,1),1./marginalF./f_orco.fs)
        xlabel('f (hz)');ylabel('duration (s) before transition');
        xlim(ylims);
        title(['Baseline dur: ' num2str(round(1./baselineRate(i))./f_orco.fs,2)])
        suptitle(['Marginal of firing rate'])

        %%%
        figure(currFig+3);
        subplot(3,2,k)
        plot(XX(1,:),marginaldF)
        xlabel('df (spk/s^2)');ylabel('# of transitions/ # of data points');
        xlim(xlims);
        title(['Baseline rate: ' num2str(round(baselineRate(i),2))])
        subplot(3,2,k+1)
        plot(XX(1,:),1./marginaldF./f_orco.fs)
        xlabel('df (spk/s^2)');ylabel('duration (s) before transition');
        xlim(xlims);
        title(['Baseline dur: ' num2str(round(1./baselineRate(i))./f_orco.fs,2)])
        suptitle(['Marginal of change in firing rate'])
    end
    k = k+2;
end

end

















