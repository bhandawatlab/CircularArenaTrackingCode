function [Vq,newState,XX,YY] = interpretLandscape(vals,XX,YY,pdN,plotFig)
% set any nan values to the closest non-nan value
    [noDatr, noDatc] = find(isnan(vals));
    [datr, datc] = find(~isnan(vals));
    D = pdist2([noDatr, noDatc],[datr, datc]);
    [~,nearestPtNdx]=min(D,[],2);
    for i = 1:numel(noDatr)
        vals(noDatr(i), noDatc(i)) = vals(datr(nearestPtNdx(i)), datc(nearestPtNdx(i)));
    end
    
    % smooth the landscape
    half_pdN = pdN./2;
    vals_pad = padarray(vals,[pdN pdN],'replicate','both');
    newState = conv2(vals_pad, ones(half_pdN,1)/half_pdN, 'same');
    newState = newState(pdN+1:end-pdN,pdN+1:end-pdN);
    Vq = @(X,Y) interp2(XX,YY,newState,X,Y);
    [X,Y] = meshgrid([-5:0.1:5],[0:0.5:45]);
    
    if plotFig
        figure;set(gcf,'Position',[6 39 409 736])
        subplot(2,1,1);surf(XX(:,:,1),YY(:,:,1),vals)
        subplot(2,1,2);surf(X,Y,Vq(X,Y))
    end
end