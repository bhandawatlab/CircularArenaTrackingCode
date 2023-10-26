function [g,tt,params,err,errF,errP] = getModelParams(val,XX,YY,lambda,plotFig)

badNdx = isnan(squeeze(nanmean(val,[1,2])));
val(:,:,badNdx) = [];
XX(:,:,badNdx) = [];
YY(:,:,badNdx) = [];
xlims = [min(XX(:)) max(XX(:))];
ylims = [min(YY(:)) max(YY(:))];

tt = (0:size(val,3)-1).*300./30;

zmin = floor(min(val(:)));
zmax = ceil(max(val(:)));

% function
g = @(c1,c2,y,phi,a1,a2,b1,b2,freq,dFreq) ((y-freq).*(dFreq.^2)./1000.*(phi)+1).*...
        (c1./(1+exp(-(freq-a1)./b1)))+(c2./(1+exp((freq-a2)./b2)));
    
%     g = @(c1,c2,y,phi,a1,a2,b1,b2,alpha,freq,dFreq) (phi.*(y-freq))./(1+(dFreq./alpha).^2)...
%         .*(c1./(1+exp(-(freq-a1)./b1)))+(c2./(1+exp((freq-a2)./b2)));

% set optimization bounds and initial conditions
% lb = [-5,-5, 0,-1, 0, 0, 0, 0];
% ub = [5,5, 50, 1,40,15,20,10];
lb = [-5,-5,-50,-100,-40,-50,-20,-1000];
ub = [ 5, 5, 50, 100, 40, 50, 20,0];
x0 = [1,2,1,25,-1,18,5,5];
options = optimoptions('fmincon','Display','off');


% optimize all parameters
xAll = zeros(size(val,3),numel(x0));
for t = 1:size(val,3)
    testVal = val(:,:,t);
    if t == 1
        optimFunc = @(P) sqrt(nansum((g(P(1),P(2),P(3),P(4),P(5),P(6),P(7),P(8),...
            YY(:,:,t),XX(:,:,t))-testVal).^2,'all'));
    else
        optimFunc = @(P) sqrt(nansum((g(P(1),P(2),P(3),P(4),P(5),P(6),P(7),P(8),...
            YY(:,:,t),XX(:,:,t))-testVal).^2,'all')+lambda.*nansum(abs(P-x0)));
    end
    x2 = fmincon(optimFunc,x0,[],[],[],[],lb,ub,[],options);
    tmp = g(x2(1),x2(2),x2(3),x2(4),x2(5),x2(6),x2(7),x2(8),YY(:,:,t),XX(:,:,t));
    tmp(isnan(testVal)) = nan;
    
    x0 = x2;
    xAll(t,:) = x2;
end

% set c1,c2,y, and phi to the average value
c1 = nanmean(xAll(:,1));
c2 = nanmean(xAll(:,2));
x0 = [1,2,1,25,-1,18,5,5];


if plotFig
    figure;set(gcf,'Position',[2 42 789 954])
end
% reoptimize the rest of the parameters based on constant c1,c2,y, and phi
xAll = zeros(1,numel(x0));err = zeros(1,size(val,3));
for t = 1:size(val,3)
    testVal = val(:,:,t);
    
    if t == 1
        optimFunc = @(P) sqrt(nansum((g(c1,c2,P(3),P(4),P(5),P(6),P(7),P(8),...
        YY(:,:,t),XX(:,:,t))-testVal).^2,'all'));
    
        optimFuncF = @(P) nansum((g(c1,c2,P(3),P(4),P(5),P(6),P(7),P(8),...
            YY(:,:,t),XX(:,:,t))-testVal).^2,'all');
        optimFuncP = @(P) nansum(abs(P-P));
    else
        optimFunc = @(P) sqrt(nansum((g(c1,c2,P(3),P(4),P(5),P(6),P(7),P(8),...
        YY(:,:,t),XX(:,:,t))-testVal).^2,'all')+lambda.*nansum(abs(P-x0)));
    
        optimFuncF = @(P) nansum((g(c1,c2,P(3),P(4),P(5),P(6),P(7),P(8),...
            YY(:,:,t),XX(:,:,t))-testVal).^2,'all');
        optimFuncP = @(P) nansum(abs(P-x0));
    end
    
    x2 = fmincon(optimFunc,x0,[],[],[],[],lb,ub,[],options);
    tmp = g(c1,c2,x2(3),x2(4),x2(5),x2(6),x2(7),x2(8),YY(:,:,t),XX(:,:,t));

    tmp2 = tmp;
    tmp(isnan(testVal)) = nan;
    
    x0 = x2;
    xAll(t,:) = x2;
    err(t) = optimFunc(x2);
    
    errF(t) = optimFuncF(x2);
    errP(t) = optimFuncP(x2);
    
    % plot the 3d surf plots 
    if plotFig
        subplot(6,4,(t-1).*2+1);
        surf(XX(:,:,1),YY(:,:,1),testVal, 'edgecolor','none');view(45,45)%view(2);
        xlim(xlims);ylim(ylims);zlim([zmin zmax]);title(['raw, t=' num2str(tt(t))])
        caxis([round(min(val(:)),1), round(max(val(:)),1)])
        subplot(6,4,t.*2);hold on
        surf(XX(:,:,1),YY(:,:,1),tmp2, 'edgecolor','none', 'FaceAlpha',0.2);
        surf(XX(:,:,1),YY(:,:,1),tmp, 'edgecolor','none');hold off;
        view(45,45)%view(2);
        xlim(xlims);ylim(ylims);zlim([zmin zmax]);title(['model, t=' num2str(tt(t))])
        caxis([round(min(val(:)),1), round(max(val(:)),1)])
        
        if t == 1
            subplot(6,4,(t-1).*2+1);
            xlabel('df');ylabel('f');
            subplot(6,4,t.*2);
            xlabel('df');ylabel('f');
        end
    end
end

valList = {'c1','c2','y','phi','a1','a2','b1','b2'};

xAll(:,1) = c1;
xAll(:,2) = c2;

%x = (0:size(xAll,1)-1).*300./30;
x2 = (0:0.2:size(val,3)-1).*300./30;
if plotFig
    figure;set(gcf,'Position',[1 40 585 600])
    for i = 1:size(xAll,2)
        yy=interp1(tt,xAll(:,i)',x2,'spline');
        
        subplot(4,2,i);
        plot(tt,xAll(:,i));hold on
        plot(x2,yy);
        xlabel('time (s)')
        title(valList{i})
        if i-4 == 1
            legend({'param fit','spline'})
        end
    end
end

params.c1 = c1;
params.c2 = c2;
params.y = xAll(:,3);
params.phi = xAll(:,4);
params.a1 = xAll(:,5);
params.a2 = xAll(:,6);
params.b1 = xAll(:,7);
params.b2 = xAll(:,8);

end