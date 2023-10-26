function [self] = getMeanBound(self,fe)
boundNdx = find(cellfun(@(c)strcmp(c,'boundary'),self.states.key,'UniformOutput',true));
boundaryCond = self.states.ndx==boundNdx;

phi = self.calcPhi;

durTmp = [];phiTmp = [];
durTmpBef = [];phiTmpBef = [];
for fly = 1:self.nFly
    % get the start and end index of each boundary track
    [startNdx,endNdx,type] = startEndSeq(boundaryCond(fly,:));
    startNdx(type==0) = [];
    endNdx(type==0) = [];
    type(type==0) = [];
    
    for track = 1:numel(type)
        % duration and angle travels around the border
        tmpDur = endNdx(track)-startNdx(track)+1;
        tmpPhi = phi(fly,startNdx(track):endNdx(track));
        
        % assign to either before or after first entry
        if startNdx(track)<fe(fly)
            durTmpBef = [durTmpBef tmpDur];
            phiTmpBef = [phiTmpBef mean(tmpPhi)];
        else
            durTmp = [durTmp tmpDur];
            phiTmp = [phiTmp mean(tmpPhi)];
        end
    end
end
phiTmp = abs(phiTmp);
phiTmpBef = abs(phiTmpBef);


CorrMat1 = corrcoef([phiTmpBef;durTmpBef]');
qBounds1 = [quantile(phiTmpBef,.05)  quantile(phiTmpBef,.95);
    quantile(durTmpBef,.05)  quantile(durTmpBef,.95)];
pHat1 = [lognfit(abs(phiTmpBef)+eps); lognfit(durTmpBef+eps)];

CorrMat2 = corrcoef([phiTmp;durTmp]');
qBounds2 = [quantile(phiTmp,.05)  quantile(phiTmp,.95);
    quantile(durTmp,.05)  quantile(durTmp,.95)];
pHat2 = [lognfit(abs(phiTmp)+eps); lognfit(durTmp+eps)];

self.model.boundary.state.state = 'boundary';
self.model.boundary.corrBaseline = CorrMat1;
self.model.boundary.corr = CorrMat2;
self.model.boundary.qBoundsBaseline = qBounds1;
self.model.boundary.qBounds = qBounds2;
self.model.boundary.baseline = pHat1;
self.model.boundary.during = pHat2;


end