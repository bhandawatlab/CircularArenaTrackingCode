% dateVal = "02-Dec-2022";
% getBadFlies
function [fileslist2, badFlies] = getBadFlies(dateVal)
C = readcell(['C:\Users\sbm94\Desktop\4ChamberBehaviorTracking_2023.xlsx']);

CStr = string(C(:, 1));
dateValIdx = find(CStr == dateVal);
badFlies2 = C(dateValIdx, 4)';
fileslist2 = C(dateValIdx, 5)';

missingIdx = zeros(1, numel(badFlies2));
for ww = 1:numel(badFlies2)
    try
    tempVal = sum(badFlies2{ww});
    missingIdx(1, ww) = ww;
%     ismissing(badFlies2{ww})
    catch
    if ismissing(badFlies2{ww})
       % do nothing
    else
        missingIdx(1, ww) = ww;
    end
    end
end
missingIdx = missingIdx(1, find(missingIdx(1, :) > 0));
% missingIdx = missingIdx(1, find(missingIdx(1, :) == 0));
badFlies = badFlies2;
%fileslist2 = fileslist2(1, missingIdx);

standard_chambers = [1,2,3,4];
for ww = 1:numel(badFlies)
    tempVal = badFlies{ww};
    
    goodfiles = str2num(string(tempVal));
    setdiff(standard_chambers,goodfiles);
    badFlies{1,ww} = setdiff(standard_chambers,goodfiles); %str2num(string(tempVal));
end % ww

% for ww = 1:numel(badFlies)
%     tempTot = [1 2 3 4];
%     tempVal = badFlies{ww};
%     badFlies{ww} = str2num(string(tempVal));
% end % ww


end %



