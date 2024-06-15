 %%find alpha parameters using fooof
%outputs: 
% CFAll: Centre frequncy, PHAll: peak height, BWAll: band width and
%APAll: alpha power above aperidic activity
%APNew : full alpha power including aperiodic activity 
function [CFAll, PHAll, BWAll, APAll, APNew] = getAlphaParameters(fooof_results,freqVals,SpecPower, medianFlag)

if ~exist('medianFlag','var');   medianFlag=1;   end

alphaRange = [7 12];
 [~,~,~,electrodeGroupList] = electrodePositionOnGrid(64,'EEG');

 if length(fooof_results)~=64
     for j = length(fooof_results)+1:64
     fooof_results{j}=[];
     end
 end
 fooofOcc = fooof_results(electrodeGroupList{1});
 specOcc = SpecPower(electrodeGroupList{1},:);
 alphaPos = intersect(find(freqVals>=alphaRange(1)),find(freqVals<=alphaRange(2)));
 freqRes = freqVals(2)-freqVals(1);
 for i=1:length(fooofOcc)
     if ~isempty(fooofOcc{i})
         peakParams = fooofOcc{i}(1).peak_params;
         cFNum = find(peakParams(:,1)>alphaRange(1) & peakParams(:,1)<alphaRange(2));
         if length(cFNum)>1
             d = find(peakParams(:,2)==max(peakParams(cFNum,2)));
             
%               CF(i) = peakParams(d,1);  %centre frequency
%               PH(i) = peakParams(d,2);   %peak height 
         else
             d = cFNum;
        
         end
         CF(i) = nanmean(peakParams(d,1));  %centre frequency
            PH(i) = nanmean(peakParams(d,2));   %peak height
         
            BW(i) = nansum(peakParams(cFNum,3));    %bandwidth
         AP(i) = 1.064467*PH(i)*BW(i);  %alpha power gen from fooof params 
         AP2(i) = nansum(specOcc(i,alphaPos))*freqRes; %alpha power from old theory
     else
         CF(i) = nan; PH(i) = nan; BW(i) = nan; AP(i) = nan; AP2(i) = nan;
     end
 end

 if medianFlag
     CFAll = nanmedian(CF); PHAll = nanmedian(PH); BWAll = nanmedian(BW);   APAll = nanmedian(AP);  APNew = nanmedian(AP2);
 else
      CFAll = nanmean(CF); PHAll = nanmean(PH); BWAll = nanmean(BW);   APAll = nanmean(AP); APNew = nanmedian(AP2);
 end
end