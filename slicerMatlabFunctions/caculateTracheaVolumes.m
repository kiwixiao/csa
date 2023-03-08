function [ volume ] = caculateTracheaVolumes( firstTrachPos, arcLength, area )
%Calculate Tracheal Volumes from 1st trach ring to end
%   Detailed explanation goes here

arcLength = arcLength(firstTrachPos:end);
area = area(firstTrachPos:end);

diffArc = diff(arcLength);
diffArc(end + 1) = mean(diffArc);

volume =  sum(diffArc .* area);


end

