function [ lengthOfConstriction65, allLengthOfConstriction ] = CalculateConstrictionLengths( firstTrachPos, area, arcLength )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Length of Constriction
constrictionFactor = 1.05;

% Crop area and arclength vectors to start at 1st tracheal ring
area = area(firstTrachPos:end);
arcLength = arcLength(firstTrachPos:end);
arcLength = arcLength - arcLength(1);

for j = 1:20
    constrictionFactor = constrictionFactor - 0.05;

    constrictedCoords = find(area < area(1) * constrictionFactor);

    % Check there's only one constriction
    if sum(diff(constrictedCoords) == 1) / (length(constrictedCoords) - 1) == 1
        lengthOfConstriction(j) = (arcLength(constrictedCoords(end)) - arcLength(constrictedCoords(1))) / arcLength(end);
    elseif size(constrictedCoords,1) > 1
       gaps = diff(constrictedCoords);

        t = diff([false;gaps==1;false]);
        p = find(t==1);
        q = find(t==-1);
        [maxlen,ix] = max(q-p);
        first = p(ix);
        last = q(ix)-1;
        %display({'i = ', num2str(i), ', j = ', num2str(j)});
        lengthOfConstriction(j) = (arcLength(constrictedCoords(last)) - arcLength(constrictedCoords(first))) / arcLength(end);
    elseif size(constrictedCoords,1) == 0
        % No constriction at this size
        lengthOfConstriction(j) = 0;    
    elseif size(constrictedCoords,1) == 1
        % 1 plane constriction at this size
        lengthOfConstriction(j) = ((arcLength(constrictedCoords+1) - arcLength(constrictedCoords))/2) / arcLength(end);
    else
        warning('How did I get here?');
    end
end

allLengthOfConstriction = lengthOfConstriction .* 100;

lengthOfConstriction65 = allLengthOfConstriction(8);

% ftab = figure;
% cnames = {'1.00','0.95','0.90','0.85','0.80','0.75','0.70','0.65','0.60','0.55','0.50','0.45','0.40','0.35','0.30','0.25','0.20','0.15'};
% rnames = {'N (A)', 'G1 (B)','G2 (C)','G3','G4','G5 (D)','G6','G7','G8 (E)', 'G11', 'G12', 'G13'};
% uitable('Parent', ftab, 'Data', lengthOfConstriction,'ColumnName',cnames,'RowName',rnames);


end

