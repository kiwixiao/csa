function [ deviation ] = calculateDeviation( firstTrachPos, pos )
% Calculate the deviation of the tracheal centreline from a straight line
% running from the centreline at the first tracheal ring to the centreline
% at the end of the trachea
%   Detailed explanation goes here

    % Crop pos to start at the first tracheal ring
    %pos = pos(firstTrachPos:end,:);

    firstPos = pos(firstTrachPos,:);
    endPos = pos(end,:);

    for i = 1:size(pos,1)
        if i <= firstTrachPos
            deviation(i) = 0;
        else

            xA = firstPos(1);
            xB = endPos(1);
            xC = pos(i,1);

            yA = firstPos(2);
            yB = endPos(2);
            yC = pos(i,2);

            zA = firstPos(3);
            zB = endPos(3);
            zC = pos(i,3);

            % Area calculated based on http://en.wikipedia.org/wiki/Triangle#Using_coordinates
            a = [xA, xB, xC; yA, yB, yC; 1, 1, 1 ];
            b = [yA, yB, yC; zA, zB, zC; 1, 1, 1 ];
            c = [zA, zB, zC; xA, xB, xC; 1, 1, 1 ];

            twiceArea = sqrt((abs(det(a))^2) + (abs(det(b))^2) + (abs(det(c))^2));

            base = sqrt ( ((xB - xA) ^2) + ((yB - yA) ^2) + ((zB - zA) ^2));

            deviation(i) = twiceArea / base;
        end        
        

    end
    %[maxDev, maxDevPos] = max(deviation);
end

