function [ diameter ] = DiameterCalc( Fsurface, Vsurface, perimeterEdges )
% Calculate the diameter of the trachea, called by PerimeterhydDiamCalc


%% perimeterEdges is a list of nodes on the edge of the STL. Make a unique 
%  list of these. This is most probably all the nodes in it, but just to be
%  sure all nodes are on the edge.

%[a,b] = size(perimeterEdges);
%perimeterNodes = reshape(perimeterEdges,a*b,1);

%perimeterNodes = unique(perimeterNodes);

%perimeterCoords = Vsurface(perimeterNodes,:);


%% Put the nodes in order around the edge

orderedNodes = perimeterEdges(1,:)';
count = 0;
%while size(orderedNodes,1) < size(perimeterNodes,1)

%% Check for non manifold edges
uniqueNodes = unique(perimeterEdges);
% noUniqueNodes = size(uniqueNodes,1);
% nodeCount = zeros(noUniqueNodes,1);
% for i = 1:noUniqueNodes
%    nodeCount(i) = sum(reshape(perimeterEdges,[a*b,1])==uniqueNodes(i));
% end
% if sum(nodeCount~=2)>0
%     display('Non-manifold edge found');
% end
%%

while true
    count = count + 1;
    [rows,~] = find(perimeterEdges == orderedNodes(end,1));
    possibleNodes = perimeterEdges(rows,:);
    
    if ~isempty(find(possibleNodes==orderedNodes(1))) && count > 1
       if size(orderedNodes,1) < size(uniqueNodes,1) && size(orderedNodes,1) < 4
           % We have found a small subloop round the perimeter and not the main loop
           % - remove this from the calculation and start the while loop again.
           perimeterEdges(orderedNodes,:) = [];
           orderedNodes = perimeterEdges(1,:)';
           count = 0;
           possibleNodes = [];
           display('Missing nodes from diameter due to non manifold surface - remove these and try again!');
           continue
       else
           %display('Complete loop found!!');
           break
       end
    end
    
    foundNewNodes = setdiff(possibleNodes,orderedNodes);
    if size(foundNewNodes) == 1
        orderedNodes(end+1,1) = foundNewNodes;
    elseif size(foundNewNodes,1) > 1
        orderedNodes(end+1,1) = foundNewNodes(1);
        alternates= foundNewNodes(2:end);
    elseif isempty(foundNewNodes) && exist('alternates','var') && ~isempty(alternates)
        % This is a bit of a hack
        alternates = setdiff(alternates,orderedNodes);
        if ~isempty(alternates)
            orderedNodes(end+1,1) = alternates(1);
            alternates(1) = [];
        else
            warning('Couldnt make a complete loop - no diameter calculated')
            diameter = [NaN NaN NaN];
            return
        end
    else
        warning('Couldnt make a complete loop - no diameter calculated')
        diameter = [NaN NaN NaN];
        return
    end
end
% orderedNodes = perimeterEdges(1,:);
% columnToSearch = 1;
% columnToTake = 2;
% direction = 'first';
% for i = 2:size(perimeterNodes,1)  
%     [row,~] = find(perimeterEdges(:,columnToSearch) == orderedNodes(end),1,direction);
%     if isempty(row)
%         columnToSearch = setdiff([1 2],columnToSearch);
%         columnToTake = setdiff([1 2],columnToSearch);
%         if strcmp(direction,'first')
%             direction = 'last';
%         elseif strcmp(direction,'last')
%             direction = 'first';
%         else
%            display('How did I get here?'); 
%         end
%         continue
%     else
%         orderedNodes(end+1) = perimeterEdges(row,columnToTake);
%     end
%     
%     
% end

perimeterCoords = Vsurface(orderedNodes,:);

%% Find the distance between each pair of nodes in the list

for i=1:size(perimeterCoords,1)
    for j=1:size(perimeterCoords,1)
        if i>j
            d(i,j) = sqrt(((perimeterCoords(i,1) - perimeterCoords(j,1)).^2) + ((perimeterCoords(i,2) - perimeterCoords(j,2)).^2) + ((perimeterCoords(i,3) - perimeterCoords(j,3)).^2));
        end
    end
end

%% Find the major axis - the biggest distance between two nodes

bigDist = max(max(d));
[bigDistNodes(1),bigDistNodes(2)]  = find(d == bigDist,1); %If mulitple distances the same size are found, arbitrarily take the first one



%% Now find the two nodes closest to 90 degrees from that - the minor axes

bigDist3D = perimeterCoords(bigDistNodes(1),:) - perimeterCoords(bigDistNodes(2),:);

midpoint = perimeterCoords(bigDistNodes(1),:) - (bigDist3D .* 0.5);

% For each point, find the vector from that point to the midpoint

% midVecs = bsxfun(@minus,perimeterCoords, midpoint);
midVecs = bsxfun(@minus,Vsurface, midpoint);
for i= 1:size(midVecs,1)
    distFromPlane(i) = dot(midVecs(i,:), bigDist3D) / (norm(midVecs(i,:),2) .* norm(bigDist3D,2));
end

% Split the nodes into two halves separated by the major axes
% bigNodeOneOrdered = find(orderedNodes == bigDistNodes(1));
% bigNodeTwoOrdered = find(orderedNodes == bigDistNodes(2));

bigNodeOneOrdered = bigDistNodes(1);
bigNodeTwoOrdered = bigDistNodes(2);

if bigNodeOneOrdered > bigNodeTwoOrdered
    nodeListOne = orderedNodes(bigNodeTwoOrdered:bigNodeOneOrdered);
    nodeListTwo = [orderedNodes(bigNodeOneOrdered:end);orderedNodes(1:bigNodeTwoOrdered)];
else
    nodeListOne = orderedNodes(bigNodeOneOrdered:bigNodeTwoOrdered);
    nodeListTwo = [orderedNodes(bigNodeTwoOrdered:end);orderedNodes(1:bigNodeOneOrdered)];
end


% if bigNodeOneOrdered > bigNodeTwoOrdered
%     nodeListOne = orderedNodes(bigNodeTwoOrdered:bigNodeOneOrdered);
%     nodeListTwo = [orderedNodes(bigNodeOneOrdered:end);orderedNodes(1:bigNodeTwoOrdered)];
% else
%     nodeListOne = orderedNodes(bigNodeOneOrdered:bigNodeTwoOrdered);
%     nodeListTwo = [orderedNodes(bigNodeTwoOrdered:end);orderedNodes(1:bigNodeOneOrdered)];
% end
    
% Find min ABSOLUTE distance in node list one and two, they should be the
% minor axis nodes!

% minOne = min(abs(distFromPlane(orderedNodes(nodeListOne))));
% minTwo = min(abs(distFromPlane(orderedNodes(nodeListTwo))));
minOne = min(abs(distFromPlane(nodeListOne)));
minTwo = min(abs(distFromPlane(nodeListTwo)));

if isempty(minOne) || isempty(minTwo)
    display('Unable to calculate diameters due to non-manifold edges or free surfaces not on the edge.');
    display('Returning NaN for diameters');
    diameter = [NaN NaN NaN];
    return
end

minNodeTemp = find(abs(distFromPlane) == minOne);

if size(minNodeTemp,2) == 1
    minNode(1) = minNodeTemp;
elseif size(minNodeTemp,2) > 1
    if size(intersect(minNodeTemp, nodeListOne),1) == 1
        try
            minNode(1) = intersect(minNodeTemp, nodeListOne);
        catch
            display('what now')
        end
    else
        %% This is a hack - should find a proper way to determine this
         minNode(1) = minNodeTemp(1);
         warning('Arbitrarily choosing first value as minimum. There are matching minimum distance values from plane');
    end
else
    warning('No minimum detected');
end

minNodeTemp = find(abs(distFromPlane) == minTwo);
if size(minNodeTemp,2) == 1
    minNode(2) = minNodeTemp;
elseif size(minNodeTemp,2) > 1
    if size(intersect(minNodeTemp, nodeListOne),1) == 1
        minNode(2) = intersect(minNodeTemp, nodeListOne);
    else
        %% This is a hack - should find a proper way to determine this
         minNode(2) = minNodeTemp(1);
         warning('Arbitrarily choosing first value as minimum. There are matching minimum distance values from plane');
    end
else
    warning('No minimum detected');
    diameter = [NaN NaN NaN];
    return
end


% % % %close all
% % % figure;
% % % % plot3(Vsurface(:,1),Vsurface(:,2),Vsurface(:,3),'rx')
% % % hold on
% % % patch('Vertices',Vsurface,'Faces',Fsurface,'facecolor', 'b','edgecolor','none')
% % % %plot3(Vsurface(orderedNodes(bigDistNodes),1),Vsurface(orderedNodes(bigDistNodes),2),Vsurface(orderedNodes(bigDistNodes),3),'k-o')
% % % plot3(Vsurface(orderedNodes(bigDistNodes),1),Vsurface(orderedNodes(bigDistNodes),2),Vsurface(orderedNodes(bigDistNodes),3),'k-')
% % % %plot3(Vsurface(bigDistNodes,1),Vsurface(bigDistNodes,2),Vsurface(bigDistNodes,3),'k-o')
% % % plot3(midpoint(:,1),midpoint(:,2),midpoint(:,3),'bx')
% % % %plot3(Vsurface(orderedNodes(minOrderedNode),1),Vsurface(orderedNodes(minOrderedNode),2),Vsurface(orderedNodes(minOrderedNode),3),'gx');
% % % %plot3(Vsurface(minNode,1),Vsurface(minNode,2),Vsurface(minNode,3),'g-o');
% % % plot3(Vsurface(minNode,1),Vsurface(minNode,2),Vsurface(minNode,3),'g-');
% % % daspect([1 1 1])


% figure
% hold on
% for i = 1:size(perimeterNodes,1) 
%     x = Vsurface(perimeterEdges(i,:),1);
%     y = Vsurface(perimeterEdges(i,:),2);
%     z = Vsurface(perimeterEdges(i,:),3);
%     plot3(x,y,z,'x-');
% end

%% Find angle between major axes vector (projected onto xy plane) and Y axis (Ap - I think X is Rl)

% Project major axis onto XY plane = just ignore the Z component
bigDist3Dproj = [bigDist3D(1) bigDist3D(2) 0 ];

% Find the angle between our major diameter vector and [0 1 ]
yax = [ 0 1 0];
majorDiAngleFromAPdeg = atan2d(norm(cross(bigDist3Dproj,yax)),dot(bigDist3Dproj,yax));

minDist = sqrt(((Vsurface(minNode(1),1) - Vsurface(minNode(2),1)).^2) + ((Vsurface(minNode(1),2) - Vsurface(minNode(2),2)).^2) + ((Vsurface(minNode(1),3) - Vsurface(minNode(2),3)).^2));

%text(midpoint(:,1),midpoint(:,2),midpoint(:,3) + 2,num2str(minDist/bigDist));

diameter = [bigDist minDist majorDiAngleFromAPdeg];

end