function [area, centroid, growingnodes, tri] = sublobeMulti_DG_AB4(PlaneNormal, PlanePoint, Vsurface, Fsurface, planenumber)

    %clc, clear all
    %figure;
    %hold on
    
    % SplitElements - each row is a list of 3 indexes of nodes making up
    % one triangle. The nodes are stored in SplitNodes. The fourth column
    % is a number representing individual separate loops.
   

%     %vertices, triangles and normals
%     [Vsurface,Fsurface,normal] = STL_Import('TestCylinder.stl');
% 
%     %Point and normal of plane
%     PlaneNormal=[0,0,1];
%     PlanePoint=[0,0,0.505];

    %normalise
    PlaneNormal=PlaneNormal/norm(PlaneNormal);
    
    [SplitElements, SplitNodes, WithNodeNumber] = SplitLobes(Fsurface, Vsurface, PlanePoint, PlaneNormal);
    if isempty(SplitElements)
        % This means the current plane does not split the mesh, so return
        % -1 to make this clear.
        warning(['Plane number ',num2str(planenumber), ' not split at this point']);
        area{1} = -1;
        centroid = PlanePoint;
        growingnodes = cell(3,1);
        tri = cell(3,1);
        return;
    else
        % Calculate centroid of plane to use as a starting point for future
        % Planes. Currently this is not a centroid, but the mean of points.
        % Close enough for reasonably uniformly distributed points?
        %%%SplitElements = DetectConnectedElements(SplitElements);
        SplitElements = DetectTwoNodeConnectedElements(SplitElements);
    end
    NoOfCuts = max(SplitElements(:,4));
    %SplitElements(:,4) = 1;
    %NoOfCuts = 1;
    
    
    tri=cell(3,1);
    area=cell(3,1);
    centroid=cell(3,1);
    growingnodes=cell(3,1);
    numberOfNodes=cell(3,1);
    allNodeList=cell(3,1);
    allCentroidPos=cell(3,1);
    splitPos=cell(3,NoOfCuts);
    %allCentroidPos=cell(3,1);
    %notCheckedNodes = zeros(3,1);
    
    L = 0;
    R = 0;
    B = 0;
    
    for i = 1:NoOfCuts
        ConnectedElements = SplitElements(:,4) == i;
        ConnectedSplitElements = SplitElements(ConnectedElements,1:3);
        %ConnectedSplitElements = OrderSplitElements(ConnectedSplitElements,Vsurface);
        [ConnectedSplitElements, Ordered] = NewOrderSplitElements(ConnectedSplitElements, planenumber,Vsurface);
         
        if ~Ordered 
            warning('Unable to split plane at this point');
            area{1} = -1;
            centroid = -1;
            growingnodes = -1;
            tri = -1;
            return;
        end


        [I] = FindPoints(WithNodeNumber, ConnectedSplitElements, SplitNodes, PlanePoint, PlaneNormal);
        NewVect = DeleteDuplicateNodes(I);
        centroidPos(i,:) = mean(NewVect(:,1:3),1);
        
        %% TEST DEBUG ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Test to see if what we've cut is near where we meant it to be this is only for tracheas...
        if abs(norm(centroidPos(i,:) - PlanePoint)) > 10 %0.1
            badPlane = true;
        else
            badPlane = false;
        end
        badPlane = false;
        if ~badPlane
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if NoOfCuts == 1 % I assume if this is missed here due to badPlane it is dealt with after doughnut tidy up.
                % Assume we're in the post septum joined up region
                side = 3;
                B = B + 1;
                j = B;
            elseif centroidPos(i,1) > PlanePoint(1)
                % Assume this is part of the left nostril
                side = 1;
                L = L + 1;
                j = L;
            else
                % Assume this is part of the right nostril
                side = 2;
                R = R + 1;
                j = R;
            end

            % Store list of all nodes so we can check for doughnuts
            allNodeList{side, j} = NewVect;
            allCentroidPos{side,j} = centroidPos(i,:);
        end
    end
    
    % Check for "doughnuts" i.e. planes within other planes that should be
    % holes
    doughnutsPlanes = CheckForDoughnuts(allNodeList);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Rearrange node lists to take doughnuts into account
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for side = 1:3
%         for otherSide = 1:3
     for side = 3:-1:1
         for otherSide = 3:-1:1

            if ~isempty(doughnutsPlanes{side,otherSide})
                for i = size(doughnutsPlanes{side,otherSide},1):-1:1
                    if otherSide == 1
                        if isempty(splitPos{side,doughnutsPlanes{side}(i,1)})
                            splitPos{side,doughnutsPlanes{side,otherSide}(i,1)}(end + 1) = size(allNodeList{side,doughnutsPlanes{side,otherSide}(i,2)},1) ;
                        else
                            splitPos{side,doughnutsPlanes{side,otherSide}(i,1)}(end + 1) = size(allNodeList{side,doughnutsPlanes{side,otherSide}(i,2)},1) + splitPos{side,doughnutsPlanes{side,otherSide}(i,1)}(end);
                        end
                    else
                        if isempty(splitPos{otherSide-1,doughnutsPlanes{side,otherSide}(i,1)}) % This is equivalent to the j before - this could all be combined.
                            splitPos{otherSide-1,doughnutsPlanes{side,otherSide}(i,1)}(end + 1) = size(allNodeList{side,doughnutsPlanes{side,otherSide}(i,2)},1) ;
                        else
                            splitPos{otherSide-1,doughnutsPlanes{side,otherSide}(i,1)}(end + 1) = size(allNodeList{side,doughnutsPlanes{side,otherSide}(i,2)},1) + splitPos{otherSide-1,doughnutsPlanes{side,otherSide}(i,1)}(end);
                        end
                    end
                end
            end
        end
    end
    
    for side = 1:3
        for otherSide = 1:3
            if ~isempty(doughnutsPlanes{side,otherSide})
                for i = 1:size(doughnutsPlanes{side,otherSide},1)

                
    %                 if size(allNodeList{side,doughnutsPlanes{side}(i,2)},1) > size(allNodeList{side,doughnutsPlanes{side}(i,1)},1)
    %                     largerCut = 2;
    %                 else
    %                     largerCut = 1;
    %                 end

    %                 if i == 1
                        %splitPos{side,doughnutsPlanes{side}(i,1)}(end + 1) = size(allNodeList{side,doughnutsPlanes{side}(i,2)},1) ;
    %                 else
    %                     splitPos{side,doughnutsPlanes{side}(i,1)}(end + 1) = size(allNodeList{side,doughnutsPlanes{side}(i,largerCut)},1) + splitPos{side,doughnutsPlanes{side}(i,1)}(end);
    %                 end


                    if otherSide == 1
                        allCentroidPos{side,doughnutsPlanes{side}(i,1)} = ((centroidPos(doughnutsPlanes{side}(i,1),:) * size(allNodeList{side,doughnutsPlanes{side}(i,1)},1))...
                                         + (centroidPos(doughnutsPlanes{side}(i,2),:) * size(allNodeList{side,doughnutsPlanes{side}(i,2)},1)))...
                                         / (size(allNodeList{side,doughnutsPlanes{side}(i,1)},1) + size(allNodeList{side,doughnutsPlanes{side}(i,2)},1));
                                     
                        allCentroidPos{side,doughnutsPlanes{side}(i,2)} = [];
                        
                        allNodeList{side,doughnutsPlanes{side}(i,1)} = [allNodeList{side,doughnutsPlanes{side}(i,2)};allNodeList{side,doughnutsPlanes{side}(i,1)}];
                        
                        allNodeList{side,doughnutsPlanes{side}(i,2)} = [];


                    else
                        % This means we've found doughnuts with the inside
                        % and outside on different sides (will be 1 and 2
                        % or 2 and 1)

                        
                        allCentroidPos{otherSide-1,doughnutsPlanes{side,otherSide}(i,1)} = ((centroidPos(doughnutsPlanes{side,otherSide}(i,1),:) * size(allNodeList{otherSide-1,doughnutsPlanes{side,otherSide}(i,1)},1))...
                                         + (centroidPos(doughnutsPlanes{side,otherSide}(i,2),:) * size(allNodeList{side,doughnutsPlanes{side,otherSide}(i,2)},1)))...
                                         / (size(allNodeList{otherSide-1,doughnutsPlanes{side,otherSide}(i,1)},1) + size(allNodeList{side,doughnutsPlanes{side,otherSide}(i,2)},1));
                        
                        allCentroidPos{side,doughnutsPlanes{side,otherSide}(i,2)} = [];
                        
                        allNodeList{otherSide-1,doughnutsPlanes{side,otherSide}(i,1)} = [allNodeList{side,doughnutsPlanes{side,otherSide}(i,2)};allNodeList{otherSide-1,doughnutsPlanes{side,otherSide}(i,1)}];
                        
                        allNodeList{side,doughnutsPlanes{side,otherSide}(i,2)} = [];
                        
                        
                        % TODO
                        % Change splitpos setup to handle this
                        % Check what happens if I delete a plane from the
                        % middle of a list on one side or the other
                        % Check going to fewer cuts still works below
                        % Check I've done this the right way round, so
                        % things end up on the correct side (if not going
                        % to side 3)
                    end
                    
                    NoOfCuts = NoOfCuts - 1;
                    if NoOfCuts == 1
                        % The Doughnut detection has made an area that had two
                        % cuts into an area with one cut.
                        for k = 1:3
                            if ~isempty(allNodeList{k})
                                remove = k;
                            end
                        end
                        allNodeList{3,1} = allNodeList{remove,1};
                        allNodeList{remove,1} = [];
                        splitPos{3,1} = splitPos{remove,1};
                        splitPos{remove,1} = [];

                        
                        allCentroidPos{3,1} = allCentroidPos{remove,1}; %% HAS BEEN CHANGED by Valerian
                        allCentroidPos{remove,1} = [];
                    end
                end
            end
        end
    end
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % AJB addition on 8/Dec/2020 - it seems we can remove
    % the first element on a "side" while leaving later
    % elements. Make each row not have any gaps in it.
    
    
    side1Noderow = allNodeList(1,~cellfun('isempty',allNodeList(1,:)));
    side2Noderow = allNodeList(2,~cellfun('isempty',allNodeList(2,:)));
    side3Noderow = allNodeList(3,~cellfun('isempty',allNodeList(3,:)));
    side1Centrow = allCentroidPos(1,~cellfun('isempty',allCentroidPos(1,:)));
    side2Centrow = allCentroidPos(2,~cellfun('isempty',allCentroidPos(2,:)));
    side3Centrow = allCentroidPos(3,~cellfun('isempty',allCentroidPos(3,:)));
    
    while size(side1Noderow,2) < NoOfCuts
        side1Noderow{1,end+1} = [];
        side1Centrow{1,end+1} = [];
    end
    while size(side2Noderow,2) < NoOfCuts
        side2Noderow{1,end+1} = [];
        side2Centrow{1,end+1} = [];
    end
    while size(side3Noderow,2) < NoOfCuts
        side3Noderow{1,end+1} = [];
        side3Centrow{1,end+1} = [];
    end
    
    allNodeList = [side1Noderow; side2Noderow; side3Noderow];
    allCentroidPos = [side1Centrow; side2Centrow; side3Centrow];

                        
 %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    
    L = 0;
    R = 0;
    B = 0;
%     for i = 1:NoOfCuts
%         
%         if NoOfCuts == 1
%             % Assume we're in the post septum joined up region
%             side = 3;
%             B = B + 1;
%             j = B;
%         elseif centroidPos(i,1) > xSplit
%             % Assume this is part of the left nostril
%             side = 1;
%             L = L + 1;
%             j = L;
%         else
%             % Assume this is part of the right nostril
%             side = 2;
%             R = R + 1;
%             j = R;
%         end 

    for side = 1:3
        for j = 1:size(allNodeList,2)
            if ~isempty(allNodeList{side,j})
                NewVect = allNodeList{side,j};
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [n1,n2] = size(NewVect);
                if n1 < 3 || n2 < 3
                    'hello'
                    save('debugInfo.mat');
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%


                %[area{side}(j), nodes, triangles] = MeshPlane(NewVect,splitPos{side,j});
                
                
                %% DEBUG TEST ONLY
                if (planenumber < 0)
                    [area{side}(j), nodes, triangles] = NewMeshPlane(NewVect,splitPos{side,j},PlaneNormal);
                    
                    triangles = triangles + size(growingnodes{side},1);
                    growingnodes{side} = [ growingnodes{side}; nodes];
                    tri{side} = [tri{side}; triangles];
                    numberOfNodes{side}(j) = size(NewVect,1);
                    centroid{side}(j,:) = allCentroidPos{side,j} * numberOfNodes{side}(j);
                elseif j == 1
                    [area{side}(j), nodes, triangles] = NewMeshPlane(NewVect,splitPos{side,j},PlaneNormal);
                    triangles = triangles + size(growingnodes{side},1);
                    growingnodes{side} = [ growingnodes{side}; nodes];
                    tri{side} = [tri{side}; triangles];
                    numberOfNodes{side}(j) = size(NewVect,1);
                    centroid{side}(j,:) = allCentroidPos{side,j} * numberOfNodes{side}(j);
                else
                    [areaTemp, nodes, triangles] = NewMeshPlane(NewVect,splitPos{side,j},PlaneNormal);
                    %% HACK on 13th January 2017 as I keep losing areas I need due to curvature of centreline
                    %%% Triple commented on the left is removal of this
                    %%% hack in Spetember 2018 as I need to ignore further
                    %%% points now...
                    %if abs(norm(mean(nodes) - PlanePoint)) < abs(norm(mean(growingnodes{side}) - PlanePoint))
                    if abs(norm(mean(nodes) - PlanePoint)) < abs(norm(mean(growingnodes{side}) - PlanePoint))
                        if size(nodes,1) > size(growingnodes{side},1)
                            display(['*** planenumber = ' num2str(planenumber) ', Losing one area - choosing closer one']);
                            area{side}(1) = areaTemp;
                            growingnodes{side} = nodes;
                            tri{side} = triangles;
                            %%% 26th October 2016
                            %centroid{side} = mean(nodes);
                            %centroid{side} = allCentroidPos{side,2} * size(growingnodes{side},1);
                            centroid{side} = allCentroidPos{side,j} * size(growingnodes{side},1);
                            numberOfNodes{side} = size(growingnodes{side},1);
                        elseif size(nodes,1) > 0.1 * size(growingnodes{side},1)
                            display('#######################################################');
                            display('#### New area is closer but bigger than 0.1 of the old plane - using the new plane');
                            display(['#### planenumber = ' num2str(planenumber) ]);
                            display('#######################################################');
                            area{side}(1) = areaTemp;
                            growingnodes{side} = nodes;
                            tri{side} = triangles;
                            %%% 26th October 2016
                            %centroid{side} = mean(nodes);
                            %centroid{side} = allCentroidPos{side,2} * size(growingnodes{side},1);
                            centroid{side} = allCentroidPos{side,j} * size(growingnodes{side},1);
                            numberOfNodes{side} = size(growingnodes{side},1);
                        else
                            display('#######################################################');
                            display('#### New area is closer but smaller - ignore it for now');
                            display(['#### planenumber = ' num2str(planenumber) ]);
                            display('#######################################################');
                        end
                         %%%
% % %                     elseif size(nodes,1) > size(growingnodes{side},1)
% % %                         display('#######################################################');
% % %                         display('#### New area is further but bigger - we will use it for now');
% % %                         display('#######################################################');
% % %                                                     area{side}(1) = areaTemp;
% % %                             growingnodes{side} = nodes;
% % %                             tri{side} = triangles;
% % %                             %%% 26th October 2016
% % %                             %centroid{side} = mean(nodes);
% % %                             centroid{side} = allCentroidPos{side,2} * size(growingnodes{side},1);
% % %                             numberOfNodes{side} = size(growingnodes{side},1);
                    else
                        if size(nodes,1) > size(growingnodes{side},1)
                            display('#######################################################');
                            display('#### New area is further but LARGER');
                            display('#######################################################');
                            
                            if areaTemp > 6 * area{side}(1) % New area is 5 times bigger or more, we prefer
                                display('#######################################################');
                                display('#### New plane more than 10x larger, using the larger area');
                                display('#### This is the case where we prefer further plane when the plane is in the nose');
                                display(['#### planenumber = ' num2str(planenumber) ]);
                                display('#######################################################');
                                

                                display(['*** planenumber = ' num2str(planenumber) ', Losing one area - choosing closer one']);
                                area{side}(1) = areaTemp;
                                growingnodes{side} = nodes;
                                tri{side} = triangles;
                                %%% 26th October 2016
                                %centroid{side} = mean(nodes);
                                centroid{side} = allCentroidPos{side,2} * size(growingnodes{side},1);
                                numberOfNodes{side} = size(growingnodes{side},1);
                            else
                                display('#######################################################');
                                display('#### New plane less than 10x larger, using smaller area');
                                display(['#### planenumber = ' num2str(planenumber) ]);
                                display('#######################################################');
                            end
                        else
                            display('#######################################################');
                            display('#### New area is further and smaller - ignore it for now');
                            display('#######################################################');
                        end
                        
                        
                        % 05-28-2019 Changed to correct for the wrong plane chosen and not correcting the centroid calculation
                      %% %% %%  centroid{side} = allCentroidPos{side,1} * size(growingnodes{side},1);%%%%%%%%%%%%%%%%%%%check!!!
                      %% %% %%  numberOfNodes{side} = size(growingnodes{side},1);
                   %   display(['*** planenumber = ' num2str(planenumber) ', Ignoring one area - too far away']);
                    end
                    
                    
%                     %% TEST TURNED OFF ONLY %%%%%%%%%%%%%
%                     if isempty(area{side}) || areaTemp > area{side}(1)
%                         display(['*** planenumber = ' num2str(planenumber) ', Losing one area']);
%                         area{side}(1) = areaTemp;
%                         growingnodes{side} = nodes;
%                         tri{side} = triangles;
%                     end
%                     %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                      if abs(norm(allCentroidPos{side,j} - PlanePoint)) > 
%                          
%                      end

                    
                end
            end
        end
    end
    
    
    
% % %     for side = 1:3
% % %         if notCheckedNodes(side) > 0 && size(growingnodes{side},1) > notCheckedNodes(side)
% % %             in = inpolygon(growingnodes{side}(1:notCheckedNodes(side),1),growingnodes{side}(1:notCheckedNodes(side),2),growingnodes{side}(notCheckedNodes(side)+1:end,1),growingnodes{side}(notCheckedNodes(side)+1:end,2));
% % %             if sum(in) > 0
% % %                 % We've found doughnuts.
% % %                 
% % %                 figure;
% % %                 patch('Vertices',growingnodes{side},'Faces',tri{side},'facecolor', 'g')
% % %                 hold on
% % %                 %plot(growingnodes{side}(1:notCheckedNodes(side),1),growingnodes{side}(1:notCheckedNodes(side),2), 'bx-')
% % %                 %plot(growingnodes{side}(notCheckedNodes(side)+1:end,1),growingnodes{side}(notCheckedNodes(side)+1:end,2),'rx-')
% % %                 plot(growingnodes{side}(1:notCheckedNodes(side)-1,1),growingnodes{side}(1:notCheckedNodes(side)-1,2), 'cx-')
% % %                 plot(growingnodes{side}(notCheckedNodes(side):end,1),growingnodes{side}(notCheckedNodes(side):end,2),'mx-')
% % %                 warning('Mmmm, doughnuts');
% % %             end
% % %         end
% % %     end

        
    for side = 1:3
        centroid{side} = sum(centroid{side},1)/sum(numberOfNodes{side});
    end
%     plot3(centroid(1),centroid(2),centroid(3),'k+')
%     xlabel('X');
%     ylabel('Y');
%     zlabel('Z');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Split lobes into 3 different sections: upper, lower and split %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SplitElements, SplitNodes, WithNodeNumber] = SplitLobes(Fsurface, Vsurface, PlanePoint, PlaneNormal)
    %Uncomment if you want to see the original stl
    %patch('Vertices',Vsurface,'Faces',Fsurface,'facecolor', 'r')

    SplitElements=[];
    SplitNodes=[];
    WithNodeNumber=[];
    %BottomSurface=[];
    %TopSurface=[];
%     maxSize = length(Fsurface);
%     SplitElements=zeros(maxSize,3);
%     SplitNodes=zeros(maxSize,3);
%     WithNodeNumber=zeros(maxSize,4);
    
%     
%     % Need to set up a point on the plane and the normal
%     d=1*(PlanePoint(1)*PlaneNormal(1)+PlanePoint(2)*PlaneNormal(2)+PlanePoint(3)*PlaneNormal(3));
%     [X,Y] = meshgrid(-2:0.1:2);
% 
%     Z = (d-PlaneNormal(1)*X - PlaneNormal(2)*Y)/PlaneNormal(3);

    FsurfaceSize = length(Fsurface);
    
    for i=1:FsurfaceSize

        % determine whether the dot product of each vertex of a face is
        % negative or positive
        a = Vsurface(Fsurface(i,1),:) -PlanePoint;
        b = Vsurface(Fsurface(i,2),:) -PlanePoint;
        c = Vsurface(Fsurface(i,3),:)-PlanePoint;
        
        m=[PlaneNormal;PlaneNormal;PlaneNormal];
        v=[a;b;c]';
        t=m*v;
        
        NodeTest1=t(1,1);
        NodeTest2=t(1,2);
        NodeTest3=t(1,3);
        %NodeTest1=dot(PlaneNormal,a);
        %NodeTest2=dot(PlaneNormal,b);
        %NodeTest3=dot(PlaneNormal,c);

        %Calculate distance between point and plane
        D1=abs(PlaneNormal(1)*(PlanePoint(1)-Vsurface(Fsurface(i,1),1))+PlaneNormal(2)*(PlanePoint(2)-Vsurface(Fsurface(i,1),2))+PlaneNormal(3)*(PlanePoint(3)-Vsurface(Fsurface(i,1),3)));
        D2=abs(PlaneNormal(1)*(PlanePoint(1)-Vsurface(Fsurface(i,2),1))+PlaneNormal(2)*(PlanePoint(2)-Vsurface(Fsurface(i,2),2))+PlaneNormal(3)*(PlanePoint(3)-Vsurface(Fsurface(i,2),3)));
        D3=abs(PlaneNormal(1)*(PlanePoint(1)-Vsurface(Fsurface(i,3),1))+PlaneNormal(2)*(PlanePoint(2)-Vsurface(Fsurface(i,3),2))+PlaneNormal(3)*(PlanePoint(3)-Vsurface(Fsurface(i,3),3)));

        
%         %TEST ONLY        
%         a = 1e-16;
%         if D1 < a || D2 < a || D3 < a
%            display('here'); 
%         end
        
        % All possible scenarios of where points lie. Note: if the distance is
        % less than 1e-16 then it is assumed to be on the same side as all 
        % other nodes. If we hit the "else" then we have located an element
        % that is split by the plane and it is added to the "split element
        % array".
        
        % Change Oct 2014
        % Triangles were missing where the cut hit a node directly,
        % therefore in this case we now return the triangle below this node
        % (arbitrarily) so we have the correct number of triangles for
        % ordering etc.
        
        if(NodeTest1 < 0 && NodeTest2 < 0 && NodeTest3 < 0)
            %BottomSurface(end+1,:)=Fsurface(i,:);
        elseif(NodeTest1 > 0 && NodeTest2 > 0 && NodeTest3 > 0)
            %TopSurface(end+1,:)=Fsurface(i,:);
        elseif (D1 < 1e-16 && (NodeTest2>0 && NodeTest3 >0))
            %TopSurface(end+1,:)=Fsurface(i,:);
        elseif (D2 < 1e-16 && (NodeTest1>0 && NodeTest3 >0))
            %TopSurface(end+1,:)=Fsurface(i,:);
        elseif (D3 < 1e-16 && (NodeTest1>0 && NodeTest2 >0))
            %TopSurface(end+1,:)=Fsurface(i,:);
% % % %         elseif (D1 < 1e-16 && (NodeTest2<0 && NodeTest3 <0))
% % % %             %BottomSurface(end+1,:)=Fsurface(i,:);
% % % %         elseif (D2 < 1e-16 && (NodeTest1<0 && NodeTest3 <0))
% % % %             %BottomSurface(end+1,:)=Fsurface(i,:);
% % % %         elseif (D3 < 1e-16 && (NodeTest1<0 && NodeTest2 <0))
            %BottomSurface(end+1,:)=Fsurface(i,:);
        else

            locateFirstNode=find(SplitElements==Fsurface(i,1));

            if (isempty(locateFirstNode))
                FillMatrixValues=[Fsurface(i,1) Vsurface(Fsurface(i,1),:)];
                SplitNodes(end+1,:) = Vsurface(Fsurface(i,1),:);
                WithNodeNumber(end+1,:)=FillMatrixValues;
%                 SplitNodes(i,:) = Vsurface(Fsurface(i,1),:);
%                 WithNodeNumber(i,:)=FillMatrixValues;
            end

            locateSecondNode=find(SplitElements==Fsurface(i,2));
            if (isempty(locateSecondNode))
                FillMatrixValues=[Fsurface(i,2) Vsurface(Fsurface(i,2),:)];
                SplitNodes(end+1,:) = Vsurface(Fsurface(i,2),:);
                WithNodeNumber(end+1,:)=FillMatrixValues;
%                 SplitNodes(i,:) = Vsurface(Fsurface(i,2),:);
%                 WithNodeNumber(i,:)=FillMatrixValues;
            end

            locateThirdNode=find(SplitElements==Fsurface(i,3));
            if (isempty(locateThirdNode))
                FillMatrixValues=[Fsurface(i,3) Vsurface(Fsurface(i,3),:)];
                SplitNodes(end+1,:) = Vsurface(Fsurface(i,3),:);
                WithNodeNumber(end+1,:)=FillMatrixValues;
%                 SplitNodes(i,:) = Vsurface(Fsurface(i,3),:);
%                 WithNodeNumber(i,:)=FillMatrixValues;
            end
            SplitElements(end+1,:)=Fsurface(i,:);
            %SplitElements(i,:)=Fsurface(i,:);
        end
    end
    
    %Remove blank Rows from SplitElements, SplitNodes and WithNodeNumber
    %SplitElements = SplitElements(~all(SplitElements==0,2),:);
    %SplitNodes = SplitNodes(~all(SplitNodes==0,2),:);
    %WithNodeNumber = WithNodeNumber(~all(WithNodeNumber==0,2),:);
    
 %   patch('Vertices',Vsurface,'Faces',TopSurface,'facecolor', 'b') 
 %   patch('Vertices',Vsurface,'Faces',BottomSurface,'facecolor', 'g')
 %figure;
 %hold on
 %    patch('Vertices',Vsurface,'Faces',SplitElements,'facecolor', 'g')
end

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%%%%%%%% This section detects connected elements %%%%%%%%%%%%%%%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function SplitElements = DetectConnectedElements(SplitElements)
     [m,n] = find(SplitElements(1,1)==SplitElements(:,1:3));
     
     CurrentConnection = 1;
     ConnectedElements = [];
    %  SplitElements(1,4) = 1;
     SplitElements(m,4) = CurrentConnection;
     for i = 1:size(m)
         if n(i) ~= 1
            ConnectedElements(end+1) = SplitElements(m(i),1);
         end
         if n(i) ~= 2
            ConnectedElements(end+1) = SplitElements(m(i),2);
         end
         if n(i) ~= 3
            ConnectedElements(end+1) = SplitElements(m(i),3);
         end
     end
     while find(SplitElements(:,4)==0)
        while size(ConnectedElements > 0)
            [m,n] = find(ConnectedElements(end)==SplitElements(:,1:3));
            ConnectedElements = ConnectedElements(1:end - 1);
            for i = 1:size(m)
                if SplitElements(m(i),4) == 0
                    SplitElements(m(i),4) = CurrentConnection;
                     if n(i) ~= 1
                        ConnectedElements(end+1) = SplitElements(m(i),1);
                     end
                     if n(i) ~= 2
                        ConnectedElements(end+1) = SplitElements(m(i),2);
                     end
                     if n(i) ~= 3
                        ConnectedElements(end+1) = SplitElements(m(i),3);
                     end
                end
            end
        end
        CurrentConnection = CurrentConnection + 1;
        ConnectedElements = SplitElements(find(SplitElements(:,4)==0,1),1);
     end
 end


 function SplitElements = DetectTwoNodeConnectedElements(SplitElements)
  
    CurrentConnection = 1;
    rows = 1;
    ConnectedRows = [];
    
    SplitElements(1,4) = 1;

    while find(SplitElements(:,4)==0)
        while ~isempty(rows)
            [rows,~] = find(sum(ismember(SplitElements(:,1:3), SplitElements(rows,1:3)),2) >= 2);
            rows = rows(~ismember(rows,ConnectedRows));
            SplitElements(rows,4) = CurrentConnection;
            ConnectedRows = [ConnectedRows; rows];    
        end
        CurrentConnection = CurrentConnection + 1;
        [rows,~] = find(SplitElements(:,4)==0,1);
    end
    



  
% % %      [m,n] = find(SplitElements(1,1)==SplitElements(:,1:3));
% % %      [m2,n2] = find(SplitElements(1,2)==SplitElements(:,1:3));
% % %      [m3,n3] = find(SplitElements(1,3)==SplitElements(:,1:3));
% % %      
% % %      m= m(ismember(m,m2) | ismember(m,m3));
% % %      n= n(ismember(m,m2) | ismember(m,m3));
% % %      
% % %      
% % %      CurrentConnection = 1;
% % %      ConnectedElements = [];
% % %      FoundElements = [];
% % %     %  SplitElements(1,4) = 1;
% % %      SplitElements(m,4) = CurrentConnection;
% % %      for i = 1:size(m)
% % %          if n(i) ~= 1
% % %             ConnectedElements(end+1) = SplitElements(m(i),1);
% % %             FoundElements(end+1) = SplitElements(m(i),1); 
% % %          end
% % %          if n(i) ~= 2
% % %             ConnectedElements(end+1) = SplitElements(m(i),2);
% % %             FoundElements(end+1) = SplitElements(m(i),2);
% % %          end
% % %          if n(i) ~= 3
% % %             ConnectedElements(end+1) = SplitElements(m(i),3);
% % %             FoundElements(end+1) = SplitElements(m(i),3);
% % %          end
% % %      end
% % %      while find(SplitElements(:,4)==0)
% % %         while size(ConnectedElements > 0)
% % %             [m,n] = find(ConnectedElements(end)==SplitElements(:,1:3));
% % %             % Check that each connection connects to at least two points
% % %             if size (FoundElements,2) > 1
% % %                 m = m(logical(sum(ismember(SplitElements(m,1:3),FoundElements(1:end-1)),2)));
% % %                 n = n(logical(sum(ismember(SplitElements(m,1:3),FoundElements(1:end-1)),2)));
% % %             else
% % %                 [m2,n2] = find(SplitElements(ConnectedElements,2)==SplitElements(:,1:3));
% % %                 [m3,n3] = find(SplitElements(ConnectedElements,3)==SplitElements(:,1:3));
% % %                 m = m(ismember(m,m2) | ismember(m,m3));
% % %                 n = n(ismember(m,m2) | ismember(m,m3));
% % %             end
% % %             ConnectedElements = ConnectedElements(1:end - 1);
% % %             for i = 1:size(m)
% % %                 if SplitElements(m(i),4) == 0
% % %                     SplitElements(m(i),4) = CurrentConnection;
% % %                      if n(i) ~= 1
% % %                         ConnectedElements(end+1) = SplitElements(m(i),1);
% % %                         FoundElements(end+1) = SplitElements(m(i),1);
% % %                      end
% % %                      if n(i) ~= 2
% % %                         ConnectedElements(end+1) = SplitElements(m(i),2);
% % %                         FoundElements(end+1) = SplitElements(m(i),2);
% % %                      end
% % %                      if n(i) ~= 3
% % %                         ConnectedElements(end+1) = SplitElements(m(i),3);
% % %                         FoundElements(end+1) = SplitElements(m(i),3);
% % %                      end
% % %                 end
% % %             end
% % %             ConnectedElements = unique(ConnectedElements,'stable');
% % %             FoundElements = unique(FoundElements,'stable');
% % %         end
% % %         CurrentConnection = CurrentConnection + 1;
% % %         %ConnectedElements = SplitElements(find(SplitElements(:,4)==0,1),1);
% % %         [row,~] = find(SplitElements(:,4)==0,1);
% % %         ConnectedElements = SplitElements(row,1:3);
% % %         FoundElements = ConnectedElements;
% % %      end
 end


 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Check for dougnuts - Holes within planes                 %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [doughnuts] = CheckForDoughnuts(allNodes)
% Check for doughnuts
doughnuts = cell(3,3);

% doughnuts cell
% 1st column is matches that have already been detected on the same side.
% 2nd column is matches where the outside of the doughnut is side 1 and the
%   inside is side 2
% 3rd column is matches where the outside of the doughnut is side 2 and the
%   inside is side 1



for side = 1:3
    for otherSide = 1:3
        for i = 1:size(allNodes,2)
            if ~isempty(allNodes{side,i})
                for j = 1:size(allNodes,2)
                    if i ~= j || side ~= otherSide % Don't test the same cut with itself
                        if ~isempty(allNodes{otherSide,j})
                            
%                             figure
%                             plot(allNodes{side,i}(:,1),allNodes{side,i}(:,2),'bx')
%                             hold on
%                             plot(allNodes{otherSide,j}(:,1),allNodes{otherSide,j}(:,2),'rx')
                            
                            in = inpolygon(allNodes{side,i}(:,1),allNodes{side,i}(:,2),allNodes{otherSide,j}(:,1),allNodes{otherSide,j}(:,2));
                            if sum(in) > 0
                                
%                             figure
%                             plot(allNodes{side,i}(:,1),allNodes{side,i}(:,2),'bx')
%                             hold on
%                             plot(allNodes{otherSide,j}(:,1),allNodes{otherSide,j}(:,2),'rx')
                                
                                % We've found doughnuts.
                                % 
    %                             figure;
    %                             plot(allNodes{side,i}(:,1),allNodes{side,i}(:,2),'bx')
    %                             hold on
    %                             plot(allNodes{side,j}(:,1),allNodes{side,j}(:,2),'rx')
                                %warning('Mmmm, doughnuts');
                                if side == otherSide
                                    if i < j
                                        doughnuts{side}(end+1,:) = [i,j];
                                    else
                                        doughnuts{side}(end+1,:) = [j,i];
                                    end
                                else
                                    %Doughnuts of different sides...
                                    if i < j
                                        doughnuts{side,otherSide + 1}(end+1,:) = [i,j];
                                    elseif j < i
                                        doughnuts{side,otherSide + 1}(end+1,:) = [j,i];
                                    else
                                        doughnuts{side,otherSide + 1}(end+1,:) = [i,j]; % Does i,j order matter here
                                    end
                                    
                                    
%                                     figure
%                                     plot(allNodes{side,i}(:,1),allNodes{side,i}(:,2),'bx')
%                                     hold on
%                                     plot(allNodes{otherSide,j}(:,1),allNodes{otherSide,j}(:,2),'rx')
                                end
%                             else
%                                 figure
%                                 plot(allNodes{side,i}(:,1),allNodes{side,i}(:,2),'gx')
%                                 hold on
%                                 plot(allNodes{otherSide,j}(:,1),allNodes{otherSide,j}(:,2),'kx')
                            end
                        end
                    end
                end
            end
        end
    end
end
    
%     for i = 1:max(SplitElements(:,4))
%         TestElements = SplitElements(:,4) == i;
%         TestElements = SplitElements(TestElements,1:3);
%         for j = 1:max(SplitElements(:,4))
%             if i ~= j
%                 OtherElements = SplitElements(:,4) == j;
%                 OtherElements = SplitElements(OtherElements,1:3);
%                 in = inpolygon(Vsurface(TestElements,1),Vsurface(TestElements,2),Vsurface(OtherElements,1),Vsurface(OtherElements,2));
%                 if sum(in) > 0
%                     % We've found doughnuts.
% 
%                     plot(Vsurface(TestElements,1),Vsurface(TestElements,2),'bx')
%                     hold on
%                     plot(Vsurface(OtherElements,1),Vsurface(OtherElements,2),'rx')
%                     warning('Mmmm, doughnuts');                    
%                 end
%             end
%         end
%     end
% end


end
 
 
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%% Find point of line which defines a split triangle %%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [I] = FindPoints(WithNodeNumber, SplitElements, SplitNodes, PlanePoint, PlaneNormal)

     %This function determines the points at which the plane intersects each
     %triangle in the array "SplitElement" and returns thes points in the array
     %"I"
     I=[];
     NewSplitElements=[];
     
     if isempty(WithNodeNumber) || isempty(SplitElements)
         return
     end
     
     for i=1:size(SplitElements,1)

         %Match element nodes with coordinates
         Element2Node1=find(WithNodeNumber(:,1)==SplitElements(i,1));
         Element2Node2=find(WithNodeNumber(:,1)==SplitElements(i,2));
         Element2Node3=find(WithNodeNumber(:,1)==SplitElements(i,3));

         NewSplitElements(end+1,:)=[Element2Node1 Element2Node2 Element2Node3 SplitElements(i,1) SplitElements(i,2) SplitElements(i,3)];

         NodeTest1=sign(dot(PlaneNormal,SplitNodes(Element2Node1,:)-PlanePoint));
         NodeTest2=sign(dot(PlaneNormal,SplitNodes(Element2Node2,:)-PlanePoint));
         NodeTest3=sign(dot(PlaneNormal,SplitNodes(Element2Node3,:)-PlanePoint));
         
         c=PlanePoint;
         
         if(NodeTest1 ~= NodeTest2)
             A1=SplitNodes(Element2Node1,:);
             A2=SplitNodes(Element2Node2,:);
             Iba = A2-A1;

             t=(dot(c-A1,PlaneNormal))/dot(Iba,PlaneNormal);
             MidNode= A1+ t.*Iba;
             Dist1=sqrt(power(MidNode(1,1)-A1(1,1),2)+power(MidNode(1,2)-A1(1,2),2)+power(MidNode(1,3)-A1(1,3),2));
             Dist2=sqrt(power(MidNode(1,1)-A2(1,1),2)+power(MidNode(1,2)-A2(1,2),2)+power(MidNode(1,3)-A2(1,3),2));
             I(end+1,:) = [MidNode Dist1 Dist2 Element2Node1 Element2Node2];
         end
         if(sign(NodeTest2) ~= sign(NodeTest3))
             A1=SplitNodes(Element2Node2,:);
             A2=SplitNodes(Element2Node3,:);
             Iba = A2-A1;
             t=(dot(c-A1,PlaneNormal))/dot(Iba,PlaneNormal);
             MidNode= A1+ t.*Iba;
             Dist1=sqrt(power(MidNode(1,1)-A1(1,1),2)+power(MidNode(1,2)-A1(1,2),2)+power(MidNode(1,3)-A1(1,3),2));
             Dist2=sqrt(power(MidNode(1,1)-A2(1,1),2)+power(MidNode(1,2)-A2(1,2),2)+power(MidNode(1,3)-A2(1,3),2));
             I(end+1,:) = [MidNode Dist1 Dist2 Element2Node2 Element2Node3];
         end
         if(sign(NodeTest1) ~= sign(NodeTest3))
             A1=SplitNodes(Element2Node1,:);
             A2=SplitNodes(Element2Node3,:);
             Iba = A2-A1;
             t=(dot(c-A1,PlaneNormal))/dot(Iba,PlaneNormal);
             MidNode= A1+ t.*Iba;
             Dist1=sqrt(power(MidNode(1,1)-A1(1,1),2)+power(MidNode(1,2)-A1(1,2),2)+power(MidNode(1,3)-A1(1,3),2));
             Dist2=sqrt(power(MidNode(1,1)-A2(1,1),2)+power(MidNode(1,2)-A2(1,2),2)+power(MidNode(1,3)-A2(1,3),2));
             I(end+1,:) = [MidNode Dist1 Dist2 Element2Node1 Element2Node3];
         end
     end
 end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%% This section deletes duplicated nodes %%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [NewVect] = DeleteDuplicateNodes(I)
     %Since each triangle is connected to another triangle the points of
     %intersection are not unique, the following removes duplicated nodes so
     %the list of intersection points is unique.
     %I;
%      numberToSearch = I(1,7);
%      nodes = find(numberToSearch==I(:,6:7));
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     write=1;
     NewVect=[];
%      OrderVec=[];
    ISize = size(I,1);
    ISize1 = ISize - 1;

     for i=1:ISize1
         TestVector=I(i,:);
         tempindex=i+1;
         for j=tempindex:ISize
             TestFirst=TestVector(1,6)-I(j,7);
             TestSecond=TestVector(1,7)-I(j,6);
             TestThird=TestVector(1,6)-I(j,6);
             TestForth=TestVector(1,7)-I(j,7);
             if(TestFirst==0 && TestSecond==0 || TestThird==0 && TestForth==0)
                 write=0;
% % %                  if j == length(I)
% % %                      write = 1; %Otherwise last element missed out
% % %                  end
                 % If we're here it means we've matched two elements record
                 % what they're attached to.
                 %OrderVec(end+1,:) = [i,j];
             end
         end
         if (write~=0)
             NewVect(end+1,:)=TestVector;
         end
         write=1;

     end

     NewVect(end+1,:)=I(end,:);
     %plot3(NewVect(:,1),NewVect(:,2),NewVect(:,3),'or')
     
     %%% TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      for i=1:length(OrderVec)
%         testVect(i,:) = I(OrderVec(i),:);
%      end
%      NewVect = testVect;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 end
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%% This section meshes the plane %%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [area, rotNodes, tri] = NewMeshPlane(NewVect,SplitPos, PlaneNormal)
 % As MeshPlane, but maps the plane nodes into 2D using plane normal,
 % meshes, then rotates the resulting triangulation back.
 
nodes = NewVect(:,1:3);
nodes = unique(nodes,'rows','stable');


 
 %Calculate angle to rotate and the vector to rotate about

 zaxis = [ 0 0 1];
 if sum(zaxis == PlaneNormal) ~= 3
    theta = acosd(dot(PlaneNormal, zaxis)/(norm(PlaneNormal)*norm(zaxis)));
    trans = cross(PlaneNormal, zaxis);
 
    if theta ~= 0 && theta ~= 180
        [rotNodes, R, t] = AxelRot(nodes', theta, trans, [0 0 0]);
        rotNodes = rotNodes';
    else
        rotNodes = nodes;
    end
    
 else
     rotNodes = nodes;
 end
 
 %Check we have a 2D plane (z axis should be equal)
 if max(rotNodes(:,3))-min(rotNodes(:,3)) > 1e-10
     warning('Plane no longer planar when rotated');
 end
 zcoord = mean(rotNodes(:,3));
 
constrainedEdges = [];
lastStop = 1;
%We need to define contrained edges so for cases involving doughnuts

 
 for i = 1:size(rotNodes,1)
    if ~ismember(i,SplitPos) && i < size(rotNodes,1)
        constrainedEdges(end+1,1) = i;
        constrainedEdges(end,2) = i + 1;
    else
        constrainedEdges(end+1,1) = i;
        constrainedEdges(end,2) = lastStop;
        lastStop = i + 1;
    end
end
 
 
 
% tri=DelaunayTri(rotNodes(:,1),rotNodes(:,2), constrainedEdges);
meshed = false;
while ~meshed
    try
        tri = delaunayTriangulation(rotNodes(:,1),rotNodes(:,2), constrainedEdges);
        meshed = true;
    catch
        %warning('Delaunay failed: Trying again');
        difs = diff(rotNodes);
        sumDifs = sum(abs(difs),2);
        [minDif, minPos] = min(sumDifs);
        rotNodes(minPos+1, :) = [];

        %TODO: Hack constrained edges for now - should be properly handled
        constrainedEdges(end,:) = [];
        constrainedEdges(end,2) = 1;
    end
end
% TODO : Handle cases when this returns an unusable triangulation (area =
% 0) and isInterior sums to 0

 %rotNodes  = tri.X;
 rotNodes = tri.Points;
 rotNodes(:,3) = zcoord;

%isInside = tri.inOutStatus();
isInside = tri.isInterior;

if sum(isInside) == 0
   warning('Triangluation may be wrong') ;
   tri = tri(:,:);
   % I should really work out why this happens
else
    tri = tri(isInside,:); 
end

 % Calculate area with cross product of two vectors from X node to
 % other nodes
 area = 0.0;
 for i = 1:size(tri,1)
     XY = rotNodes(tri(i,2),:) - rotNodes(tri(i,1),:);
     XZ = rotNodes(tri(i,3),:) - rotNodes(tri(i,1),:);
     area = area + (norm(cross(XZ,XY)) / 2);
 end
 
% Rotate the planes back to the original coordinate system so we can save
% them etc.
if sum(zaxis == PlaneNormal) ~= 3
  [rotNodes, R, t] = AxelRot(rotNodes', theta, -1 * trans, [0 0 0]);
  rotNodes = rotNodes';
end

 end 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%% This section meshes the plane %%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [area, growingnodes, trifaces] = MeshPlane(NewVect,SplitPos)
 %warning('off','MATLAB:DelaunayTri:ConsConsSplitWarnId')

growingnodes = NewVect(:,1:3);
growingnodes = unique(growingnodes,'rows','stable');

if size(growingnodes,1) < 3
    area = 0.0;
    growingnodes= [];
    trifaces=[];
    return
end

% We mesh in 2D - so need make sure the plane we choose in 3D doesn't collapse to
% a line in 2D.
difX = max(growingnodes(:,1))-min(growingnodes(:,1));
difY = max(growingnodes(:,2))-min(growingnodes(:,2));
difZ = max(growingnodes(:,3))-min(growingnodes(:,3));


constrainedEdges = [];
lastStop = 1;
%We need to define contrained edges so for cases involving doughnuts


for i = 1:size(growingnodes,1)
    if ~ismember(i,SplitPos) && i < size(growingnodes,1)
        constrainedEdges(end+1,1) = i;
        constrainedEdges(end,2) = i + 1;
    else
        constrainedEdges(end+1,1) = i;
        constrainedEdges(end,2) = lastStop;
        lastStop = i + 1;
    end
end


 
 noX = false;
 noY = false;
 noZ = false;
 
 if difY ~= 0
    tri=DelaunayTri(growingnodes(:,1),growingnodes(:,2), constrainedEdges);
 else
     noY = true;
 end
 
 if difZ ~= 0
    tri2=DelaunayTri(growingnodes(:,1),growingnodes(:,3), constrainedEdges);
 else
    noZ = true; 
 end
 
 if difX ~= 0
    % This should always be true...
 else
    noX = true;
    warning('No X variation - this is not currently handled');
 end
 
 if ~noY && ~noZ
    growingnodes = tri.X;
    temp = tri2.X;
    %DEBUG
    if size(growingnodes,1) ~= size(temp,1)
        for i = 1:size(growingnodes,1)
            %for j = 1:size(temp,1)
            if i <= size(temp,1)
                if growingnodes(i,1) == temp(i,1)
                   growingnodes(i,3) = temp(i,2);
                else
                    warning('Unmatched Nodes');
                    if size(growingnodes,1) > size(temp,1)
                        growingnodes(i,:) = [];
                    elseif size(growingnodes,1) < size(temp,1)
                        temp(i,:) = [];
                    else
                        growingnodes(i,:) = [];
                        temp(i,:) = [];
                    end
                    i = i - 1;
                   %break;
                end
            else
                if i < size(growingnodes,1) 
                    warning('Unmatched Nodes');
                end
            end
            %end
        end
    else
        growingnodes(:,3) = temp(:,2);
    end
 elseif noY
    y = growingnodes(1,2);
    growingnodes(:,1) = tri2.X(:,1);
    growingnodes(:,3) = tri2.X(:,2);
    growingnodes(:,2) = y;
 elseif noZ
    z = growingnodes(1,3);
    growingnodes = tri.X;
    growingnodes(:,3) = z;
 else
    area = 0.0;
    growingnodes= [];
    trifaces=[];
    return
 end
 
 
 isInside = tri.inOutStatus();
 tri = tri(isInside,:); 
 
 
 
 
 trifaces=tri;


     % Calculate area with cross product of two vectors from X node to
     % other nodes
     area = 0.0;
     for i = 1:size(tri,1)
         XY = growingnodes(tri(i,2),:) - growingnodes(tri(i,1),:);
         XZ = growingnodes(tri(i,3),:) - growingnodes(tri(i,1),:);
         area = area + (norm(cross(XZ,XY)) / 2);
     end
     
     warning('on','MATLAB:DelaunayTri:ConsConsSplitWarnId')
 end
 
 
 function [OrderedVec, Ordered] = NewOrderSplitElements(SplitElements, planenumber,Vsurface)
 %Vsurface only for debug
    AlternativeRoute={};
    BadMeshRoute={};
    OrderedVec = [];
    MatchedRows = [];
    ForkElements = [];
    %Choose first triangle
    RowToMatch = 1;
    Columns = [1 2 3];
    ColumnToMatch = 1;
    iteration = 0;
    restart = false;
    restartNo = 0;
    debugCount = 0;
    Ordered = true;
    loopNumber = 0;

    while true
    
        debugCount = debugCount + 1;
%         display(num2str(debugCount));
        
        iteration = iteration + 1;
        ColumnToFind = Columns(Columns ~= ColumnToMatch);
        %Add triangle to "used list"
        OrderedVec = [OrderedVec; SplitElements(RowToMatch,:)];
        MatchedRows = [MatchedRows; RowToMatch];

        %Find triangles that share at least two of these nodes
        [r1,c1] = find(SplitElements == SplitElements(RowToMatch, ColumnToMatch));
        [r2,c2] = find(SplitElements == SplitElements(RowToMatch, ColumnToFind(1)));
        [r3,c3] = find(SplitElements == SplitElements(RowToMatch, ColumnToFind(2)));
        
%         r1r2 = ismember(r1,r2);
%         r1r3 = ismember(r1,r3);
%         sharedRow = r1r2 | r1r3;
% 
%         PossibleRow = r1(sharedRow);
%         rowTest = ismember(PossibleRow,MatchedRows);
%         PossibleRow = PossibleRow(~rowTest);

        r1r2 = fastintersect(r1,r2);
        r1r3 = fastintersect(r1,r3);
        PossibleRow = [r1r2; r1r3];
        PossibleRow = setdiff_ns(PossibleRow, MatchedRows);
        
        
        
        if length(PossibleRow) == 2 
            % We have alternative Routes we could take
            % I assume from each triangle there are only at most 2 alternative
            % routes (you can't go back on yourself. We'll store
            % the second option, and use the first option for now.
            AlternativeRoute{end + 1, 1} = PossibleRow(2);
            %AlternativeRoute{end, 2} = find(~ismember(SplitElements(PossibleRow(2),:),OrderedVec(end,:)));
            AlternativeRoute{end, 3} = MatchedRows;
            AlternativeRoute{end, 4} = OrderedVec;
            ForkElements(end + 1, :) = RowToMatch;
            ForkElements = unique(ForkElements);
        elseif length(PossibleRow) == 3 && iteration == 1
            % The first cell could possibly be joined to three others.
            % Subsequent ones can't as they can't go back to where they
            % came from, so only have two alternate paths.
            % Store first alternate route
            AlternativeRoute{end + 1, 1} = PossibleRow(2);
            %AlternativeRoute{end, 2} = find(~ismember(SplitElements(PossibleRow(2),:),OrderedVec(end,:)));
            AlternativeRoute{end, 3} = MatchedRows;
            AlternativeRoute{end, 4} = OrderedVec;
            % Store second alternate route
            AlternativeRoute{end + 1, 1} = PossibleRow(3);
            %AlternativeRoute{end, 2} = find(~ismember(SplitElements(PossibleRow(2),:),OrderedVec(end,:)));
            AlternativeRoute{end, 3} = MatchedRows;
            AlternativeRoute{end, 4} = OrderedVec;
        elseif length(PossibleRow) > 2
            %If here, the surface is non manifold. A triangle can only
            %possibly join on to 3 others (where it came from plus two
            %alternate directions.) The case of the first triangle handled
            %above.
            warning('Too many possible Routes found -- non manifold surface?');
%             figure;
            %patch('Vertices',Vsurface,'Faces',SplitElements,'facecolor', 'g');
            return;
        elseif isempty(PossibleRow)
            % We have run out of places to go - the end of the current loop
            % attempt. Check to see if this route matches the two finishing
            % criteria
            if size(OrderedVec,1) == size(SplitElements,1)
                % All elements used in loop, now check second criteria that
                % the last element connects to the first.
                [rt1,ct1] = find(SplitElements(1,:) == SplitElements(RowToMatch, ColumnToMatch));
                [rt2,ct2] = find(SplitElements(1,:) == SplitElements(RowToMatch, ColumnToFind(1)));
                [rt3,ct3] = find(SplitElements(1,:) == SplitElements(RowToMatch, ColumnToFind(2)));
                testSize = length(rt1) + length(rt2) + length (rt3);
                if testSize == 2
                    % We've found a complete route round the loop
                    break
                else
                    loopNumber = loopNumber + 1;
                    display(['Plane number ' num2str(planenumber) ', Loop number ' num2str(loopNumber) ' Complete Loop Found, but cant return to start, continuing']);
                    if loopNumber >= 10
                        AlternativeRoute = {};
                        BadMeshRoute = {};
                        loopNumber = 0;
                    end
                end
            end
            if ~isempty(AlternativeRoute)
                % The loop we were on has failed one of the end criteria. Start
                % searching again from the last saved Alternative route.
                PossibleRow = AlternativeRoute{end,1};
                %PossibleColumn = AlternativeRoute{end,2};
                MatchedRows = AlternativeRoute{end,3};
                OrderedVec =  AlternativeRoute{end,4};
                % Store this alternative route as a point to come back to
                % if things really don't work - ie. a bad mesh...
                BadMeshRoute(end+1,:) = AlternativeRoute(end, :);
                % Remove this option we're trying so we don't go down this
                % route again, until we're in the bad mesh handling stage
                AlternativeRoute(end, :) = [];
            elseif ~isempty(BadMeshRoute)
                %display('Unable to find path round loop using each triangle only once');
                %figure;
                %display(iteration);
                %patch('Vertices',Vsurface,'Faces',SplitElements,'facecolor', 'g');
                %patch('Vertices',Vsurface,'Faces',SplitElements(ForkElements,:),'facecolor', 'r');
                %We need to handle these poor meshes
                PossibleRow = BadMeshRoute{end,1};
                MatchedRows = BadMeshRoute{end,3};
                OrderedVec =  BadMeshRoute{end,4};
                % Remove fork elements from the list, so we may be able
                % to use them as crossing points.
                MatchedRows = MatchedRows(~ismember(MatchedRows,ForkElements));
                BadMeshRoute(end, :) = [];
                %return;
            elseif size(OrderedVec,1) > size(SplitElements,1)
                display('Returning Loop with repeated elements');
                return;
            else
                restartNo = restartNo + 1;
%                 if restartNo == 11
                if restartNo == 1
%                     figure;
%                     patch('Vertices',Vsurface,'Faces',SplitElements,'facecolor', 'g');
%                     patch('Vertices',Vsurface,'Faces',SplitElements(ForkElements,:),'facecolor', 'r');
                    warning('Unable to order elements - Returned area may not be correct');
                    Ordered = false;
                    return
                end
                display(['Plane number ' num2str(planenumber) ' Unable to find loop - trying different starting point ' num2str(restartNo)]);
                %display('panic');
                %figure;
                %patch('Vertices',Vsurface,'Faces',SplitElements,'facecolor', 'g');
                %patch('Vertices',Vsurface,'Faces',SplitElements(ForkElements,:),'facecolor', 'r');
                %%%Maybe try starting from the other end of split elements
                %%%instead?
                AlternativeRoute={};
                BadMeshRoute={};
                OrderedVec = [];
                MatchedRows = [];
                ForkElements = [];
                %Choose first triangle
                RowToMatch = floor((restartNo/10) *size(SplitElements,1));
                Columns = [1 2 3];
                ColumnToMatch = 1;
                iteration = 0;
                restart = true;

                
                %return
            end
% % % % % %         elseif length(PossibleRow) == 1 && isempty(AlternativeRoute) && iteration == 1
% % % % % %             display([ 'here' num2str(planenumber) ' ' num2str(debugCount)]);
        else
            if debugCount > 100000
                %display('Here');%Shouldn't get here, debug only
                AlternativeRoute = {};
                BadMeshRoute = {};
                debugCount = 0;
            end
        end

        if ~restart
            RowToMatch = PossibleRow(1);
            rows = SplitElements(PossibleRow(1),:);
            lastRow = OrderedVec(end,:);
            %ColumnToMatch = find(~ismember(rows,lastRow));%Commented out
            %as I think line below does the same thing faster
            %ColumnToMatch = find(setdiff_ns(rows,lastRow));
            
            %DEBUG
%             ColumnToMatch = find(rows == setdiff_ns(rows,lastRow));
            rowToCheck = setdiff_ns(rows,lastRow);
            if size(rowToCheck,2) == 1 && size(rows,2) == 3
                ColumnToMatch = find(rows == rowToCheck);
            else
                display('here');
            end
            
%             if ColumnToMatch ~= tempTest
%                 display('');
%             end
%             tempTest2 = fastintersect(rows,lastRow);

        end
        restart = false;
    end
    
    %Find a route until you can't go anymore. Everytime you get to a split
    %store the current location and alternative direction. When your route
    %can't go anymore, check 1) has it used all elements 2) From the end
    %point can get back to the beginning. If these both pass, the route is
    %OK, otherwise go back and try each of your alternate directions. There
    %should be 2 possible routes (going left and right out the first
    %element...)
 
 end