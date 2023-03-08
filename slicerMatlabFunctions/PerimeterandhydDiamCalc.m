function [perimeter, hydDiam, diameter] = PerimeterandhydDiamCalc(n)

% Perimeter calc
% Loads stl planes from slicer tool and calcs the perimeter of each

%     if nargin == 2
%         doHydDiameter = false;
%     end
    
        doHydDiameter = true;
%     planeFilename = [filenameStart '-Planes-'];

 disp('you picked to calculate the left side diameter first')
 disp('now we need to double check both left side and right side have the same amount time points data')
 left_numFiles = length(dir(['./LeftNoseDecendingSlicedSTLs'  '/out_*-All.stl']));
 right_numFiles = length(dir(['./RightNoseSlicedSTLs' '/out_*-All.stl']));
        if left_numFiles == right_numFiles
            disp('left and right side have the same time points data')
            fileNumberCheck = true;
        else
            disp('left and right time points are not the same')
            fileNumberCheck = false;
            return
        end
              
FolderName = {'LeftNoseDecendingSlicedSTLs','RightNoseSlicedSTLs'};

    %n = str2num(n)
    fileFolderName = FolderName{n};
    disp(['calculating the diameter of ', fileFolderName])
    
    tpFile = dir(['./' fileFolderName '/out_*-All.stl']);
    noOfTimePoints = length(tpFile);
    t = linspace(0,(noOfTimePoints-1)*100,noOfTimePoints);
        
 for k = 1:noOfTimePoints
     
     file = dir(['./' fileFolderName '/'])
     filenames = {file.name}
     timePoints_idx = find(contains(filenames, [sprintf('%06d', t(k)) '-Planes']));
     timePoints_planeName = filenames(timePoints_idx)
%      planeFilename = [filenameStart sprintf('%06d', t(k)) '-Planes-']
     noOfPlanes = (length(timePoints_planeName)-1)
     perimeter = zeros(noOfPlanes,1);
     hydDiam =  zeros(noOfPlanes,1);
   for i = 1:noOfPlanes
        % Read in stl files.
        display(['Reading plane no ' num2str(i)]);
%         filename = [planeFilename num2str(i,'%03d') '.stl'];
        filename = ['./' fileFolderName '/' timePoints_planeName{i}]
        %filename = [planeFilename num2str(i,'%d') '.stl'];
        [Vsurface,Fsurface,~] = STL_Import(filename);
        
        % Find edges
        % Edges are unique pairs of nodes
        
        edge=[];
        
        % From node list, list all edges
        for j = 1:size(Fsurface,1)
                
            % Clear any faces that have only two edges - i.e. one edge = 0
            if sum(Fsurface == 0) > 0
                display('Removing non-manifold edge from surface');
                Fsurface(j,:) = [];
                j = j-1;
                continue
            end
            
            % Handle non manifold edges/vertices
            if Fsurface(j,1) == Fsurface(j,2) || Fsurface(j,1) == Fsurface(j,3) || Fsurface(j,2) == Fsurface(j,3)
                %This shouldn't happen, but does when the surface of the
                %trachea touches itself - like a hole stuck to the wall at
                %a point - should be handled in splitedges in sublobemulti,
                %but isn't as yet.
                display('Non-manifold edges found');


                %If there is a legitimate edge on the triangle use that.
                %i.e. if the nodes are [a, b, b], assume a-b is a good edge
                nonManEdges = unique(Fsurface(j,:));
                if size(nonManEdges,2) == 2
                   edge(end+1,:) = nonManEdges; 
                end
                
                continue
            end
            node(1) = min(Fsurface(j,:));
            node(3) = max(Fsurface(j,:));
            minNode = find(Fsurface(j,:) == node(1));
            maxNode = find(Fsurface(j,:) == node(3));
            midNode = setdiff( [ 1 2 3],[minNode maxNode]);
            node(2) = Fsurface(j,midNode);
            
            edge(end+1,:) = [node(1) node(2)];
            edge(end+1,:) = [node(2) node(3)];
            edge(end+1,:) = [node(1) node(3)];
           
            
        end
        
        %perimeterEdges = unique(edge,'rows');
        [a,c,b] = unique(edge,'rows');
        % Count the number of times each edge appears, take the ones which
        % only appear once.
        n = accumarray(b,1);
        [row, col] = find(n == 1);
        perimeterEdges = a(row,:);
        
%         %% Test Only
%         figure
%         pTest = patch('Vertices', Vsurface, 'Faces', Fsurface);
%         set(pTest,'facecolor','none','edgecolor','black');
%         hold on
%         for z = 1:size(perimeterEdges,1)
%             plot3([Vsurface(perimeterEdges(z,1),1), Vsurface(perimeterEdges(z,2),1)],...
%                   [Vsurface(perimeterEdges(z,1),2), Vsurface(perimeterEdges(z,2),2)],...
%                   [Vsurface(perimeterEdges(z,1),3), Vsurface(perimeterEdges(z,2),3)],'r');
%         end
        

        %% Calcualate perimeter
        
        for j = 1:size(perimeterEdges,1)
           perimeter(i) = perimeter(i) +  sqrt((Vsurface(perimeterEdges(j,1),1) - Vsurface(perimeterEdges(j,2),1))^2 + (Vsurface(perimeterEdges(j,1),2) - Vsurface(perimeterEdges(j,2),2))^2 + (Vsurface(perimeterEdges(j,1),3) - Vsurface(perimeterEdges(j,2),3))^2);
        end
        
        %% Run the code to also calculate the Major Diameter, the minor diameter and the angle between the major diameter and Y axis (AP)
        diameter(i,:)  = DiameterCalc( Fsurface, Vsurface, perimeterEdges );
   end
    save(['out_' sprintf('%06d', t(k)) '-Perimeters-' fileFolderName(1) '.mat'], 'perimeter');
    save(['out_' sprintf('%06d', t(k)) '-Diameters-' fileFolderName(1) '.mat'], 'diameter');
    
  
   
%    save([filenameStart sprintf('%06d', t(j)) '-Perimeters.mat'], 'perimeter');
%    save([filenameStart sprintf('%06d', t(j)) '-Diameters.mat'], 'diameter');
  
    if doHydDiameter
        %read in area data, assuming it starts filenamestart + Data.csv
        areaFilename = ['./' fileFolderName '/out_' sprintf('%06d', t(k)) '-Data.csv'];
        area = csvread(areaFilename);
        area = area(end:-1:1,:);
               
        if size(area,1) == size(perimeter,1)
            hydDiam = 4 .* area(:,2) ./ (perimeter);
%             fDh = figure;
%             plot(area(:,1),hydDiam)
%             xlabel('Position (mm)');
%             ylabel('Hydraulic Diameter (mm)');
%             savefig(fDh, [filenameStart 'HydraulicDiameters.fig']);
%             saveas(fDh, [filenameStart 'HydraulicDiameters.png'],'png');
            
        else
            display('Cant calculate hydraulic diameters, have a different number of areas and perimeters');
        end
        
    end
    save(['out_' sprintf('%06d', t(k)) '-Hydraulicdiameters-' fileFolderName(1) '.mat'], 'hydDiam');
    
% %     figure
% %     plot(diameter(:,2)./diameter(:,1))
    
 end
 
 

 
end
    
    
     