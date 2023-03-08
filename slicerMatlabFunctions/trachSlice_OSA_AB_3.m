function [arcLength, area, pos, centroid, goodpos, missedPlane] = trachSlice_OSA_AB_3(STLfilename,VTKfilename,znormal)

if nargin == 2
    znormal = false;
end

%Read in stl
disp('Reading in STL');
[Vsurface,Fsurface,~] = STL_Import(STLfilename);
disp('Done reading in STL')

%Read in VTK file
[pos, ~] = read_vtk(VTKfilename);

% Gives unique vertices without sorting it in a particlar order
pos = unique(pos','rows','stable');


%% Commented out based on Alister's suggestion to just do this once. 1/8/21
% zaxis = pos1(3,:);
% result = [];
% for i = 1:numel(zaxis)-1
%     %result = [result, zaxis(i) < zaxis(i+1)];
%     result = [result, abs(zaxis(i) - zaxis(i+1)) > 5];
% end
% condtru = find(result==1);
%
%     % Test for multiple condtru
%     Find_num = numel(condtru_pre);
%     condtru = condtru_pre(Find_num);

% if numel(pos) == numel(pos1)
%     pos = pos';
%     pos = pos(end:-1:1,:);
% else
%     pos = unique(pos','rows','stable')';
%     trach_bronchi_1 = pos(:,1:condtru);
%     justbronchi_2 = pos(:,condtru+1:end);
%     
%     % Test: compare all elements of justbronch_2 with trach_bronchi not
%     % just one to one
%     
%     %test
%     %        trach_bronchi_1 = flip(pos(:,1:condtru),2);
%     %        justbronchi_2 = flip(pos(:,condtru+1:end),2);
%     
%     sizebronch2 = numel(justbronchi_2(1,:));
%     
%     % discard = abs(zaxis(1:sizebronch2) - justbronchi2(1:end));
%     %discard = (trach_bronchi_1(:,1:sizebronch2) - justbronchi_2(:,1:end));
%     discard = (trach_bronchi_1(:,1:sizebronch2)) - (justbronchi_2(:,1:end));
%     %testwa = abs(trach_bronchi_1(:,1:sizebronch2)) - abs(justbronchi_2(:,1:end));
%     % discard1 = (discard(1,:)); WORKED WELL TILL 08/06/2019
%     discard1 = (discard(1,:)); % Trying something new
%     
%     if all(discard1<0)
%         discard = abs(trach_bronchi_1(:,1:sizebronch2)) - abs(justbronchi_2(:,1:end));
%         if all(discard(1,:)<0)
%             discard = abs(abs(trach_bronchi_1(:,1:sizebronch2)) - abs(justbronchi_2(:,1:end)));
%         end
%     else
%         discard = (trach_bronchi_1(:,1:sizebronch2)) - (justbronchi_2(:,1:end));
%     end
%     discard1 = (discard(1,:));
%     threshold = 1.5; %Threshold for avoiding overlap of normals
%     discard1(discard1(1,:)<threshold) = [];
%     h = numel(discard1);
%     final_justbronchi_2 = justbronchi_2(:,1:h);
%     %     pos = horzcat(trach_bronchi_1,final_justbronchi_2); #Firstoneworks
%     pos = horzcat(final_justbronchi_2, trach_bronchi_1);
%     
%     pos = unique(pos','rows','stable')';
%     
%     pos = pos';
%     
%     pos = pos(end:-1:1,:); %correct one
% end

%%

if znormal == false
    % Add different normals for Trachea-bronch1 and bronchi2
    normals = pos(2,:) - pos(1,:);
    
    for i = 2:size(pos,1) - 1
        
        normals(i,:) = pos(i+1,:) - pos(i-1,:);
        
    end
    
    normals(end+1,:) = pos(end,:) - pos(end-1,:);
else
    for i = 1:size(pos,1)
        normals(i,:) = [ 0 0 1]; % This has been changed for the balloon experiment - should be [ 0 0 1]
    end
end

nodes = [];
tri = [];
j = 1;
missedPlane = [];
deep = [];
normals = smoothdata(normals);
%deepu = [];


for i = 1:size(pos,1)
    
    display(['Testing Plane number ', num2str(i), ', of ', num2str(size(pos,1))]);
    [currentArea, currentCentroid, nodeCell, triCell] = sublobeMulti_DG_AB4(normals(i,:), pos(i,:), Vsurface, Fsurface, i);
    
    %to see which cell the centroid is in
    if currentArea{1} == -1
        
        % We haven't successfully cut the geom here. we need to
        % handle this without crashing the code! Set answers to area,
        % centroid etc. that won't crash the code...
        
        %We are trying to NOT use missedPlane here, rather set values
        %that are needed to NaN.
        
        area(j) = NaN;
        centroid(j,:) = [pos(i,1), pos(i,2), pos(i,3)];% We don't have a real centroid, so just assume it was cut at the centerline pos
        goodpos(j,:) = pos(i,:);
        nodesSeparate{j} = NaN;
        triSeparate{j} = NaN;
        j = j+1;
        
        
        
        %             missedPlane = [missedPlane i];
        continue;
    end
    deep = [deep currentCentroid];
    goodpos(j,:) = pos(i,:);
    % If multiple areas detected, use area with centroid closest to centreline point
    % rather than sum.
    if isempty(currentCentroid{3})
        if isempty(currentCentroid{1})
            if isempty(currentCentroid{2})
                warning('No area detected');
                continue
            else
                sideToUse = 2;
            end
        else
            if isempty(currentCentroid{2})
                sideToUse = 1; % decide which side of the nose to use?
            else
                if abs(norm(goodpos(j,:) - currentCentroid{1}) - norm(goodpos(j,:) - currentCentroid{2})) < 2
                    if currentArea{1,1} > currentArea{2,1}    
                        sideToUse = 1;
                    else
                        sideToUse = 2;
                    end
                elseif norm(goodpos(j,:) - currentCentroid{1}) < norm(goodpos(j,:) - currentCentroid{2})
                    sideToUse = 1;
                else 
                    sideToUse = 2;
                end
            end
        end
        
        currentCentroid{3} = currentCentroid{sideToUse};
        currentArea{3} = currentArea{sideToUse};
        nodeCell{3} = nodeCell{sideToUse};
        triCell{3} = triCell{sideToUse};
        
        %threshold_centroid = 2; %mm
        %away = norm(goodpos(j,:) - currentCentroid{sideToUse});
        %deepu = [deepu away];
        
    end
    %change here
    
    %          threshold_centroid = 1; %mm
    %          if norm(goodpos(j,:) - currentCentroid{3}) > threshold_centroid
    %              goodpos(j,:) = [];
    %          else
    %              continue;
    %          end
    
    centroid(j,:) = currentCentroid{3};
    area(j) = currentArea{3};
    
    % test
    %          percentdiff(1) = abs(1-goodpos(j,1)/currentCentroid{3}(1))*100;
    %          percentdiff(2) = abs(1-goodpos(j,2)/currentCentroid{3}(2))*100;
    %          percentdiff(3) = abs(1-goodpos(j,3)/currentCentroid{3}(3))*100;
    % test end
    threshold_centroid1 = 100;
    %if abs(norm(goodpos(j,:) - currentCentroid{3})) > threshold_centroid1
    X = abs(abs(goodpos(j,:)) - abs(currentCentroid{3}));
    for g = 1:3
        %if percentdiff(g) > 80 %X(g) > threshold_centroid1
        if X(g) > threshold_centroid1
            goodpos(j,:) = NaN;
            area(j) = NaN;
        else
            continue;
        end
    end
    
    nodesSeparate{j} = nodeCell{3};
    triSeparate{j} = triCell{3};
    triCell{3} = triCell{3} + size(nodes,1);
    nodes = [nodes; nodeCell{3}];
    tri = [tri; triCell{3}];
    
    j = j + 1;
end

%    away = goodpos-centroid;
%     t = size(away,1);
%     for x = 1:t
%         threshold = 1;
%         if abs(away(x,:)) > threshold
%             away(x,:) = [];
%         else
%             continue;
%         end
%     end
%     for x = 1:size(pos,1)
%          threshold_centroid = 1; %mm
%          if abs(norm(goodpos(x,:) - centroid(x,:))) > threshold_centroid
%              goodpos(x,:) = [];
%          else
%              continue;
%          end
%     end
%Difference_pos = abs(abs(pos) - abs(centroid));
%Difference_goodpos = abs(abs(goodpos) - abs(centroid));
arcLength(size(goodpos,1)) = 0.0;
for i = 2:size(goodpos,1)
    arcLength(i) = arcLength(i-1) + norm(goodpos(i,:) - goodpos(i-1,:));
    %arcLength(i) = arcLength(i-1) + norm(pos(i,:) - pos(i-1,:));
end

if max(arcLength) < 1
    % Data is in m, convert to mm
    arcLength = arcLength .* 1000;
    area = area .* 1e6; %   m^2 to mm^2
end

dataToWrite(:,1) = arcLength(end:-1:1)';
%dataToWrite(:,2) = area(end:-1:1);
dataToWrite(:,2) = area(end:-1:1);
csvwrite([STLfilename(1:end-4) '-Data' '.csv'], dataToWrite);


stlwrite(['./' STLfilename(1:end-4) '-Planes-All' '.stl'], tri, nodes);
for i = 1:size(goodpos,1)
    if size(triSeparate{i},2) == 3 && size(nodesSeparate{i},2) == 3 %~isnan(triSeparate{i}) && ~isnan(nodesSeparate{i})
        stlwrite([STLfilename(1:end-4) '-Planes-' num2str(i, '%03d') '.stl'], triSeparate{i}, nodesSeparate{i});
    end
end
%% AB debug
%     figure
%     subplot(2,2,1)
%     plot(area,'x-')
%     xlabel('Plane number')
%     ylabel('Area')
%
%     subplot(2,2,2)
%     plot(arcLength, area,'x-')
%     xlabel('Arc length')
%     ylabel('Area')
%
%     subplot(2,2,3)
%     plot3(centroid(:,1),centroid(:,2),centroid(:,3),'x-')
%     xlabel('X')
%     ylabel('Y')
%     zlabel('Z')
%     daspect([1,1,1])
%
%     subplot(2,2,4)
%     plot(arcLength,'x-')
%     ylabel('Arc length')
%ylabel('Area')

end
