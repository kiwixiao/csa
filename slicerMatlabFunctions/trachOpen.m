function [ output_args ] = trachOpen( STLfilename, axis )
% Cut the ends of the STL open, leaving planar openings.
%

    if nargin == 1
        axis = 3;
    end

    display('Cutting STL');
    
    %Read in stl
    display('Reading in STL');
    [Vsurface,Fsurface,normals] = STL_Import(STLfilename);
    display('Done reading in STL');
    
    stlwrite([STLfilename(1:end-4) '-UnCut' '.stl'], Fsurface, Vsurface);

    % Find the extreme positions in the axis of interest
    maxPos = max(Vsurface(:,axis));
    minPos = min(Vsurface(:,axis));

    height = maxPos - minPos;
    
    cutHeight = 0.05 * height; % 1 per cent
    maxCut = maxPos - cutHeight;
    minCut = minPos + cutHeight;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Based on Lulu's segmentations, hard code these values instead
    maxCut = -626.5;
    minCut = -782.5;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    VtoDelete = find(Vsurface(:,axis) > maxCut);
    VtoDelete = [VtoDelete; find(Vsurface(:,axis) < minCut)];
    
    
    
    %% Map vertices into an array that doesn't include the deleted ones
%     Vmap = zeros(length(Vsurface),1);
%     for i = 1:length(VtoDelete)
%         
%         Vmap(VtoDelete(i):end) = Vmap(VtoDelete(i):end) - 1;
%         
%     end
%     
%     Vsurface(VtoDelete,:) = [];
%     Vmap(VtoDelete) = inf;

    %%
    FtoDelete = ismember(Fsurface,VtoDelete) ;
    
    normalsToDelete = abs(normals(:,3)) > 0.999;
    
    FRowToDelete = sum(FtoDelete,2) > 0 & normalsToDelete;
    
    Fsurface = Fsurface(~FRowToDelete,:);
    
    %%
    nodeList = 1:length(Vsurface); 
    Vunused = ismember(nodeList,Fsurface);
    %Vunused = sum(Vunused,2) == 0;
    
    Vunused = nodeList(Vunused)';
    
    for i = 1:length(Vsurface)
        for j = 1:3
            replacement = find(Vunused == i);
%             if length(replacement) ~= 1
%                 display('')
%             end
            Fsurface(Fsurface == i) = replacement;   
        end
    end
    
    VtoReallyDelete = setdiff(nodeList,Vunused);
    
    Vsurface(VtoReallyDelete,:) = [];
            
%     Vmap = zeros(length(Vsurface),1);
%          
%     for i = 1:length(Vunused)
%         
%         Vmap(Vunused(i):end) = Vmap(Vunused(i):end) - 1;
%         
%     end
%     
%     Vsurface(Vunused,:) = [];
%     Vmap(Vunused) = inf;
%     
%     for i = 1:length(Vmap)
%         for j = 1:3
%             Fsurface(Fsurface(:,j) == i,j) = Fsurface(Fsurface(:,j) == i,j) + Vmap(i);
%         end
%     end
    %%  
    stlwrite([STLfilename(1:end-4) '.stl'], Fsurface, Vsurface);

end

