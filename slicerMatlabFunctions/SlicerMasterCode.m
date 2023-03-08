% Code to run TrachSlice for OSA in a loop for all surfaces
function [] = SlicerMasterCode(subject, path)


disp(['this is for subject', subject]);
disp(['slicing partition: ',path]);

%path = 'LeftNoseDecending';
%subject = 'DYMOSA801';

vtkfiles = ['../',path,'/FFD/vtk/*.vtk'];
stlfiles = ['../',path,'/FFD/stl/*.stl'];

tic
VTKdir = dir(vtkfiles); %Name of the directory with centerlines
STLdir = dir(stlfiles); %Name of the directory with STLs

folder = ['./',path,'SlicedSTLs'];

if ~exist(folder,'dir')
    mkdir(folder);
%      cd('SlicedSTLs');
elseif exist(folder,'dir')
    rmdir (folder, 's')
    mkdir(folder);
%      cd('SlicedSTLs');
end

addpath(VTKdir(1).folder);
addpath(STLdir(1).folder);


arcLength_test = {}';
area_test = {}';
missedPlane_test = {}';
goodpos_test = {}';

condtru_test = {}';

tn1_co = {}';
n2_co = {}';

postest = {}';
trachnose1_Area = {}';
nose2_Area = {}';
trachnose1_Arclength = {}';
nose2_Arclength = {}';

for i = 1:numel(VTKdir)
    STLfilename = STLdir(i).name;
    VTKfilename = VTKdir(i).name;
    
    [arcLength, area, pos, centroid, goodpos, missedPlane] = trachSlice_OSA_AB_3(STLfilename, VTKfilename);
    
    goodpos3 = goodpos(:,3);
    result = [];
    
    if i == 1
        for p = 1:numel(goodpos3)-1
            result = [result, abs(goodpos3(p) - goodpos3(p+1)) > 5];
        end
        condtru = find(result==1);
    end
    
    TN1_co = pos(1:condtru,:);
    N2_co = pos(condtru+1:end,:);
    
    trachnose1_area = area(1:condtru);
    nose2_area = area(condtru+1:end);
    
    trachnose1_arclength = arcLength(1:condtru);
    nose2_arclength = arcLength(condtru+1:end);
    
    condtru_test{i,1} = condtru;
    
    arcLength_test{i,1} = arcLength; %#ok<*NASGU>
    area_test{i,1} = area;
    missedPlane_test{i,1} = missedPlane;
    goodpos_test{i,1} = goodpos;
    
    % Flipping the order for plotting correctly
    
    T = diff(trachnose1_arclength);
    T1 = flip(T);
    a = 0;
    T1_final = horzcat(a,T1);
    T1_final = cumsum(T1_final);
    trachnose1_arclength_final = T1_final;
    
    D = diff(nose2_arclength);
    D1 = flip(D);
    b = 0;
    D1_final = horzcat(b,D1);
    D1_final = cumsum(D1_final);
    nose2_arclength_final = D1_final;
    
    % saving variables from each run in a cell for later use in
    % plotting/debugging
    
    postest{i,1} = pos;
    tn1_co{i,1} = TN1_co;
    n2_co{i,1} = N2_co;
    trachnose1_Area{i,1} = flip(trachnose1_area);
    nose2_Area{i,1} = flip(nose2_area);
    trachnose1_Arclength{i,1} = trachnose1_arclength_final;
    nose2_Arclength{i,1} = nose2_arclength_final;
    fclose('all');
    movefile('out_*Planes*.stl', folder)
    movefile('out_*.csv', folder)
    
end

save([subject,'_',path,'.mat']) %Change the name of the file accordingly. 
toc

end



