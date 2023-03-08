function [] = ChangeDisplacementTableToVelocities(CSVin)
% Read in a table of displacement's from andreas and convert to a velocity
% table

    timeBetweenImages = 0.3;
    CFDdt = 0.05;
    
    velocityCSVout = [CSVin(1:end-4) '_Velocities.csv'];
    displacementCSVout = [CSVin(1:end-4) '_displacementPerTimestep.csv'];

    % Read CSV titles
    fid = fopen(CSVin, 'r');
    headers = textscan(fid,'%s',1);
    headersCell = textscan(char(headers{1,1}),'%s','Delimiter',',');
    displacementHeaders = headersCell{1,1};
    fclose(fid);
    
    velocityHeaders = strrep(displacementHeaders, 'D', 'V');
    
    rawData = csvread(CSVin, 1, 0);
    
    positions = rawData(:,1:3) ./ 1000; %scale to m
    displacements = rawData(:,4:end);

    clear rawData;
    
    velocities = displacements ./ (timeBetweenImages * 1000); % Convert mm displacement to m/s velocity
    
    displacementPerDT = displacements * CFDdt ./ (timeBetweenImages * 1000); % Convert mm displacement to m displacement per CFD timestep
    
    
    csvwrite_with_headers(velocityCSVout,[positions, velocities],velocityHeaders);
    csvwrite_with_headers(displacementCSVout,[positions, displacementPerDT],displacementHeaders);
    
end