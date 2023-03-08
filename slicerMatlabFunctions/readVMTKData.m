function [ pos, MISR, curvature ] = readVMTKData( filename )
%Read in the data from Andrew's VMTK script and return 
%   Detailed explanation goes here


%X	Y	Z	MaximumInscribedSphereRadius	Curvature	Torsion	FrenetTangent0	FrenetTangent1	FrenetTangent2	FrenetNormal0	FrenetNormal1	FrenetNormal2	FrenetBinormal0	FrenetBinormal1	FrenetBinormal2	Abscissas	ParallelTransportNormals0	ParallelTransportNormals1	ParallelTransportNormals2

data = dlmread(filename, ' ', 1, 0);

%The data is reversed compared to the trachslice data

data = data(end:-1:1,:);

pos = data(:,1:3);
MISR = data(:,4);
curvature = data(:,5);

end

