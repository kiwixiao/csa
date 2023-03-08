function [colours] = ChangingColours(noOfStates)

R = zeros(noOfStates,1);
G = zeros(noOfStates,1);
B = zeros(noOfStates,1);

Q1 = 1 : round(noOfStates / 4);
Q2 = Q1(end) + 1 : round(noOfStates / 2);
Q3 = Q2(end) + 1 : Q1(end) + Q2(end);
Q4 = Q3(end) : noOfStates - 1;

g1 = 0: 1 / length(Q1) : 1;
b = 1: -1 / length(Q2) : 0;
r = 0: 1 / length(Q3) : 1;
g2 = 1: -1 / length(Q4) : 0;

R(Q1) = 0;
G(Q1) = g1(1:end-1);
B(Q1) = 1;

R(Q2) = 0;
G(Q2) = 1; 
B(Q2) = b(1:end-1);

R(Q3) = r(1:end-1);
G(Q3) = 1;
B(Q3) = 0;

R(Q4) = 1;
G(Q4) = g2(1:end-1);
B(Q4) = 0;

R(end) = 1;
G(end) = 0;
B(end) = 0;

colours = [R, G, B];

end
