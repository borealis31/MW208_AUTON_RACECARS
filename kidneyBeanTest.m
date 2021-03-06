function [scenario, testName, roadWidth] = kidneyBeanTest()
% createDrivingScenario Returns the drivingScenario defined in the Designer

% Generated by MATLAB(R) 9.10 (R2021a) and Automated Driving Toolbox 3.3 (R2021a).
% Generated on: 03-Jun-2021 21:27:11

testName = 'Kidney Bean Test';

% Construct a drivingScenario object.
scenario = drivingScenario('SampleTime',0.05,'StopTime',45);

% Add all road segments
roadCenters = [ 5.6400   13.6000 0;
               -9.3600   14.1800 0;
               -8.7800  -19.6200 0;
               24.0000  -19.7800 0;
               24.1600   13.9000 0;
               14.4000   13.6000 0;
               13.9600  -12.3400 0;
                4.6200  -12.0400 0;
                5.6400   13.6000 0];

roadWidth = 3;
road(scenario, roadCenters, roadWidth, 'Name', 'Road');

