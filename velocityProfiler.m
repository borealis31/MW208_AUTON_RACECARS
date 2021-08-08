function [vProfile, vProfileByS, bestLapSet] = velocityProfiler(timeStep, waypointsIn, lapsQty, ftMax, fnMax, vMax, mass)
%VELOCITYPROFILER Velocity Profiler uses corrected point-mass kinematics from the below thesis to generate a time-velocity profile that is optimized for performance
% -------------------------------------------------------------------------
% 
% Based on:
% Friedl, Tyler, "Optimized Trajectory Generation for Car-Like Robots on a Closed Loop Track" (2017).
% Electronic Theses and Dissertations. 1370.
% https://digitalcommons.du.edu/etd/1370
%
% -------------------------------------------------------------------------
% 
% INPUTS: timeStep — Resolution of time by which to determine velocity
%                     profile
%         waypointsIn — Optimized curvature path
%         lapsQty — Number of laps to calculate velocty profile for
%         ftMax — Maximum tangential force constraint
%         fnMax — Maximum normal (centripetal) force constraint
%         mass — Mass of the point mass
%
% OUTPUT: vProfile — 8-Column matrix detailing the velocity profile of the
%                     point-mass [velocity, time, cumulative distance, X,
%                     Y, tangential force, normal force, lap number]
%         vProfileByS — 3-Column matrix detailing the velocity of the mass
%                        in terms of distance traveled as well as the
%                        critical velocities at each point for reference
%                        [distance traveled, velocity, critical velocity]
%         bestLapSet — 4-Column matrix detailing the velocity profile of
%                       the point-mass during the middle-most lap
%                       [velocity, time, X, Y]
% 
% -------------------------------------------------------------------------
%

% Initialize constants used in kinematic equations to reduce calculations
ftFactor = ftMax^2;
ftfnFactor = (ftMax/fnMax)^2;

% If we want multiple laps, concatenate the lap points excluding
% overlapping points
for idxLap = 1:lapsQty
    if idxLap == 1
        waypointsTracked = waypointsIn;
    else
        waypointsTracked = [waypointsTracked; waypointsIn(2:end,:)];
    end
end

% Initialize number of points and number of profiles
numPoints = length(waypointsTracked);
numProfiles = 1;

% Calculate distance traveled between each point as well as distance
% traveled cumulatively at each point
sDiff = ([0; sqrt(sum((waypointsTracked(2:end,:)-waypointsTracked(1:end-1,:)).^2,2))]);
sAtPt = cumsum(sDiff);

% Normalize waypoints and distances to ensure that points along path are
% evenly spaced
normFactor = floor(sAtPt(end))/sAtPt(end);
sNormalized = sAtPt.*normFactor;
waypointsTrackedNormalized = waypointsTracked.*normFactor;
waypointsTracked = [spline(sNormalized, waypointsTrackedNormalized(:,1),0:0.01:sNormalized(end))',...
                    spline(sNormalized, waypointsTrackedNormalized(:,2),0:0.01:sNormalized(end))']./normFactor;
sDiff = ([0; sqrt(sum((waypointsTracked(2:end,:)-waypointsTracked(1:end-1,:)).^2,2))]);
sAtPt = cumsum(sDiff);

% Calculate a ds value for use to calculate finer X, Y, DDX, and DDY
% values across the path later on
ds = sAtPt(end)/(10*numPoints);

X = spline(sAtPt,waypointsTracked(:,1));
Y = spline(sAtPt,waypointsTracked(:,2));

DDX = fnder(X,2);
DDY = fnder(Y,2);

% Create a high-definition (HD) set of cumulative distance-travelled values
sHD = (0:ds:sAtPt(end))';

% Pass cumulative distances as well as DDX and DDY splines to the local
% curvature function to calculate radii of curvature across the track. Use
% these values to create a spline for radius of curvature values. See the
% local function for more information
rTrack = 1./kAtPoint(sHD, DDX, DDY);
rInterp = spline(0:ds:sAtPt(end), rTrack);

% Calculate maximum velocities along the path by radius of curvature and
% force constraints
vCrit = sqrt(ppval(rInterp,sHD).*(fnMax/mass));

% Find local minima on the vCrit array and return a 2-column matrix
% [cumulative distance traveled, critical velocity at cumulative distance
% traveled]. Output the number of profiles to be generated
vCritMins = [sHD(islocalmin(vCrit,'FlatSelection','center','MaxNumExtrema',floor(sAtPt(end)/ds))), vCrit(islocalmin(vCrit,'FlatSelection','center','MaxNumExtrema',floor(sAtPt(end)/ds)))];
fprintf([num2str(size(vCritMins,1),'%.0f'), ' Profiles Identified', newline]);

% Initialize number of profiles needed to create cumulative profile
profiles = cell(length(vCritMins)+2, 1);

% For each minima, generate a velocity profile based on kinematics
% equations. Additionally, generate a velocity profile for initial and final
% constraints (vo = 0, vf = 0). See local function for more information
for idxMin = 1:length(vCritMins)+2
    if idxMin == length(vCritMins)+1
        profiles{idxMin} = generateVByS(ds, [0,0], numProfiles, sAtPt, rInterp, ftFactor, ftfnFactor, vMax, mass);
    elseif idxMin == length(vCritMins)+2
        profiles{idxMin} = generateVByS(ds, [sAtPt(end), 0], numProfiles, sAtPt, rInterp, ftFactor, ftfnFactor, vMax, mass);
    else
        profiles{idxMin} = generateVByS(ds, vCritMins(idxMin,:), numProfiles, sAtPt, rInterp, ftFactor, ftfnFactor, vMax, mass);
    end
    
    if numProfiles == 1
        vProfComposite = profiles{numProfiles};
    else
        vProfComposite = min([vProfComposite, profiles{numProfiles}], [], 2);
    end
    
    numProfiles = numProfiles + 1;
end

% Initialize composite array of velcoity profiles and select the minimum
% velocity values for each distance traveled to create optimized velocity
% profile
vProfileByS = [(0:ds:sAtPt(end))', vProfComposite, vCrit];

% Initalize time-velocity profile values. NOTE: the initial starting
% position is ds rather than 0, as starting at 0 results in a long pause at
% the beginning of the time-velocity profile
v = 0;
vIdx = 1;
t = 0;
sTraveled = ds;
lap = 1;
ft = 0;
fn = 0;

% Create a velocity spline in terms of distance traveled
vInterp = spline(vProfileByS(:,1), vProfileByS(:,2));

% Start timer for time-velocity calculations
vTime = tic;
fprintf(['Beginning Time-Velocty Calculations...' newline])

% While the distance traveled is less than the total distance, run state
% updaters and velocity calculations. Store data at each time point in a
% variety of arrays
while (sAtPt(end) - sTraveled(end)) > 0.001
    sTraveled(vIdx+1,1) = sTraveled(vIdx,1) + v(vIdx)*timeStep;
    vIdx = vIdx + 1;
    v(vIdx,1) = ppval(vInterp, sTraveled(end));
    t(vIdx,1) = t(vIdx-1,1) + timeStep;
    ft(vIdx,1) = abs(mass*(v(vIdx,1)-v(vIdx-1,1))/timeStep);
    fn(vIdx,1) = abs(mass*(v(vIdx,1)^2)/ppval(rInterp,sTraveled(vIdx,1)));
    lap(vIdx,1) = ceil((sTraveled(end)/sAtPt(end))*lapsQty);
end
fprintf(['Time-Velocity Calculations Complete: %f sec' newline], toc(vTime))

% Compile arrays into velocity profile output
vProfile = [v, t, sTraveled, ppval(X, sTraveled), ppval(Y, sTraveled), ft, fn, lap];

% Find the logic map of all values that exist during the middle-most lap
% and record them in the bestLapSet output variable
bestLapMap = lap == ceil(lapsQty/2);
bestLapSet = [vProfile(bestLapMap,1), vProfile(bestLapMap,2), vProfile(bestLapMap,4), vProfile(bestLapMap,5)];

end

function prof = generateVByS(ds, critMins, numProfiles, sIn, rInterp, ftFactor, ftfnFactor, vMax, mass)
%LOCALGENERATEVBYS This local function generates a profile for velocity in terms of distance traveled
% -------------------------------------------------------------------------
%
% INPUTS: ds — Distance steps by which to calculate the profile
%         critMins — The point of the critical minimum [S, vCrit]
%         numProfiles — Current profile number
%         sIn — Distance traveled array
%         rInterp — Radius of curvature spline in terms of distance
%                   traveled
%         ftFactor — Constant for calculations based on known
%                     characteristics
%         ftfnFactor — Constant for calculations based on known
%                       characteristics
%         mass — Mass of the point-mass
%
% OUTPUTS: prof — Velocity profile in terms of distance traveled [S, V]
%
% -------------------------------------------------------------------------

% Initialize arrays for acceleration and deceleration profiles
vDecel = critMins(2);
vAccel = critMins(2);
decelS = (critMins(1):-ds:0)';
accelS = (critMins(1):ds:sIn(end))';

% Eliminate accelaration and deceleration S value arrays if any are above
% maximum or below minimum
if any(accelS > sIn(end))
    accelS = [];
end

if any(decelS < 0)
    decelS = [];
end

% Set previous velocity value
vPrev = critMins(2);

% For each ds point between the critical point's S and 0, calculate the
% velocity value subject to kinematic constraints
rDecel = ppval(rInterp,decelS);
for idxDecel = 2:length(decelS)
    vNew = real(sqrt((vPrev^2)+((2*abs(decelS(idxDecel)-decelS(idxDecel-1))/mass)*sqrt(ftFactor - ftfnFactor*((mass*vPrev^2)/rDecel(idxDecel-1))^2))));
    vPrev = vNew;
    vDecel(end+1,1) = vPrev;
end

% Reset previous velocty value
vPrev = critMins(2);

% For each ds point between the critical point's S and max S, calculate the
% velocity value dubject to kinematic constraints
rAccel = ppval(rInterp,accelS);
for idxAccel = 2:length(accelS)
    vNew = real(sqrt((vPrev^2)+((2*abs(accelS(idxAccel)-accelS(idxAccel-1))/mass)*sqrt(ftFactor - ftfnFactor*((mass*vPrev^2)/rAccel(idxAccel))^2))));
    vPrev = vNew;
    vAccel(end+1,1) = vPrev;
end


% Compile results into a new matrix for outputting velocity profile in
% terms of S
prof = [[flip(vDecel,1); vAccel(2:end)], [flip(decelS,1); accelS(2:end,1)]];
if mod(numProfiles, 25) == 0
    fprintf([num2str(numProfiles, '%.0f'), ' Profiles Generated', newline]);
end

profVMaxTemp = prof(:,1);
profVMaxTemp(profVMaxTemp > vMax) = vMax;
prof(:,1) = profVMaxTemp;
end

function kSection = kAtPoint(sIn, DDX, DDY)
%LOCALKATPOINT This local function finds and returns the curvature at a given point in terms of total distance traveled
% -------------------------------------------------------------------------
%
% INPUT: sIn — Set of cumulative distances traveled to determine 
%               curvatures for
%        DDX — DDX spline in terms of cumulative distance traveled
%        DDY — DDY spline in terms of cumulative distance traveled
%
% OUTPUT: kSetion — Set of curvatures returned for input distances
%
% -------------------------------------------------------------------------

% Calculate all curavtures at once using vector inputs to ppval functions
kSection = sqrt(sum([[ppval(DDX, sIn)].^2, [ppval(DDY, sIn)].^2],2));

end