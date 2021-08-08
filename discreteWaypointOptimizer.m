function [waypointsOut, videoSetInter, videoSetPt] = discreteWaypointOptimizer(waypointsIn, widthNormalized, maxIterations, mindK, riseAllowance, bufferSize, title)
%DISCRETEWAYPOINTOPTIMIZER This script optimizes the input waypoints for curvature using an iterative patternsearch algorithm that focuses on localized curvature of a set of discrete points
% -------------------------------------------------------------------------
% 
%   INPUTS: waypointsIn — Initial normalized set of waypoints (normalized 
%                          center-points of track)
%           widthNormalized — Width of track multiplied by normalizing 
%                              factor from previous script
%           maxIteration — Maximum number of iterations the optimizer will
%                           run through before exiting
%           mindK — Minimum reduction in K that we can experience across
%                    iterations before exiting
%           riseAllowance — Maximum amount the track's curvature can
%                            increase across iterations without exiting
%                            (precautionary)
%           bufferSize — Percentage of track we reserve for a buffer region
%                         to ensure that clipping does not occur between 
%                         optimized track and original bounds
%
%   OUTPUTS: waypointsOut — Normalized, curvature optimized set of
%                            waypoints
%            videoSetInter — Cell array of normalized, curvature optimized
%                             sets of waypoints in order of iteration, also
%                             includes information on each set's iteration 
%                             and curvature
%            videoSetPt — Cell array of normalized, curvature optimized
%                          sets of waypoints in order of POINTS for EACH 
%                          iteration, also includes information on each 
%                          set's iteration, curvature, and point position
%                          after optimization
%            sFit — The length of the new raceline following optimization
% 
% -------------------------------------------------------------------------
%

% Initialize variables for optimization script
iteration = 1;
pt = 1;
dK = 1e2;
kPrev = 1e3;
minK = 1e4;
rowCount = size(waypointsIn,1);
halfTrackWidthAllowed = (widthNormalized/2)*(1-bufferSize);

% Move waypointsIn value to waypointsOut to prevent overwrite
waypointsOut = waypointsIn;

% Initialize Video Set cell containers for future output efficiency
videoSetInter = cell(maxIterations+1,4);
videoSetInter{1,1} = waypointsIn;
videoSetInter{1,2} = 0;
videoSetInter{1,3} = 999;

videoSetPt = cell((maxIterations)*length(waypointsIn)+1,3);
videoSetPt{1,1} = waypointsIn;
videoSetPt{1,2} = 0;
videoSetPt{1,3} = 999;
videoSetPt{1,4} = waypointsOut(1,:);

% Begin optimization script and execute until a criteria fails
while (iteration <= maxIterations) && abs(dK) > mindK && kPrev - minK < riseAllowance
    % Repeat optimization for every point (excluding the LAST point, see
    % below)
    for idxWP = 1:rowCount-1
        % Select the current point as well as the 2 before (M1, M2) and 2
        % after (P1, P2)
        switch idxWP
            % SPECIAL CASES: The point being optimized is the 1st, 2nd,
            % 3rd-to-last, or 2nd-to-last point. The LAST point is not
            % used, as it is equal to the 1st point
            case 1
                ptM2 = waypointsOut(end - 2, :);
                ptM1 = waypointsOut(end - 1, :);
                ptOp = waypointsOut(idxWP, :);
                ptP1 = waypointsOut(idxWP + 1, :);
                ptP2 = waypointsOut(idxWP + 2, :);
            case 2
                ptM2 = waypointsOut(end - 1, :);
                ptM1 = waypointsOut(idxWP - 1, :);
                ptOp = waypointsOut(idxWP, :);
                ptP1 = waypointsOut(idxWP + 1, :);
                ptP2 = waypointsOut(idxWP + 2, :);
            case rowCount - 1
                ptM2 = waypointsOut(idxWP - 2, :);
                ptM1 = waypointsOut(idxWP - 1, :);
                ptOp = waypointsOut(idxWP, :);
                ptP1 = waypointsOut(1, :);
                ptP2 = waypointsOut(2, :);
            case rowCount - 2
                ptM2 = waypointsOut(idxWP - 2, :);
                ptM1 = waypointsOut(idxWP - 1, :);
                ptOp = waypointsOut(idxWP, :);
                ptP1 = waypointsOut(idxWP + 1, :);
                ptP2 = waypointsOut(1, :);
            otherwise
                ptM2 = waypointsOut(idxWP - 2, :);
                ptM1 = waypointsOut(idxWP - 1, :);
                ptOp = waypointsOut(idxWP, :);
                ptP1 = waypointsOut(idxWP + 1, :);
                ptP2 = waypointsOut(idxWP + 2, :);
        end
        
        % Call the optimizing script for the given set of points — see the
        % local function declaration for more details
        points = optim([ptP2; ptP1; ptOp; ptM1; ptM2], waypointsIn, idxWP, halfTrackWidthAllowed, iteration);
        
        % Update the point's values in the waypoint set for future
        % iterations
        waypointsOut(idxWP, 1) = points(1);
        waypointsOut(idxWP, 2) = points(2);
        waypointsOut(end,:) = waypointsOut(1,:);
        
        % Output Video Set point data 
        videoSetPt{pt+1,1} = waypointsOut;
        videoSetPt{pt+1,2} = iteration;
        videoSetPt{pt+1,3} = kPrev;
        videoSetPt{pt+1,4} = [waypointsOut(idxWP, 1), waypointsOut(idxWP, 2)];
        
        % Progress point
        pt = pt + 1;
    end

    % Recalculate distance traveled at a point. Use these values as well as
    % the curvature-optimized waypoints to generate splines so that we can 
    % obtain DDX and DDY values at a point along the path using 
    % fnder and ppval
    sAtPt = cumsum([0; sqrt(sum((waypointsOut(2:end,:)-waypointsOut(1:end-1,:)).^2,2))]);
    X = spline(sAtPt,waypointsOut(:,1));
    Y = spline(sAtPt,waypointsOut(:,2));
    DDX = fnder(X,2);
    DDY = fnder(Y,2);

    % Define an equation for curvature in terms of distance traveled and
    % integrate along the distance traveled to get the total curvature of
    % the path
    curvature = @(s) norm([ppval(DDX, s), ppval(DDY, s)]);
    kCurrent = integral(@(sIn) curvature(sIn), 0, sAtPt(end),'ArrayValued',true);
    
    % Reset while loop comparison values for next iteration conditions
    dK = minK - kCurrent;
    kPrev = kCurrent;
    
    % Print out optimization progress based on iteration and curvature
    fprintf(['Iteration: ', num2str(iteration, '%.0f'), ' | k: ', num2str(kPrev, '%.8f'), newline]);
    
    % Reset the vale for minK if the k of the new path is less than it
    if kPrev < minK
        bestPoints = waypointsOut;
        bestIt = iteration;
        minK = kPrev;
        bestPtFn = pt;
    end
    
    % Output Video Set iteration data
    videoSetInter{iteration+1,1} = waypointsOut;
    videoSetInter{iteration+1,2} = iteration;
    videoSetInter{iteration+1,3} = kPrev;
    
    % Progress iteration
    iteration = iteration + 1;
    
    % If the while loops exits due to exceeding rise allowance, reset the
    % waypoints and relevant things to the best recorded state
    if ~(kPrev - minK < riseAllowance)
        fprintf(['Rise allowance exceeded; returning to best iteration: ' num2str(bestIt, '%.0f')]);
        waypointsOut = bestPoints;
        videoSetInter(bestIt+1:end,:) = [];
        videoSetPt(bestPtFn+1:end,:) = [];
        iteration = bestIt;
    end
end


% Once iterations complete, remove empty cells from video sets and reshape
% to original column size
videoSetInter = reshape(videoSetInter(~cellfun('isempty',videoSetInter)),[],3);
videoSetPt = reshape(videoSetPt(~cellfun('isempty',videoSetPt)),[],4);

end

function optimPtOut = optim(pointsIn, waypointsIn, ptIdx, halfTrackWidthAllowed, iteration)
%LOCALOPTIMPTOUT This local function uses discrete methods of calculating curvature to find new point positions via patternsearch optimization
%
% INPUTS: pointsIn — Set of local points in which the position of the
%                     middle-most point is being optimized
%         waypointsIn — Original set of waypoints used to determine the
%                        chord along which a given point can be optimized
%         ptIdx — The point number at which the script is currently
%                  optimizing
%         halfTrackWifthAllowed — Range that the point can shift along the
%                                  chord based on normalized distance and
%                                  allawble buffer from previous section
%         iteration — The number of the current iteration of curvature
%                      optimization
%
% OUTPUTS: optimPtOut — The new X,Y position of the optimized waypoint
% -------------------------------------------------------------------------
%

% Calculate the heading of the middle point so that the chord can be
% created for optimization, accounting for the 1st and final point overlap
if ptIdx == size(waypointsIn,1)
    ang = atan2(waypointsIn(1,2)-waypointsIn(ptIdx,2), waypointsIn(1,1)-waypointsIn(ptIdx,1));
else
    ang = atan2(waypointsIn(ptIdx+1,2)-waypointsIn(ptIdx,2), waypointsIn(ptIdx+1,1)-waypointsIn(ptIdx,1));
end

% Rotate and center all five points so that the chord for optimization is
% centered at the y-axis
rotPointsIn = [(pointsIn(:,1) - waypointsIn(ptIdx,1))*cos(ang) + (pointsIn(:,2) - waypointsIn(ptIdx,2))*sin(ang), ...
    -(pointsIn(:,1) - waypointsIn(ptIdx,1))*sin(ang) + (pointsIn(:,2) - waypointsIn(ptIdx,2))*cos(ang)];

% Create expressions for tangent vectors (1st derivatives) in terms of y
% (dR/ds == T)
T2 = (rotPointsIn(1,:)-rotPointsIn(2,:))/sqrt(sum((rotPointsIn(1,:)-rotPointsIn(2,:)).^2,2));
T1 = @(y) (rotPointsIn(2,:)-[0, y])/sqrt(sum((rotPointsIn(2,:)-[0, y]).^2,2));
Tn1 = @(y) ([0, y]-rotPointsIn(4,:))/sqrt(sum(([0, y]-rotPointsIn(4,:)).^2,2));
Tn2 = (rotPointsIn(4,:)-rotPointsIn(5,:))/sqrt(sum((rotPointsIn(4,:)-rotPointsIn(5,:)).^2,2));

% Use tangent vectors to create expressions for normal vectors (2nd
% derivatives) in terms of y (dT/ds == N)
N1 = @(y) T2-T1(y);
N0 = @(y) T1(y)-Tn1(y);
Nn1 = @(y) Tn1(y)-Tn2;

% Use dT/ds curvature equation to create expressions for curvature in terms
% of y (K == magnitude(dT/ds))
K1 = @(y) norm(N1(y));
K0 = @(y) norm(N0(y));
Kn1 = @(y) norm(Nn1(y));

% Use the three curvature expressions to create a fitness function for the
% patternsearch algorithm, weighting the curren point's curvature slightly
% more than the immediate previous and future curvatures
fitnessFunction = @(y) (0.5*Kn1(y) + K0(y) + 0.5*K1(y));

% If this is the first iteration, make the initial y value for the
% patternseach a random number within track constraints
if iteration == 1
    rotPointsIn(3,2) = rand(1)*halfTrackWidthAllowed;
end

% Implement patternsearch algorithm to return a y value that optimizes the
% localized curvature and suppresses its own output
[newYValue, ~, ~, ~] = patternsearch(@(y) fitnessFunction(y), rotPointsIn(3,2), [], [], [], [], -halfTrackWidthAllowed, halfTrackWidthAllowed, [], optimoptions(@patternsearch, 'MeshRotate', 'On', 'AccelerateMesh', true, 'Display', 'off'));

% Unrotate and retranslate new point to be stored in place of previous 
% point being optimized 
optimPtOut = [0, newYValue]*[cos(-ang) -sin(-ang); sin(-ang) cos(-ang)]+[waypointsIn(ptIdx,1), waypointsIn(ptIdx,2)];

end