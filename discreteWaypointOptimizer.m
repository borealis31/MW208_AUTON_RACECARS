function [waypointsOut, videoSetInter, videoSetPt, sFit] = discreteWaypointOptimizer(waypointsIn, widthNormalized, maxIterations, mindK, riseAllowance, bufferSize, title)
%DISCRETEWAYPOINTOPTIMIZER Summary of this function goes here
%   Detailed explanation goes here
iteration = 1;
pt = 1;
kPrev = 1e3;
minK = 1e4;
dK = 1e2;

rowCount = size(waypointsIn,1);

waypointsOut = waypointsIn;

videoSetInter = cell(maxIterations+1,4);
videoSetInter{1,1} = waypointsIn;
videoSetInter{1,2} = 0;
videoSetInter{1,3} = 999;

videoSetPt = cell((maxIterations)*length(waypointsIn)+1,3);
videoSetPt{1,1} = waypointsIn;
videoSetPt{1,2} = 0;
videoSetPt{1,3} = 999;
videoSetPt{1,4} = waypointsOut(1,:);

while (iteration <= maxIterations) && abs(dK) > mindK && kPrev - minK < riseAllowance
    for idxWP = 1:rowCount-1
        switch idxWP
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
        points = optim([ptP2; ptP1; ptOp; ptM1; ptM2], widthNormalized, waypointsIn, idxWP, bufferSize, iteration);
%         points = optim2(waypointsOut, waypointsIn, widthNormalized, idxWP, bufferSize)
        waypointsOut(idxWP, 1) = points(1);
        waypointsOut(idxWP, 2) = points(2);
        waypointsOut(end,:) = waypointsOut(1,:);
        
        videoSetPt{pt+1,1} = waypointsOut;
        videoSetPt{pt+1,2} = iteration;
        videoSetPt{pt+1,3} = kPrev;
        videoSetPt{pt+1,4} = [waypointsOut(idxWP, 1), waypointsOut(idxWP, 2)];
        
        pt = pt + 1;
    end

    arcLenAtPt = cumsum([0; sqrt(sum((waypointsOut(2:end,:)-waypointsOut(1:end-1,:)).^2,2))]);
    X = spline(arcLenAtPt,waypointsOut(:,1));
    Y = spline(arcLenAtPt,waypointsOut(:,2));
    DX = fnder(X,1);
    DDX = fnder(X,2);
    DY = fnder(Y,1);
    DDY = fnder(Y,2);

    curvature = @(s) (norm(cross([ppval(DX, s), ppval(DY, s), 0],[ppval(DDX, s), ppval(DDY, s), 0])))/(norm([ppval(DX, s), ppval(DY, s), 0]))^3;
    kCurrent = integral(@(sIn) curvature(sIn), 0, arcLenAtPt(end),'ArrayValued',true);
    dK = kPrev - kCurrent;
    kPrev = kCurrent;
    fprintf(['Iteration: ', num2str(iteration, '%.0f'), ' | k: ', num2str(kPrev, '%.8f'), newline]);
    
    if kPrev < minK
        minK = kPrev;
    end
    
    sFit = arcLenAtPt(end);
    
    videoSetInter{iteration+1,1} = waypointsOut;
    videoSetInter{iteration+1,2} = iteration;
    videoSetInter{iteration+1,3} = kPrev;
    
    iteration = iteration + 1;
end
videoSetInter = reshape(videoSetInter(~cellfun('isempty',videoSetInter)),[],3);
videoSetPt = reshape(videoSetPt(~cellfun('isempty',videoSetPt)),[],4);
end

function optimPtOut = optim(pointsIn, widthNormalized, waypointsIn, ptIdx, bufferSize, it)

    if ptIdx == size(waypointsIn,1)
        ang = atan2(waypointsIn(1,2)-waypointsIn(ptIdx,2), waypointsIn(1,1)-waypointsIn(ptIdx,1));   
    else
        ang = atan2(waypointsIn(ptIdx+1,2)-waypointsIn(ptIdx,2), waypointsIn(ptIdx+1,1)-waypointsIn(ptIdx,1));   
    end
    
    rotPointsIn = [(pointsIn(:,1) - waypointsIn(ptIdx,1))*cos(ang) + (pointsIn(:,2) - waypointsIn(ptIdx,2))*sin(ang), ...
                   -(pointsIn(:,1) - waypointsIn(ptIdx,1))*sin(ang) + (pointsIn(:,2) - waypointsIn(ptIdx,2))*cos(ang)];

    T2 = (rotPointsIn(1,:)-rotPointsIn(2,:))/sqrt(sum((rotPointsIn(1,:)-rotPointsIn(2,:)).^2,2));
    T1 = @(y) (rotPointsIn(2,:)-[0, y])/sqrt(sum((rotPointsIn(2,:)-[0, y]).^2,2));
    Tn1 = @(y) ([0, y]-rotPointsIn(4,:))/sqrt(sum(([0, y]-rotPointsIn(4,:)).^2,2));
    Tn2 = (rotPointsIn(4,:)-rotPointsIn(5,:))/sqrt(sum((rotPointsIn(4,:)-rotPointsIn(5,:)).^2,2));
    
    N1 = @(y) T2-T1(y);
    N0 = @(y) T1(y)-Tn1(y);
    Nn1 = @(y) Tn1(y)-Tn2;
    
    K1 = @(y) norm(N1(y));
    K0 = @(y) norm(N0(y));
    Kn1 = @(y) norm(Nn1(y));
    
    fitnessFunction = @(y) (0.75*Kn1(y) + K0(y) + 0.75*K1(y));
    
    if it == 1
        rotPointsIn(3,2) = rand(1)*(widthNormalized/2)*(1-bufferSize);
    end
    
    [newYValue, ~, ~, ~] = patternsearch(@(y) fitnessFunction(y), rotPointsIn(3,2), [], [], [], [], (-widthNormalized/2)*(1-bufferSize), (widthNormalized/2)*(1-bufferSize), [], optimoptions(@patternsearch, 'MeshRotate', 'On', 'AccelerateMesh', true, 'Display', 'off'));
    
    optimPtOutRot = [0, newYValue];
    
    optimPtOut = optimPtOutRot*[cos(-ang) -sin(-ang); sin(-ang) cos(-ang)]+[waypointsIn(ptIdx,1), waypointsIn(ptIdx,2)];
        
end

function optimPtOut = optim2(pointsIn, waypointsIn, widthNormalized, ptIdx, bufferSize)

if ptIdx == size(waypointsIn,1)
    ang = atan2(waypointsIn(1,2)-waypointsIn(ptIdx,2), waypointsIn(1,1)-waypointsIn(ptIdx,1));
else
    ang = atan2(waypointsIn(ptIdx+1,2)-waypointsIn(ptIdx,2), waypointsIn(ptIdx+1,1)-waypointsIn(ptIdx,1));
end

rotPointsIn = [(pointsIn(:,1) - waypointsIn(ptIdx,1))*cos(ang) + (pointsIn(:,2) - waypointsIn(ptIdx,2))*sin(ang), ...
    -(pointsIn(:,1) - waypointsIn(ptIdx,1))*sin(ang) + (pointsIn(:,2) - waypointsIn(ptIdx,2))*cos(ang)];

ptCount = length(pointsIn);

if ptIdx == 1
    rotPointsInFn = @(y) [0,y; rotPointsIn(2:end,:)];
elseif ptIdx == ptCount
    rotPointsInFn = @(y) [rotPointsIn(1:end-1,:); 0,y];
else
    rotPointsInFn = @(y) [rotPointsIn(1:ptIdx-1,:); 0,y; rotPointsIn(ptIdx+1:end,:)];
end

arcLenDiff = @(y) [0; sqrt(sum((subsref(rotPointsInFn(y), struct('type', '()', 'subs', {{2:ptCount,':'}})) - subsref(rotPointsInFn(y), struct('type', '()', 'subs', {{1:ptCount-1,':'}}))).^2,2))];
arcLenAtPtFn = @(y) cumsum(arcLenDiff(y));

%     arcLenAtPt = @(y) cumsum([0; sqrt(sum((waypointsOut(2:end,:)-waypointsOut(1:end-1,:)).^2,2))]);
X = @(y) spline(arcLenAtPtFn(y),subsref(rotPointsInFn(y), struct('type', '()', 'subs', {{':',1}})));
Y = @(y) spline(arcLenAtPtFn(y),subsref(rotPointsInFn(y), struct('type', '()', 'subs', {{':',2}})));
DX = @(y) fnder(X(y),1);
DDX = @(y) fnder(X(y),2);
DY = @(y) fnder(Y(y),1);
DDY = @(y) fnder(Y(y),2);

curvature = @(y, sIn) norm([ppval(DDX(y), sIn), ppval(DDY(y), sIn)])/norm([ppval(DX(y), sIn), ppval(DDY(y), sIn)]);
fitnessFunction = @(y) integral(@(sIn) 0.001*(curvature(y, sIn))^2, 0, subsref(arcLenAtPtFn(y), struct('type', '()', 'subs', {{ptCount}})), 'WayPoints', arcLenAtPtFn(y), 'ArrayValued',true,'AbsTol',1e-4,'RelTol',1e-2);

%[newYValue, ~, ~, out] = ga(@(y) fitnessFunction(y), 1, [], [], [], [], (-widthNormalized/2)*(1-bufferSize), (widthNormalized/2)*(1-bufferSize), [], optimoptions(@ga, 'MaxTime',20));
[newYValue, ~, ~, out] = patternsearch(@(y) fitnessFunction(y), rotPointsIn(ptIdx,2), [], [], [], [], (-widthNormalized/2)*(1-bufferSize), (widthNormalized/2)*(1-bufferSize), [], optimoptions(@patternsearch,'Cache','off'));

optimPtOutRot = [0, newYValue];

optimPtOut = optimPtOutRot*[cos(-ang) -sin(-ang); sin(-ang) cos(-ang)]+[waypointsIn(ptIdx,1), waypointsIn(ptIdx,2)];

end
















