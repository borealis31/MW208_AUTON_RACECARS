function [waypointsOut, videoSetInter, videoSetPt] = discreteWaypointOptimizer(waypointsIn, widthNormalized, maxIterations, mindK, bufferSize, title)
%DISCRETEWAYPOINTOPTIMIZER Summary of this function goes here
%   Detailed explanation goes here
iteration = 1;
pt = 1;
k_prev = 999;
dK = 999;

rowCount = size(waypointsIn,1);

waypointsOut = waypointsIn;

videoSetInter = cell(maxIterations+1,3);
videoSetInter{1,1} = waypointsIn;
videoSetInter{1,2} = 0;
videoSetInter{1,3} = 999;

videoSetPt = cell((maxIterations)*length(waypointsIn)+1,3);
videoSetPt{1,1} = waypointsIn;
videoSetPt{1,2} = 0;
videoSetPt{1,3} = 999;

while (iteration <= maxIterations) && dK > mindK
    fprintf('Iteration: %f/n',iteration);
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
        waypointsOut(idxWP, 1) = points(1);
        waypointsOut(idxWP, 2) = points(2);
        waypointsOut(end,:) = waypointsOut(1,:);
        
        videoSetPt{pt+1,1} = waypointsOut;
        videoSetPt{pt+1,2} = iteration;
        videoSetPt{pt+1,3} = k_prev;
        
        pt = pt + 1;
    end

    X = spline(0:rowCount-1,[waypointsOut(:,1)]);
    Y = spline(0:rowCount-1,[waypointsOut(:,2)]);
    DX = fnder(X,1);
    DDX = fnder(X,2);
    DY = fnder(Y,1);
    DDY = fnder(Y,2);

    arcLenAtPt = cumsum([0; sqrt(sum((waypointsOut(2:end,:)-waypointsOut(1:end-1,:)).^2,2))]);
    curvature = @(s) norm(cross([ppval(DX, s), ppval(DY, s), 0],[ppval(DDX, s), ppval(DDY, s), 0]));
    k_current = integral(@(sIn) curvature(sIn), 0, arcLenAtPt(end),'ArrayValued',true);
    dK = k_prev - k_current;
    k_prev = k_current;

    videoSetInter{iteration+1,1} = waypointsOut;
    videoSetInter{iteration+1,2} = iteration;
    videoSetInter{iteration+1,3} = k_prev;
    
    iteration = iteration + 1;
end
videoSetInter = reshape(videoSetInter(~cellfun('isempty',videoSetInter)),[],3);
videoSetPt = reshape(videoSetPt(~cellfun('isempty',videoSetPt)),[],3);
end

function optimPtOut = optim(pointsIn, widthNormalized, waypointsIn, s, bufferSize, it)

    if s == size(waypointsIn,1)
        ang = atan2(waypointsIn(1,2)-waypointsIn(s,2), waypointsIn(1,1)-waypointsIn(s,1));   
    else
        ang = atan2(waypointsIn(s+1,2)-waypointsIn(s,2), waypointsIn(s+1,1)-waypointsIn(s,1));   
    end
    
    rotPointsIn = [(pointsIn(:,1) - waypointsIn(s,1))*cos(ang) + (pointsIn(:,2) - waypointsIn(s,2))*sin(ang), ...
                   -(pointsIn(:,1) - waypointsIn(s,1))*sin(ang) + (pointsIn(:,2) - waypointsIn(s,2))*cos(ang)];               
%     arcLen = @(y) cumsum([0; ...
%                           sqrt(sum((rotPointsIn(2,:)-rotPointsIn(1,:)).^2,2)); ...
%                           sqrt(sum(([0, y]-rotPointsIn(2,:)).^2,2)); ...
%                           sqrt(sum((rotPointsIn(4,:)-[0,y]).^2,2)); ...
%                           sqrt(sum((rotPointsIn(5,:)-rotPointsIn(4,:)).^2,2))]);
%     
%     X = @(y) spline(arcLen(y), [rotPointsIn(1:2,1); 0; rotPointsIn(4:5,1)]);
%     Y = @(y) spline(arcLen(y), [rotPointsIn(1:2,2); y; rotPointsIn(4:5,2)]);
%     
%     DX = @(y) fnder(X(y),1);
%     DDX = @(y) fnder(X(y),2);
%     DY = @(y) fnder(Y(y),1);
%     DDY = @(y) fnder(Y(y),2);
%     
%     curvature1 = @(y) sum(norm([ppval(DDX(y), arcLen(y)), ppval(DDY(y), arcLen(y))]));
%     curvature2 = @(y) sum(norm(cross([sum(ppval(DX(y), arcLen(y))), sum(ppval(DY(y), arcLen(y))), 0],[sum(ppval(DDX(y), arcLen(y))), sum(ppval(DDY(y), arcLen(y))), 0])));
    
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
    
%     fitnessFunction1 = @(y) ((1/3)*Kn1(y) + K0(y) + (1/3)*K1(y))*100;
    fitnessFunction2 = @(y) (0.75*Kn1(y) + K0(y) + 0.75*K1(y));
    if it == 1
        rotPointsIn(3,2) = rand(1)*(widthNormalized/2)*(1-bufferSize);
    end
%     newYValue = ga(@(y) fitnessFunction(y), 1, [], [], [], [], -widthNormalized/2, widthNormalized/2, [], optimoptions(@ga, 'FunctionTolerance', 1e-12));
    [newYValue, ~, ~, out] = patternsearch(@(y) fitnessFunction2(y), rotPointsIn(3,2), [], [], [], [], (-widthNormalized/2)*(1-bufferSize), (widthNormalized/2)*(1-bufferSize), [], optimoptions(@patternsearch, 'MeshRotate', 'On', 'AccelerateMesh', true));
%     fprintf('%f ',out.iterations)
    
    optimPtOutRot = [0, newYValue];
    
    optimPtOut = optimPtOutRot*[cos(-ang) -sin(-ang); sin(-ang) cos(-ang)]+[waypointsIn(s,1), waypointsIn(s,2)];

        
end


















