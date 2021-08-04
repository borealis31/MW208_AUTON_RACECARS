function [vProfile, vProfByS, vCrit ,minTLap] = velocityProfilerV2(timeStep, waypointsIn, lapsQty, sLap, ftMax, fnMax, mass)
%UNTITLED Summary of this function goes here
%   Based on:
% Friedl, Tyler, "Optimized Trajectory Generation for Car-Like Robots on a Closed Loop Track" (2017).
% Electronic Theses and Dissertations. 1370.
% https://digitalcommons.du.edu/etd/1370

ftFactor = ftMax^2;
ftfnFactor = (ftMax/fnMax)^2;

for idxLap = 1:lapsQty
    if idxLap == 1
        waypointsTracked = waypointsIn;
    else
        waypointsTracked = [waypointsTracked; waypointsIn(2:end,:)];
    end
end

numPoints = length(waypointsTracked);
numProfiles = 1;

sDiff = ([0; sqrt(sum((waypointsTracked(2:end,:)-waypointsTracked(1:end-1,:)).^2,2))]);
sAtPt = cumsum(sDiff);


% normFactor = floor(sAtPt(end))/sAtPt(end);
% sNormalized = sAtPt.*normFactor;
% waypointsTrackedNormalized = waypointsTracked.*normFactor;
% waypointsTracked = [spline(sNormalized, waypointsTrackedNormalized(:,1),0:sNormalized(end))',...
%                     spline(sNormalized, waypointsTrackedNormalized(:,2),0:sNormalized(end))']./normFactor;
% sDiff = ([0; sqrt(sum((waypointsTracked(2:end,:)-waypointsTracked(1:end-1,:)).^2,2))]);
% sAtPt = cumsum(sDiff);





ds = sAtPt(end)/(10*numPoints);

X = spline(sAtPt,waypointsTracked(:,1));
Y = spline(sAtPt,waypointsTracked(:,2));

DX = fnder(X,1);
DDX = fnder(X,2);
DY = fnder(Y,1);
DDY = fnder(Y,2);

sHD = ([0:ds:sAtPt(end)])';

rTrack = 1./kAtPoint(sHD, DX, DDX, DY, DDY);
rInterp = spline(0:ds:sAtPt(end), rTrack);

vCrit = fnMax*sqrt(ppval(rInterp,sHD)./mass);

vCritMins = [sHD(islocalmin(vCrit,'FlatSelection','all','MaxNumExtrema',floor(sAtPt(end)))), vCrit(islocalmin(vCrit,'FlatSelection','first','MaxNumExtrema',floor(sAtPt(end))))];

profiles = cell(length(vCritMins)+1, 1);

for idxMin = 1:length(vCritMins)+2
    if idxMin == length(vCritMins)+1
        profiles{idxMin} = generateVByS(ds, [0,0], numProfiles, sAtPt, rInterp, ftFactor, ftfnFactor, mass);
    elseif idxMin == length(vCritMins)+2
        profiles{idxMin} = generateVByS(ds, [sAtPt(end), 0], numProfiles, sAtPt, rInterp, ftFactor, ftfnFactor, mass);
    else
        profiles{idxMin} = generateVByS(ds, vCritMins(idxMin,:), numProfiles, sAtPt, rInterp, ftFactor, ftfnFactor, mass);
    end
    numProfiles = numProfiles + 1;
end

vProfComposite = [];
for idxProf = 1:numel(profiles)
    vProfComposite = [vProfComposite, profiles{idxProf}];
end

vProfByS = min(vProfComposite,[],2);

v = [0];
vIdx = 1;
t = [0];
sTraveled = [ds];
lap = 1;
ft = 0;
fn = 0;

vInterp = spline(0:ds:sAtPt(end)-ds, vProfByS);


vTime = tic;
fprintf(['Beginning Time-Velocty Calculations...' newline])
while sTraveled(end) < sAtPt(end)
    sTraveled(end+1,1) = stateRefresher(sTraveled(vIdx,1), v(vIdx,1), timeStep);
    vIdx = vIdx + 1;
    v(vIdx,1) = ppval(vInterp, sTraveled(end));
    t(vIdx,1) = t(vIdx-1,1) + timeStep;
    lap(vIdx,1) = ceil(sTraveled(end)/sLap);
end
fprintf(['Time-Velocity Calculations Complete: %f sec' newline], toc(vTime))
vVsT = spline(t, v);

ft = mass*(ppval(fnder(vVsT,1),t));
fn = mass*((ppval(vVsT,t)).^2)./(ppval(rInterp,sTraveled));

vCrit = [(0:ds:sAtPt(end))', vCrit];
vProfile = [v, t, sTraveled, ppval(X, sTraveled), ppval(Y, sTraveled), ft, fn, lap];
tBestLapSet = t(lap == ceil(lapsQty/2));
minTLap = tBestLapSet(end) - tBestLapSet(1);

end

function prof = generateVByS(ds, critMins, numProfiles, sIn, rInterp, ftFactor, ftfnFactor, mass)

vDecel = [critMins(2)];
vAccel = [critMins(2)];
decelS = (critMins(1):-ds:0)';
accelS = (critMins(1):ds:sIn(end)-ds)';

if accelS >= sIn(end)
    accelS = [];
end

vPrev = critMins(2);

for idxDecel = 2:length(decelS)
    vNew = real(sqrt(((-2*(decelS(idxDecel)-decelS(idxDecel-1))*sqrt(ftFactor-(mass/(ppval(rInterp,decelS(idxDecel-1))))*ftfnFactor*(vPrev^2)))/(mass))+vPrev^2));
    vPrev = vNew;
    vDecel(end+1,1) = vPrev;
end

vPrev = critMins(2);

for idxAccel = 2:length(accelS)
    vNew = real(sqrt(((2*(accelS(idxAccel)-accelS(idxAccel-1))*sqrt(ftFactor-(mass/(ppval(rInterp,accelS(idxAccel))))*ftfnFactor*(vPrev^2)))/(mass))+vPrev^2));
    vPrev = vNew;
    vAccel(end+1,1) = vPrev;
end



profProc = [flip(vDecel,1); vAccel(2:end)];
% profProcSplineNormalizer = spline([flip(decelS,1); accelS(2:end,1)], real(profProc));
% prof = [ppval(profProcSplineNormalizer, (0:0.01:sIn(end))'), (0:0.01:sIn(end))'];
prof = [profProc, [flip(decelS,1); accelS(2:end,1)]];
fprintf(['Profile ', num2str(numProfiles, '%.0f'), ' Generated...' newline]);
if numel(accelS) == 0 || numel(decelS) == 0
    prof = prof(1:end-1,:);
end

end

function sOut = stateRefresher(sIn, vIn, tStep)
sOut = sIn + vIn*tStep;
end

function kSection = kAtPoint(sIn, DX, DDX, DY, DDY)
kSection = zeros(length(sIn), 1);
for idxS = 1:length(sIn)
    %     kSection(idxS) = norm(cross([ppval(DX, sIn(idxS)), ppval(DY, sIn(idxS)), 0],[ppval(DDX, sIn(idxS)), ppval(DDY, sIn(idxS)), 0]))/norm([ppval(DX, sIn(idxS)), ppval(DY, sIn(idxS))]);
    kSection(idxS) = norm([ppval(DDX, sIn(idxS)), ppval(DDY, sIn(idxS)), 0]);
end
end