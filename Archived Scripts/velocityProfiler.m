function [vProfile, minTLap] = velocityProfiler(timeStep, waypointsIn, sPath, ftMax, fnMax, mass)
%UNTITLED Summary of this function goes here
%   Based on:
% Friedl, Tyler, "Optimized Trajectory Generation for Car-Like Robots on a Closed Loop Track" (2017).
% Electronic Theses and Dissertations. 1370.
% https://digitalcommons.du.edu/etd/1370

v = [1e-6];
t = [0];
vIdx = 1;
sTraveled = [0];
lap = [1];
ft = [0];
fn = [0];

sDiff = ([0; sqrt(sum((waypointsIn(2:end,:)-waypointsIn(1:end-1,:)).^2,2))]);
sAtPt = cumsum(sDiff);

X = spline(sAtPt,waypointsIn(:,1));
Y = spline(sAtPt,waypointsIn(:,2));
DX = fnder(X,1);
DDX = fnder(X,2);
DY = fnder(Y,1);
DDY = fnder(Y,2);

vCrit = [fnMax*sqrt(1/(mass*kAtPoint(0, DX, DDX, DY, DDY)))];

while sTraveled(end) < sAtPt(end)
    sTraveled(end+1,1) = stateRefresher(sTraveled(vIdx,1), v(vIdx,1), timeStep);
    vIdx = vIdx + 1;
    vCrit(vIdx,1) = fnMax*sqrt(1/(mass*kAtPoint(sTraveled(end), DX, DDX, DY, DDY)));
    v(vIdx,1) = velocityGenerator(v(vIdx-1,1), sTraveled, ftMax, fnMax, mass, DX, DDX, DY, DDY);
    t(vIdx,1) = t(vIdx-1,1) + timeStep;
    lap(vIdx,1) = ceil(sTraveled(end)/sPath);
    ft(vIdx, 1) = mass*(v(vIdx,1)-v(vIdx-1,1))/timeStep;
    fn(vIdx, 1) = mass*(v(vIdx,1)^2)*kAtPoint(sTraveled(vIdx,1), DX, DDX, DY, DDY);
end

vProfile = [v, t, sTraveled, ppval(X, sTraveled), ppval(Y, sTraveled), ft, fn, vCrit];

tLap2 = t(lap == 2);

minTLap = t(end) - tLap2(end);

end

function vNew = velocityGenerator(vPrev, sIn, ftMax, fnMax, mass, DXIn, DDXIn, DYIn, DDYIn)

ftFactor = ftMax^2;
ftfnFactor = (ftMax/fnMax)^2;

tempV = sqrt(((2*(sIn(end,1)-sIn(end-1,1))*sqrt(ftFactor-mass*kAtPoint(sIn(end-1,1), DXIn, DDXIn, DYIn, DDYIn)*ftfnFactor*vPrev^2))/(mass))+vPrev^2);

if imag(tempV) == 0
    vNew = tempV;
else
    vNew = sqrt(((-2*(sIn(end,1)-sIn(end-1,1))*sqrt(ftFactor+mass*kAtPoint(sIn(end-1,1), DXIn, DDXIn, DYIn, DDYIn)*ftfnFactor*vPrev^2))/(mass))+vPrev^2);
end

if vNew > 40
    vNew = 40;
end

end

function sOut = stateRefresher(sIn, vIn, tStep)
sOut = sIn + vIn*tStep;
end

function kSection = kAtPoint(sIn, DX, DDX, DY, DDY)
kSection = norm(cross([ppval(DX, sIn), ppval(DY, sIn), 0],[ppval(DDX, sIn), ppval(DDY, sIn), 0]));
end