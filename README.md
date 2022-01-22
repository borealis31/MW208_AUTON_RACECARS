# MW208_AUTON_RACECARS
An attempt at [Project 208](https://github.com/mathworks/MathWorks-Excellence-in-Innovation/tree/main/projects/Path%20Planning%20for%20Autonomous%20Race%20Cars)
from the [MathWorks Excellence in Innovation program](https://github.com/mathworks/MathWorks-Excellence-in-Innovation).<br>
For an in-depth overview and demonstration of the project, see the [YouTube Overview](https://youtu.be/WzsFHxG3lDw).

This optimizer utilizes a curvature-optimization approach that is coupled with a point-mass velocity-profile optimization method subject to force constraints and mass.

The optimizer begins by creating an initial path determined by the centerline between the outer and inner bounds. Following this, the points are passed to the discrete path optimizer. These points are then optimized along a "chord", which is a radial path from the inner to outer bound crossing through the point, using a patternsearch optimizer from MATLAB's Global Optimization Toolbox. This process is iterated until a maximum number of iterations is met or until the total curvature decreases less than the minimum decrement. See the below GIFs for a visualization of this process (discreteWaypointOptimizer.m).

Kidney Bean Track            |  Oblong Track
:-------------------------:|:-------------------------:
![Optimization Kidney Bean Visualization GIF](https://github.com/borealis31/MW208_AUTON_RACECARS/blob/main/KidneyBeanTest.gif)  |  ![Optimization Oblong Visualization GIF](https://github.com/borealis31/MW208_AUTON_RACECARS/blob/main/OblongTest.gif)

Once these paths are optimized according to the provided track constraints, a velocity profile is generated based on point-mass kinematics constraints. Using maximum/critical velocity calculations, we can identify minimum velocity locations and generate acceleration and deceleration profiles for after and before each minimum velocity. These profiles are then all compared against each other, and for each distance on every profile, the minimum velcoty is selected in order to ensure controlability along the track while minimizing lap time. All of this data is returned and can be used to generate a video that showcases the position of the mass in real-time. Best lap values are also returned based on the middle-most lap and plotted using a colored scatter plot to show the velocity along the path. See the below pictures for a visualization of each track's maximum velocity gradient (velocityProfilerV2.m). We can also plot each velocity and critical velocity by its distance along the track.

Kidney Bean Track            |  Oblong Track
:-------------------------:|:-------------------------:
![Velocity Kidney Bean Visualization PNG](https://github.com/borealis31/MW208_AUTON_RACECARS/blob/main/velocityGradientBestLapKidneyBeanTestEXAMPLE.png)  |  ![Velocity Oblong Visualization PNG](https://github.com/borealis31/MW208_AUTON_RACECARS/blob/main/velocityGradientBestLapOblongTestEXAMPLE.png)
![VvsS Kidney Bean Visualization PNG](https://github.com/borealis31/MW208_AUTON_RACECARS/blob/main/velocityByDistanceKidneyBeanTestEXAMPLE.png) | ![VvsS Oblong Visualization PNG](https://github.com/borealis31/MW208_AUTON_RACECARS/blob/main/velocityByDistanceOblongTestEXAMPLE.png) 
Maximum Tangential Force: 1500N | Maximum Tangential Force: 2400N
Maximum Normal Force: 1500N | Maximum Normal Force: 900N
Mass: 300kg | Mass: 100kg
Maximum Velocity: 100m/s | Maximum Velocity: 175m/s
