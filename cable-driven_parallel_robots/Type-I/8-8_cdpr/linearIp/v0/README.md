The initial version of the codes.

Execute the file 'runfor.m' in MATLAB.
Upon successful completion, it will create the following files for post-processing the results.
1. statevar.dat -- an array of positions and velocities of joints.
2. timevar.dat -- an array of time instances at which the positions and velocities are recorded in the above file.
3. posvelacc.dat -- an array of time, positions, velocities, and accelerations of the joints.
4. jtorque.dat -- an array of generalised forces due to spring and damper elements, and the generalised constraint forces. 
5. lambda.dat -- an array of Lagrange multipliers, i.e., the reaction forces at the anchoring points of the cables.
6. elLen.dat -- an array of the instantaneous lengths of the elements of cables.

