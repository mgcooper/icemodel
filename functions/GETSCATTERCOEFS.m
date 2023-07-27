function [g,qext,ss_coalb,wavelength] = GETSCATTERCOEFS(opts,scattercoefs)
%GETSCATTERCOEFS get the scattering coefficients
%
% read in mie.dat and store g, qext, omega, and lambda

g = scattercoefs(1:opts.nradii,1:opts.nwavl);
qext = scattercoefs(opts.nradii+1:2*opts.nradii,1:opts.nwavl);
ss_coalb = scattercoefs(2*opts.nradii+1:3*opts.nradii,1:opts.nwavl);
wavelength = scattercoefs(3*opts.nradii+1:4*opts.nradii,1:opts.nwavl);