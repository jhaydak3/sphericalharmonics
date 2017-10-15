# sphericalharmonics
Matlab code for producing stl files/3d graphs of Atomic orbitals and Spherical Harmonics. I made these a while ago, so some of these descriptions might be slightly inaccurate. Most of the code is commented.

The function AO(n,l,m) produces a 3d plot of the corresponding atomic (hydrogen) electron orbitals and produces an stl file ("Hnlmorbital.stl")of the orbital. Because these often consist of several separated surfaces, modification for stability will need to be done with a separate program.
ac
sharmonic(l,m) will produce the corresponding spherical harmonic plot.

orbitalSTLmaker is the cool function. This makes a stl file of the spherical harmonic. There are some options for improving the stability of the stl file when printing by trying to make an inner structure. This is somewhat successful, but best way to improve stability is probably to just edit it with some other CAD software. Some examples of what this function can make can be found at 


These programs make use of the publicly available stlwrite.m and surf2solid.m :
https://www.mathworks.com/matlabcentral/fileexchange/20922-stlwrite-filename--varargin-
https://www.mathworks.com/matlabcentral/fileexchange/42876-surf2solid-make-a-solid-volume-from-a-surface-for-3d-printing
