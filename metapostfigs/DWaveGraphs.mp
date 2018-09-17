% DWaveGraphs.mp -- drawing the Pegasus graph.
% Szilárd Szalay, Wigner Research Centre for Physics,
%
% compiling:
% the required components mpost and a2ping are contained 
% in the standard TexLive distribution of Ubuntu
% simple compiling: 
%     mpost DWaveGraphs.mp
% if there are LaTeX labels:
%     TEX=latex mpost DWaveGraphs.mp
% The script a2ping is for converting .eps figures to .pdf figures,
% while not screwing up the bounding box.
% It drops an error message, but it works fine...
% One can use the bash script for multiple conversions:
%     for f in fig*.eps; do a2ping $f; done
%

% ** declarations for metaposting **********************************************
prologues := 3;	% needed if there are LaTeX labels and also for converting to pdf (as far as I remember...)


% ** color definitions *********************************************************
color cblack;	cblack := ( 0*16+ 0, 0*16+ 0, 0*16+ 0)/256;	% hsv:   ?,  0,  0
color cgrey;	cgrey  := ( 7*16+ 0, 7*16+ 0, 7*16+ 0)/256;	% hsv:   
color cred;	cred   := (15*16+15, 3*16+ 6, 3*16+ 6)/256;	% hsv:   0, 79,100
color cblue;	cblue  := ( 8*16+ 4, 8*16+ 4,15*16+15)/256;	% hsv: 240, 22, 67
color cgreeno;	cgreeno:= ( 9*16+10,15*16+15, 9*16+10)/256;	% hsv: 120, 40,100

color cgreen;	cgreen := ( 0*16+ 0, 9*16+10, 4*16+ 2)/256;	% hsv: 146,100, 60
color corange;	corange:= (15*16+15, 9*16+15, 0*16+ 0)/256;	% hsv:  37,100,100

color textgrey; textgrey := (35,31,32)/256;	% for beamer poster/slide



% ** units *********************************************************************

% base units (for unit cells)
numeric ut; ut := 28bp;		% base unit for the whole figure
numeric uh; uh := ut;		% horizontal
numeric uv; uv := ut; %3/4*ut;	% vertical

% base units (for full graph)
numeric vt; vt := 3ut; 		% base unit
numeric vh; vh := vt;		% horizontal
numeric vv; vv := 4/3*vt;	% vertical

% tilting factors (for the unit cells)
numeric th; th :=-0.3;		% horizontal
numeric tv; tv :=+0.1;		% vertical




% ** drawing routines **********************************************************

% default vertex
vardef grvertex(expr A) =
	draw A withpen pencircle scaled .15ut;
enddef;

% default edge for the graph
vardef gredge(expr A,B,clr) =
	draw A--B withpen pencircle scaled .7bp withcolor clr;
enddef;




% ** definition of the vertex positions ***************************************

% calculates che coordinates of vertices, based on the .pdf dated Sept 8, 2018
vardef PV(expr c,i,j,u,p,q) =
	c*(vh,vv)		% c: label of the Chimera graphs (horizontal index)
	- i*(2vh,0)		% i: horizontal index of the unit cells
	- j*(0,2vv)		% j: vertical index of the unit cells
	+ u*(uh,0)		% u: horizontal index of the vertices inside the unit cells
	- (2*p+q)*(0,uv)	% 2*p+q: vertical index of the vertices inside the unit cells
	+ ( (2*p+q)*th*uh, (-1)**u*tv*uv )	% tilting of the unit cells
%	+ ( (2*p+q)*th*uh, u*tv*uv )	% tilting of the unit cells
enddef;

% only five selected unit cells
numeric selected[][];
selected[0][0] := 1;	% c
selected[0][1] := 1;	% i
selected[0][2] := 1;	% j ...
selected[1][0] := 2;
selected[1][1] := 1;
selected[1][2] := 1;
selected[2][0] := 2;
selected[2][1] := 1;
selected[2][2] := 2;
selected[3][0] := 2;
selected[3][1] := 2;
selected[3][2] := 1;
selected[4][0] := 2;
selected[4][1] := 2;
selected[4][2] := 2;




% ** Pegasus graph ************************************************************



% three Chimera graphs with modified K_{4,4} unit cells
picture ChimerasX;
ChimerasX = image(
	numeric c,i,j;
	% equation (8): K_{4,4} edges in unit cells
	for mi = 0 upto 4:
	c := selected[mi][0];
		i := selected[mi][1];
			j := selected[mi][2];
					for p = 0 upto 1:
					for pp = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(c,i,j,0,p,q), PV(c,i,j,1,pp,qp), cgrey );
						endfor;
						endfor;
					endfor;
					endfor;
	endfor;
	% equation (1): additional edges in unit cells
	for mi = 0 upto 4:
	c := selected[mi][0];
		i := selected[mi][1];
			j := selected[mi][2];
				for u = 0 upto 1:
					for p = 0 upto 1:
							gredge( PV(c,i,j,u,p,0), PV(c,i,j,u,p,1), cblack );
					endfor;
				endfor;
	endfor;
	% equation (9): vertical couplings
	for mi = 1, 3:
	c := selected[mi][0];
		i := selected[mi][1];
			j := selected[mi][2];
					for p = 0 upto 1:
						for q = 0 upto 1:
							gredge( PV(c,i,j,0,p,q), PV(c,i,j+1,0,p,q), cgrey );
						endfor;
					endfor;
	endfor;
	% equation (10):  horizontal couplings
	for mi = 1, 2:
	c := selected[mi][0];
		i := selected[mi][1];
			j := selected[mi][2];
					for p = 0 upto 1:
						for q = 0 upto 1:
							gredge( PV(c,i,j,1,p,q), PV(c,i+1,j,1,p,q), cgrey );
						endfor;
					endfor;
	endfor;
); % of ChimeraX




% all the dots for the vertices
picture Dots;
Dots = image(
	numeric c,i,j;
	for mi = 0 upto 4:
	c := selected[mi][0];
		i := selected[mi][1];
			j := selected[mi][2];
				for u = 0 upto 1:
					for p = 0 upto 1:
						for q = 0 upto 1:
							grvertex( PV(c,i,j,u,p,q) );
						endfor;
					endfor;
				endfor;
	endfor;
); % of Dots




% extra edges for the Pegasus
picture PegasusEE;
PegasusEE = image(
	% ** equation (2) for u=0 ******
					for pp = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(1,1,1,0,0,q), PV(2,1,1,1,pp,qp), cblue );
						endfor;
						endfor;
					endfor;
	% ** equation (2) for u=1 ******
					for pp = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(1,1,1,1,0,q), PV(2,1,1,0,pp,qp), cred );
						endfor;
						endfor;
					endfor;
	% ** equation (3) **************
					for pp = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(1,1,1,1,1,q), PV(2,1,2,0,pp,qp), corange );
						endfor;
						endfor;
					endfor;
	% ** equation (4) **************
					for pp = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(1,1,1,0,1,q), PV(2,2,1,1,pp,qp), cgreen );
						endfor;
						endfor;
					endfor;

); % of PegasusEE



% output
outputtemplate := "fig_Pegasus.eps";
beginfig(0)
	draw ChimerasX;
	draw PegasusEE;
	draw Dots;
endfig;



end;
