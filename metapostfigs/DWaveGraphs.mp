% DWaveGraphs.mp -- drawing the Pegasus graph.
% Szilárd Szalay, Wigner Research Centre for Physics,
% 2018.09.11.
%
% compiling:
% the required components mpost and a2ping are contained 
% in the standard TexLive distribution of Ubuntu
% simple compiling: 
%     mpost DWaveGraphs.mp
% if there are LaTeX labels:
%     TEX=latex mpost DWaveGraphs.mp
% The script ./a2ping.pl is for converting .eps figures to .pdf figures,
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
numeric uv; uv := 3/4*ut;	% vertical

% base units (for full graph)
numeric vt; vt := 3ut; 		% base unit
numeric vh; vh := vt;		% horizontal
numeric vv; vv := 4/3*vt;	% vertical

% tilting factors (for the unit cells)
numeric th; th :=-.2;		% horizontal
numeric tv; tv := .2;		% vertical

% tilting factors (for the full graph)
%numeric tgh; tgh :=-.2;	% horizontal
%numeric tgv; tgv := .2;	% vertical




% ** drawing routines **********************************************************

% default vertex
vardef grvertex(expr A) =
	draw A withpen pencircle scaled .15ut;
enddef;

% default edge for the graph
vardef gredge(expr A,B,clr) =
	draw A--B withpen pencircle scaled .7bp withcolor clr;
enddef;

vardef gredged(expr A,B,clr) =
	draw A--B withpen pencircle scaled .7bp withcolor clr dashed evenly;
enddef;




% ** definition of the vertex positions ***************************************

numeric L; L := 3;		% size of the Chimera blocks, $L\times L$

% calculates che coordinates of vertices, based on the .pdf dated Sept 8, 2018
vardef PV(expr c,i,j,u,p,q) =
	c*L*(vh,0)		% c: label of the Chimera graphs (horizontal index)
	+ i*(vh,0)		% i: horizontal index of the unit cells
	+ j*(0,vv)		% j: vertical index of the unit cells
	+ u*(uh,0)		% u: horizontal index of the vertices inside the unit cells
	+ (2*p+q)*(0,uv)	% 2*p+q: vertical index of the vertices inside the unit cells
	+ ( (2*p+q)*th*uh, (-1)**u*tv*uv )	% tilting of the unit cells
enddef;




% ** Pegasus graph ************************************************************



% three Chimera graphs with modified K_{4,4} unit cells
picture ChimerasX;
ChimerasX = image(
	% equation (8): K_{4,4} edges in unit cells
	for c = 0 upto 2:
		for i = 0 upto L-1:
			for j = 0 upto L-1:
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
		endfor;
	endfor;
	% equation (1): additional edges in unit cells
	for c = 0 upto 2:
		for i = 0 upto L-1:
			for j = 0 upto L-1:
				for u = 0 upto 1:
					for p = 0 upto 1:
							gredge( PV(c,i,j,u,p,0), PV(c,i,j,u,p,1), cgrey );
					endfor;
				endfor;
			endfor;
		endfor;
	endfor;
	% equation (9): vertical couplings
	for c = 0 upto 2:
		for i = 0 upto L-1:
			for j = 0 upto L-2:	% !
					for p = 0 upto 1:
						for q = 0 upto 1:
							gredge( PV(c,i,j,0,p,q), PV(c,i,j+1,0,p,q), cgrey );
						endfor;
					endfor;
			endfor;
		endfor;
	endfor;
	% equation (10):  horizontal couplings
	for c = 0 upto 2:
		for i = 0 upto L-2:	% !
			for j = 0 upto L-1:
					for p = 0 upto 1:
						for q = 0 upto 1:
							gredge( PV(c,i,j,1,p,q), PV(c,i+1,j,1,p,q), cgrey );
						endfor;
					endfor;
			endfor;
		endfor;
	endfor;
); % of ChimeraX




% all the dots for the vertices
picture Dots;
Dots = image(
	for c = 0 upto 2:
		for i = 0 upto L-1:
			for j = 0 upto L-1:
				for u = 0 upto 1:
					for p = 0 upto 1:
						for q = 0 upto 1:
							grvertex( PV(c,i,j,u,p,q) );
						endfor;
					endfor;
				endfor;
			endfor;
		endfor;
	endfor;
); % of Dots




% extra edges for the Pegasus
picture PegasusEE;
PegasusEE = image(
	numeric c,i,j,u,p,q;
	% ** equation (2) for u=0 ******
%	for c = 0 upto 2:
	c := 1;
%		for i = 0 upto L-1:
		i := 1;
%			for j = 0 upto L-1:
			j := 1;
				u := 0;
					for pp = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(c,i,j,u,0,q), PV((c+1)mod3,i,j,(u+1)mod2,pp,qp), cblue );
						endfor;
						endfor;
					endfor;
%			endfor;
%		endfor;
%	endfor;
	% ** equation (2) for u=1 ******
%	for c = 0 upto 2:
	c := 1;
%		for i = 0 upto L-1:
		i := 1;
%			for j = 0 upto L-1:
			j := 1;
				u :=1;
					for pp = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(c,i,j,u,0,q), PV((c+1)mod3,i,j,(u+1)mod2,pp,qp), cred );
						endfor;
						endfor;
					endfor;
%			endfor;
%		endfor;
%	endfor;
	% ** equation (3) **************
%	for c = 0 upto 2:
	c := 1;
%		for i = 0 upto L-1:
		i := 1;
%			for j = 0 upto L-2:	% !
			j := 1;
					for pp = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(c,i,j,1,1,q), PV((c+1)mod3,i,j+1,0,pp,qp), corange );
						endfor;
						endfor;
					endfor;
%			endfor;
%		endfor;
%	endfor;
	% ** equation (4) **************
%	for c = 0 upto 2:
	c := 1;
%		for i = 0 upto L-2:	% !
		i := 1;
%			for j = 0 upto L-1:
			j := 1;
					for pp = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
%							gredge( PV(c,i,j,u,p,q), PV(c,i,j,u,p,q), cblue );
							gredge( PV(c,i,j,0,1,q), PV((c+1)mod3,i+1,j,1,pp,qp), cgreen );
						endfor;
						endfor;
					endfor;
%			endfor;
%		endfor;
%	endfor;
	% ** equation (5) for u=0 ****** 
%	for c = 0 upto 2:
	c := 1;
%		for i = 0 upto L-1:
		i := 1;
%			for j = 0 upto L-1:
			j := 1;
				u := 0;
					for p = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredged( PV(c,i,j,u,p,q), PV((c-1)mod3,i,j,(u+1)mod2,0,qp), cred );
						endfor;
						endfor;
					endfor;
%			endfor;
%		endfor;
%	endfor;
	% ** equation (5) for u=1 ******
%	for c = 0 upto 2:
	c := 1;
%		for i = 0 upto L-1:
		i := 1;
%			for j = 0 upto L-1:
			j := 1;
				u := 1;
					for p = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredged( PV(c,i,j,u,p,q), PV((c-1)mod3,i,j,(u+1)mod2,0,qp), cblue );
						endfor;
						endfor;
					endfor;
%			endfor;
%		endfor;
%	endfor;
	% ** equation (6) **************
%	for c = 0 upto 2:
	c := 1;
%		for i = 0 upto L-1:
		i := 1;
%			for j = 1 upto L-1:	% !
			j := 1;
					for p = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredged( PV(c,i,j,0,p,q), PV((c-1)mod3,i,j-1,1,1,qp), corange );
						endfor;
						endfor;
					endfor;
%			endfor;
%		endfor;
%	endfor;
	% ** equation (7) **************
%	for c = 0 upto 2:
	c := 1;
%		for i = 1 upto L-1:	% !
		i := 1;
%			for j = 0 upto L-1:
			j := 1;
					for p = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredged( PV(c,i,j,1,p,q), PV((c-1)mod3,i-1,j,0,1,qp), cgreen );
						endfor;
						endfor;
					endfor;
%			endfor;
%		endfor;
%	endfor;

); % of PegasusEE



% output
outputtemplate := "fig_Pegasus.eps";
beginfig(0)
	draw ChimerasX;
	draw PegasusEE;
	draw Dots;
endfig;



end;
