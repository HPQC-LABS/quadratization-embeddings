% DWaveGraphs.mp -- drawing the Chimera and the Pegasus graph.
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
numeric ut; ut := 32bp;		% base unit for the whole figure
numeric uh; uh := ut;		% horizontal
numeric uv; uv := 0.4*ut; %ut;	% vertical

% base units (for full graph)
%numeric vt; vt := 2.7ut; 	% base unit
numeric vt; vt := 2ut; 		% base unit
numeric vh; vh := vt;		% horizontal
%numeric vv; vv := 4/3*vt;	% vertical
numeric vv; vv := 0.8vt;	% vertical

numeric vs; vs := 1.1;		% scale factor for the L-by-L graph
%numeric vs; vs := 1.2;		% scale factor for the L-by-L graph

% tilting factors (for the unit cells)
numeric th; th :=-0.4;		% horizontal
numeric tv; tv :=-0.1;		% vertical




% ** drawing routines **********************************************************

% default vertex
vardef grvertex(expr A) =
	draw A withpen pencircle scaled .1ut;
enddef;

% default edge for the graph
vardef gredge(expr A,B,clr) =
	draw A--B withpen pencircle scaled .7bp withcolor clr;
enddef;




% ** definition of the vertex positions ***************************************

% calculates che coordinates of vertices, based on the .pdf dated Sept 8, 2018
vardef PV(expr c,i,j,u,p,q) =

	if drawingmode = 1:	% this is the 4+1 unit cell mode
		numeric pq;
		if ((i+j)mod2=0) and (u=0) and (c=2):
			pq := 2*p+q;
		else:
			pq := 2*(1-p)+(1-q);	
		fi;
		% position of the unit cell:
		c*(vh,vv)			% c: index of the Chimera layers
		- i*(2vh,0)			% i: horizontal index of the unit cells
		- j*(0,2vv)			% j: vertical index of the unit cells
		% position of the vertex inside the unit cell:
		+ (u-0.5)*(uh,0)		% u: horizontal index of the vertices inside the unit cells
		+ (pq-1.5)*(0,uv)		% 2*p+q: vertical index of the vertices inside the unit cells
		+ (c-1)*( (-1)**(i+j)*(pq-1.5)*th*uh, (-1)**u*tv*uv )	% tilting of the unit cells
	elseif drawingmode = 2:	% this is the two-layer L-by-L unit cell mode
		% position of the unit cell:
		c*vs*(.75vh,1.23vv)		% c: index of the Chimera layers
%		- i*vs*(2vh,0)			% i: horizontal index of the unit cells
		- i*vs*(1.7vh,0.5vv)		% i: horizontal index of the unit cells
		- j*vs*(0,2vv)			% j: vertical index of the unit cells
		% position of the vertex inside the unit cell:
		+ (u-0.5)*(.5uh,0)		% u: horizontal index of the vertices inside the unit cells
		- (2*p+q-1.5)*(0,uv)		% 2*p+q: vertical index of the vertices inside the unit cells
		+ ( -(2*p+q-1.5)*th*.5*uh, (-1)**u*tv*uv )				% tilting of the unit cells
	elseif drawingmode = 3:	% this is the two-layer L-by-L unit cell mode with diamond shaped unit cells
		pair ucl;
%		if (u<>p):
%			ucl := ((-1)**u)* (2*p+q-1.5)*(u,(1-u))*ut*.4 rotated -90;
%		else:
%			ucl := ((-1)**(1-u))* (2*(1-p)+(1-q)-1.5)*(1-u,u)*ut*.4 rotated -90;
%		fi;
		% position of the unit cell:
		c*vs*(0.75vt,0.75vt)		% c: index of the Chimera layers
		- i*vs*(1.5vt,0)		% i: horizontal index of the unit cells
		- j*vs*(0,1.5vt)		% j: vertical index of the unit cells
		% position of the vertex inside the unit cell:
		+ ((-1)**u)* (2*p+q-1.5)*(u,(1-u))*ut*.4 rotated +90		% 2*p+q:  u: vertical or horizontal
%		+ ucl
	fi
enddef;

% only five selected unit cells for the 4+1 mode
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


numeric L; L := 3;	% size of the Chimera blocks, $L\times L$



% ** Pegasus graph ************************************************************

% drawing edges defined by the equations, for given unit cell coordinates c,i,j
vardef edgeeq(expr eq, c,i,j) =
	if eq = "8":		% equation (8): K_{4,4} edges in unit cells
					for p = 0 upto 1:
					for pp = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(c,i,j,0,p,q), PV(c,i,j,1,pp,qp), cgrey );
						endfor;
						endfor;
					endfor;
					endfor;
	elseif eq = "9":	% equation (9): vertical couplings
					for p = 0 upto 1:
						for q = 0 upto 1:
							gredge( PV(c,i,j,0,p,q), PV(c,i,j+1,0,p,q), cgrey );
						endfor;
					endfor;
	elseif eq = "10":	% equation (10): vertical couplings
					for p = 0 upto 1:
						for q = 0 upto 1:
							gredge( PV(c,i,j,1,p,q), PV(c,i+1,j,1,p,q), cgrey );
						endfor;
					endfor;
	elseif eq = "1":	% equation (1): additional edges in unit cells
				for u = 0 upto 1:
					for p = 0 upto 1:
							gredge( PV(c,i,j,u,p,0), PV(c,i,j,u,p,1), cblack );
					endfor;
				endfor;
	elseif eq = "2a":	% equation (2) for u=0
					for pp = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(c,i,j,0,0,q), PV((c+1)mod3,i,j,1,pp,qp), cblue );
						endfor;
						endfor;
					endfor;
	elseif eq = "2b":	% equation (2) for u=1
					for pp = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(c,i,j,1,0,q), PV((c+1)mod3,i,j,0,pp,qp), cred );
						endfor;
						endfor;
					endfor;
	elseif eq = "3":	% equation (3)
					for pp = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(c,i,j,1,1,q), PV((c+1)mod3,i,j+1,0,pp,qp), cgreen ); %corange );
						endfor;
						endfor;
					endfor;
	elseif eq = "4":	% equation (4)
					for pp = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(c,i,j,0,1,q), PV((c+1)mod3,i+1,j,1,pp,qp), cgreen );
						endfor;
						endfor;
					endfor;
	fi;
enddef;


% drawing vertices for unit cell of coordinates c,i,j
vardef ucvertex(expr c,i,j) =
				for u = 0 upto 1:
					for p = 0 upto 1:
						for q = 0 upto 1:
							grvertex( PV(c,i,j,u,p,q) );
%							label(decimal(u)&decimal(i)&decimal(j), PV(c,i,j,u,p,q) );
						endfor;
					endfor;
				endfor;
enddef;



% setting drawingmode to the 2x2+1 layout
numeric drawingmode; drawingmode := 1;

picture ChimeraTT;
ChimeraTT = image(
	numeric mi,c,i,j;
	for mi = 1 upto 4:
	c := selected[mi][0];
		i := selected[mi][1];
			j := selected[mi][2];
				edgeeq("8",c,i,j);
	endfor;
	for mi = 1, 3:
	c := selected[mi][0];
		i := selected[mi][1];
			j := selected[mi][2];
				edgeeq("9",c,i,j);
	endfor;
	for mi = 1, 2:
	c := selected[mi][0];
		i := selected[mi][1];
			j := selected[mi][2];
				edgeeq("10",c,i,j);
	endfor;
	for mi = 1 upto 4:
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
); % of ChimeraTT



picture PegasusTT;
PegasusTT = image(
	draw ChimeraTT;
	numeric mi,c,i,j;
	mi := 0;
	c := selected[mi][0];
		i := selected[mi][1];
			j := selected[mi][2];
				edgeeq("8",c,i,j);
	for mi = 0 upto 4:
	c := selected[mi][0];
		i := selected[mi][1];
			j := selected[mi][2];
				edgeeq("1",c,i,j);
	endfor;
	mi := 0;
	c := selected[mi][0];
		i := selected[mi][1];
			j := selected[mi][2];
				edgeeq("2a",c,i,j);
				edgeeq("2b",c,i,j);
				edgeeq("3",c,i,j);
				edgeeq("4",c,i,j);
); % of PegasusTT


picture ChimeraTTDots;
ChimeraTTDots = image(
	numeric mi,c,i,j;
	for mi = 1 upto 4:
	c := selected[mi][0];
		i := selected[mi][1];
			j := selected[mi][2];
				ucvertex(c,i,j);
	endfor;
); % of ChimeraTTDots

picture PegasusTTDots;
PegasusTTDots = image(
	draw ChimeraTTDots;
	numeric mi,c,i,j;
	mi := 0;
	c := selected[mi][0];
		i := selected[mi][1];
			j := selected[mi][2];
				ucvertex(c,i,j);
); % of PegasusTTDots


% output
outputtemplate := "fig_ChimeraTT.eps";
beginfig(0)
	draw ChimeraTT;
	draw ChimeraTTDots;
endfig;


outputtemplate := "fig_PegasusTT.eps";
beginfig(0)
	draw PegasusTT;
	draw PegasusTTDots;
endfig;




% setting drawingmode to the 2xLxL layout
drawingmode := 2;

picture ChimeraLL;
ChimeraLL := image(
	numeric c;
	c := 2;
		for i = 0 upto L-1:
			for j = 0 upto L-1:
				edgeeq("8",c,i,j);
			endfor;
		endfor;
		for i = 0 upto L-1:
			for j = 0 upto L-2:
				edgeeq("9",c,i,j);
			endfor;
		endfor;
		for i = 0 upto L-2:
			for j = 0 upto L-1:
				edgeeq("10",c,i,j);
			endfor;
		endfor;


); % of ChimeraLL



picture PegasusLL;
PegasusLL := image(
	numeric c;
	for c = 1 upto 2:
		for i = 0 upto L-1:
			for j = 0 upto L-1:
				edgeeq("8",c,i,j);
				edgeeq("1",c,i,j);
			endfor;
		endfor;
	endfor;
	for c = 1 upto 2:
		for i = 0 upto L-1:
			for j = 0 upto L-2:
				edgeeq("9",c,i,j);
			endfor;
		endfor;
	endfor;
	for c = 1 upto 2:
		for i = 0 upto L-2:
			for j = 0 upto L-1:
				edgeeq("10",c,i,j);
			endfor;
		endfor;
	endfor;
	c := 1;
		for i = 0 upto L-1:
			for j = 0 upto L-1:
				edgeeq("2a",c,i,j);
				edgeeq("2b",c,i,j);
			endfor;
		endfor;
		for i = 0 upto L-1:
			for j = 0 upto L-2:
				edgeeq("3",c,i,j);
			endfor;
		endfor;
		for i = 0 upto L-2:
			for j = 0 upto L-1:
				edgeeq("4",c,i,j);
			endfor;
		endfor;
); % of PegasusLL


picture ChimeraLLDots;
ChimeraLLDots := image(
	numeric c;
	c := 2;
		for i = 0 upto L-1:
			for j = 0 upto L-1:
				ucvertex(c,i,j);
			endfor;
		endfor;
); % of ChimeraLLDots



picture PegasusLLDots;
PegasusLLDots := image(
	for c = 1 upto 2:
		for i = 0 upto L-1:
			for j = 0 upto L-1:
				ucvertex(c,i,j);
			endfor;
		endfor;
	endfor;
); % of PegasusLLDots


outputtemplate := "fig_ChimeraLL.eps";
beginfig(0)
	draw ChimeraLL;
	draw ChimeraLLDots;
endfig;

outputtemplate := "fig_PegasusLL.eps";
beginfig(0)
	draw PegasusLL;
	draw PegasusLLDots;
endfig;












drawingmode := 3;


picture PegasusLL;
PegasusLL := image(
	numeric c;
	for c = 1 upto 2:
		for i = 0 upto L-1:
			for j = 0 upto L-1:
				edgeeq("8",c,i,j);
				edgeeq("1",c,i,j);
			endfor;
		endfor;
	endfor;
	for c = 1 upto 2:
		for i = 0 upto L-1:
			for j = 0 upto L-2:
				edgeeq("9",c,i,j);
			endfor;
		endfor;
	endfor;
	for c = 1 upto 2:
		for i = 0 upto L-2:
			for j = 0 upto L-1:
				edgeeq("10",c,i,j);
			endfor;
		endfor;
	endfor;
	c := 1;
		for i = 0 upto L-1:
			for j = 0 upto L-1:
				edgeeq("2a",c,i,j);
				edgeeq("2b",c,i,j);
			endfor;
		endfor;
		for i = 0 upto L-1:
			for j = 0 upto L-2:
				edgeeq("3",c,i,j);
			endfor;
		endfor;
		for i = 0 upto L-2:
			for j = 0 upto L-1:
				edgeeq("4",c,i,j);
			endfor;
		endfor;
); % of PegasusLL

picture PegasusLLDots;
PegasusLLDots := image(
	for c = 1 upto 2:
		for i = 0 upto L-1:
			for j = 0 upto L-1:
				ucvertex(c,i,j);
			endfor;
		endfor;
	endfor;
); % of PegasusLLDots

outputtemplate := "fig_PegasusLLd.eps";
beginfig(0)
	draw PegasusLL;
	draw PegasusLLDots;
endfig;

end;
