% DWaveGraphs.mp -- drawing the Chimera and the Pegasus graphs.
% Szilárd Szalay, Wigner Research Centre for Physics, Budapest.
%
% compiling:
% simple compiling: 
%     mpost DWaveGraphs.mp
% if there are LaTeX labels:
%     TEX=latex mpost DWaveGraphs.mp
% The script a2ping is for converting .eps figures to .pdf figures,
% while not screwing up the bounding box.
% It gives an error message, but it works fine...
% One can use the bash script for multiple conversions:
%     for f in fig*.eps; do a2ping $f; done
% For cropping the figures:
%     pdfcrop --margins '0 0 -160 -160' fig_ChimeraLL.pdf fig_ChimeraLL_crop.pdf
%     pdfcrop --margins '0 0 -125 -125' fig_PegasusLL.pdf fig_PegasusLL_crop.pdf
%     pdfcrop --margins '0 0 -200 -200' fig_PegasusLLfull.pdf fig_PegasusLLfull_crop.pdf

% ** declarations for metaposting **********************************************
prologues := 3;	% needed if there are LaTeX labels and also for converting to pdf

% for writing out the indices:
%input TEX;
%TEXPRE("%&latex" & char(10) & "\documentclass{article}\begin{document}");
%TEXPOST("\end{document}");


% ** color definitions *********************************************************
color cblack;	  cblack     := ( 0*16+ 0, 0*16+ 0, 0*16+ 0)/256;	% hsv:   ?,  0,  0
color cgrey;	  cgrey      := (12*16+ 0,12*16+ 0,12*16+ 0)/256;	% hsv:   
color cred;	  cred       := (15*16+15, 3*16+ 6, 3*16+ 6)/256;	% hsv:   0, 79,100
color clightblue; clightblue := ( 180    , 206    , 256    )/256;	% hsv: 
color cblue;	  cblue      := ( 8*16+ 4, 8*16+ 4,15*16+15)/256;	% hsv: 240, 22, 67
color cdarkblue;  cdarkblue  := ( 45     , 98     ,     255)/256;	% hsv: 
color cgreeno;	  cgreeno    := ( 9*16+10,15*16+15, 9*16+10)/256;	% hsv: 120, 40,100

color cgreen;	  cgreen     := ( 0*16+ 0, 9*16+10, 4*16+ 2)/256;	% hsv: 146,100, 60
color corange;	  corange    := (15*16+15, 9*16+15, 0*16+ 0)/256;	% hsv:  37,100,100



% ** units *********************************************************************

% base unit (general, and for unit cells)
numeric ut; ut := 16bp;

% base unit (for graph)
numeric vt; vt := 4.5ut;

numeric urt; urt := 00; % rotation of the unit cell


% ** drawing routines **********************************************************

% default vertex
vardef grvertex(expr A) =
	draw A withpen pencircle scaled .2ut;
enddef;

% default edge for the graph
vardef gredge(expr A,B,clr,emode) =
	if emode = "continuous":
		draw A--B withpen pencircle scaled .7bp withcolor clr;
	elseif emode = "dashed":
		draw A--B withpen pencircle scaled .7bp dashed evenly withcolor clr;
	elseif emode = "dotted":
		draw A--B withpen pencircle scaled .7bp dashed withdots withcolor clr;
	fi;
enddef;



% ** definition of the vertex positions ***************************************

% calculates che coordinates of vertices, based on the .pdf dated Sept 8, 2018
vardef PV(expr c,i,j,u,p,q) =
% this is for the 2-layer 
%	c*(vt,vt)		% c: index of the Chimera layers
%	- i*(2vt,0)		% i: horizontal index of the unit cells
%	- j*(0,2vt)		% j: vertical index of the unit cells
%	+ ((-1)**u)*(2*p+q-1.5)*(u,(1-u))*ut rotated +90		% 2*p+q:  u: vertical or horizontal
	- c*(vt,vt)	% c: index of the Chimera layers
	+ i*(3.25vt,0)		% i: horizontal index of the unit cells
	+ j*(0,3.25vt)		% j: vertical index of the unit cells
	+ ((-1)**u)*(2*p+q-1.5)*(u,(1-u))*ut rotated (-90+urt)		% 2*p+q:  u: vertical or horizontal
enddef;

% only five selected unit cells for the 4+1 mode
numeric selected[][];
selected[0][0] := 0;	% c
selected[0][1] := 1;	% i
selected[0][2] := 1;	% j ...
selected[1][0] := 1;
selected[1][1] := 1;
selected[1][2] := 1;
selected[2][0] := 1;
selected[2][1] := 1;
selected[2][2] := 2;
selected[3][0] := 1;
selected[3][1] := 2;
selected[3][2] := 1;
selected[4][0] := 1;
selected[4][1] := 2;
selected[4][2] := 2;

numeric L; L := 3;	% size of the Chimera blocks, $L\times L$




% ** Graph definitions ********************************************************

% drawing edges defined by the equations, for given unit cell coordinates c,i,j
vardef edgeeq(expr eq, c,i,j) =
	if eq = "8":		% equation (8): K_{4,4} edges in unit cells
					for p = 0 upto 1:
					for pp = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(c,i,j,0,p,q), PV(c,i,j,1,pp,qp), cgrey, "continuous" );
						endfor;
						endfor;
					endfor;
					endfor;
	elseif eq = "9":	% equation (9): vertical couplings
					for p = 0 upto 1:
						for q = 0 upto 1:
							gredge( PV(c,i,j,0,p,q), PV(c,i,j+1,0,p,q), cgrey, "continuous" );
						endfor;
					endfor;
	elseif eq = "10":	% equation (10): horizontal couplings
					for p = 0 upto 1:
						for q = 0 upto 1:
							gredge( PV(c,i,j,1,p,q), PV(c,i+1,j,1,p,q), cgrey, "continuous" );
						endfor;
					endfor;
	elseif eq = "1":	% equation (1): additional edges in unit cells
				for u = 0 upto 1:
					for p = 0 upto 1:
							gredge( PV(c,i,j,u,p,0), PV(c,i,j,u,p,1), cblack, "continuous" );
					endfor;
				endfor;
	elseif eq = "2a":	% equation (2) for u=0
					for pp = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(c,i,j,0,0,q), PV((c+1)mod3,i,j,1,pp,qp), cdarkblue, "continuous" );
						endfor;
						endfor;
					endfor;
	elseif eq = "2b":	% equation (2) for u=1
					for pp = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(c,i,j,1,0,q), PV((c+1)mod3,i,j,0,pp,qp), cdarkblue, "continuous" );
						endfor;
						endfor;
					endfor;
	elseif eq = "3":	% equation (3)
					for pp = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(c,i,j,1,1,q), PV((c+1)mod3,i,j+1,0,pp,qp), cdarkblue, "continuous" ); %corange %corange );
						endfor;
						endfor;
					endfor;
	elseif eq = "4":	% equation (4)
					for pp = 0 upto 1:
						for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(c,i,j,0,1,q), PV((c+1)mod3,i+1,j,1,pp,qp), cdarkblue, "continuous" );
						endfor;
						endfor;
					endfor;
	elseif eq = "5a":	% equation (5) for u=0
					for p = 0 upto 1:
					for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(c,i,j,0,p,q), PV((c-1)mod3,i,j,1,0,qp), clightblue, "continuous" );
						endfor;
					endfor;
					endfor;
	elseif eq = "5b":	% equation (5) for u=1
					for p = 0 upto 1:
					for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(c,i,j,1,p,q), PV((c-1)mod3,i,j,0,0,qp), clightblue, "continuous" );
						endfor;
					endfor;
					endfor;
	elseif eq = "6":	% equation (6)
					for p = 0 upto 1:
					for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(c,i,j,0,p,q), PV((c-1)mod3,i,j-1,1,1,qp), cdarkblue, "continuous" ); % green %corange );
						endfor;
					endfor;
					endfor;
	elseif eq = "7":	% equation (7)
					for p = 0 upto 1:
					for q = 0 upto 1:
						for qp = 0 upto 1:
							gredge( PV(c,i,j,1,p,q), PV((c-1)mod3,i-1,j,0,1,qp), cdarkblue, "continuous" );
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
% uncomment the following line to write out the indices
%							label.top(TEX("\tiny $"&decimal(c)&decimal(i)&decimal(j)&decimal(u)&decimal(p)&decimal(q)&"$") rotated 20, PV(c,i,j,u,p,q) );
						endfor;
					endfor;
				endfor;
enddef;



% ** 2x2+1 layout *************************************************************

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




% ** 3xLxL layout *************************************************************

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
	for c = 0 upto 2:
		for i = 0 upto L-1:
			for j = 0 upto L-1:
				edgeeq("8",c,i,j);
				edgeeq("1",c,i,j);
			endfor;
		endfor;
	endfor;
	for c = 0 upto 2:
		for i = 0 upto L-1:
			for j = 0 upto L-2:
				edgeeq("9",c,i,j);
			endfor;
		endfor;
	endfor;
	for c = 0 upto 2:
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

		for i = 0 upto L-1:
			for j = 0 upto L-1:
				edgeeq("5a",c,i,j);
				edgeeq("5b",c,i,j);
			endfor;
		endfor;
		for i = 0 upto L-1:
			for j = 1 upto L-1:
				edgeeq("6",c,i,j);
			endfor;
		endfor;
		for i = 1 upto L-1:
			for j = 0 upto L-1:
				edgeeq("7",c,i,j);
			endfor;
		endfor;
); % of PegasusLL


picture PegasusLLfull;
PegasusLLfull := image(
	numeric c;
	for c = 0 upto 2:
		for i = 0 upto L-1:
			for j = 0 upto L-1:
				edgeeq("8",c,i,j);
				edgeeq("1",c,i,j);
			endfor;
		endfor;
	endfor;
	for c = 0 upto 2:
		for i = 0 upto L-1:
			for j = 0 upto L-2:
				edgeeq("9",c,i,j);
			endfor;
		endfor;
	endfor;
	for c = 0 upto 2:
		for i = 0 upto L-2:
			for j = 0 upto L-1:
				edgeeq("10",c,i,j);
			endfor;
		endfor;
	endfor;
	for c = 0 upto 2:
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
	endfor;
); % of PegasusLLfull


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
	for c = 0 upto 2:
		for i = 0 upto L-1:
			for j = 0 upto L-1:
				ucvertex(c,i,j);
			endfor;
		endfor;
	endfor;
); % of PegasusLLDots


% output
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

outputtemplate := "fig_PegasusLLfull.eps";
beginfig(0)
	draw PegasusLLfull;
	draw PegasusLLDots;
endfig;







end;
