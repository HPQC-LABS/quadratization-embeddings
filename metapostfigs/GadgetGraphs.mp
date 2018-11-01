% GadgetGraphs.mp -- drawing the graphs of gadgets.
% Szilárd Szalay, Wigner Research Centre for Physics, Budapest.
%
% compiling:
% simple compiling: 
%     mpost GadgetGraphs.mp
% if there are LaTeX labels:
%     TEX=pdflatex mpost GadgetGraphs.mp
% The script a2ping is for converting .eps figures to .pdf figures,
% while not screwing up the bounding box.
% It gives an error message, but it works fine...
% One can use the bash script for multiple conversions:
%     for f in fig_Gadget*.eps; do ./a2ping.pl $f; done
%     << or >>
%     for f in fig_Gadget*.eps; do epspdf $f; done

% setting the sizes (also in DWaveGraphs.mp) are needed

% ** declarations for metaposting **********************************************
prologues := 3;	% needed if there are LaTeX labels and also for converting to pdf

% ** units *********************************************************************

numeric ut; ut := 24bp; 

numeric dt; dt := .2ut;			% size of a vertex
numeric it; it := 1.0bp; % 0.7bp	% width of an edge

% ** drawing routines **********************************************************

% default vertex for data qubit
vardef grvertexdata(expr A) =
	draw A withpen pencircle scaled dt withcolor black;
enddef;

% default vertex for aux qubit
vardef grvertexaux(expr A) =
	draw A withpen pencircle scaled dt withcolor black;
%	draw A withpen pencircle scaled (dt-2it) withcolor red;
	draw A withpen pencircle scaled (dt-it) withcolor red;
enddef;

% default edge for the graph
vardef gredge(expr A,B) =
	draw A--B withpen pencircle scaled it withcolor black;
enddef;



% ** 3+1 gadgets ***************************************************************

outputtemplate := "fig_Gadget31_AllToAux.eps";
beginfig(0)
pair A,B;
	A := (0,0);
	for i = 0 upto 2:
		B := ut*dir(i*360/3+90);
		gredge(A,B);
		grvertexdata(B);
	endfor;
	grvertexaux(A);
endfig;


outputtemplate := "fig_Gadget31_AllToAuxT.eps";
beginfig(0)
pair A,B;
	A := (0,0);
	for i = 1 upto 3:
		B := ut*dir(i*360/4+90);
		gredge(A,B);
		grvertexdata(B);
	endfor;
	grvertexaux(A);
endfig;



% ** 4+1 gadgets ***************************************************************

outputtemplate := "fig_Gadget41_AllToAux.eps";
beginfig(0)
pair A,B;
	A := (0,0);
	for i = 0 upto 3:
		B := ut*dir(i*360/4);
		gredge(A,B);
		grvertexdata(B);
	endfor;
	grvertexaux(A);
endfig;

outputtemplate := "fig_Gadget41_AllToAuxv.eps";
beginfig(0)
pair A,B;
	A := (0,0);
	for i = 0 upto 3:
		B := ut*dir(i*360/4+45);
		gredge(A,B);
		grvertexdata(B);
	endfor;
	grvertexaux(A);
endfig;

outputtemplate := "fig_Gadget41_Tripod.eps";
beginfig(0)
pair A,B;
	A := (0,0);
	for i = 0 upto 2:
		B := ut*dir(i*360/4);
		gredge(A,B);
		grvertexdata(B);
	endfor;
	B := ut*dir(3*360/4);
	gredge(A,B);
	grvertexaux(B);
	grvertexdata(A);
endfig;

outputtemplate := "fig_Gadget41_K5.eps";
beginfig(0)
pair A,B;
	A := (0,0);
	for i = 0 upto 3:
		B := ut*dir(i*360/4+45);
		gredge(A,B);
		gredge(B,ut*dir((i+1)*360/4+45));
	endfor;
	for i = 0 upto 3:
		B := ut*dir(i*360/4+45);
		grvertexdata(B);
	endfor;
	grvertexaux(A);
endfig;



% ** 4+2 gadgets ***************************************************************

pair P[];
for i = 0 upto 1:
	P[3i] := .25ut*dir(i*360/2);
	for j = 1 upto 2:
		P[3i+j] := P[3i] + ut*dir((2i+j)*360/4-135);
	endfor;
endfor;

outputtemplate := "fig_Gadget42_AllToAuxp1.eps";
beginfig(0)
	gredge(P[5],P[1]);
	for i = 0 upto 1:
		for j = 1 upto 2:
			gredge(P[3i],P[3i+j]);
			gredge(P[3(1-i)],P[3i+j]);
		endfor;
	endfor;
	for i = 0 upto 1:
		for j = 1 upto 2:
			grvertexdata(P[3i+j]);
		endfor;
		grvertexaux(P[3i]);
	endfor;
endfig;

outputtemplate := "fig_Gadget42_K6m1.eps";
beginfig(0)
	gredge(P[1],P[2]);
	gredge(P[2],P[4]);
	gredge(P[4],P[5]);
	gredge(P[5],P[1]);
	for i = 0 upto 1:
		for j = 1 upto 2:
			gredge(P[3i],P[3i+j]);
			gredge(P[3(1-i)],P[3i+j]);
		endfor;
	endfor;
	for i = 0 upto 1:
		for j = 1 upto 2:
			grvertexdata(P[3i+j]);
		endfor;
		grvertexaux(P[3i]);
	endfor;
endfig;

outputtemplate := "fig_Gadget42_K6.eps";
beginfig(0)
	gredge(P[1],P[2]);
	gredge(P[2],P[4]);
	gredge(P[4],P[5]);
	gredge(P[5],P[1]);
	gredge(P[0],P[3]);
	for i = 0 upto 1:
		for j = 1 upto 2:
			gredge(P[3i],P[3i+j]);
			gredge(P[3(1-i)],P[3i+j]);
		endfor;
	endfor;
	for i = 0 upto 1:
		for j = 1 upto 2:
			grvertexdata(P[3i+j]);
		endfor;
		grvertexaux(P[3i]);
	endfor;
endfig;


% ** 3+3 gadgets ***************************************************************

for i = 0 upto 2:
	P[i] := ut*dir(i*360/3+90);
	P[3+i] := .5P[i];
endfor;

outputtemplate := "fig_Gadget33_K5p1.eps"; % "bad name"
beginfig(0)
	for i = 0 upto 2:
		gredge(P[i],P[(i+1)mod3]);
		gredge(P[3+i],P[3+(i+1)mod3]);
		gredge(P[4],P[i]);
		gredge(P[5],P[i]);
	endfor;
	for i = 0 upto 2:
		grvertexdata(P[i]);
		grvertexaux(P[3+i]);
	endfor;
endfig;

% ** 4+3 gadgets ***************************************************************

for i = 0 upto 1:
	for j = 0 upto 2:
		P[3i+j] := .5ut*dir(i*360/2)  + .5ut*dir(j*360/4-90 + i*360/2);
	endfor;
endfor;
P[6] := (0,0);

outputtemplate := "fig_Gadget43_DoubleK4.eps"; % "bad name"
beginfig(0)
	for i = 0 upto 1:
		for j = 0 upto 2:
			gredge(P[3i+j],P[6]);
		endfor;
		gredge(P[3i],P[3i+1]);
		gredge(P[3i+1],P[3i+2]);
	endfor;
	for ij = 0, 5, 6:
		grvertexaux(P[ij]);
	endfor;
	for ij = 1 upto 4:
		grvertexdata(P[ij]);
	endfor;
endfig;

end;
