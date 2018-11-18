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

% for writing out the indices:
numeric showindices;
showindices := 0;       % hide the indices
%showindices := 1;      % show the indices

if showindices=1:
        input TEX;
        TEXPRE("%&latex" & char(10) & "\documentclass{article}\begin{document}");
        TEXPOST("\end{document}");
fi;


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


% ** gadget definitions ********************************************************

string gdlabel[];
numeric gdNvert[];
pair gdvert[][];
numeric gdverttype[][];		% 0: data, 1: aux
numeric gdNedge[];
numeric gdedge[][][];

picture gadget[];

numeric gd,gdN;
gd := 0;

def gdadd(expr gdl)=
	gd := gd+1;
	gdlabel[gd] := gdl;
	gdNedge[gd] := 0;
	gdNvert[gd] := 0;
enddef;

def gdaddvert(expr A,t)=
	gdvert[gd][ gdNvert[gd] ] := A;
	gdverttype[gd][ gdNvert[gd] ] := t;
	gdNvert[gd] := gdNvert[gd]+1;
	
enddef;

def gdaddedge(expr ei,ej)=
	gdedge[gd][ gdNedge[gd] ][0] := ei;
	gdedge[gd][ gdNedge[gd] ][1] := ej;
	gdNedge[gd] := gdNedge[gd]+1;
enddef;


% ** 3+1 gadgets ***************************************************************


gdadd("31_AllToAux");
gdaddvert((0,0),1);
for i = 0 upto 2:
	gdaddvert(ut*dir(i*360/3+90),0);
	gdaddedge(0,i+1);
endfor;


gdadd("31_AllToAuxT");
gdaddvert((0,0),1);
for i = 0 upto 2:
	gdaddvert(ut*dir(i*360/4+180),0);
	gdaddedge(0,i+1);
endfor;


gdadd("31_T");
gdaddvert((0,0),0);
for i = 0 upto 2:
	gdaddvert(ut*dir(i*360/4+180),if i = 1: 1 else: 0 fi);
	gdaddedge(0,i+1);
endfor;


gdadd("31_Tv");
gdaddvert((0,0),0);
for i = 0 upto 2:
	gdaddvert(ut*dir((i+1)*360/3+90),if i = 2: 1 else: 0 fi);
	gdaddedge(0,i+1);
endfor;


gdadd("31_Coathanger");
gdaddvert((0,0),1);
for i = 0 upto 2:
	gdaddvert(ut*dir(i*360/3+90),0);
	gdaddedge(0,i+1);
endfor;
gdaddedge(2,3);


gdadd("31_K4m1");
gdaddvert((0,0),1);
for i = 0 upto 2:
	gdaddvert(ut*dir(i*360/3+90),0);
	gdaddedge(0,i+1);
endfor;
gdaddedge(1,2);
gdaddedge(1,3);


gdadd("31_K4");
gdaddvert((0,0),1);
for i = 0 upto 2:
	gdaddvert(ut*dir(i*360/3+90),0);
	gdaddedge(0,i+1);
	gdaddedge(i+1,(i+1) mod 3 +1);
endfor;


% ** 4+1 gadgets ***************************************************************


gdadd("41_AllToAux");
gdaddvert((0,0),1);
for i = 0 upto 3:
	gdaddvert(ut*dir(i*360/4),0);
	gdaddedge(0,i+1);
endfor;


gdadd("41_AllToAuxv");
gdaddvert((0,0),1);
for i = 0 upto 3:
	gdaddvert(ut*dir(i*360/4+45),0);
	gdaddedge(0,i+1);
endfor;


gdadd("41_Tripod");
gdaddvert((0,0),0);
for i = 0 upto 3:
	gdaddvert(if i mod 2 = 1: (1/sqrt(2))* fi ut*dir(i*360/4 if i mod 2 = 0: + (i-1)*360/8 fi),if i = 3: 1 else: 0 fi);
	gdaddedge(0,i+1);
endfor;


gdadd("41_Tripodv");
gdaddvert((0,0),0);
for i = 0 upto 3:
	gdaddvert(ut*dir(i*360/4),if i = 3: 1 else: 0 fi);
	gdaddedge(0,i+1);
endfor;


gdadd("41_K5");
gdaddvert((0,0),1);
for i = 0 upto 3:
	gdaddvert(ut*dir(i*360/4+45),0);
	gdaddedge(0,i+1);
	gdaddedge(i+1,(i+1) mod 4 +1);
endfor;


% ** 4+2 gadgets ***************************************************************


gdadd("42_AllToAuxp1");
for i = 0 upto 1:
	gdaddvert(.25ut*dir(i*360/2),1);
	for j = 1 upto 2:
		gdaddvert(.25ut*dir(i*360/2) + ut*dir((2i+j)*360/4-135),0);
		gdaddedge(3i,3i+j);
		gdaddedge(3(1-i),3i+j);
	endfor;
endfor;
gdaddedge(5,1);


gdadd("42_K6m1");
for i = 0 upto 1:
	gdaddvert(.25ut*dir(i*360/2),1);
	for j = 1 upto 2:
		gdaddvert(.25ut*dir(i*360/2) + ut*dir((2i+j)*360/4-135),0);
		gdaddedge(3i,3i+j);
		gdaddedge(3(1-i),3i+j);
	endfor;
endfor;
gdaddedge(1,2);
gdaddedge(2,4);
gdaddedge(4,5);
gdaddedge(5,1);


gdadd("42_K6");
for i = 0 upto 1:
	gdaddvert(.25ut*dir(i*360/2),1);
	for j = 1 upto 2:
		gdaddvert(.25ut*dir(i*360/2) + ut*dir((2i+j)*360/4-135),0);
		gdaddedge(3i,3i+j);
		gdaddedge(3(1-i),3i+j);
	endfor;
endfor;
gdaddedge(1,2);
gdaddedge(2,4);
gdaddedge(4,5);
gdaddedge(5,1);
gdaddedge(0,3);


% ** 3+3 gadgets ***************************************************************


gdadd("33_K5p1");
for i = 0 upto 2:
	gdaddvert(ut*dir(i*360/3+90),0);
	gdaddvert(.5ut*dir(i*360/3+90),1);
	gdaddedge(2i,2*((i+1)mod3));
	gdaddedge(2i+1,2*((i+1)mod3)+1);
	gdaddedge(3,2*i);
	gdaddedge(5,2*i);
endfor;


% ** 4+3 gadgets ***************************************************************

gdadd("43_DoubleK4v");
for i = 0 upto 1:
	for j = 0 upto 2:
		gdaddvert(.5ut*dir(i*360/2) + .5ut*dir(j*360/4-90 + i*360/2), if j=1: 1 else: 0 fi);
		gdaddedge(3i+j,6);
		gdaddedge(3i+j,3i);
	endfor;
	gdaddedge(3i+1,3i+2);
endfor;
gdaddvert((0,0),1);



gdadd("43_DoubleK4");
for i = 0 upto 1:
	gdaddvert(.5ut*dir(i*360/2),1);
	for j = 0 upto 2:
		if j<>0: gdaddvert(.5ut*dir(i*360/2)+ut*dir(j*360/4-135+i*360/2),0) fi;
		gdaddedge(3i+j,6);
		gdaddedge(3i+j,3i);
	endfor;
	gdaddedge(3i+1,3i+2);
endfor;
gdaddvert((0,0),1);


% ** drawing the figures and output ********************************************

gdN := gd;

for gd = 1 upto gdN:
gadget[gd]:= image(
	for ie = 0 upto gdNedge[gd]-1:
		gredge( gdvert[gd][gdedge[gd][ie][0]], gdvert[gd][gdedge[gd][ie][1]] );
	endfor;
	for iv = 0 upto gdNvert[gd]-1:
		if gdverttype[gd][iv] = 0:
			grvertexdata(gdvert[gd][iv]);
		elseif gdverttype[gd][iv] = 1:
			grvertexaux(gdvert[gd][iv]);
		fi;
		if showindices=1: label.urt(TEX("\tiny $"&decimal(iv)&"$"), gdvert[gd][iv] ); fi;
	endfor;
);
endfor;


for gd = 1 upto gdN:
	show gdlabel[gd];
	outputtemplate := "fig_Gadget"&gdlabel[gd]&".eps";
	beginfig(0) draw gadget[gd]; endfig;
endfor;

end;

