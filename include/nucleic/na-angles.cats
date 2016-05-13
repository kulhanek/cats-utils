// init master topology
var mt = new Topology("md0.parm7");
mt.printInfo();

// init trajectory pool
//var mp = new Trajectory(mt,"prod001.traj");
var mp = new TrajPool(mt);
mp.addTrajListFrom(1,2,".");
mp.printInfo();

// what should be analyzed
var ms = new RSelection(mt,":1-30");

// time
var time_unit = "ns";
var stfac = 100;  // number of snapshots in time_unit

// -----------------------------------------------------------------------------
// directory for data
var path = "data";
CATs.system(sprintf("mkdir %s",path));

// helper environment
var ss = new Snapshot(mt);

// -----------------------------------------------------------------------------
// prepare internal structures

var cvs = new OFile("cvs.in");

cvs.printf("{MAIN}\n");
cvs.printf("{PMFLIB}\n");
cvs.printf("[units]\n");
cvs.printf("angle deg\n\n");
cvs.printf("{CVS}\n");

var residues = [];

for(var i=0; i < ms.getNumOfSelectedResidues(); i++) {
        var residue = {};
        residue.topres = ms.getResidue(i);
        residue.name = sprintf("%03d%s",residue.topres.getIndex()+1,residue.topres.getName().trim());
        residue.resid = residue.topres.getIndex();

        residue.type = 0;  // middle
        if( residue.topres.getName().indexOf("5") > -1 ){
            residue.type = 1;
        }
        if( residue.topres.getName().indexOf("3") > -1 ){
            residue.type = 2;
        }        

        residue.angles = [];
        
        if( residue.type != 1 ){
            var angle = {};            
            angle.name = "alpha";
            cvs.printf("[DIH]\n");
            cvs.printf("name %s\n",residue.name.concat(angle.name));
            cvs.printf("group_a :%d@O3'\n",residue.resid);
            cvs.printf("group_b :%d@P\n",residue.resid+1);
            cvs.printf("group_c :%d@O5'\n",residue.resid+1);
            cvs.printf("group_d :%d@C5'\n",residue.resid+1);
            residue.angles.push(angle);
            
            var angle = {};
            angle.name = "beta";
            cvs.printf("[DIH]\n");
            cvs.printf("name %s\n",residue.name.concat(angle.name));
            cvs.printf("group_a :%d@P\n",residue.resid+1);
            cvs.printf("group_b :%d@O5'\n",residue.resid+1);
            cvs.printf("group_c :%d@C5'\n",residue.resid+1);
            cvs.printf("group_d :%d@C4'\n",residue.resid+1);
            residue.angles.push(angle);             
        }
        
        var angle = {};
        angle.name = "gamma";
        cvs.printf("[DIH]\n");
        cvs.printf("name %s\n",residue.name.concat(angle.name));
        cvs.printf("group_a :%d@O5'\n",residue.resid+1);
        cvs.printf("group_b :%d@C5'\n",residue.resid+1);
        cvs.printf("group_c :%d@C4'\n",residue.resid+1);
        cvs.printf("group_d :%d@C3'\n",residue.resid+1);        
        residue.angles.push(angle);
        
        var angle = {};
        angle.name = "delta";
        cvs.printf("[DIH]\n");
        cvs.printf("name %s\n",residue.name.concat(angle.name));
        cvs.printf("group_a :%d@C5'\n",residue.resid+1);
        cvs.printf("group_b :%d@C4'\n",residue.resid+1);
        cvs.printf("group_c :%d@C3'\n",residue.resid+1); 
        cvs.printf("group_d :%d@O3'\n",residue.resid+1);        
        residue.angles.push(angle);

        if( residue.type != 2 ){  
            var angle = {};
            angle.name = "epsilon";
            cvs.printf("[DIH]\n");
            cvs.printf("name %s\n",residue.name.concat(angle.name));
            cvs.printf("group_a :%d@C4'\n",residue.resid+1);
            cvs.printf("group_b :%d@C3'\n",residue.resid+1); 
            cvs.printf("group_c :%d@O3'\n",residue.resid+1);
            cvs.printf("group_d :%d@P\n",residue.resid+2);        
            residue.angles.push(angle);
            
            var angle = {};
            angle.name = "zeta";
            cvs.printf("[DIH]\n");
            cvs.printf("name %s\n",residue.name.concat(angle.name));
            cvs.printf("group_a :%d@C3'\n",residue.resid+1); 
            cvs.printf("group_b :%d@O3'\n",residue.resid+1);
            cvs.printf("group_c :%d@P\n",residue.resid+2);
            cvs.printf("group_d :%d@O5'\n",residue.resid+2);            
            residue.angles.push(angle);            
        }

	var angle = {};
        angle.name = "chi";
        cvs.printf("[DIH]\n");
        cvs.printf("name %s\n",residue.name.concat(angle.name));

        if ((residue.name.indexOf('G') === -1) && (residue.name.indexOf('A') === -1)){
            cvs.printf("group_a :%d@O4'\n",residue.resid+1); 
            cvs.printf("group_b :%d@C1'\n",residue.resid+1);
            cvs.printf("group_c :%d@N1\n",residue.resid+1);
            cvs.printf("group_d :%d@C2\n",residue.resid+1);   
	} else {
	    cvs.printf("group_a :%d@O4'\n",residue.resid+1); 
            cvs.printf("group_b :%d@C1'\n",residue.resid+1);
            cvs.printf("group_c :%d@N9\n",residue.resid+1);
            cvs.printf("group_d :%d@C4\n",residue.resid+1);   
	}          
        residue.angles.push(angle);
        
        for(var j = 0; j < residue.angles.length; j++) {
            var name = residue.name.concat(residue.angles[j].name);
            residue.angles[j].ofs = new OFile(path.concat("/").concat(name.concat(".out")));
            residue.angles[j].ofs.printf("# snapshot %s\n",name);
            residue.angles[j].ofs.printf("# -------- ----------\n");
            residue.angles[j].hst = new Histogram();
            residue.angles[j].hst.setName(name);
            residue.angles[j].hst.setMinValue(-180.0);            
            residue.angles[j].hst.setMaxValue(180.0);
            residue.angles[j].hst.setNumOfBins(180); 
	    residue.angles[j].last = 0;           
        }
        
        if( residue.type != -1 ){
            residues.push(residue);
        }
}

cvs.close();

PMFLib.init(ss,"cvs.in");
printf("\n");

// -----------------------------------------------------------------------------
// prepare internal structures

printf("# Analyzing ...\n");
printf("# ------------------------------------------------------------------------------\n");

var snap = 0;
while( mp.read(ss) ){
    snap++;    
    PMFLib.setCoordinates(ss);
    var dihid = 0;
    for(var i = 0; i < residues.length; i++) {
        for(var j = 0; j < residues[i].angles.length; j++) {
            var dih = PMFLib.getCVValue(dihid);
	  
	  if (Math.abs(residues[i].angles[j].last - dih) > 180)
	  {
	      if (dih > 0)
	      {
	          residues[i].angles[j].ofs.printf("%10d %10.3f\n\n",snap-0.5,-180);
	          residues[i].angles[j].ofs.printf("%10d %10.3f\n",snap-0.5,180);
  	      }
	      else
	      {
		  residues[i].angles[j].ofs.printf("%10d %10.3f\n\n",snap-0.5,180);
	          residues[i].angles[j].ofs.printf("%10d %10.3f\n",snap-0.5,-180);
	      }
            }

            residues[i].angles[j].ofs.printf("%10d %10.3f\n",snap,dih);
            residues[i].angles[j].hst.addSample(dih);
	    residues[i].angles[j].last = dih;
            dihid++;
        }
    }
    mp.printProgress();
}

// -----------------------------------------------------------------------------

printf("\n");
printf("# Generating figures ...\n");
printf("# ------------------------------------------------------------------------------\n");

var maxhist = 0.0;

for(var i = 0; i < residues.length; i++) {
    for(var j = 0; j < residues[i].angles.length; j++) {
        var angle = residues[i].angles[j];
        angle.ofs.close();
        angle.hst.normalize();
        if( maxhist < angle.hst.getMaxOccupancy() ){
            maxhist = angle.hst.getMaxOccupancy();
        }
        var name = residues[i].name.concat(angle.name);
        angle.hst.save(path.concat("/").concat(name).concat(".hist"));
    }
}

var gnuplot = new OFile("gnuplot.in");

gnuplot.printf("#!/usr/bin/env gnuplot\n");

//create a placeholder image (for angles that are not available)
gnuplot.printf("set term png size 1917,717 font \"arial,45\"\n");
gnuplot.printf("set output \"%s\"\n",path.concat("/placeholder.png"));
gnuplot.printf("unset xtics\n");
gnuplot.printf("unset ytics\n");
gnuplot.printf("set label 10 at 0,0 \"This torsion angle is not present in the structure.\" front center\n");
gnuplot.printf("set yrange [-1:1]\n");
gnuplot.printf("set key off\n");
gnuplot.printf("plot 2 with lines\n");

gnuplot.printf("unset label 10\n# -----------------------------------------------------------------\n");

gnuplot.printf("set terminal postscript enhanced eps color size 6.40,2.40 font \"Arial\" 14 dashlength 3\n");
gnuplot.printf("set border lw 2\n");
gnuplot.printf("set bmargin 3.5\n");
gnuplot.printf("set rmargin 2.0\n");
gnuplot.printf("set nokey\n");

for(var i = 0; i < residues.length; i++) {
    for(var j = 0; j < residues[i].angles.length; j++) {
        var angle = residues[i].angles[j];
        var name = path.concat("/").concat(residues[i].name).concat(angle.name);
        
        gnuplot.printf("# ------------------------------------------------------------------------------\n");
        gnuplot.printf("set output \"%s\"\n",name.concat(".eps"));
        gnuplot.printf("\n");
        gnuplot.printf("set multiplot layout 1,2 title \"%s - %s\"\n",residues[i].name,angle.name);
        gnuplot.printf("\n");    
        gnuplot.printf("set size 0.6,0.95\n");
        gnuplot.printf("\n");
        gnuplot.printf("set xlabel \"t [%s]\"\n",time_unit);
        gnuplot.printf("set xrange[0:%f]\n",snap/stfac);
        gnuplot.printf("set xtics autofreq\n");
        gnuplot.printf("set format x \"%%.0f\"\n");
        gnuplot.printf("\n");
        gnuplot.printf("set ylabel \"%s [{/Symbol \\260}]\"\n",angle.name);        
        gnuplot.printf("set yrange[-180:180]\n");
        gnuplot.printf("set ytics 60\n");        
        gnuplot.printf("set format y \"%%.0f\"\n"); 
        gnuplot.printf("\n");
        gnuplot.printf("plot '%s' using ($1/%f):2 with lines\n",name.concat(".out"),stfac);          
        gnuplot.printf("\n");      
        gnuplot.printf("set origin 0.6,0.0\n");
        gnuplot.printf("set size 0.4,0.95\n");
        gnuplot.printf("\n");
        gnuplot.printf("set xlabel \"%s [{/Symbol \\260}]\"\n",angle.name);          
        gnuplot.printf("set xrange[-180:180]\n");
        gnuplot.printf("set xtics 60\n");
        gnuplot.printf("set format x \"%%.0f\"\n");
        gnuplot.printf("\n");
        gnuplot.printf("set ylabel \"occupancy\"\n"); 
        gnuplot.printf("set yrange[0:%f]\n",maxhist);
        gnuplot.printf("set ytics autofreq\n");        
        gnuplot.printf("set format y \"%%.3f\"\n");        
        gnuplot.printf("\n");
        gnuplot.printf("plot '%s' with boxes\n",name.concat(".hist"));
        gnuplot.printf("\n");
        gnuplot.printf("unset multiplot\n");
        gnuplot.printf("\n");        
    }
}

var dtypes = [ "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi" ];

for(var t=0; t < dtypes.length; t++ ){
    gnuplot.printf("# ------------------------------------------------------------------------------\n");
    gnuplot.printf("set output \"%s\"\n",path.concat("/").concat(dtypes[t]).concat(".eps"));

    gnuplot.printf("set xlabel \"%s [{/Symbol \\260}]\"\n",dtypes[t]);          
    gnuplot.printf("set xrange[-180:180]\n");
    gnuplot.printf("set xtics 60\n");
    gnuplot.printf("set format x \"%%.0f\"\n");
    gnuplot.printf("\n");
    gnuplot.printf("set ylabel \"occupancy\"\n"); 
    gnuplot.printf("set yrange[0:%f]\n",maxhist);
    gnuplot.printf("set ytics autofreq\n");        
    gnuplot.printf("set format y \"%%.3f\"\n");        
    gnuplot.printf("\n");
    gnuplot.printf("plot ");

    var prev = false;
    for(var i = 0; i < residues.length; i++) {
        for(var j = 0; j < residues[i].angles.length; j++) {
            var angle = residues[i].angles[j];
            if( angle.name == dtypes[t] ) {
                var name = path.concat("/").concat(residues[i].name).concat(angle.name);
                if( prev == true ){
                    gnuplot.printf(",\\\n");
                }
                gnuplot.printf("     '%s' with lines",name.concat(".hist"));
                prev = true;
            }
        }
    }
    gnuplot.printf("\n");
    gnuplot.printf("\n");
}


gnuplot.close();

CATs.system("gnuplot gnuplot.in");

for(var t=0; t < dtypes.length; t++ ){
    printf("converting %s ...\n",path.concat("/").concat(dtypes[t]).concat(".eps"));
    CATs.system(sprintf("convert  -density 300x300 %s -units PixelsPerInch -density 300 -background white -flatten %s",
                        path.concat("/").concat(dtypes[t]).concat(".eps"),path.concat("/").concat(dtypes[t]).concat(".png")));
}

for(var i = 0; i < residues.length; i++) {
    for(var j = 0; j < residues[i].angles.length; j++) {
        var angle = residues[i].angles[j];
        var name = path.concat("/").concat(residues[i].name).concat(angle.name);
        printf("converting %s ...\n",name.concat(".eps"));
        CATs.system(sprintf("convert  -density 300x300 %s -units PixelsPerInch -density 300 -background white -flatten %s",name.concat(".eps"),name.concat(".png")))
    }
}

// -----------------------------------------------------------------------------

printf("\n");
printf("# Creating the final report (this may take a few minutes)...\n");
printf("# ------------------------------------------------------------------------------\n");

var report = new OFile("report.tex");

report.printf("\\documentclass[a4paper,oneside]{book}\n\\usepackage{graphicx}\n\\usepackage{fancyhdr}\n\\usepackage[explicit]{titlesec}\n\\usepackage{hyperref}\n\\usepackage{caption}\n\\renewcommand{\\chaptername}{}\n\\renewcommand{\\thechapter}{}\n\\usepackage[left=0.1in,right=0.1in,top=1.0in,bottom=1.0in,includefoot]{geometry}\n\\newcommand{\\cchapter}[1]{\\chapter[#1]{\\centering #1}}\n\\usepackage{tocloft}\n\\renewcommand{\\cfttoctitlefont}{\\hfill\\huge\\bfseries}\n\\renewcommand{\\cftaftertoctitle}{\\hfill}\n\n\\begin{document}\n\\pagestyle{fancy}\n\\title{Dihedral Angle Analysis}\n\\author{CATs}\n\\maketitle\n\\tableofcontents\n\n");

for(var i = 0; i < residues.length; i++) {

    var fileName1 = path.concat("/placeholder.png");
    var fileName2 = path.concat("/placeholder.png");
    var fileName3 = path.concat("/placeholder.png");
    var fileName4 = path.concat("/placeholder.png");
    var fileName5 = path.concat("/placeholder.png");
    var fileName6 = path.concat("/placeholder.png");
    var fileName7 = path.concat("/placeholder.png");
    
    var counter = 0;

if (residues[i].angles.length > counter)
    if (residues[i].angles[counter].name == "alpha")
    {
	fileName1 = path.concat("/").concat(residues[i].name).concat(residues[i].angles[counter].name).concat(".png");
	counter++;
    }

if (residues[i].angles.length > counter)
    if (residues[i].angles[counter].name == "beta")
    {
	fileName2 = path.concat("/").concat(residues[i].name).concat(residues[i].angles[counter].name).concat(".png");
	counter++;
    }

if (residues[i].angles.length > counter)
    if (residues[i].angles[counter].name == "gamma")
    {
	fileName3 = path.concat("/").concat(residues[i].name).concat(residues[i].angles[counter].name).concat(".png");
	counter++;
    }

if (residues[i].angles.length > counter)
    if (residues[i].angles[counter].name == "delta")
    {
	fileName4 = path.concat("/").concat(residues[i].name).concat(residues[i].angles[counter].name).concat(".png");
	counter++;
    }

if (residues[i].angles.length > counter)
    if (residues[i].angles[counter].name == "epsilon")
    {
	fileName5 = path.concat("/").concat(residues[i].name).concat(residues[i].angles[counter].name).concat(".png");
	counter++;
    }

if (residues[i].angles.length > counter)
    if (residues[i].angles[counter].name == "zeta")
    {
	fileName6 = path.concat("/").concat(residues[i].name).concat(residues[i].angles[counter].name).concat(".png");
	counter++;
    }

if (residues[i].angles.length > counter)
    if (residues[i].angles[counter].name == "chi")
    {
	fileName7 = path.concat("/").concat(residues[i].name).concat(residues[i].angles[counter].name).concat(".png");
    }

    report.printf("\\cchapter{Residue %s}\n\n",residues[i].name);
    
    report.printf("\\vspace{11pt}\n");
    report.printf("\\begin{figure}[h!]\n");
    report.printf("\\centering\n");
    report.printf("\\includegraphics[width=0.45\\linewidth]{%s}\n",fileName1);
    report.printf("$\\hspace{0.3cm}$\n");
    report.printf("\\includegraphics[width=0.45\\linewidth]{%s}\n",fileName2);
    report.printf("$\\vspace{0.75cm}$\n");
    report.printf("\\includegraphics[width=0.45\\linewidth]{%s}\n",fileName3);
    report.printf("$\\hspace{0.3cm}$\n");
    report.printf("\\includegraphics[width=0.45\\linewidth]{%s}\n",fileName4);
    report.printf("$\\vspace{0.75cm}$\n");
    report.printf("\\includegraphics[width=0.45\\linewidth]{%s}\n",fileName5);
    report.printf("$\\hspace{0.3cm}$\n");
    report.printf("\\includegraphics[width=0.45\\linewidth]{%s}\n",fileName6);
    report.printf("$\\vspace{0.75cm}$\n");
    report.printf("\\includegraphics[width=0.45\\linewidth]{%s}\n",fileName7);
    report.printf("\\caption*{Torsion angles for residue %s.}\n",residues[i].name);
    report.printf("\\end{figure}\n\n");

    report.printf("{\\let\\clearpage\\relax}\n\n");
} 


report.printf("\\cchapter{Aggregate results (Part 1)}\n\n");
    
report.printf("\\vspace{11pt}\n");
report.printf("\\begin{figure}[h!]\n");
report.printf("\\centering\n");
report.printf("\\includegraphics[width=0.60\\linewidth]{%s}\n",path.concat("/alpha.png"));
report.printf("$\\vspace{0.5cm}$\n");
report.printf("\\includegraphics[width=0.60\\linewidth]{%s}\n",path.concat("/beta.png"));
report.printf("$\\vspace{0.5cm}$\n");
report.printf("\\includegraphics[width=0.60\\linewidth]{%s}\n",path.concat("/gamma.png"));
report.printf("\\caption*{Aggregate results.}\n");
report.printf("\\end{figure}\n\n");

report.printf("{\\let\\clearpage\\relax}\n\n");

report.printf("\\cchapter{Aggregate results (Part 2)}\n\n");

report.printf("\\vspace{11pt}\n");
report.printf("\\begin{figure}[h!]\n");
report.printf("\\centering\n");
report.printf("\\includegraphics[width=0.60\\linewidth]{%s}\n",path.concat("/delta.png"));
report.printf("$\\vspace{0.5cm}$\n");
report.printf("\\includegraphics[width=0.60\\linewidth]{%s}\n",path.concat("/epsilon.png"));
report.printf("$\\vspace{0.5cm}$\n");
report.printf("\\includegraphics[width=0.60\\linewidth]{%s}\n",path.concat("/zeta.png"));
report.printf("\\caption*{Aggregate results.}\n");
report.printf("\\end{figure}\n\n");

report.printf("{\\let\\clearpage\\relax}\n\n");

report.printf("\\cchapter{Aggregate results (Part 3)}\n\n");

report.printf("\\vspace{11pt}\n");
report.printf("\\begin{figure}[h!]\n");
report.printf("\\centering\n");
report.printf("\\includegraphics[width=0.60\\linewidth]{%s}\n",path.concat("/chi.png"));
report.printf("\\caption*{Aggregate results.}\n");
report.printf("\\end{figure}\n\n");

report.printf("\\end{document}\n");

report.close();

CATs.system("pdflatex -interaction=nonstopmode report.tex");
CATs.system("pdflatex -interaction=nonstopmode report.tex");
