README

***1. SYSTEM REQUIREMENTS***

All software dependencies and operating systems (including version numbers)：
Microsoft Windows 10

Software version that has been tested on：
R version 4.3.1 (2023-06-16)

Any required non-standard hardware：
Not applicable.


***2. INSTALLATION GUIDE***

Not applicable.


***3 INSTRUCTIONS***

Before running, make sure "LATEst.R" and "simulation-data-func.R" files are in "src/" folder.

Run "make-figures.R" from "src/" folder.

Before running, set the root path (line 6).

Run line 13 to import "LATEst.func" function.

Run line 14 to import "simulation.data.func" function.

Lines 16-81 generate Figure 1b.

Lines 82-195 generate source data for Figure 1c.

Lines 196-322 generate Figure 2.

Lines 323-367 generate Figure 3a.

Lines 368-420 generate Figure 4b and Figure 4e.

Lines 421-471 generate Figure 5b and Figure 5e.

Results are saved in "./output" folder.


***4. POSTSCRIPT***

Our R code implementation generates fully reproducible results, because random seeds are set to fixed values to guarantee identical randomization.

Total expected run time on a 'normal' desktop computer:
1 hours
