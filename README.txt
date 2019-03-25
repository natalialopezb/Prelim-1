Prelim 1

The following folder contains all files require to solve Prelim 1. 
A written response (mostly by-hand) can be found in the document "Soluntion.pdf"

The solution to each specific problem can be found as follows:
	
	Problem 1:
		a. "Solution.pdf" pages 1 and 2
		b. "Solution.pdf" page 2
		
	Problem 2:
		a. The figures can be obtained by running the MatLab code entitled "P2a.m". This should generate the figures shown in files "mRNA steady-state.png", "Proteins steady-state.png", "mRNA inducer.png" and "Proteins inducer.png". Comments on the elaboration of the code are found in pages 3 and 4 of "Solution.pdf".
		b. The scaled sensitivity coefficients matrixes are found in the files "BigS1.pdf", "BigS2.pdf" and "BigS3.pdf", where BigS1 corresponds to the matrix of the window without inducer, BigS2 to the matrix of the window at early-times of the inducer and BigS3 to the matrix of the window at late-times of the inducer. The matrixes can be reproduced by running first the MatLab code found in "P2b.m" followed by "P2c.m". Matrixes are stored in variables with same names.
		c. The U arrays for the solution are found in the files "U1.pdf", "U2.pdf" and "U3.pdf", where U1 corresponds to the window without inducer, U2 to the window of early-times with inducer and U3 to the window of late-times with inducer. A comment on the ranking is found in page 5 of "Solution.pdf". To reproduce the arrays, run "P2b.m" followed by "P2c.m".
		
	Problem 3:
		a. Comments on the construction of the stochiometric matrix and the matrix itself are found in pages 6 and 7 of "Solution.pdf"
		b. The required figure is found inside the folder named "P3b" under the name "Protein.pdf". To reproduce the figure, import the folder into your favorite Julia compiler and run "include("main.jl")". This should both show and save the figure file.
