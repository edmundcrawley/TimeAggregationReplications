This is a readme that details how to run the files to replicate the main results for the BPS Paper. 

** NOTE **
The replication attempt with Time aggregation is currently on the replication_attempt branch on git!
Make sure you clone appropriately.

Step 1: 
Clone the project from git. 
git clone <hash>  

Step 2: 
Open shell.do. Change the path appropriately. Change the global variable "do" on line 9. Set this to point to ../Papers/CIFLS 
in the cloned directory. Run up to line 

Step 3: 
open bootstrap.do. Set the globals for table and column. Baseline results are table = 4, column = 1. 

Step 4: 
Open the paramfiles for the table/column you set in bootstrap.do. For the example of table 4 column 1, open ../Papers/paramfiles/T4C1.do.
This file sets all the parameters for the bootstrap. Read through the comments and set the parameters as you see fit. The final output files will contain the 
parameter names you set here.  

Step 5:
Run bootstrap.do 

Step 6: 
Find the output of the run you just did in ../Papers/CIFLS/T#C#_bs. There will be a csv starting with "estimates" with all the parameter values used in the file name as well. 
