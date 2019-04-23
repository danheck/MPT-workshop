
################################################################################
#
# This R script serves as instruction and test to install TreeBUGS.
#
# For questions, please contact Daniel Heck (heck@uni-mannheim.de)
#
################################################################################


################## Preliminary steps

# 1a. Install R:
#     => http://cran.r-project.org/

# 1b. Install RStudio:
#     => https://www.rstudio.com/products/rstudio/download/

# 1c. Install JAGS
#     => https://sourceforge.net/projects/mcmc-jags/files/



################## Install TreeBUGS

# 2a. Open this R script within RStudio
#     => e.g., via "File -> Open File"

# 2b. Install TreeBUGS
#     => Select the following line and click on the "Run" button (upper right)
install.packages("TreeBUGS")




################## Check TreeBUGS installation

# 3a. Check whether TreeBUGS can be loaded
#     => Select the following line and click on the "Run" button (upper right)
library("TreeBUGS")

# 3b. If you get no errors, everything is fine :-)

# 3c. Seldomly, Windows users get the following error:
#     "Failed to locate any version of JAGS version 4"
#           => this means that R does not find the install folder of JAGS
#           => to solve this, follow these steps (Windows 10):
#    1. Click on the windows button
#    2. Search for "edit the system environment variables"
#       (German Windows: "Systemumgebungsvariablen bearbeiten")
#    3. Click on "Environment Variables" (German: "Umgebungsvariablen")
#    4. In the lower field, search for "Path" (German: "Pfad")
#    5. Click on "Edit" and add the following entry:
#       C:\Program Files\JAGS\JAGS-4.3.0\x64\bin
#       (check that the JAGS version matches the one you installed!)


# 3d. Run the following lines (again using the "Run" button)
#     => Everything works fine if you get something like
#        "MCMC sampling started at  2018-08-......"

library("TreeBUGS")
EQNfile <- system.file("MPTmodels/2htm.eqn", package="TreeBUGS")
genDat <- genTraitMPT(N = 20, numItems = c(Target=150, Lure=150),
                      eqnfile = EQNfile,
                      mean = c(Do=.7, Dn=.7, g=.5),
                      sigma =   c(Do=.3, Dn=.3, g=.15))
fit <- traitMPT(EQNfile, genDat$data, n.thin=5, n.iter = 2000, n.burnin = 200,
                restrictions=list("Do=Dn"))
summary(fit)

