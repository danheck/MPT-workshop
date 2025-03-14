
# get newest TreeBUGS version from GitHub (requires build tools)
# devtools::install_github("danheck/TreeBUGS")

library(TreeBUGS)

# Decision Inertia model
# (without independence assumption: separate parameters c1 and c2)

di <- "
optimal        opt_stay       d*c1
optimal        opt_stay       d*(1-c1)
optimal        opt_stay       (1-d)*c2
optimal        opt_stay       (1-d)*(1-c2)*g
optimal        opt_switch     (1-d)*(1-c2)*(1-g)
suboptimal     sub_switch     d*c1
suboptimal     sub_stay       d*(1-c1)
suboptimal     sub_switch     (1-d)*c2
suboptimal     sub_stay       (1-d)*(1-c2)*g
suboptimal     sub_switch     (1-d)*(1-c2)*(1-g)
"

# Extend the MPT model equations for the 2x2 design
# --> before copying to multiTree, row numbers and the column "EQN" must be removed
# --> the newest TreeBUGS version from Gtihub does this automatically

withinSubjectEQN(eqnfile = di,
                 labels = c("free80", "free60", "fixed80", "fixed60"),
                 constant = "g")
