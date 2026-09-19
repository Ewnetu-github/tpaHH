# tpaHH
This is a source code for a hybrid hazard-based models using two-piece distributions. 
We proposed parametric, semi-parametric and non-parametric regression functions for the time-scale and relative hazard-scale change. 
Three parametric baseline distirbutions with  a log-link function are included in the code: TPA normal, TPA logistic and TPA Laplace distributions. 
# The code implements **hybrid hazard-based survival regression models** built
from **two-piece asymmetric (TPA)** baseline distributions. In these models,
a covariate can act on the survival time on two different scales at once:

* a **time-scale** (acceleration) effect, which stretches or compresses time
  itself (as in an Accelerated Failure Time, AFT, model), and
* a **hazard-scale** (relative-risk) effect, which multiplies the hazard up
  or down (as in a Proportional Hazards, PH, model).

Combining both gives the **Hybrid Hazard (HH)** model; letting only one of
the two operate recovers the classical **AFT**, **PH**, **Accelerated
Hazards (AH)**, or **General Hazards (GH)** models as special cases. The
covariate effect on either scale can be modelled:

