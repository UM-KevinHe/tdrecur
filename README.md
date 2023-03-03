# tdrecur
## Overview
tdrecur is a package dealing with the recurrent events model with time-dependent covariates. Users are allowed to choose time-independent version if their data does not contain time-dependent covariates, which would be faster. Parallel computation based on OpenMP is also provided for the time-dependent version when the user's computation environment supports. 

The function will estimate block effects, beta, and baseline function.

(we might want to incorporate function of t in the covariates, or time since last event, these should not be step function as the current ones like covid. maybe put these function of computing covariates values outside the program and call them when needed)
(try to make the package more general)

## Installation
```{r setup}
devtools::install_github("XuemeiDing/tdrecur")
```

## Usage
The data is saved as a large matrix, where each row specifies information for one participant. The result will be the same if a participant's record is broken into several rows due to changing block, etc. The following columns are required, but the order is not important: time of entering the study, time of leaving the study, time of having the events, time-independent covariates, block the participant belongs to. If the time-independent covariates are present, columns specifying the time points of changing values as well as the start value and values after each change of the time-independent covariates are also required. In order to make the input data a matrix, NAs need to be added to the data so that each patient will have the same number of columns for the time-independent covariates. For further details of function arguments, please refer to the help document of the function.

Besides the time-independent covariates that require users to specify their values and changing points, functions of time (since 0) or time since last event and their interactions with other covariates are also allowed.
```{r}
library(tdrecur)
load(example_data)
tdrecur(example_data, ..., tol = 1.0e-6, max.iter = 10000, time-dependent = T, omp = T)
```
