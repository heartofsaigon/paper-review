# paper-review

## Objective

We aim to review the Bayesian optimization (BO) and validate the approach by optimizing a pre-defined function using the proposed method.
The slide deck and report can be found in the following links.

1. [presentation slide](https://heartofsaigon.github.io/paper-review/02_presentation/BO-spresentation.html#/title-slide)
 
2. [Report (PDF)](https://heartofsaigon.github.io/paper-review/01_report/final-project-report.pdf)

We review two papers relating to Bayesian optimization. The first paper introduces Bayesian optimization as a derivative-free optimization algorithm, suitably used for objective functions with 
complex or unknown structure. The author well explained how BO works based on two primary steps: modelling observed data and suggesting potential data points for the next evaluation using the so-called acquisition functions. 
Pseudo-code of the algorithm is also provided. While the first paper focuses on the standard BO without constraints on the function domain, 
The second paper introduced a new acquisition function for constraint optimization.

We validate the algorithm by optimizing a defined function using both `R` built-in function `optimize` and BO algorithm and compare the results.
As the structure of the defined function is known, we can evaluate whether BO finds optimal point effectively compared to that obtained using the `R` built-in function.  

## Structure of repository

- `01_report`: the final report--PDF and quarto files.

- `02_presentation`: The presentation slide-- html and quarto files.

- `03_papers`: The two main papers reviewed.

- `04_code`: The only file `runBO_v2.R` includes codes for defining the object function and optimizing the function using both `optimize` function and BO algorithm.
The function optimized is defined as follows:

$$
f(x) = -6\times\mathrm \phi(x|4,1) - 7.5\times\mathrm \phi(x|7,1) + 1.3,
$$

where $\phi$ is the PDF of normal distribution. This code is only used for BO with EI acquisition function. We also assume the error-free model, i.e., 
All points in the domain can be evaluated without measurement error. Extensions, such as CBO or additional hyperparameters, can be achieved by modifying this code.  



