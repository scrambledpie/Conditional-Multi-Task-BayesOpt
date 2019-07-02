# Continuous Multi-Task Bayesian Optimization with Correlation

Toy code for reproducing the simpler experiments from the paper:

[Continuous Multi-Task Bayesian Optimization with Correlation](https://www.sciencedirect.com/science/article/abs/pii/S0377221718302261)

The code is in R and is all self-contained in a single file. It requires R packages which are installed by running these command in the R console:

```
install.packages("FastGP")
install.packages("MASS")
install.packages("Rcpp")
```


The algortihm takes a function of two inputs and finds all the optima of one input conditional on the other input. In other words, optimize every vertical slice in the picture below, each slice is it's own Bayesian optiimzation problem, we just assume there is more than one, or even a continuum of problems.


```
@article{pearce2018continuous,
  title={Continuous multi-task Bayesian Optimisation with correlation},
  author={Pearce, Michael and Branke, Juergen},
  journal={European Journal of Operational Research},
  volume={270},
  number={3},
  pages={1074--1085},
  year={2018},
  publisher={Elsevier}
}
```

![alt text](https://github.com/scrambledpie/Conditional-Multi-Task-BayesOpt/blob/master/MTKGREVI_sparse.gif)

![alt text](https://raw.githubusercontent.com/scrambledpie/Conditional-Multi-Task-BayesOpt/master/Screen%20Shot%202019-05-21%20at%205.23.01%20PM.png)
