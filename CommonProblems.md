# Common Problems
Here are some common problems and their solutions.

## `Mu optimization returned NaN objective: restricting search space`
This error indicates that when optimizing the mu value for a guide, the log likelihood function returned NaN. There are three usual reasons for this:
1) Data that has really low coverage in the "input" and no pseudocount added (e.g. if A is 0% of the library, the likelihood of >1 reads in bin A is NaN)
2) Data where one "guide" has really high abundance. (e.g. if a "guide" makes up 89% of your library, a Z score of 1 might mean that you now expect OVER 100% to end up in bin F)
3) incorrect bin parameters (e.g. the bin covers 0% of the distribution [accidentally], but we saw reads there).

This warning does not necessarily mean there is a problem, however if this is followed by the subsequent error:
```
Error in optim(0, fn = function(mu) { :
  L-BFGS-B needs finite values of 'fn'
```
This could indicate that the mu resulted in NaN likelihoods for all tested vaues of mu, and parameters will have to be adjusted.
To solve (1), ensure that you add a pseudocount to the data.
To solve (2), usually MAUDE automatically restricting the search space will solve this.
To solve (3), double check your bin parameters. Make sure you remove any null bins from your bin data.frame (e.g. all bins must have some non-zero bin fraction (i.e. binEndQ-binStartQ > 0).

## Other problems
Please submit an Issue, or contact the authors for any other problems you encounter so that we can fix/clarify things as necessary.
