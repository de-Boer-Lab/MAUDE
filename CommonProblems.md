# Common Problems
Here are some common problems and their solutions.

## `Mu optimization returned NaN objective: restricting search space`
This warning indicates that when optimizing the mu value for a guide, the log likelihood function returned NaN. 
This warning does not necessarily mean there is a problem.

There are three usual reasons for this:
### 1) Data that has really low coverage in the "input" and no pseudocount added (e.g. if A is 0% of the library, the likelihood of >1 reads in bin A is NaN)
To resolve the warning, please ensure that you are including a pseudocount to the data.

### 2) Data where one "guide" has really high abundance. (e.g. if a "guide" makes up 89% of your library, a Z score of 1 might mean that you now expect OVER 100% to end up in bin F)
Usually, MAUDE automatically restricting the search space will resolve this. The warning will stay, but you don't have to worry about it.

### 3) The mu search limits are too wide. (an extreme mu is so unlikely that the log likelihood is infinite, leading to this error)
If you only get this message a few times, there's no need to adjust anything; MAUDE automatically adjusts the search space. If many warnings are issued, it is recommended that you adjust the `limits` parameter. It defaults to c(-4,4). Adjusting it closer to 0 will eliminate the warning message. 

## Other problems
Please submit an Issue, or contact the authors for any other problems you encounter so that we can fix/clarify things as necessary.
