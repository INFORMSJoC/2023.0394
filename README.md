# 2023.0394

This archive consists of the codes of three examples (Example 1, Example 2 and Example 3) in main text and one example (jump diffusion model) in Appendix. In each example, we upload the codes for different methods.

1. The codes for our proposed method is in the file "our method" while those for the EL procedure proposed by Lan et al. (2010) is in the file "Lan".

2. For our proposed method, to get a CI of CVaR, run file "get_results.m". The parameters including time horizon and confidence level are problem dependent.
Results are stored in the folder "results". For coverage probabilities and average CIs of multiple replications, run "get_results.m", and then run "gather_results.m".

3. For the EL procedure, the codes for parameter finding are in the folder "findPara". After finding the optimal parameter setting among 9 settings, update parameters in "get_results.m" in the "run" folder.
In the folder "findPara", the codes for the i-th (i=1,...,9) parameter setting are named "group i". Set parameters in "get_results.m" and run it to obtain the CI of CVaR. For coverage probabilities and average CIs of multiple replications, run "get_results.m", and then run "gather_results.m".

4. Specifically, for Example 1, the codes for the case of out-of-sample (without sample recycling) are in the folder "\Example 1\Out-of-sample".
For Example 1, the codes for different basis functions (standard polynomials up to fourth order and Laguerre polynomials up to second order) are stored in the folder "\Example 1\Standard polynomials up to fourth order" and "\Example 1\Laguerre polynomials up to second order".

5. The "results" folder contains subfolders named after sample sizes (e.g., "100000"). For example, if a subfolder is named "100000", it contains files storing the end-points of a CI with simulation sample size 100000 for one replication.