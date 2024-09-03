### ![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)

# Simulating Confidence Intervals for Conditional Value-at-Risk via Least-Squares Metamodels

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data that were used in the research reported on in the paper [Simulating Confidence Intervals for Conditional Value-at-Risk via Least-Squares Metamodels](https://doi.org/10.1287/ijoc.2023.0394) by Qidong Lai, Guangwu Liu, Bingfeng Zhang and Kun Zhang.

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

http://doi.org/10.1287/ijoc.2023.0394

http://doi.org/10.1287/ijoc.2023.0394.cd

Below is the BibTex for citing this snapshot of the repository.

```
@article{Lai2024CICVaR,
  author =        {Lai, Qidong and Liu, Guangwu and Zhang, Bingfeng and Zhang, Kun},
  publisher =     {INFORMS Journal on Computing},
  title =         {Github repository: Simulating Confidence Intervals for Conditional Value-at-Risk via Least-Squares Metamodels.},
  year =          {2024},
  doi =           {10.1287/ijoc.2023.0394.cd},
  note =          {Available for download at https://github.com/INFORMSJoC/2023.0394},
}  
```

## Description

The goal of this software is to demonstrate the effect of constructing confidence intervals for CVaR proposed in "Simulating Confidence Intervals for Conditional Value-at-Risk via Least-Squares Metamodels" by Qidong Lai, Guangwu Liu, Bingfeng Zhang and Kun Zhang.

Specifically, this archive consists of the codes of three examples (in folders [Example 1](https://github.com/KennethKZH/2023.0394/tree/main/Example1), [Example 2](https://github.com/KennethKZH/2023.0394/tree/main/Example2) and [Example 3](https://github.com/KennethKZH/2023.0394/tree/main/Example3), respectively) in main text and one example (in the folder [jump diffusion model](https://github.com/KennethKZH/2023.0394/tree/main/Jumpdiffusionmodel)) in Appendix. Each example folder contains the codes for different methods. For instance, in the folder [Example 1](https://github.com/KennethKZH/2023.0394/tree/main/Example1), the codes for our proposed method are in the subfolder [our method](https://github.com/KennethKZH/2023.0394/tree/main/Example1/Our_method) while those for the EL procedure proposed by Lan et al. (2010) are in the subfolder [Lan](https://github.com/KennethKZH/2023.0394/tree/main/Example1/Lan).

## Building

1. Clone this repository

2. Install Matlab (https://www.mathworks.com/products/matlab.html). Matlab is a commercial mathematical software.

3. The codes are entirely written in Matlab m-code. To install, unzip the downloaded files in any directory and add the subdirectory to the Matlab path.

## Results

All detailed results are available in the [Results](https://github.com/KennethKZH/2023.0394/tree/main/Results) folder.

## Replicating

Using Example 1 as a reference, the following steps outline the process to replicate the results: 

1. For our proposed method, to get a CI of CVaR, run `get_results.m`. Parameters can be adjusted to suit your purpose, such as the number of replications `replications`, the time horizon `T`, the confidence level of the confidence interval `z_beta`, etc. Results will be stored in the subfolder [results](https://github.com/KennethKZH/2023.0394/tree/main/Example1/Our_method/results). For coverage probabilities and average CIs of multiple replications, run `get_results.m`, and then run `gather_results.m`.

2. For the EL procedure, the codes for parameter finding are in the subfolder [findPara](https://github.com/KennethKZH/2023.0394/tree/main/Example1/Lan/findPara). After finding the optimal parameter setting among 9 settings, update parameters in `get_results.m` in the [run](https://github.com/KennethKZH/2023.0394/tree/main/Example1/Lan/run) subfolder. In the subfolder [findPara](https://github.com/KennethKZH/2023.0394/tree/main/Example1/Lan/findPara), the codes for the i-th (i=1,...,9) parameter setting are named "group i". Set parameters in `get_results.m` and run it to obtain the CI of CVaR. For coverage probabilities and average CIs of multiple replications, run "get_results.m", and then run `gather_results.m`.

3. Specifically, for Example 1, the codes for the case of out-of-sample (without sample recycling) are in the subfolder [Out-of-sample](https://github.com/KennethKZH/2023.0394/tree/main/Example1/Out-of-sample). For Example 1, the codes for different basis functions (standard polynomials up to fourth order and Laguerre polynomials up to second order) are stored in the subfolder [Standard polynomials up to fourth order](https://github.com/KennethKZH/2023.0394/tree/main/Example1/Standard_polynomials_up_to_fourth_order) and [Laguerre polynomials up to second order](https://github.com/KennethKZH/2023.0394/tree/main/Example1/Laguerre_polynomials_up_to_second_order), respectively.

4. After running `get_results.m`, the [results](https://github.com/KennethKZH/2023.0394/tree/main/Example1/Our_method/results) subfolder will contain subfolders named after sample sizes (e.g., "100000"). For example, if a subfolder is named "100000", it contains files storing the end-points of a CI with simulation sample size 100000 for one replication.

## Contact

If you have any questions about the codes, please contact Kun Zhang at kunzhang@ruc.edu.cn.



