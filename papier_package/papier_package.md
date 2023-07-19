---
title: "FunQuant"
output: md_document
tags:
  - R
  - Quantization
  - Statistics
  - Metamodel
  - Importance Sampling
authors:
  - name: Charlie Sire
    orcid: 0000-0002-4432-7940
    affiliation: "1, 2, 3" # (Multiple affiliations must be quoted)
    corresponding: true
  - name: Yann Richet
    orcid: 0000-0002-5677-8458
    affiliation: 3
  - name: Rodolphe Le Riche
    orcid: 0000-0002-3518-2110
    affiliation: 1
  - name: Didier Rullière
    affiliation: 1
  - name: Jérémy Rohmer
    orcid: 0000-0001-9083-5965
    affiliation: 2
  - name: Lucie Pheulpin
    affiliation: 3


affiliations:
 - name: Mines Saint-Etienne, Univ. Clermont Auvergne, CNRS, UMR 6158 LIMOS
   index: 1
 - name: BRGM
   index: 2
 - name: IRSN
   index: 3
bibliography: biblio.bib
---


# Summary

Quantization helps understand continuous distributions by providing a discrete approximation [@Pages]. Among the widely adopted methods for data quantization is the K-Means algorithm, which partitions the space into Voronoï cells, that can be seen as clusters, and constructs a discrete distribution based on their centroids and probabilistic masses. K-Means investigates the optimal centroids in a minimal expected distance sense [@Bock], but this approach poses significant challenges in scenarios where data evaluation is costly, and relates to a rare event that accumulates the majority of the probabilistic mass in a single cluster. In this context, a metamodel is required and adapted sampling methods are relevant to increase the precision of the computations on the rare clusters.

# Statement of need

FunQuant is an R package that has been specifically developed for carrying out quantization in the realm of rare events. While numerous cutting-edge packages facilitate straightforward implementation of the K-Means algorithm, they lack the incorporation of any probabilistic factors, treating all data points equally in terms of weighting. Conversely, FunQuant employs Importance Sampling estimators [@Paananen] instead of traditional Monte Carlo approach for calculating the centroids. To be more precise, when data $Y$ depends on probabilistic inputs $X$, the centroid of a cluster $C$ is estimated by the following formula: 

$$\frac{\frac{1}{n} \sum^{n}_{k=1} Y(\tilde{X}_{k})\mathbb{1}_{Y(\tilde{X}_{k})\in C}\frac{f_{X}(\tilde{X}^k)}{g(\tilde{X}_{k})}}{\frac{1}{n} \sum^{n}_{k=1} \mathbb{1}_{Y(\tilde{X}^k)\in C} \frac{f_{X}(\tilde{X}_k)}{g(\tilde{X}_{k})}}$$
where $f_{X}$ is the known density function of the inputs $X$, and $(\tilde{X}_k)^{n}_{k=1}$ i.i.d. random variables of density function $g$.
Importance Sampling is employed with the aim of reducing the variance of the estimators of the centroids when compared to classical Monte Carlo methods. FunQuant provides various approaches for implementing these estimators, depending on the sampling density denoted as $g$. The simplest method involves using the same function $g$ for each iteration and every cluster, which is straightforward to work with and still yields significant variance reduction. More advanced implementations enable the adaptation of the sampling density for each cluster at every iteration.

In addition, FunQuant is designed to mitigate the computational burden associated with the evaluation of costly data. While users have the flexibility to utilize their own metamodels to generate additional data, FunQuant offers several functions tailored specifically for a metamodel dedicated to spatial outputs such as maps. This metamodel relies on Functional Principal Component Analysis and Gaussian Processes, based on the work of [@Perrin], adapted with the rlibkriging R package [@rlib]. FunQuant assists in the fine-tuning of its hyperparameters for a quantization task, with different performance metrics involved.

Additional theoretical information can be found in [@sire]. The paper provides a comprehensive exploration of the application of FunQuant to the quantization of flooding maps.
# Illustrative example

We consider a random variable $X = (R cos(\Theta), R sin(\Theta)) \in \mathbb{R}^{2}$ with $R$ and $\Theta$ 2 independant random variables defined by the following probability density functions:
$$\left\{
    \begin{array}{ll}
        f_{R}(r) = 0.99 \delta_{0} + 0.01\times 2(1-r)\mathbb{1}_{]0,1]}(r)\\
        f_{\Theta}(\theta) = \frac{1}{2\pi}\mathbb{1}_{[0,2\pi]}(\theta)
    \end{array}
\right.$$

The density function of $X$, denoted $f_{X}$, is represented in the following figure.
<p align="center">
  <img src="fX.jpg" width="350" title="hover text">

$99\%$ of the probability mass is concentrated at $(0,0)$.

We want to perform quantization on $X$ here (or $Y(X)$ with $Y$ the identity function).
 If the classical K-Means algorithm is performed with a budget of $1000$ points, it leads to the following outcome, with only a few sampled points different of $(0,0)$. Then, the centroids of the Voronoi cells that do not contain $(0,0)$ are computed with a very small number of points, leading to a very high variance.

 <p align="center">
  <img src="kmeans_quanti.jpg" width="350" title="K-Means centroids and Voronoi cells">

The FunQuant package to adapt the sampling and consider the probabilistic weights of each sample, with is the ratio $\frac{f_{X}}{g}$ where $g$ is density function associated to the adapted sampling. 

A possible function $g$ is $g = g_{R}\times g_{\Theta}$ 
$$\left\{
    \begin{array}{ll}
        g_{R}(r) = 0.05 \delta_{0} + 0.95\mathbb{1}_{]0,1]}(r)\\
        g_{\Theta}(\theta) = \frac{1}{2\pi}\mathbb{1}_{[0,2\pi]}(\theta)
    \end{array}
\right.$$


    sample_g = function(n){
      u = runif(n)
      vec_r = c(rep(0, sum(u>0.8)),runif(sum(u<=0.8)))
      vec_theta = runif(n, 0, pi/2)
      return(cbind(vec_r*cos(vec_theta), vec_r*sin(vec_theta)))
    }

    g = function(x){
        r = sqrt(x[1]^2+x[2]^2)
      if(r==0){return(0.2)}
      else if(r<=1){return(0.8/(2*pi))}
      else(return(0))
    }

    inputs = sample_g(1000)
    density_ratio = compute_density_ratio(f = fX, g = g, inputs = inputs)
    quantization = find_prototypes(nb_cells = 5,multistart = 3,data = t(inputs),density_ratio = density_ratio)


Figure 3 shows the sampled points, their associated probabilistic weights, and the obtained prototypes. It clearly appears that this sampling brings more information about each Voronoi cells. 

 <p align="center">
  <img src="is_quanti.jpg" width="350" title="K-Means centroids and Voronoi cells"> 


FunQuant allows to estimate the standard deviations of the estimators of the centroids for each Voronoi cell, highlighting the variance reduction obtained with the adapted sampling for the cells that do not contain $(0,0)$.


    large_sample = sample_fX(10^5)

    std_centroid_kmeans = std_centroid(data = t(large_sample), prototypes_list = list(protos_kmeans), density_ratio = rep(1, nrow(large_sample)), cells = 1:5, nv = 1000)

    std_centroid_kmeans #the cells are ordered by increasing "x" coordinate of their centroid

    ## [[1]]
    ## [[1]][[1]]
    ## [1] 7.786362e-05 7.366411e-05
    ## 
    ## [[1]][[2]]
    ## [1] 0.03561213 0.05429926
    ## 
    ## [[1]][[3]]
    ## [1] 0.08347199 0.10650128
    ## 
    ## [[1]][[4]]
    ## [1] 0.05506697 0.04682364
    ## 
    ## [[1]][[5]]
    ## [1] 0.09770415 0.12095001

    large_sample_g = sample_g(10^5)

    std_centroid_funquant = std_centroid(data = t(large_sample_g), prototypes_list = list(protos_funquant), density_ratio = rep(1, nrow(large_sample)), cells = 1:5, nv = 1000)

    std_centroid_funquant #the cells are ordered by increasing "x" coordinate of their centroid

    ## [[1]]
    ## [[1]][[1]]
    ## [1] 0.002392146 0.002499710
    ## 
    ## [[1]][[2]]
    ## [1] 0.005661370 0.005629974
    ## 
    ## [[1]][[3]]
    ## [1] 0.007004364 0.012310445
    ## 
    ## [[1]][[4]]
    ## [1] 0.014033466 0.006112396
    ## 
    ## [[1]][[5]]
    ## [1] 0.01082570 0.01121463

# Acknowledgments

This research was conducted with the support of the consortium in
Applied Mathematics CIROQUO, gathering partners in technological and
academia in the development of advanced methods for Computer
Experiments.

# References

Bock, Hans-Hermann. 2008. “Origins and Extensions of the k-Means
Algorithm in Cluster Analysis.” *Journal Électronique d’Histoire Des
Probabilités Et de La Statistique \[Electronic Only\]* 4 (January):
Article 14.

Havé, Pascal, Yann Richet, Yves Deville, Conrad Sanderson, Colin Fang,
Ciyou Zhu, Richard Byrd, Jorge Nocedal, and Jose Luis Morales. 2022.
“Rlibkriging: Kriging Models Using the ’libKriging’ Library.” *GitHub
Repository*. GitHub. <https://github.com/libKriging/rlibkriging>.

Paananen, Topi, Juho Piironen, Paul-Christian Bürkner, and Aki Vehtari.
2021. “Implicitly Adaptive Importance Sampling.” *Statistics and
Computing* 31 (2). <https://doi.org/10.1007/s11222-020-09982-2>.

Pagès, Gilles. 2014. “<span class="nocase">Introduction to optimal
vector quantization and its applications for numerics</span>.” LPMA.
<https://hal.archives-ouvertes.fr/hal-01034196>.

Perrin, T. V. E., O. Roustant, J. Rohmer, Olivier Alata, J. P. Naulin,
D. Idier, R. Pedreros, D. Moncoulon, and P. Tinard. 2021. “Functional
Principal Component Analysis for Global Sensitivity Analysis of Model
with Spatial Output.” *Reliability Engineering and System Safety* 211
(July): 107522. <https://doi.org/10.1016/j.ress.2021.107522>.

Sire, Charlie, Rodolphe Le Riche, Didier Rullière, Jérémy Rohmer, Lucie
Pheulpin, and Yann Richet. 2023. “Quantizing Rare Random Maps:
Application to Flooding Visualization.” *Journal of Computational and
Graphical Statistics* 0: 1–31.
<https://doi.org/10.1080/10618600.2023.2203764>.
