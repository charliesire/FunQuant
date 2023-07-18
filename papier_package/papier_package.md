---
title: "FunQuant"
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
date: 12 July 2023
bibliography: biblio.bib
---


# R Markdown

Quantization helps understand continuous distributions by providing a discrete approximation [@Pages]. Among the widely adopted methods for data quantization is the K-Means algorithm, which partitions the space into Voronoï cells, that can be seen as clusters, and constructs a discrete distribution based on their centroids and probabilistic masses. K-Means investigates the optimal centroids in a minimal expected distance sense [@Bock], but this approach poses significant challenges in scenarios where data evaluation is costly, and relates to a rare event that accumulates the majority of the probabilistic mass in a single cluster. In this context, a metamodel is required and adapted sampling methods are relevant to increase the precision of the computations on the rare clusters.

# Statement of need

FunQuant is an R package that has been specifically developed for carrying out quantization in the realm of rare events. While numerous cutting-edge packages facilitate straightforward implementation of the K-Means algorithm, they lack the incorporation of any probabilistic factors, treating all data points equally in terms of weighting. Conversely, FunQuant employs Importance Sampling estimators [@Paananen] instead of traditional Monte Carlo approach for calculating the centroids. To be more precise, when data $Y$ depends on probabilistic inputs $X$, the centroid of a cluster $C$ is estimated by the following formula: 

$$\frac{\frac{1}{n} \sum^{n}_{k=1} Y(\tilde{X}_{k})\mathbb{1}_{Y(\tilde{X}_{k})\in C}\frac{f_{X}(\tilde{X}^k)}{g(\tilde{X}_{k})}}{\frac{1}{n} \sum^{n}_{k=1} \mathbb{1}_{Y(\tilde{X}^k)\in C} \frac{f_{X}(\tilde{X}_k)}{g(\tilde{X}_{k})}}$$
where $f_{X}$ is the known density function of the inputs $X$, and $(\tilde{X}_k)^{n}_{k=1}$ i.i.d. random variables of density function $g$.
Importance Sampling is employed with the aim of reducing the variance of the estimators of the centroids when compared to classical Monte Carlo methods. FunQuant provides various approaches for implementing these estimators, depending on the sampling density denoted as $g$. The simplest method involves using the same function $g$ for each iteration and every cluster, which is straightforward to work with and still yields significant variance reduction. More advanced implementations enable the adaptation of the sampling density for each cluster at every iteration.

In addition, FunQuant is designed to mitigate the computational burden associated with the evaluation of costly data. While users have the flexibility to utilize their own metamodels to generate additional data, FunQuant offers several functions tailored specifically for a metamodel dedicated to spatial outputs such as maps. This metamodel relies on Functional Principal Component Analysis and Gaussian Processes, based on the work of @Perrin, adapted with the rlibkriging R package [@rlib]. FunQuant assists in the fine-tuning of its hyperparameters for a quantization task, with different performance metrics involved.

Additional theoretical information can be found in @sire. The paper provides a comprehensive exploration of the application of FunQuant to the quantization of flooding maps.

# Illustrative example

![Screenshot](fX.jpg)

# Acknowledgments

This research was conducted with the support of the consortium in
Applied Mathematics CIROQUO, gathering partners in technological and academia in the development of advanced methods for Computer Experiments.

# References


