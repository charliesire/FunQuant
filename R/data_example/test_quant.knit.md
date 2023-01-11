
<!-- rnb-text-begin -->

---
title: "R Notebook"
output: html_notebook
---

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*. 

<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuc291cmNlLmFsbChwYXRoID1cIkM6L1VzZXJzL3NpcmUtY2hhL0RvY3VtZW50cy9GdW5RdWFudC9SXCIpXG5gYGAifQ== -->

```r
source.all(path ="C:/Users/sire-cha/Documents/FunQuant/R")
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiKioqIHNvdXJjZSAxZF90b19tYXRyaXguUiBcbioqKiBzb3VyY2UgY29tcHV0ZV9jZW50cm9pZHMuUiBcbioqKiBzb3VyY2UgY3YuUiBcbioqKiBzb3VyY2UgZGVuc2l0eV9yYXRpby5SIFxuKioqIHNvdXJjZSBkaXN0YW5jZTJnYW1tYS5SIFxuKioqIHNvdXJjZSBlc3RpbV9kZW5vbV9jZW50cm9pZC5SIFxuKioqIHNvdXJjZSBlc3RpbV9udW1fY2VudHJvaWQuUiBcbioqKiBzb3VyY2UgZ2V0X251bS5SIFxuKioqIHNvdXJjZSBnZXRfcHJvYmFzLlIgXG4qKiogc291cmNlIGhlbGxvLlIgXG4qKiogc291cmNlIHNvcnRfZ2FtbWEuUiBcbiJ9 -->

```
*** source 1d_to_matrix.R 
*** source compute_centroids.R 
*** source cv.R 
*** source density_ratio.R 
*** source distance2gamma.R 
*** source estim_denom_centroid.R 
*** source estim_num_centroid.R 
*** source get_num.R 
*** source get_probas.R 
*** source hello.R 
*** source sort_gamma.R 
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubGlicmFyeShHR2FsbHkpXG5saWJyYXJ5KHNwKVxubGlicmFyeShnZ3Bsb3QyKVxubGlicmFyeShndG9vbHMpXG5saWJyYXJ5KGRwbHlyKVxubGlicmFyeShQT1QpXG5saWJyYXJ5KHByYWNtYSlcbmxpYnJhcnkoRGljZURlc2lnbilcbmxpYnJhcnkocmFuZHRvb2xib3gpXG5saWJyYXJ5KGFiaW5kKVxubGlicmFyeShtaWNlYWRkcylcbmxpYnJhcnkod2F2ZXNsaW0pXG5saWJyYXJ5KGZvcmVhY2gpXG5saWJyYXJ5KERpY2VLcmlnaW5nKVxubGlicmFyeShvcnRob2dvbmFsc3BsaW5lYmFzaXMpXG5saWJyYXJ5KHZpcmlkaXMpXG5saWJyYXJ5KGNvbnRvdXJlUilcbmxpYnJhcnkoTUFTUylcbmxpYnJhcnkocmVzaGFwZTIpXG5saWJyYXJ5KFJDb2xvckJyZXdlcilcbmxpYnJhcnkocmFuZG9tRm9yZXN0KVxubGlicmFyeShwbmcpXG5saWJyYXJ5KENsaW1Qcm9qRGlhZ3MpXG5saWJyYXJ5KGRpc21vKVxuYGBgIn0= -->

```r
library(GGally)
library(sp)
library(ggplot2)
library(gtools)
library(dplyr)
library(POT)
library(pracma)
library(DiceDesign)
library(randtoolbox)
library(abind)
library(miceadds)
library(waveslim)
library(foreach)
library(DiceKriging)
library(orthogonalsplinebasis)
library(viridis)
library(contoureR)
library(MASS)
library(reshape2)
library(RColorBrewer)
library(randomForest)
library(png)
library(ClimProjDiags)
library(dismo)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZ2FtbWEgPSBnYW1tYVxuZGlzdGFuY2VfZnVuYyA9IGRpc3RfTDJcbm5jb2VmZl92ZWMgPSBjKDUwMCwxNTAwLDI1MDAsMzAwMClcbm5wY192ZWM9Mjo1XG5mb3JtdWxhID0gfjFcbmRlc2lnbiA9IGRvZV9mdWxsW3F1aV9uYiwxOjVdXG5yZXNwb25zZSA9IGZwXG5cbmBgYCJ9 -->

```r
gamma = gamma
distance_func = dist_L2
ncoeff_vec = c(500,1500,2500,3000)
npc_vec=2:5
formula = ~1
design = doe_full[qui_nb,1:5]
response = fp

```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiRXJyb3I6IG9iamVjdCAnZnAnIG5vdCBmb3VuZFxuIn0= -->

```
Error: object 'fp' not found
```



<!-- rnb-output-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuc291cmNlLmFsbChwYXRoID1cIkM6L1VzZXJzL3NpcmUtY2hhL0RvY3VtZW50cy9GdW5RdWFudC9SXCIpXG5cbmBgYCJ9 -->

```r
source.all(path ="C:/Users/sire-cha/Documents/FunQuant/R")

```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiKioqIHNvdXJjZSAxZF90b19tYXRyaXguUiBcbioqKiBzb3VyY2UgY29tcHV0ZV9jZW50cm9pZHMuUiBcbioqKiBzb3VyY2UgY3YuUiBcbioqKiBzb3VyY2UgZGVuc2l0eV9yYXRpby5SIFxuKioqIHNvdXJjZSBkaXN0YW5jZTJnYW1tYS5SIFxuKioqIHNvdXJjZSBlc3RpbV9kZW5vbV9jZW50cm9pZC5SIFxuKioqIHNvdXJjZSBlc3RpbV9udW1fY2VudHJvaWQuUiBcbioqKiBzb3VyY2UgZ2V0X251bS5SIFxuKioqIHNvdXJjZSBnZXRfcHJvYmFzLlIgXG4qKiogc291cmNlIGhlbGxvLlIgXG4qKiogc291cmNlIHNvcnRfZ2FtbWEuUiBcbiJ9 -->

```
*** source 1d_to_matrix.R 
*** source compute_centroids.R 
*** source cv.R 
*** source density_ratio.R 
*** source distance2gamma.R 
*** source estim_denom_centroid.R 
*** source estim_num_centroid.R 
*** source get_num.R 
*** source get_probas.R 
*** source hello.R 
*** source sort_gamma.R 
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxucHAgPSBwcm9iYXNfbG9vKG91dHB1dHMgPSBvdXRwdXRzWywsIHF1aV9uYl0sIGRlbnNpdHlfcmF0aW8gPSBkZW5zaXR5X3JhdGlvW3F1aV9uYl0sZ2FtbWEgPSBnYW1tYV90dCwgZGlzdGFuY2VfZnVuYyA9IGRpc3RfTDIsIG5jb2VmZl92ZWMgPSBjKDUwMCwxNTAwLDI1MDAsMzAwMCksIG5wY192ZWM9Mjo1LCAgZm9ybXVsYSA9IH4xLCBkZXNpZ24gPSBkb2VfZnVsbFtxdWlfbmIsMTo1XSwgY292dHlwZSA9IFwibWF0ZXJuNV8yXCIsIGNvZWYudHJlbmQgPSBOVUxMLCBjb2VmLnZhciA9IE5VTEwsIGNvZWYuY292ID0gTlVMTCwgY29udHJvbCA9IGxpc3QodHJhY2UgPSAwKSwgbnVnZ2V0ID0gMC4xLG11bHRpc3RhcnQ9NSwgbG93ZXIgPSByZXAoMC4xLDUpLCBtb2RlbF90dW5pbmcgPSBtb2RkLCBiaWFzID0gYyhjc3RlX2JpYXMsIHJlcCgwLDQpKSlcbmBgYCJ9 -->

```r
pp = probas_loo(outputs = outputs[,, qui_nb], density_ratio = density_ratio[qui_nb],gamma = gamma_tt, distance_func = dist_L2, ncoeff_vec = c(500,1500,2500,3000), npc_vec=2:5,  formula = ~1, design = doe_full[qui_nb,1:5], covtype = "matern5_2", coef.trend = NULL, coef.var = NULL, coef.cov = NULL, control = list(trace = 0), nugget = 0.1,multistart=5, lower = rep(0.1,5), model_tuning = modd, bias = c(cste_bias, rep(0,4)))
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiWzFdIDFcblsxXSAyXG5bMV0gM1xuWzFdIDRcblsxXSA1XG5bMV0gNlxuWzFdIDdcblsxXSA4XG5bMV0gOVxuWzFdIDEwXG5bMV0gMTFcblsxXSAxMlxuWzFdIDEzXG5bMV0gMTRcblsxXSAxNVxuWzFdIDE2XG4ifQ== -->

```
[1] 1
[1] 2
[1] 3
[1] 4
[1] 5
[1] 6
[1] 7
[1] 8
[1] 9
[1] 10
[1] 11
[1] 12
[1] 13
[1] 14
[1] 15
[1] 16
```



<!-- rnb-output-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxub3V0cHV0cyA9IG91dHB1dHNbLCwgcXVpX25iXVxubmJfZm9sZHMgPSA3XG5kZW5zaXR5X3JhdGlvID0gZGVuc2l0eV9yYXRpb1txdWlfbmJdXG5nYW1tYSA9IGdhbW1hX3R0XG5kaXN0YW5jZV9mdW5jID0gZGlzdF9MMlxubmNvZWZmX3ZlYyA9IGMoNTAwLDE1MDAsMjUwMCwzMDAwKVxubnBjX3ZlYz0yOjVcbmZvcm11bGEgPSB+MVxuZGVzaWduID0gZG9lX2Z1bGxbcXVpX25iLDE6NV1cbmNvdnR5cGUgPSBcIm1hdGVybjVfMlwiXG5jb2VmLnRyZW5kID0gTlVMTFxuY29lZi52YXIgPSBOVUxMXG5jb2VmLmNvdiA9IE5VTExcbmNvbnRyb2wgPSBsaXN0KHRyYWNlID0gMClcbm51Z2dldCA9IDAuMVxubXVsdGlzdGFydD01XG5sb3dlciA9IHJlcCgwLjEsNSlcbm1vZGVsX3R1bmluZyA9IG1vZGRcbmJpYXMgPSBjKGNzdGVfYmlhcywgcmVwKDAsNCkpXG5gYGAifQ== -->

```r
outputs = outputs[,, qui_nb]
nb_folds = 7
density_ratio = density_ratio[qui_nb]
gamma = gamma_tt
distance_func = dist_L2
ncoeff_vec = c(500,1500,2500,3000)
npc_vec=2:5
formula = ~1
design = doe_full[qui_nb,1:5]
covtype = "matern5_2"
coef.trend = NULL
coef.var = NULL
coef.cov = NULL
control = list(trace = 0)
nugget = 0.1
multistart=5
lower = rep(0.1,5)
model_tuning = modd
bias = c(cste_bias, rep(0,4))
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxucHBfZm9sZCA9IHByb2Jhc19rX2ZvbGQob3V0cHV0cyA9IG91dHB1dHMsIG5iX2ZvbGRzID0gNywgZGVuc2l0eV9yYXRpbyA9IGRlbnNpdHlfcmF0aW8sZ2FtbWEgPSBnYW1tYV90dCwgZGlzdGFuY2VfZnVuYyA9IGRpc3RfTDIsIG5jb2VmZl92ZWMgPSBjKDUwMCwxNTAwLDI1MDAsMzAwMCksIG5wY192ZWM9Mjo1LCAgZm9ybXVsYSA9IH4xLCBkZXNpZ24gPSBkb2VfZnVsbFtxdWlfbmIsMTo1XSwgY292dHlwZSA9IFwibWF0ZXJuNV8yXCIsIGNvZWYudHJlbmQgPSBOVUxMLCBjb2VmLnZhciA9IE5VTEwsIGNvZWYuY292ID0gTlVMTCwgY29udHJvbCA9IGxpc3QodHJhY2UgPSAwKSwgbnVnZ2V0ID0gMC4xLG11bHRpc3RhcnQ9NSwgbG93ZXIgPSByZXAoMC4xLDUpLCBtb2RlbF90dW5pbmcgPSBtb2RkLCBiaWFzID0gYyhjc3RlX2JpYXMsIHJlcCgwLDQpKSlcbmBgYCJ9 -->

```r
pp_fold = probas_k_fold(outputs = outputs, nb_folds = 7, density_ratio = density_ratio,gamma = gamma_tt, distance_func = dist_L2, ncoeff_vec = c(500,1500,2500,3000), npc_vec=2:5,  formula = ~1, design = doe_full[qui_nb,1:5], covtype = "matern5_2", coef.trend = NULL, coef.var = NULL, coef.cov = NULL, control = list(trace = 0), nugget = 0.1,multistart=5, lower = rep(0.1,5), model_tuning = modd, bias = c(cste_bias, rep(0,4)))
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiWzFdIDFcblsxXSAyXG5bMV0gM1xuWzFdIDRcblsxXSAxXG5bMV0gMlxuWzFdIDNcblsxXSA0XG5bMV0gMVxuWzFdIDJcblsxXSAzXG5bMV0gNFxuWzFdIDFcblsxXSAyXG5bMV0gM1xuWzFdIDRcblsxXSAxXG5bMV0gMlxuWzFdIDNcblsxXSA0XG5bMV0gMVxuWzFdIDJcblsxXSAzXG5bMV0gNFxuWzFdIDFcblsxXSAyXG5bMV0gM1xuWzFdIDRcbiJ9 -->

```
[1] 1
[1] 2
[1] 3
[1] 4
[1] 1
[1] 2
[1] 3
[1] 4
[1] 1
[1] 2
[1] 3
[1] 4
[1] 1
[1] 2
[1] 3
[1] 4
[1] 1
[1] 2
[1] 3
[1] 4
[1] 1
[1] 2
[1] 3
[1] 4
[1] 1
[1] 2
[1] 3
[1] 4
```



<!-- rnb-output-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Ctrl+Alt+I*.

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Ctrl+Shift+K* to preview the HTML file).

The preview shows you a rendered HTML copy of the contents of the editor. Consequently, unlike *Knit*, *Preview* does not run any R code chunks. Instead, the output of the chunk when it was last run in the editor is displayed.

<!-- rnb-text-end -->

