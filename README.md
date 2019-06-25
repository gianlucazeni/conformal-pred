# conformal-pred

The code generalizes the R package [conformalInference](https://github.com/ryantibs/conformal), precisely when dealing with the normalization of nonconformity scores.

### Details

Normalized nonconformity scores are computed as:

![](https://latex.codecogs.com/gif.latex?\dpi{100}&space;R_i&space;=&space;\frac{\lvert&space;\hat{y}_i&space;-&space;y_i&space;\rvert}{\sigma_i},)

where the absolute error concerning the *i*th example is scaled using the expected accuracy $\sigma_i$ of the underlying model; see, e.g., [1,2].

$\sigma_i$ is an estimate of the difficulty of predicting the label $y_i$. 
There is a wide choice of estimates of the accuracy available in the literature. 
A common practice is to train another model to predict errors. 
More in details, once the underlying algorithm $f$ has been trained and the residual errors computed, a different model $g$ is fit using the object $x_1, \dots, x_n$ and the residuals. Then, $\sigma_i$ could be set equal to $g(x_i) + \beta$, where $\beta$ is a sensitivity parameter that has to be set in a way to give more or less weight to the impact of normalization. Increasing $\beta$ results in a less sensitive nonconformity measure with respect to changes of $g(x)$, and vice versa. As an alternative, $\sigma_i = \exp (\beta g(x_i))$. 

### References

[1]. Papadopoulos, H., Vovk, V., and Gammerman, A. (2011). Regression conformal prediction with nearest neighbours. *Journal of Artificial Intelligence Research*, 40:815-840        
[2]. Papadopoulos, H. and Haralambous, H. (2011). Reliable prediction intervals with regression neural networks. *Neural Networks*, 24(8):842-851   
