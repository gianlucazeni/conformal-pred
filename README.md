# conformal-pred
Conformal prediction R functions.

The code generalizes the package conformalInference, precisely when dealing with the normalization of nonconformity scores.

Normalized nonconformity scores are computed as:

    R_i = \frac{\lvert \hat{y}_i - y_i \rvert}{\sigma_i},

where the absolute error concerning the $i$th example is scaled using the expected accuracy $\sigma_i$ of the underlying model; see, e.g., \cite{papadopoulos2011reliable}, and \cite{papadopoulos2011regression}.

$\sigma_i$ is an estimate of the difficulty of predicting the label $y_i$. 
There is a wide choice of estimates of the accuracy available in the literature. 
A common practice is to train another model to predict errors. 
More in details, once the underlying algorithm $f$ has been trained and the residual errors computed, a different model $g$ is fit using the object $x_1, \dots, x_n$ and the residuals. Then, $\sigma_i$ could be set equal to $g(x_i) + \beta$, where $\beta$ is a sensitivity parameter that has to be set in a way to give more or less weight to the impact of normalization. Increasing $\beta$ results in a less sensitive nonconformity measure with respect to changes of $g(x)$, and vice versa. As an alternative, $\sigma_i = \exp\,(\beta g(x_i))$. 

