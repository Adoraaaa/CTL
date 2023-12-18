This paper was presented at the 22nd China Automation Conference：[paper](https://ieeexplore.ieee.org/abstract/document/10055437)
## BACKGROUND
* Sparse representation models includes two categories: Dictionary learning (DL) and Transform learning (TL).
In traditional DL and TL, the complete signal is separated
into multiple overlapping signal blocks, which are sparsely
processed respectively, and then the obtained sparse features
are combined into a complete sparse representation.
*  However, this method destroys the internal structure of the whole signal
and is not effective in dealing with large data sets.
* the convolutional transform learning (CTL) is
proposed by the researchers to overcome the problem of destroying the complete structure of signals in the TL.
* At present, among the convolution kernel updating methods
of CTL, an orthonormality constraint is used to avoid generating trivial solutions. However, for the convolution kernel
learning with an orthogonality constraint, the traditional methods
analyzes the convolution kernel in Euclidian space with the
orthonormal constraint. Such a method may make it difficult
to capture the inherent structure of the data and requires
additional operations to meet orthonormality constraint.
## CONTRIBUTIONS
* To capture the inherent structure of the data, we employ a manifold optimization method to update the convolution kernel,
which not only make learned convolution kernels capture the inherent structure of the signal but also simplify convolution kernels update stage.
* In order to extract the deep semantic features of the
complex signals, we propose a multi-layer CTL model, which
is the expansion of single-layer CTL model, to extract higher
quality and more discriminative sparse features from input
signals, developing a multi-layer CTL algorithm based on
manifold optimization method.
## EXPERIMENTAL RESULT
We used public YALE face dataset and divided it into a training set and a test set.
* the manifold optimization method can be used for improving the sparse feature extraction effect of the CTL model
 compared to the proximal gradient descent method
* multi-layer CTL-MF-log can extract deeper semantic information compared to the single-layer CTL-MF-log.
