This is an expiriment as to whetehr we can get useable marginal effects from ICE curves. 

See http://www.tandfonline.com/doi/abs/10.1080/10618600.2014.907095 or https://arxiv.org/abs/1309.6392 for reference on ICE curves.

The general hypothesis is that we can take the mean of the ICE derivatives across observations and the range of `x` to get an estimate of the marginal effect `x` has on `y`. In the linear and log linear cases, this should be proveably true. Then, we can extend this to less-transparent models such as random forest and neural networks.

If this expiriment seems (reasonably) successful, I will turn the underlying functions into an R package to publish on CRAN.
