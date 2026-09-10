# [MatlabStan](https://github.com/brian-lau/MatlabStan)
<a href="http://mc-stan.org">
<img src="https://raw.githubusercontent.com/stan-dev/logos/master/logo.png?raw=true" width=100 alt="Stan Logo"/>
</a>

A Matlab interface to [Stan](http://mc-stan.org), a package for Bayesian inference.

For more information on Stan and its modeling language, see the Stan User's Guide and Reference Manual at [http://mc-stan.org/](http://mc-stan.org/documentation/).

## Installation
Details can be found in the [Getting started](https://github.com/brian-lau/MatlabStan/wiki/Getting-Started) page of the wiki.

## Example
The following is the classic 'eight schools' example from Section 5.5 of [Gelman et al (2003)](http://stat.columbia.edu/~gelman/book/). The output can be compared to that obtained using the [Rstan](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started#example-1-eight-schools) and [Pystan](https://github.com/stan-dev/pystan/blob/develop/README.rst#example) interfaces.
```
schools_code = {
   'data {'
   '    int<lower=0> J; // number of schools '
   '    real y[J]; // estimated treatment effects'
   '    real<lower=0> sigma[J]; // s.e. of effect estimates '
   '}'
   'parameters {'
   '    real mu; '
   '    real<lower=0> tau;'
   '    real eta[J];'
   '}'
   'transformed parameters {'
   '    real theta[J];'
   '    for (j in 1:J)'
   '    theta[j] <- mu + tau * eta[j];'
   '}'
   'model {'
   '    eta ~ normal(0, 1);'
   '    y ~ normal(theta, sigma);'
   '}'
};
  
schools_dat = struct('J',8,...
                     'y',[28 8 -3 7 -1 1 18 12],...
                     'sigma',[15 10 16 11 9 11 10 18]);

fit = stan('model_code',schools_code,'data',schools_dat);

print(fit);

eta = fit.extract('permuted',true).eta;
mean(eta)

```
## Within-chain threading
MatlabStan supports CmdStan within-chain parallelization using `threads_per_chain`.

- `threads_per_chain = 1` (default): threading off (current behavior).
- `threads_per_chain = N` (`N >= 1`): explicit threads per chain.
- `threads_per_chain = 'auto'`: uses `floor(available_cores/chains)` with minimum `1`.

Example:
```
fit = stan('model_code',model_code,'data',data,...
           'chains',4,'threads_per_chain','auto');
```

Behavior by CmdStan version:

- CmdStan `>= 2.28`: passes `num_threads=<threads_per_chain>` on the command line.
- CmdStan `< 2.28`: sets environment variable `STAN_NUM_THREADS=<threads_per_chain>` at runtime.

When `threads_per_chain > 1`, model compilation is invoked with `STAN_THREADS=true`.
To force a clean rebuild (`make clean-all`/`clean` and `make build STAN_THREADS=true`) before threaded model compilation, set:
```
fit = stan(...,'threads_per_chain',4,'rebuild_for_threads',true);
```

Note: threading only speeds up models that actually use threaded Stan constructs (for example `reduce_sum`).

## Need help?
You may be able to find a solution in the [wiki](https://github.com/brian-lau/MatlabStan/wiki/). Otherwise, open an [issue](https://github.com/brian-lau/MatlabStan/issues).

Contributions
--------------------------------
MatlabStan Copyright (c) 2017 Brian Lau [brian.lau@upmc.fr](mailto:brian.lau@upmc.fr), [BSD-3](https://github.com/brian-lau/MatlabStan/blob/master/LICENSE.txt)

[PSIS package](https://github.com/avehtari/MatlabPSIS) Copyright (c) 2015 Aki Vehtari, [GPL-3](http://www.gnu.org/licenses/gpl-3.0.en.html)

Please feel free to [fork](https://github.com/brian-lau/MatlabStan/fork) and contribute!
