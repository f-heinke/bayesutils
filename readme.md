# bayesutils
This R package offers a lightweight set of tools for working with posterior draws from Markov Chain Monte Carlo (MCMC) samplers and associated Stan code. It provides 

- flexible, interface-agnostic utilities for 
- extracting, 
- transforming, 
- comparing, 
- and visualizing parameter draws, 
- as well as performing prior and posterior predictive simulations.

The package is designed to complement existing tools such as `CmdStanR` and `rstanarm`, as well as the widely-used `rethinking`package, addressing common workflow gaps in Bayesian modeling. While still a work in progress, it aims to support both teaching and research use cases by streamlining analysis across various Stan-based interfaces.

## Installation

This package is not on CRAN. Instead you can install it directly from GitHub. Yet some additional packages are required (some you might have already installed). 

```
install.packages(
    c("Cairo",
    "cmdstanr",
    "coda",
    "jsonlite",
    "posterior",
    "rstan",
    "stringr",
    "ggplot2")	
)
devtools::install_github("f-heinke/bayesutils")
```

## Some Functionality demonstrated

To showcase some of the package features, let's consider the Palmer penguins dataset. This table holds data for 333 penguins recorded near Palmer Station, Antarctica. We shall explore how ``bayesutils`` can support building a statistical model for infering penguin log-body mass, and how MCMC draws can be analyzed.   

We start off by retrieving the data (f.e. by means of the ``palmerpenguins`` package). 


```
library(palmerpenguins)
library(dplyr)
d <- penguins %>% na.omit %>% mutate(
  log.bm.kg = log( body_mass_g / 1000 ),
  bl.cm = bill_length_mm / 10
)
```

### Structured Composing of Stan Code

Next we want to build a simple linear regression model in Stan for log-body mass prediction with bill length, sex and species as covariates. Now imagine you are in the spot where you have to teach the model step-by-step, or manage multiple modelling ideas. Here the `bayesutils` idea provides a simple workflow, where single Stan code components are composed into one functioning Stan program. 

In this example, we define the following Stan code for the model block only:

```
cdat <- "
  int N;
  // predictors
  array[N] int sex; // 1: fem, 2: male
  array[N] int spec; // 1: Adel, 2: Chin, 3: Gent
  vector[N] l; // bill length (cm)

  // outcome
  vector[N] m; // log-body mass (kg)
"
```
Imagine this beaing a code chunk in a notebook, giving detailed explanations to students on what is going on in the code.

In this fashion we continue by writing the model block:


```
cmdl <- "
  for(k in 1:3) B[k,] ~ normal(0, 1);
  a ~ normal(0, 1);
  sigma ~ exponential(1);
  for(i in 1:N) m[i] ~ normal( B[ spec[i], sex[i]] + a * l[i], sigma  );
"
```

(Imagine here yourself writing notes or talking to students.)

Finally, the parameters block:

```
cpar <- "
  matrix[ 3, 2 ] B;
  real a;
  real<lower = 0> sigma;
"
```

Then we compose the full Stan program from these blocks:

```
code <- stancode( cdat, cpar, cmdl )
```

To get the model to run, it is best to write it to a file and invoke the Stan C++ toolchain to do its magic.

```
sc.to.stan.file("peng.stan", code = code)

library("cmdstanr")
m <- cmdstan_model("peng.stan")
```

Note, that `sc.to.stan.file` calls ``rstans`` pre-compiler for syntax checking. If there are any syntax errors in your code, you'd be informed.

Next, we aggregate the needed data and start the sampler:

```
dat <- list(
  N = nrow(d),
  sex = as.numeric(d$sex),
  spec = as.numeric(d$species),
  l = d$bl.cm,
  m = d$log.bm.kg
)

fit <- m$sample(data = dat, parallel_chains = 2)
```

### Monitoring the sampling process


The HMC sampler finished quickly in the previous example. But now imagine that, given some other model to fit, sampling is slow -- very slow. This could be due to 1. simply a lot of data to process for each sampling iteration 2. non-optimal code, or 3. a mix of reasons including model misspecification and non-identifiability. The latter are most concerning, as they result in pathological chain dynamics and poor sampling. Anyway as a user you can't tell which one of the reasons might be the case until sampling is finished and you go through your post-sampling diagnostics checklist; unless you inspect the draws log files which are continuously saved while sampling. And that's exactly what the following wrapper function does, by providing traceplot visualizations of the ongoing sampling process:


``` 
fit <- sample.with.monitor( m, dat, 
	parallel_chains = 2, 
	chains = 2)
``` 



Notice that a folder and an R script is produced in your working directory. The folder holds draw files while the R script takes care of rendering continously updated traceplots that are then pushed to your webrowser--while the sampler is still going. 

Since sampling happens so quickly in this example, the monitor might not start-up and get to a refresh cycle before the sampler is already finished. Thus try the following code with increased iteration numbers, the `save_warmup` option enabled and `adapt_delta` set to an overly causious value, just for putting the brakes on artificially. This should give the monitor some time to load and to demonstrate some refresh cycles.

``` 
fit <- sample.with.monitor( m, dat, 
	parallel_chains = 2, 
	chains = 2, 
	iter_sampling = 5000,
	iter_warmup = 5000,
	save_warmup = T,
	adapt_delta = 0.999)
``` 


<div align="center">
  <img src="./demo/monitor.png" alt="Centered Image" width="800"/>
  <br/>
  <em>Browser view on parameter sample traceplots, frequently updated while the sampler continues.</em>
</div>

</br>
This monitoring mechanism still needs a lot of polishing, but it can already aid as an "early warning system".



### Producing Posterior Summaries and Visualizations

We next generate some summary statistics and visualizations from draws obtained in the previous example. Of course, one would now check diagnostics first, but this is dependend of your personal workflow.

```
fit$diagnostic_summary()
fit$summary()
fit
bayesplot::mcmc_trace(fit$draws())

```

``bayesutils`` has some summary stats generating functionality as well.

```
mcmc.summary( fit )
```

```
          mean    sd  median PI50.lwr PI50.upr PI90.lwr
lp__   692.645 2.012 692.970  691.527  694.125  688.928
B[1,1]   0.874 0.065   0.875    0.829    0.918    0.765
B[1,2]   1.028 0.071   1.030    0.980    1.076    0.909
B[2,1]   0.835 0.082   0.837    0.779    0.891    0.698
B[2,2]   0.904 0.089   0.906    0.844    0.964    0.755
B[3,1]   1.129 0.080   1.131    1.075    1.183    0.996
B[3,2]   1.252 0.086   1.254    1.193    1.311    1.109
a        0.091 0.017   0.090    0.079    0.102    0.062
sigma    0.075 0.003   0.075    0.073    0.077    0.070
       PI90.upr
lp__    695.266
B[1,1]    0.980
B[1,2]    1.143
B[2,1]    0.969
B[2,2]    1.052
B[3,1]    1.258
B[3,2]    1.394
a         0.120
sigma     0.080
```

These stats can be quickly visualized in the usual sense:


```
mcmc.forestplot(fit, vlines.at = 0)
```

<div align="center">
  <img src="./demo/forestplot.svg" alt="Centered Image" width="600"/>
</div>



Changing the variavble labels gives these inference some context:


```
altnames <- expand.grid( levels(d$sex), levels(d$species)) %>%
  apply(1, paste,simplify = F, collapse = "-") %>%
  unlist() %>%
  c("alpha", "sigma")

mcmc.forestplot(draws, vlines.at = 0,labels = altnames)
```


<div align="center">
  <img src="./demo/forestplot.alt.svg" alt="Centered Image" width="600"/>
</div>



Alternatively you can extract the samples from the fit object and use these directly instead:

```
draws <- extract_samples( fit )
```

Users of the `rethinking` package are familiar with the structure of the resulting object, which is a named list of matrices and vectors.

### Comparing Inferences

Continuing from the example, consider a second model which is similar to the first, but only stratifications happens on the level of the sexes. Thus the Stan model block is:

```
cmdl2 <- "
  B ~ normal(0, 1);
  a ~ normal(0, 1);
  sigma ~ exponential(1);
  for(i in 1:N) m[i] ~ normal( B[ sex[i]] + a * l[i], sigma  );
"
```

Notice that `B` is not a matrix anymore. This becomes apparent in the parameters block:

```
cpar2 <- "
  vector[2] B;
  real a;
  real<lower = 0> sigma;
"
```

Since the data block remains unchanged (Stan does not care about unused data), we can compose and compile the model quickly:


```
sc.to.stan.file("peng2.stan", code = code2)
m2 <- cmdstan_model("peng2.stan")

fit2 <- m2$sample(data = dat, parallel_chains = 2)
```

How do the inferences compare to the first model? A quick visualization to put them side by side:

```
mcmc.forestplot( list(fit, fit2), xlim = c(-0.25, 1.5) )

```

<div align="center">
  <img src="./demo/forestplot.compare.svg" alt="Centered Image" width="600"/>
</div>



Including labels:

```
mcmc.forestplot( list("Sex+Spec" = fit, "Sex" = fit2), xlim = c(-0.4, 1.5) )

```

<div align="center">
  <img src="./demo/forestplot.compare.lbl.svg" alt="Centered Image" width="600"/>
</div>

Focusing on parameters present in both models:

```
mcmc.forestplot( list(fit, fit2), vars = c("a", "sigma") )
```



<div align="center">
  <img src="./demo/forestplot.compare.sel.svg" alt="Centered Image" width="600"/>
</div>

Currently, the 50%, 90% and 99% credible intervals are shown. This can be quickly modified if needed -- here we tell the plotting function to use 50% and 80% credible intervals instead:

```
piplotter <- PI.boxed.plotter( p = c(0.5, 0.8),
                               col = c(bu.color(2, alpha = 1),
                                       bu.color(2, alpha = 0.5)
                                       )
)

mcmc.forestplot( list(fit, fit2), vars = c("a", "sigma"), pi.plotter = piplotter)

```

<div align="center">
  <img src="./demo/forestplot.compare.sel.col.svg" alt="Centered Image" width="600"/>
</div>



Finally, we visualize the inferences for the predictor coefficient as a density plot, highlighting the non-overlapping 90% credible intervals.

```
plot.dens( draws2$a, col = bu.color(1), xlim = c(0, 0.25))
lines.dens(draws$a, col = bu.color(2))
polyg.intv(x = draws2$a, prob = 0.90, col = acol( bu.color(1) ))
polyg.intv(x = draws$a, prob = 0.90, col = acol( bu.color(2) ))
```

<div align="center">
  <img src="./demo/dens.plot.svg" alt="Centered Image" width="600"/>
</div>

The same plot with more CI information visualized:

```
ggplot.densities(
  xs = list( "alpha (model 1)" = draws$a,
             "alpha (model 2)" = draws2$a
  ),
  col_densline = c(bu.color(2), bu.color(1)),
  pi_lvls = c(0.5, 0.90)
)
```

<div align="center">
  <img src="./demo/ggdens.plot.svg" alt="Centered Image" width="600"/>
</div>


## Checking Prior and Posterior Draw Densities for Consistency

The `bayesutils` package further allows to quickly contrast marginal posterior and prior draw distributions. This can be helpful for checking for semi- or non-identifiability of a given model. Semi-identifiabilty means that for some subset of parameters, conditioning the model on the data does not imply a change from the prior assumptions; thus $\Pr(\Theta | X) = \Pr(\Theta)$, where $\Theta$ denotes the parameter subset. The remainder set is however updated consistently with the Bayesian framework -- yet the update is conditional on the prior assumptions for $\Theta$. Sometimes such modelling scenarios with semi-identifiability are wanted and the only means for successfull inference. However, semi-identifiability can be an unwanted pathological model defect, being undetected and ruining inference when not checked for. Such a posterior-prior check could be to inspect empirical descriptive statistics of obtained draws, or to visually compare marginal prior and posterior draw distributions. This package provide visualizations and textual outputs to quickly contrast empirical draw distributions *a priori versus a posteriori*.

To illustrate, we continue with the above example. We here construct a regression model, again with log-mass being the outcome variable. Sex and bill length are predictors.

We model the log-mass likelihood by again assuming normal residual errors. The expected values are constructed through fixed-effects and varying intercepts dependent on sex. The likelihood reads as

$$
m_i \sim \text{N}(\alpha l_i + (\beta_0 + \beta_{[s_i]}), \sigma)
$$

where $l_i$ is the i-th bill length and $s_i$ is an indicator for the sex of the i-th penguin. $\beta_0$ is a latent "baseline intercept". This model is non-sensical. Without connecting $\beta_0$ to the data (f.e. by centering $\beta_0$ and $\beta$), this parameter is not identifiable. We shall explore, how we can detect the issue using `bayesutils` built-in check.

Let's built the model code first:

```
cpar <- "
  vector[2] b;
  real b0;
  real a;
  real<lower = 0> sigma;
"
```

And we implement the Stan model block accordingly. To drive the point of this example home, we only put informative priors on $\beta_0$, $\alpha$ and $\sigma$. For $\beta$ we use an "infinite" flat, non-informative prior by not stating any assumptions for $\beta$ in the model block:

```
cmdl <- "
  b0 ~ normal(0, 1);
  a ~ normal(0, 1);
  for(i in 1:N) m[i] ~ normal( (b[ sex[i] ] + b0) + a * l[i]  , sigma);
"
```

Since we want to contrast posterior draws with prior draws, we need a second model for a prior predictive simulation. We implement it by simply removing the likelihood from the model block:

```
cmdl_priorpred <- "
  b0 ~ normal(0, 1);
  a ~ normal(0, 1);
"
```

We next compose, write and compile the Stan files. The data block code is borrowed from the examples above.

```
code_priorpred <- stancode( cdat, cpar, cmdl_priorpred )
code <- stancode( cdat, cpar, cmdl )

sc.to.stan.file(filename = "m_not-identifiable_priorpred.stan", code = code_priorpred)
  sc.to.stan.file(filename = "m_not-identifiable.stan", code = code)

m_pp <- cmdstan_model("m_not-identifiable_priorpred.stan")
m <- cmdstan_model("m_not-identifiable.stan")

```

Fitting the models and extracting draws next:

```
priorpred <- m_pp$sample(data = dat, parallel_chains = 4)
fit <- m$sample(data = dat, parallel_chains = 4)

draws_prior <- extract_samples( priorpred )
draws_post  <- extract_samples( fit )
```

When running the fitting process, the HMC sampler will greet you with unpleasant divergent transition warnings. The pathological issues did not remain undetected by the sampler here. In actual practical work, we should be worried by now.

When we contrast the prior and posterior draw densities visually, we can quickly see that conditioning on the data did not change the draw distribution of $\beta_0$ by much:

```
plot.dens(  draws_prior$b0, col = bu.color(1), lty = 2, xlim = c(-2.5, 2.5))
lines.dens( draws_post$b0, col = bu.color(2))
``` 

<div align="center">
  <img src="./demo/plt_ni_densplots.svg" alt="Centered Image" width="600"/>
</div>

The marginal prior distribution is illustrated by the blue dashed line, the red line indicates the posterior density. Both distributions are essentially identical. As such inferences for the remaining four parameters is conditional on $\beta_0 \sim N(0, 1)$. These semi-identified inferences are as follows:


```
mcmc.summary(fit)

         mean    sd  median PI50.lwr PI50.upr PI90.lwr PI90.upr
lp__  472.139 1.547 472.436  471.351  473.291  469.173  474.013
b[1]    0.588 0.964   0.578   -0.063    1.235   -0.992    2.206
b[2]    0.686 0.965   0.676    0.029    1.339   -0.896    2.299
b0      0.023 0.960   0.038   -0.603    0.677   -1.585    1.594
a       0.172 0.016   0.173    0.162    0.183    0.146    0.198
sigma   0.146 0.006   0.146    0.142    0.150    0.137    0.156
```

To diagnose identifiability issues quickly, we can use the `post.prior.check` function.

```
ppcheck <- post.prior.check( draws_prior, draws_post )
```

This function computes a simple prior-posterior contrast statistic for each parameter, which we will print next:

```
print( ppcheck )

Summary statistics for normalized quantile residuals: 

   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.0017  0.3158  0.5439  0.5646  0.9298  1.0000 


Variables with smallest normalized quantile residuals: 

         b0           a        b[1]        b[2]        lp__ 
0.001699474 0.298245614 0.368421053 0.719298246 1.000000000 
      sigma 
1.000000000 
``` 

The resulting statistics can be rougly interpretted as "distribution dissimilarities". The smaller the statistic, the more similar are the marginal prior and posterior draw distributions. Since the outputs are ordered, we see that $\beta_0$ has by far the largest agreement between posterior and prior. More formally, let $q_p$ be the empirical quantile of the prior distribution for a given parameter to some level $p \in (0, 1)$. Let $D$ be the set of $N$ parameter posterior draws $D = {d_1, \ldots,, d_N}$. Then $f_p$ is the ratio of the number of draws less-equal to $q_p$ and $N$,

$$f_p = \frac{1}{N} \sum_{i=1}^N \mathbf{1}_{\{t_i \leq q_p\}}$$,

which corresponds to the empirical cumulative distribution function (ECDF).

When the posterior and prior distributions are practically identical, we have $q_p = f_p$ for any level $p$. In the package, the 10%, 20%, ..., 90% quantiles are chosen as levels. The constrast statistic is the normalized squared sum of posterior and prior ECDF residuals $|f_p - q_p|^2$. Thus, if the distributions are practically identical, this sum is close to zero. Another interpretation of this statistic is to think of it as a discrepancy between a hypothetical line in a quantile-quantile plot under the assumption of both distributions being identical and the actually obtained line.

The latter notion supports the interpretation of the posterior-prior check visualization:


```
plot( ppcheck,annotate = T )
```  

<div align="center">
  <img src="./demo/plt_ni_ppchech_2.svg" alt="Centered Image" width="600"/>
</div>

This is the aforementioned quantile-quantile plot. Each line in the diagram, from left to right, corresponds to a parameter and the contrast between its posterior and prior distribution. Opaqueness is coded through the contrast statistic. We can quickly see, that $\sigma$ is well-identified (conditionally to $\beta_0$ though). The "jumps" from 0 to 1 in case of $\alpha$ and $\beta$ are due to the flat prior collapsing into conditionally-identified regions in the posterior. As expected, $\beta_0$ is creeping along the quantile diagonal. It is thus this single plot that can aid to detect such identifiability issues in actual works. Any parameters plotted along or somewhat close to the diagonal should raise concerns and should be checked more thoroughly (i.e. by prior calibration checks).

Unchecking the annotate flag (which is `FALSE`per default)  still communicates the issue, however in a more clean way where, however, the parameter legend is hidden:

```
plot( ppcheck )
```  

<div align="center">
  <img src="./demo/plt_ni_ppchech_1.svg" alt="Centered Image" width="600"/>
</div>

Note that the number of parameters reported in the print and plot output is limitted to the 10 smallest parameters w.r.t. the contast statistic, in order to keep the outputs clean. All parameter statistics can be accessed through the object returned by the `post.prior.check` function if needed. If only a subset of parameters shall be evaluated, then a character vector can be supplied.







## Roadmap

As stated earlier, this package is far from being a polished product. Thus there are some feature implementations and refinements to be made in the future. These include:

- sampling monitoring: 
	- replace link list with dropdown menu, since the current solution can quickly clutter the screen when the number of model parameters is large.
	- add option to display histograms including/instead of traceplots. This is of interest f.e. in retrieval testing and calibration checking; as such the monitor can help to detect discrepancies between samples and ground truth early.
	- increase CSV reading speed through dedicated packages. The current base `R` solution can become rather slow with bigger models and large sample numbers.
	- add function argument to specify which parameters to monitor.
- plotting:
	- a lightweight OO design could provide more flexibilty to the user. 
	- maybe moving to ggplot2.
- Stan code handling:
	- export and import, with JSON support
	- provide (teaching) templates 
	
