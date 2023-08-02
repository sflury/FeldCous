'''
Stats functions for use in analysis of COS spectra

Dependencies:
    numpy 1.17.3
    scipy 1.3.1

Sophia Flury 2021.02.20
'''
from numpy import arange,argsort,array,concatenate,cumsum,diff,exp,interp,isfinite,log,sqrt,sum,where,zeros
from numpy.random import poisson,rand,randn,seed
from scipy.special import erf,ndtri,gammaincc,gammaln,iv
'''
Name:
    calc_source

Purpose:
    calculate the confidence intervals on a Poisson source with background
    using Feldman & Cousins if gross > background, bootstrapping if gross <
    background, or sqrt(gross+background) if source gross-background > 100.

Arguments:
    :gross (*float*): gross measured counts
    :bkg (*float*): estiamted background counts
    :bkg_err (*float*): uncertainty in background counts

Keyword Arguments:
    :p_det (*float*): threshold probability for detection. Default is 0.1587.

Returns:
    :src (*float*): source counts
    :src_lower (*float*): 15.87th percentile error on source counts
    :src_upper (*float*): 84.13th percentile error on source counts

'''
def calc_source(gross,bkg,p_det=0.1587,sigma=1.):
    # +/- 1 sigma quantiles
    quant = 0.5*(1+erf(array([-1,1])/sqrt(2)))
    # frequentist: source = gross - bkg
    src = max([gross - bkg,0])
    # if large number of gross counts, solution converges on sqrt(gross)
    if gross > 256:
        err_lo,err_hi = poiss_conf_lim(gross,quant)
        src_err_lower = err_lo-bkg
        src_err_upper = err_hi-bkg
    # Feldman Cousins method if detected
    elif calc_det_prob(gross,bkg)[0] < p_det:
        err_lo,err_hi = feldman_cousins(gross,bkg,quant)
        src_err_lower = src-max([err_lo[0],0])
        src_err_upper = err_hi[0]-src
    # bootstrap method for negative counts
    elif bkg > 0:
        err_lo,err_hi = poiss_conf_lim(bkg,quant)
        src_err_lower = err_lo-bkg
        src_err_upper = err_hi-bkg
    # failsafe
    else:
        src_err_lower,src_err_upper = 0.,0.
    return src,src_err_lower,src_err_upper
'''
Name:
    calc_conf

Purpose:
    Calculate the confidence intervals on a Poisson source with background
    using Feldman & Cousins if gross > background, bootstrapping if gross <
    background, or sqrt(gross+background) if source gross > 256.

Arguments:
    :gross (*float*): gross measured counts
    :bkg (*float*): estiamted background counts
    :p_conf (*float*): probability intervals of the desired confidence intervals
                            (e.g., 0.68 for 1-sigma confidence interval)
Keyword Arguments:
    :p_det (*float*): threshold probability for detection. Default is 0.1587.

Returns:
    :conf_lower (*ndarray*): lower confidences on source counts
    :conf_upper (*ndarray*): upper confidences on source counts

'''
def calc_conf(gross,bkg,p_conf,p_det=0.1587):
    # note actual quantiles are
    # if large gross counts, everything is just Poisson
    if gross-bkg > 256:
        conf_lower = poiss_conf_lim(gross,0.5-0.5*p_conf)-bkg
        conf_upper = poiss_conf_lim(gross,0.5+0.5*p_conf)-bkg
    # Feldman Cousins method for small counts
    elif calc_det_prob(gross,bkg)[0] < p_det :
        conf_lower,conf_upper = feldman_cousins(gross,bkg,p_conf)
    # CDF method for negative counts and non-detections
    else:
        conf_lower = poiss_conf_lim(bkg,0.5-0.5*p_conf)-bkg
        conf_upper = poiss_conf_lim(bkg,0.5+0.5*p_conf)-bkg
    return conf_lower,conf_upper
'''
Name:
    poiss_conf_lim

Purpose:
    Compute confidence limits on a measured/known background for given quantiles
    assuming the background is the mean of a Poisson distribution.

Arguments:
    :mu (*float*): known/measured counts
    :p_conf (*float* or *ndarray*): probability intervals of the desired
                                confidence intervals

Returns:
    :conf_lower (*ndarray*): lower confidence intervals
    :conf_upper (*ndarray*): upper confidence intervals
'''
def poiss_conf_lim(mu,p_conf):
    # check dtype on bkg
    if not hasattr(p_conf,'__len__'):
        p_conf = array([p_conf])
    # Poisson CDF
    k = arange(0,mu+5*sqrt(mu),1)
    cdf = gammaincc(k+1,mu)
    # confidences from quantiles
    return interp(p_conf,cdf,k)
'''
Name:
    skellam_conf_lim

Purpose:
    Compute confidence limits on a measured/known difference of Poisson variates
    for given quantiles assuming.

Arguments:
    :mu1 (*float*): known/measured counts 1 (gross signal)
    :mu1 (*float*): known/measured counts 2 (background)
    :p_conf (*float* or *ndarray*): probability intervals of the desired
                                confidence intervals

Returns:
    :conf_lower (*ndarray*): lower confidence intervals
    :conf_upper (*ndarray*): upper confidence intervals
'''
def skellam_conf_lim(mu1,mu2,p_conf):
    # Skellam CDF
    nmax = max([mu1,mu2])
    n = arange(0,nmax+5*sqrt(nmax),1)
    cdf = [skellam_cdf(k,mu1,mu2) for k in n]
    # confidences from quantiles
    return interp(p_conf,cdf,n)
'''
Name:
    calc_det_prob

Purpose:
    Calculate the probability that measured counts were sampled from the
    distribution of background counts.

Arguments:
    :sig (*float*): gross signal counts
    :bkg (*float*): background counts

Returns:
    :p_nb (*float*): probability that the signal is sampled from
                        the background
    :signif (*float*): normal significance corresponding to P(>N|B)
                         assuming a normal distribution
'''
def calc_det_prob(sig,bkg):
    # compute the Poisson CDF for signal given background
    # recalling the CDF = sum_(k=0)^n mu^k * exp(-mu) / Gamma(k+1) = Q(k+1,mu)
    # where Q is the regularized incomplete Gamma function
    p_nb = gammaincc(sig+1,bkg)
    # if greater than 1 due to machine precision issues, set to 1
    if p_nb > 1.:
        p_nb = 1.
    # significance is probit(p)
    signif = ndtri(p_nb)
    if not isfinite(signif) and p_nb > 0.99:
        signif = 8.21
    if signif < 0:
        signif = 0.
    # probability of sampling from background is 1-detection prob
    return 1.-p_nb,signif
'''
Name:
    poisson_pmf

Purpose:
    Compute the Poisson probability of getting counts n given some mean mu.

Arguments:
    :n (*float*): number of counts
    :mu (*float*): mean counts of distribution

Returns:
    :cdf (*float*): Poisson probability mass function at n
'''
def poisson_pmf(n,mu):
    # compute the Poisson probability, recalling x!=gamma(x+1)
    return exp( n*log(mu) - (mu + gammaln(n+1)) )
'''
Name:
    skellam_pmf

Purpose:
    Compute the Skellam probability of getting some difference n between counts
    given means mu1 and mu2 for two Poisson distributions.

Arguments:
    :n (*float*): number of counts
    :mu1 (*float*): mean counts of distribution 1
    :mu2 (*float*): mean counts of distribution 3

Returns:
    :pmf (*float*): Skellam probability mass function n ~ mu1 - mu2
'''
def skellam_pmf(n,mu1,mu2):
    # compute Skellam probability (difference of two Poisson variates)
    return exp(-(mu1+mu2))*(mu1/mu2)**(n/2)*iv(n,2*sqrt(mu1*mu2))
def skellam_cdf(n,mu1,mu2):
    return sum([skellam_pmf(float(k),mu1,mu2) for k in range(int(n))])
'''
Name:
    sample_source

Purpose:
    Sample deviates from the likelihood distribution of a measured signal given
    the gross signal and the background. Determines confidence intervals of
    the source signal given the gross signal, background, and a range of
    probabilities for the desired confidence ratios. Builds log likelihood
    distribution from the confidence intervals and probabilities and samples
    this likelihood distribution a given number of times.

Arguments:
    :gross (*float*): gross signal
    :bkg (*float*): background noise in signal
    :p_conf (*float* or *ndarray*): probability intervals of the desired
                            confidence intervals

Keyword Arguments:
    :n_samp (*int*): number of samples to draw from the likelihood distribution
    :p_det (*float*): threshold probability for detection. Default is 0.1587.

Returns:
    :samples (*ndarray*): 1xn_samp array of samples of the likelihood
                            distribution for source = gross - background
'''
def sample_source(gross,bkg,p_conf,n_samp=int(100000),p_det=0.1587):
    # get confidence intervals
    conf_lower,conf_upper = calc_conf(gross,bkg,p_conf,p_det=p_det)
    # quantiles corresponding to probability intervals
    quants = concatenate([0.5-0.5*p_conf[::-1],[0.5],0.5+0.5*p_conf])
    # get standard deviations of confidence intervals
    sig = ndtri(quants)
    ## log likelihoods
    ll = -0.5*sig**2
    # if source is strongly detected
    if gross-bkg > 256:
        # merge confidence intervals
        #gconf = poiss_conf_lim(gross,quants)
        #bconf = poiss_conf_lim(bkg,quants)
        # sample deviates from likelihood distribution using inverse transform
        #samples = invtrans_interp(gconf,ll,n_samp=n_samp) - \
        #            invtrans_interp(bconf,ll,n_samp=n_samp)
        #samples = invtrans_interp(conf,ll,n_samp=n_samp)
        gross_samp = poisson(gross,n_samp).astype('float64')
        bkg_samp = poisson(bkg,n_samp).astype('float64')
        samples = gross_samp - bkg_samp
    # if source is detected, include in distribution
    elif calc_det_prob(gross,bkg)[0] < p_det :
        # Feldman & Cousins confidence intervals
        conf_lower,conf_upper = feldman_cousins(gross,bkg,p_conf)
        # merge confidence intervals
        conf = concatenate([conf_lower[::-1],[gross-bkg],conf_upper])
        cdf = cumtrapz(conf,exp(ll))
        ## clip out extra zeros and bad upper limits (max machine precision)
        i_min = int(0)
        if len(where(conf<=0.01)[0])>0:
            i_min = int(where(conf<=0.01)[0][-1])
        clip = [i for i in range(i_min,len(conf)) if conf[i] < 2**50]
        sig,conf,ll = sig[clip],conf[clip],ll[clip]
        # sample deviates from likelihood distribution using inverse transform
        samples = invtrans_interp(conf,ll,n_samp=n_samp)
    # otherwise, sample background
    else:
        # merge confidence intervals
        #conf = poiss_conf_lim(bkg,quants)-bkg
        # sample deviates from likelihood distribution using inverse transform
        #samples = invtrans_interp(conf,ll,n_samp=n_samp)
        samples = poisson(bkg,n_samp).astype('float64')-bkg
    return samples#,conf,exp(ll)
'''
Name:
    cumtrapz

Purpose:
    Compute the cumulative integral of a given integrand using trapezoid rule.
    Optionally normalizes so that the integral goes to 1 over the range of x.

Arguments:
    :x (*ndarray*): integrand abscissa
    :y (*ndarray*): integrand ordinate

Keyword Arguments:
    :norm (*bool*): Optional normalization so that total integral over x is 1.

Returns:
    :cumint (*ndarray*): cumulative integral at each value of x
'''
def cumtrapz(x,y,norm=True):
    # trapezoid integral for non-uniform spacing in x
    cumint = array(list(map(lambda i: \
                            0.5*sum((y[:i]+y[1:i+1])*diff(x[:i+1])),\
                            range(len(x)))))
    # normalize
    if norm:
        cumint/=cumint[-1]
    return cumint
'''
Name:
    inv_trans_interp

Purpose:
    Use inverse transform method to sample deviates from a log likelihood
    function given only a set of nodes defining the log likelihood. CDF is
    calculated using trapezoid rule. Variable x is interpolated with respect to
    the CDF to translate uniform deviates into the desired distribution.

Arguments:
    :x (*ndarray*): nodes of variable x at which the log likelihood is known
    :ll (*ndarray*): log likelihood nodes of values x

Keyword Arguments:
    :n_samp (*int*): number of deviates to sample from the interpolated log
                            likelihood distribution

Returns:
    :samples (*ndarray*): samples of log likelihood
'''
def invtrans_interp(x,ll,n_samp=int(100000)):
    # set seed
    #seed(123)
    # use cumulative trapezoid rule to calculate CDF
    cdf = cumtrapz(x,exp(ll),norm=True)
    # sample some uniform variates
    uvar = rand(n_samp)
    # interpolate to get deviates
    return interp(uvar,cdf,x)
'''
 Name:
    feldman_cousins

Purpose:
    Compute Neyman-Pearson confidence intervals following Feldman & Cousins
    1998, Phys.Rev.D.,57,3873, as originally implemented by Worseck et al 2016,
    ApJ 825, 144, and Makan et al 2020, in C.

Arguments:
    :gross (*float*): measured gross counts
    :bkg (*float*): measured/model background counts
    :p_conf (*float*, *list*, or *np.ndarray*): probability confidence
                            interval(s) to evaluate mu

Keyword Arguments:
    :delta_mu (*float*): resolution in the simulated measurements mu. Small
                            values will improve contraints but require much more
                            computation time. Lower limit of 1e-4 and upper
                            limit of sqrt(gross) are imposed for speed and
                            sufficient sampling. Default is delta_mu = 0.01

Returns:
    :conf_lower (*np.ndarray*): 1xN array of lower confidence intervals where
                            N is the number of provided confidence values
    :conf_upper (*np.ndarray*): 1xN array of upper confidence intervals where
                            N is the number of provided confidence values
'''
def feldman_cousins(gross,bkg,p_conf,delta_mu=0.01):
    # source is gross - background counts
    src = gross - bkg
    # pseudo error is Poissonian sqrt(gross) counts
    sq_gross = sqrt(gross)
    # check on delta_mu
    if delta_mu < 1e-4 or delta_mu > sq_gross:
        print(f"Provided delta_mu={delta_mu} is out of bounds.\n"+
             "Reverting to default delta_mu=0.01.")
        delta_mu = 0.01
    # check on p_conf data type
    if not hasattr(p_conf,'__len__'):
        p_conf = [p_conf]
    # mu values to consider range over (gross-bkg) +/- 5*sqrt(gross)
    # for increments of delta_mu in mu, imposing mu > 0
    mu_init = max([src-5*sq_gross,0])//1.
    mu_init = arange(mu_init,src+5*sq_gross+delta_mu,delta_mu)
    # fix round-off errors ("mu" in C script)
    mu = array([round(mui,6) for mui in mu_init])
    # number of mu values ("n" in C script)
    n_mu = len(mu)
    # counts to consider for each mu ("m" in C script)
    n_counts = int(3*gross+101)
    # values of n ("must_test" in C script)
    n_vals = arange(0,n_counts,1.)//1.
    # max(n-bkg,0) ("mubest" in C script)
    mu_b = array([n-bkg if n-bkg > 0 else 0 for n in n_vals])
    # mu confidence values for each simulated source mu
    mu_min = zeros((n_mu,len(p_conf)))+2**52
    mu_max = zeros((n_mu,len(p_conf)))
    # for all mu values
    for i in range(n_mu):
        # Poisson log likelihood to prevent float overflow
        lnfact = gammaln(n_vals+1)
        # probability for assumed mu given n counts
        prob = exp(n_vals*log(mu[i]+bkg)-(mu[i]+bkg)-lnfact)
        # probability for test mu (mu_b or "mubest" in C) given n counts
        prob_test = exp(n_vals*log(mu_b+bkg)-(mu_b+bkg)-lnfact)
        # ratio of mu probabilities given n_vals
        # to test counts mu_b probabilities given n_vals
        ratio = prob/prob_test
        # sort indices from largest to smallest probability ratio
        inds = argsort(ratio)[::-1]
        # cumulative sum of probabilities along sorted ratio
        cprob = cumsum(prob[inds])
        # get confidences for each prescribed confidence interval
        for j,p_max in enumerate(p_conf):
            # index where sum probabilities ~ confidence interval probability
            l = len(where(cprob<p_max)[0])+1
            # max and min mu vals are max and min n vals within interval
            mu_max[i,j] = max(n_vals[inds[:l]])
            mu_min[i,j] = min(n_vals[inds[:l]])
    # confidence interval arrays ("errordown" and "errorup" in C code)
    conf_lower = zeros(len(p_conf))
    conf_upper = zeros(len(p_conf))+2**52
    # get confidences for each prescribed confidence interval
    for j,p_max in enumerate(p_conf):
        # check max mu values to get lower confidence intervals
        # where max mu corresponds to gross counts
        min_inds = where(abs(mu_max[:,j]-gross) < delta_mu)[0]
        if len(min_inds) > 0:
            conf_lower[j] = min(mu[min_inds])
        # check min mu values to get upper confidence intervals
        # where min mu corresponds to gross counts
        max_inds = where(abs(mu_min[:,j]-gross) < delta_mu)[0]
        if len(max_inds) > 0:
            conf_upper[j] = max(mu[max_inds])
    # return confidence intervals
    return conf_lower,conf_upper
