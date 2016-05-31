##' Simulate an ensemble of outbreaks based on a branching process with negative binomial offspring distribution
##'
##' Content for details - still to be written. (need to give the formula for the negative binomial offspring distribution)
##' @title Outbreak simulation
##' @param R0 Basic reproduction number of the branching process,
##' which also is the mean of the offspring distribution.
##' @param k Dispersion parameter of the offspring distribution.
##' @param serial.interval.sample Vector of a numeric`l sample of
##' generation times.
##' @param n.outbreaks Number of outbreaks to simulate. Default
##' n.outbreaks = 1000.
##' @param outbreak.size.max Maximum number of cases to
##' simulate. Default outbreak.size.max = 1000.
##' @param n.crit No idea what this is at this moment. Default n.crit
##' = 100.
##' @param incidence.interval Interval into which the incidence is
##' aggregated. Default incidence.interva = 7 (reporting weekly
##' incidence, if time is measured in days).
##' @param n.output Maximum number of cases per outbreak to output. Default n.output = 1000. 
##' @param p.report Probability that a case is reported. Default p.report = 1, needs to be 0 < p.report <= 1.
##' @return A list of various stuff.
##' @author Tini
##' @export
simulate.branching.process = function(R0, k, serial.interval.sample, n.outbreaks = 1000, outbreak.size.max = 1000, n.crit = 100, incidence.interval = 7, n.output = 1000, p.report = 1) {

    ## matrix of case id's
    mat.id = matrix(1:outbreak.size.max, nrow = outbreak.size.max, ncol = n.outbreaks)

    ## matrix of individual R's based on the R0 value, sampled from a Poisson distribution
    ## mat.R.planned = matrix(rpois(outbreak.size.max * n.outbreaks, lambda = R0), nrow = outbreak.size.max, ncol = n.outbreaks)
    mat.R.planned = matrix(rnbinom(n = outbreak.size.max * n.outbreaks, size = k, mu = R0), nrow = outbreak.size.max, ncol = n.outbreaks)
    
    ## matrix of source id's
    mat.source.id = rbind(rep(0, n.outbreaks), sapply(1:n.outbreaks, FUN = function(outbreak.id) rep(1:outbreak.size.max, times = mat.R.planned[, outbreak.id])[1:(outbreak.size.max-1)]))
    
    ## matrix of generation times
    mat.gen.time = matrix(sample(serial.interval.sample, size = outbreak.size.max * n.outbreaks, replace = TRUE), nrow = outbreak.size.max, ncol = n.outbreaks)

    ########################################
    ## ptm = proc.time()
    ## for(outbreak.id in 1:n.outbreaks) {
    ##   mat.gen.time[, outbreak.id] = mat.gen.time[order(mat.source.id[, outbreak.id], mat.gen.time[, outbreak.id]), outbreak.id]
    ## }
    ## proc.time() - ptm
    ########################################
    
    ## matrix of infection times and generations
    mat.inf.time = mat.generation = matrix(NA, nrow = outbreak.size.max, ncol = n.outbreaks)
    mat.inf.time[1, ] = 0
    mat.generation[1, ] = 1
    for(outbreak.id in 1:n.outbreaks) {
        for(case.id in 2:outbreak.size.max) {
            mat.inf.time[case.id, outbreak.id] = mat.inf.time[mat.source.id[case.id, outbreak.id], outbreak.id] + mat.gen.time[case.id, outbreak.id]
            mat.generation[case.id, outbreak.id] = mat.generation[mat.source.id[case.id, outbreak.id], outbreak.id] + 1
        }
    }
    mat.source.id[is.na(mat.inf.time)] = NA ## knocking out cases that weren't infected.
    mat.gen.time[is.na(mat.inf.time)] = NA
    mat.R.planned[is.na(mat.inf.time)] = NA

    ## mat.generation[mat.report == 0] = NA
    

    
    ## which cases have completed all infections they were planning to make?
    mat.R.actual = matrix(0, nrow = outbreak.size.max, ncol = n.outbreaks)
    for(outbreak.id in 1:n.outbreaks) {
        tt = table(mat.source.id[, outbreak.id])[-1] ## take off the 0 from the introduction
        mat.R.actual[as.numeric(names(tt)), outbreak.id] = tt
    }
    mat.R.actual[is.na(mat.inf.time)] = NA
    mat.R.complete = (mat.R.planned - mat.R.actual == 0)
    
    ## outbreak sizes & extinction probability:
    outbreak.size = colSums(!is.na(mat.inf.time))
    p.extinction = 1 - sum(outbreak.size == outbreak.size.max)/n.outbreaks

    
    ## until which case have we got a complete line list?
    time.complete = sapply(1:n.outbreaks, function(outbreak.id) min(mat.inf.time[!mat.R.complete[, outbreak.id] & !is.na(mat.R.complete[, outbreak.id]), outbreak.id]))
    n.complete = sapply(1:n.outbreaks, function(outbreak.id) sum(sort(mat.inf.time[, outbreak.id]) <= time.complete[outbreak.id]))
    
    time.crit = sapply(n.crit/p.report, function(n.crit) sapply(1:n.outbreaks, function(outbreak.id) sort(mat.inf.time[, outbreak.id])[n.crit])) ## based on infections, scaled up to account for reporting.
    
    ## theoretical extinction probability based on the generating function of the Poisson distribution which I've assumed as offspring distribution. 
    ## p.extinction.theo = ifelse(R0 > 1, uniroot(f = function(z) (log(z)/(z-1) - R0), lower = 0, upper = 1 - 1e-10)$root, 1)
    p.extinction.theo = ifelse(R0 > 1, uniroot(f = function(z) ((1 + R0/k*(1-z))^(-k) - z), lower = 0, upper = 1 - 1e-10)$root, 1)

    ## dealing with under-reporting:
    ## matrix indicating whether a case is reported or not:
    mat.report = matrix(rbinom(n = length(mat.id), size = 1, prob = p.report), nrow = outbreak.size.max, ncol = n.outbreaks)
    mat.inf.time[mat.report == 0] = NA ## knocking out unreported cases
    outbreak.size.reported = colSums(!is.na(mat.inf.time))
    time.crit.reported = sapply(n.crit, function(n.crit) sapply(1:n.outbreaks, function(outbreak.id) sort(mat.inf.time[, outbreak.id])[n.crit])) ## based on severe cases. 
    

    ## calculating incidence:
    max.week = floor((max(mat.inf.time, na.rm = TRUE)+incidence.interval)/incidence.interval)
    mat.incidence = matrix(0, nrow = max.week, ncol = n.outbreaks)
    tt = table(cut(mat.inf.time, breaks = seq(0, incidence.interval * max.week, by = incidence.interval), labels = 1:max.week, right = FALSE), rep(1:n.outbreaks, each = outbreak.size.max))
    mat.incidence[as.numeric(rownames(tt)), ] = tt

    ## sorting the outbreaks by infection time prior to output:
    mat.order = apply(mat.inf.time, 2, order)
    ## tmp = mat.id
    ## tmp2 = sapply(1:n.outbreaks, function(outbreak.id) { mat.id[, outbreak.id] = mat.id[mat.order[outbreak.id], outbreak.id] })
    for(outbreak.id in 1:n.outbreaks) {
      mat.id[, outbreak.id] = mat.id[mat.order[, outbreak.id], outbreak.id]
      mat.source.id[, outbreak.id] = mat.source.id[mat.order[, outbreak.id], outbreak.id]
      mat.inf.time[, outbreak.id] = mat.inf.time[mat.order[, outbreak.id], outbreak.id]
      mat.generation[, outbreak.id] = mat.generation[mat.order[, outbreak.id], outbreak.id]
      mat.R.planned[, outbreak.id] = mat.R.planned[mat.order[, outbreak.id], outbreak.id]
      mat.R.actual[, outbreak.id] = mat.R.actual[mat.order[, outbreak.id], outbreak.id]
      mat.R.complete[, outbreak.id] = mat.R.complete[mat.order[, outbreak.id], outbreak.id]
      
    }

    res = list(id = mat.id[1:n.output, ],
        source.id = mat.source.id[1:n.output, ],
        inf.time = mat.inf.time[1:n.output, ],
        generation = mat.generation[1:n.output, ],
        R.planned = mat.R.planned[1:n.output, ],
        R.actual = mat.R.actual[1:n.output, ],
        R.complete = mat.R.complete[1:n.output, ],
        incidence = mat.incidence,
        outbreak.size = outbreak.size,
        outbreak.size.reported = outbreak.size.reported,
        time.complete = time.complete,
        n.complete = n.complete,
        time.crit = time.crit,
        time.crit.reported = time.crit.reported, 
        p.extinction = c(obs = p.extinction, theo = p.extinction.theo))

    return(res)
}
################################################################################
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @return 
##' @author 
my.fun = function() {

}
