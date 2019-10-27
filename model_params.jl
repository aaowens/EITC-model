paramsgen = @with_kw (
 age_lower = 20,
 age_upper = 40,
 mean_eps = 0.3,
 sigma_eps = 0.2,
 sigma_cost = 0.2,
 mean_cost = 0.6,
 beta_m = 0.,
 sigma_c = 1.01,
 alpha_h = 0.5,
 sigma_h = 2.,
 kidcost_labor = 0., # 1.3, # amplifies alpha_h
 kidcost_switch = 0., #1.3, # amplifies switching cost
 beta = 0.8,
 cbar = 1.,
 eitcparams = makeeitc())


"Compute the EITC schedule parameters as a function of desired rates, max credit, and flat region size"
function makeeitc(;rate1 = 0.2, rate2 = 0.2, max_credit = 0.2, flat_size = 0.05)
    cutoff1 = max_credit / rate1
    cutoff2 = cutoff1 + flat_size
    cutoff3 = cutoff2 + max_credit / rate2
    eitcparamsgen(cutoff1 = cutoff1, cutoff2 = cutoff2, cutoff3 = cutoff3,
    rate1 = rate1, rate2 = rate2)
end

 eitcparamsgen = @with_kw (
 cutoff1 = 0.16,
  cutoff2 = 0.2,
  cutoff3 = 0.36,
  rate1 = 0.2,
  rate2 = 0.2)
