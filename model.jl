"CRRA utility function. Preference parameters depend on state"
function u(c, l, state, params)
    @unpack sigma_c, sigma_h, alpha_h, kidcost_labor = params
    @unpack haskid = state
    alpha_h = (alpha_h + alpha_h*haskid*kidcost_labor)
    c <= 0 && return -Inf
    l >= 1 && return -Inf
    c^(1 - sigma_c) / (1 - sigma_c) + alpha_h*(1 - l)^(1 - sigma_h) / (1 - sigma_h)
end

"Consumption given the choices, state"
function consumption(choices, state, params)
    @unpack l, switching = choices
    @unpack ci, haskid = state
    @unpack beta_m, cbar, kidcost_switch = params
    w = wagefun(state, params)
    w*l + tau(w*l, state, params) - ci*(1 + haskid * kidcost_switch) * switching + cbar
end

"Wage as a function of state"
function wagefun(state, params)
    @unpack eps, manager = state
    @unpack beta_m = params
    eps + eps * manager * beta_m
end

"Maximum labor supply before earned income hits the end of the EITC schedule"
function maxl(state, params)
    w = wagefun(state, params)
    # What is the max l which hits the upper EITC?
    eitc_up = params.eitcparams.cutoff3
    # w*l = eitc_up
    eitc_up / w
end

"Intratemporal objective function. Evaluates utility given choices and state"
function intra_obj(choices, state, params)
    @unpack l = choices
    c = consumption(choices, state, params)
    u(c, l, state, params)
end

```
Generate the model object as a function of parameters
Produces markov chain of income given mean and variance
Currently there is no transition between states. P = I
Produces discretized state space for switching cost, age, manager status, and
    kid status
```
function makemodel(params)
    N = 50
    d = Normal(params.mean_eps, params.sigma_eps)
    e = expectation(d, n = N)
    mc = QuantEcon.MarkovChain(I(N), exp.(e.nodes))
    epss = mc.state_values
    @unpack mean_cost, sigma_cost = params
    cis = (mean_cost - 2sigma_cost):(sigma_cost/4):(mean_cost + 2sigma_cost)
    statevec = (age = 20:1:50, manager = (false, true), eps = epss, ci = cis, haskid = (false, true))
    Model(params, mc, statevec)
end

"Transform a state index into a state vector"
function makestate(model, state_nx::CartesianIndex)
    statevec = model.states
    age = statevec.age[state_nx[1]]
    manager = statevec.manager[state_nx[2]]
    eps = statevec.eps[state_nx[3]]
    ci = statevec.ci[state_nx[4]]
    haskid = statevec.haskid[state_nx[5]]
    State(age = age, manager = manager, eps = eps, ci = ci, haskid = haskid)
    #State(age, manager, eps, ci, haskid)
end


## Discrete choice. Switch, no switch
## Then because of the non-convexity, we
# need to consider whether they will be in the EITC region
# or not

# VF is over the whole state
# The only uncertainty is in eps, which is IID
# So EVF is just an integral over the eps variable
# Hours doesn't actually influence the state, which is good
# The only choice that influences the state is switch
# So, the conditional problem is totally static
# Just choose best hours

# If in EITC, compute upper bound on hours so that
# earned income is the top cutoff
# If not, do the reverse
# So, need to solve 4 optimization problems at each step

# Case 1: No switch, in EITC
#state params

```
Inner step for the intratemporal optimization problem, conditonal
    on the switch choice and whether they are taking the EITC or not
Conditional on that, optimal hours satisfies the FOC
This solves for hours l
```
function inner_step(switching, bounds, state, params)
    function obj(l)
        choices = Choice(l = l, switching = switching)
        intra_obj(choices, state, params)
    end
    res = optimize(x -> -obj(x), bounds.l, bounds.u)
    res.converged == false && error()
    rmax = -Optim.minimum(res); rmaximizer = Optim.minimizer(res)
    (rmax = rmax, rmaximizer = rmaximizer, switching = switching, bounds = bounds)
end

```
Performs a bellman iteration for a particular state index
**Arguments**
* `VF` : The value function array
* `state_nx` : The index of the state to iterate on
* `model` : The model structure

Computes expected value as a function of the switching choice
Computes value for each intertemporal choice and EITC combination using
    the inner step function
Chooses the best of those
Returns the optimal choices and the value of the state
```
function bellman_step(VF, state_nx, model)
    states = model.states
    params = model.params
    @unpack beta = params
    EVF = expected_VF(VF, state_nx, model) # tuple of switching choice
    state = makestate(model, state_nx)
    switchranges = (false, true)
    up = clamp(maxl(state, params), 0., 1.)
    eitcranges = ((l = 0., u = up), (l = up, u = 1.))
    out = @SMatrix [inner_step(switchranges[i], eitcranges[j], state, params) for i = 1:2, j = 1:2]
    vals = @SMatrix [out[i, j].rmax + beta*EVF[1 + out[i, j].switching] for i in 1:2, j = 1:2]

    (bestval, bestind) = findmax(vals)
    bestout = out[bestind]
    bestl = bestout.rmaximizer; switching = bestout.switching
    return (choice = Choice(l = bestl, switching = switching), val = bestval)
end

"Returns 2-tuple of EVF, for switching or not switching"
function expected_VF(VF, state_nx::CartesianIndex, model)
    EVFf = EVF_sub(VF, state_nx, 1, model)
    EVFt = EVF_sub(VF, state_nx, 2, model)
    (EVFf, EVFt)
end
```
Given today's current state_nx, we need to integrate over
tomorrow's value function conditional on the state transition rule
given the choices
Age always increases by 1
Manager depends on the switching choice
Epsilon is random, so we integrate over the epsilon distribution
```
function EVF_sub(VF, state_nx, choice_nx, model)
    # Given the current epsilon state, what's the distribution
    # of tomorrow's epsilon?
    # For now, let's make epsilon fixed
    ismanager = (state_nx[2] == 2)
    choice_nx2 = ismanager ? 2 : choice_nx
    nextage = state_nx[1] + 1
    eps_nx = state_nx[3]
    i4, i5 = state_nx[4], state_nx[5]
    d = model.IncomeProcess.p
    sum(VF[nextage, choice_nx2, i, i4, i5] * d[eps_nx, i] for i in axes(VF)[3])
end

```
Solve the model using value function iteration backwards from retirement
**Arguments**
* `VF` : An uninitialized value function array
* `model` : The model structure
***Output***
* `pol1` : Policy function for hours
* `pol2` : Policy function for occupation switching
```
function solve_model(VF, model)
    params = model.params
    state_nxs = CartesianIndices(VF)
    fill!(VF, 0.)
    # Initialize the last period
    # It is the utility of eating a retirement benefit forever
    # The dummystate is needed for the preference parameters, but it is arbitrary
    dummystate = makestate(model, state_nxs[1])
    lastutil = u(1.5, 0., dummystate, params) / (1 - params.beta)
    pol1 = zeros(size(VF))
    pol2 = zeros(Bool, size(VF))
    ageaxis = axes(state_nxs)[1]
    VF[length(ageaxis), :, :, :, :] .= lastutil
    for a_nx in (length(ageaxis) - 1):-1:1
        for state_nx in state_nxs[a_nx, :, :, :, :]
            out = bellman_step(VF, state_nx, model)
            pol1[state_nx] = out.choice.l
            pol2[state_nx] = out.choice.switching
            VF[state_nx] = out.val
        end
    end
    return (pol1 = pol1, pol2 = pol2)
end
