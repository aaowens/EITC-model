include("init.jl")

params = paramsgen()
model = makemodel(params)
sizeVF = ntuple(i -> length(model.states[i]), 5)
VF = zeros(sizeVF)
state_nxs = CartesianIndices(VF)

@time pols = solve_model(VF, model)

# Clone the model, then change the EITC parameters
model2 = deepcopy(model)
using Setfield
model2 = @set model2.params.eitcparams = makeeitc(rate1 = 0.4, rate2 = 0.4, max_credit = 0.4, flat_size = 0.05)
pols2 = solve_model(VF, model2)

# Look at the policy functions
using Plots
plot(model.IncomeProcess.state_values, pols.pol1[1, 1, :, 1, :])
plot(pols.pol2[1, 1, :, :, 1])
plot(VF[2, :, 5, 1, 1])


## Draw a random set of states to be the population for simulation
peopledist = state_nxs[:]
peopledist = rand(peopledist, 10000)

## Simulate model data and save it as a Stata dta file
# Provide the model and policy from the economy with the EITC change
# It will simulate the transition

out = makedata(peopledist, pols, model; post = (pols = pols2, model = model2))

using DataFrames, RCall
function write_dta(x, s)
    R"library(haven)"
    @rput x
    rstr = "write_dta(x, \"" * s * "\")"
    reval(rstr)
end
write_dta(out, "outdf_eitc_sim2.dta")
