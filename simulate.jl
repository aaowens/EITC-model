"Transiton rule for the state index, as a function of the policy function"
function nextstate(state_nx, pols)
    reborn = false
    x = state_nx
    age_nx, man_nx, eps_nx, ci_nx, kid_nx = x[1], x[2], x[3], x[4], x[5]
    age_nx = age_nx + 1
    # If they die, they are reborn as not a manager
    if age_nx == length(axes(pols.pol1)[1]) + 1
        age_nx = 1
        man_nx = 1
        reborn = true
    end
    switching = pols.pol2[state_nx]
    if man_nx == 2
        man_nx = 2
    else
        man_nx = switching + 1 # 1 if false, 2 if true
    end
    (state_nx2 = CartesianIndex(age_nx, man_nx, eps_nx, ci_nx, kid_nx), reborn = reborn)
end

```
Given a state index, a person ID, the time, policy functions, and model
Make a named tuple of data
1. Person ID
2. The time
3. Hours worked
4. Their wage
5. Their EITC usage
6. Earned income
7. Manager status
8. Occ change
9. Age
```
function makedata_ind(state_nx, i, t, pols, model)
    @unpack params, states = model
    state = makestate(model, state_nx)
    l = pols.pol1[state_nx]
    switching = pols.pol2[state_nx]
    choice = Choice(l, switching)
    c = consumption(choice, state, params)
    w = wagefun(state, params)
    y = w*l
    eitc = tau(y, state, params)
    ismanager = (state_nx[2] == 2)
    age = state_nx[1]
    haskid = state_nx[5]
    data = (i = i, t = t, hours = l, wage = w, eitc = eitc, y = y, ismanager = ismanager,
    switching = switching, age = age, haskid = haskid)
end



```
Produces a data frame of model simulated data given an initial distribution,
    policy functions, and the model.
Optionally provide post, a second model and policy functions for that model
If post is provided, the economy will go through a transition to the new steady state,
    starting halfway through the simulation.

Algorithm: Start with initial distribution of people
Push each person through T periods, keeping track of the person ID
If the person hits retirement, remove them, replace with
a young person, not manager, same other state
A person is a state, an ID, and a time
Init some arbitrary IDs at time 1
When the person retires, reset their state, and give them
    a new ID
When someone is reborn, they become a random draw from the initial dist,
and go through life again
There is a burn-in period so that this rebirth cycle reaches a steady state
```
function makedata(dist, pols, model; post = nothing)
    burn_T = 500
    realT = 16
    T = burn_T + realT
    ID_iterator = 1:(length(dist)*T)
    ID, state = iterate(ID_iterator)
    datas = [makedata_ind(dist[1], ID, 1, pols, model)] # get the type right
    pop!(datas) # empty the array
    sizehint!(datas, length(dist)*T)
    randstate_nx = rand(dist)
    for state_nx in dist
        ID, state = iterate(ID_iterator, state)
        person_data = makedata_ind(state_nx, ID, 1, pols, model)
        push!(datas, person_data)
        for t = 2:T
            currentpols = pols
            currentmodel = model
            if post !== nothing
                if t >= burn_T + realT ÷ 2
                    currentpols = post.pols
                    currentmodel = post.model
                end
            end
            state_nx, reborn = nextstate(state_nx, currentpols)
            if reborn
                ID, state = iterate(ID_iterator, state)
                x = randstate_nx
                # reborn as young and not a manager
                state_nx = CartesianIndex(1, 1, x[3], x[4], x[5])
                randstate_nx = rand(dist)
            end
            person_data = makedata_ind(state_nx, ID, t, currentpols, currentmodel)
            push!(datas, person_data)
        end
    end
    df = DataFrame(datas)
    @where(df, :t .> burn_T)
end
