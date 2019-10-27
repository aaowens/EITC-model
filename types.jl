
Base.@kwdef struct State
    age::Int = 20
    manager::Bool = false
    eps::Float64 = 0.
    ci::Float64 = 1.
    haskid::Bool = false
end

Base.@kwdef struct Choice
    l::Float64 = 0.5
    switching::Bool = false
end

struct Model{Tp, T, Ts}
    params::Tp
    IncomeProcess::T
    states::Ts
end
