filePath = @__DIR__
(@isdefined(FlowSystem)) ? nothing : include(string(filePath,"/FlowSystem.jl"))
include("topographies.jl")
include("initial_conditions.jl")
using .FlowSystem
using DifferentialEquations
using Plots, LaTeXStrings, ColorSchemes
using CSV, DataFrames
using Dates
Plots.PyPlotBackend()

plot_font = "Computer Modern"
default(fontfamily=plot_font, grid=false, color = ColorSchemes.berlin[1])

function fluid_plot(sol, time_skip::Int64, save_every, x, s, tmax, xmax)
    dt = floor(Int64, time_skip/save_every)
    plt = plot(x, s, dpi=200, xlabel=L"Dimensionless length $x$", ylabel=L"Dimensionless free surface height $\phi$", legend=false, color="black",
            title=L"Fluid profile at $\delta t = %$time_skip$ intervals", ylims=(0, Inf), xlims=(0, xmax))
    for i = 1:floor(Int64, tmax/time_skip)+1
        plot!(plt, x, sol.u[1 + dt*(i-1)][2:end-1] + s)
    end
    plot!(plt, x, s, color="black")
    return plt
end

function fluid_anim(sol, ymax, xmax, x, s)
    anim = @animate for i = 1:size(sol.u)[1]
        timestep = sol.t[i]
        plot(x, sol.u[i][2:end-1]+s, legend=false, ylims=(0, ymax), xlims=(0, xmax), title=L"Fluid Profile $(t=%$timestep)$", 
        xlabel=L"Dimensionless length $x$", ylabel=L"Dimensionless free surface height $\phi$")
        plot!(x, s, color="black")
    end
    return anim
end

"""
Function to save fluid profile data to a CSV. CSV saved in the format of:
Timestamp | Ghost 1 | 0 | ... | Lx | Ghost 2
# Arguments
- `sol`: solution object of DiffEq solver
"""
function fluid_data(sol, path)
    CSV.write(path * "/data.csv", DataFrame(sol))
end

function save_params(path, filename="params")
    open(path*"/"*filename*".txt", "w") do f
        domain_info = "# Domain Parameters\nnx: $nx\nLx: $Lx\ndx: $dx\n"
        topo_info = "\n# Topography Parameters\n$topo\n"
        init_cond_info = "\n# Initial Condition Parameters\n$ic_obj\n"
        ode_sys_info = "\n# ODE System Parameters\nD: $D\nC: $C\nalpha: $alpha\n"
        time_info = "\n# Time Parameters\ntime_span: $tspan\ntimestep_save: $save_every"
        write(f, domain_info * topo_info * init_cond_info * ode_sys_info * time_info)
    end
end

function main(curr_dir)
    # Discretization of x-domain
    nx = 800
    Lx = 40
    dx = Lx/nx
    x = 0.0:dx:Lx

    topo = FlowSystem.Equation(flat, ())
    s = topo.f.(x, topo.params...)

    b = 0.01
    ic_center = 5
    ic_steep = 0.5
    ic_height = 1
    init_cond = FlowSystem.Equation(reflected_sigmoid, (b, ic_center, ic_steep, ic_height))
    ic = init_cond.f.(x, init_cond.params...)

    # Adding ghost points
    pushfirst!(ic, ic[1])
    push!(ic, ic[end])

    # Parameters of the system of ODEs
    α = 0
    g = 9.8 # m/s^2

    # Oil Physical Properties
    ρ = 900 # kg/m^3
    γ = 20 * 10^(-3) # kg/s^2
    μ = 50 * 10^(-6) # kg*m/s

    # Saw Forcing Properties
    A = 8 * 10^(-10) # m 
    ω = 40*pi*10^6 # 1/s
    α₁ = 2.386
    kᵢ = -0.7683 # 1/m 

    # Scaling Factors
    h_c = 200 * 10^(-6)
    x_c = ((γ*h_c)/(ρ*g))^(1/3)
    t_c = (3*μ*x_c)/(h_c^2*ρ*g)

    # Compact Dimless Params
    Ca = (h_c^2 * ρ * g)/(3 * γ)
    D = (3*Ca)^(1/3)
    C = ((1 + α₁^2)*ω^2*A^2)/(2*g*x_c)

    xsym = Symbol.(x)
    pushfirst!(xsym, Symbol("ghost1"))
    push!(xsym, Symbol("ghost2"))

    p = (D, α, dx, C, kᵢ*x_c, topo.f, topo.params)
    tspan = (0.0, 300.0); save_every = 0.5
    f = ODEFunction(system!, syms=xsym)
    prob = ODEProblem(f, ic, tspan, p, saveat=save_every)
    sol = solve(prob, alg=Rodas4())

    date = string(Dates.format(now(), "YYYYmmdd"))
    time = string(Dates.format(now(), "HHMMSS"))
    path = "$curr_dir/../runs/$date/$time"

    time_skip = 10
    plt = fluid_plot(sol, time_skip, save_every, x, s, tspan[2], Lx)
    anim = fluid_anim(sol, 5, Lx, x, s)

    mkpath(path)

    gif(anim, path * "/flow.gif", fps = 30)
    savefig(plt, path * "/plt.png")
    fluid_data(sol, path)
    #save_params(path)
end

main(filePath)