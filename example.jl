# Generation of example figures for the TestParticle.jl presentation at JuliaConCN 2022.
#
# Hongyang Zhou, hyzhou@umich.edu 11/30/2022

using TestParticle
using TestParticle: get_gc, getB_dipole, getE_dipole, sph2cart, Rₑ
using TestParticleMakie
using OrdinaryDiffEq
using StaticArrays
using GLMakie, TestParticleMakie

# Uniform static B example
begin
	# Magnetic field
	B(x) = SA[0, 0, 1e-8]
	# Electric field
	E(x) = SA[0, 0, 0]
	# Initial conditions
	x0 = [1.0, 0.0, 0.0]
	v0 = [0.0, 1.0, 0.1]
	stateinit = [x0..., v0...]
	tspan = (0, 20)
	# Assemble particle + fields
	param = prepare(E, B, species=Proton)
	prob = ODEProblem(trace!, stateinit, tspan, param)
	# Trace trajectory and save positions & velocities
	sol = solve(prob, Tsit5(); save_idxs=[1,2,3,4,5,6])
	plot(sol)
end

# Electron and proton example
let
	B(x) = SA[0, 0, 1e-11]
	E(x) = SA[0, 0, 5e-13]
	x0 = [0.0, 0.0, 0.0]
	u0 = [1.0, 0.0, 0.0]
	stateinit = [x0..., u0...]

	param_electron = prepare(E, B, species=Electron)
	tspan_electron = (0.0, 15.0)

	param_proton = prepare(E, B, species=Proton)
	tspan_proton = (0.0, 5.0)

	prob_e = ODEProblem(trace!, stateinit, tspan_electron, param_electron)
	prob_p = ODEProblem(trace!, stateinit, tspan_proton, param_proton)

	sol_e = solve(prob_e, Tsit5(); save_idxs=[1,2,3])
	sol_p = solve(prob_p, Tsit5(); save_idxs=[1,2,3])

	f = Figure()
	#gcd = f[1:2, 1] = GridLayout()
	#gd = gcd[2, 1] = GridLayout()
	ax = Axis3(f[1,1], aspect=:data, title="Particle trajectories")
    lines!(sol_e, linewidth=3, label="e")
	lines!(sol_p, linewidth=3, label="p")
	#Legend(f[1,2], ax, nbanks=1)
	#rowgap!(gd, 1)
	#colgap!(gd, 1)
	f
end

# Adiabatic example
let
	Ek = 5e7 # [eV]

	m = TestParticle.mᵢ
	q = TestParticle.qᵢ
	c = TestParticle.c
	Rₑ = TestParticle.Rₑ   
	invRE = inv(Rₑ)

	# initial velocity, [m/s]
	v₀ = sph2cart(c*sqrt(1-1/(1+Ek*q/(m*c^2))^2), 0.0, π/4)
	# initial position, [m]
	r₀ = sph2cart(2.5*Rₑ, 0.0, π/2)
	stateinit = [r₀..., v₀...]

	param = prepare(getE_dipole, getB_dipole)
	tspan = (0.0, 20.0)

	prob = ODEProblem(trace!, stateinit, tspan, param)

	sol = solve(prob, Tsit5(); save_idxs=[1,2,3])

	f = Figure()
	ax = Axis3(f[1,1], aspect=:data, title="Proton trajectory in a dipole field")
	l = plot!(sol, color=:blue3)
	scale!(l, invRE, invRE, invRE)

	for ϕ in range(0, stop=2*π, length=10)
        lines!(TestParticle.fieldline(ϕ)..., color=:tomato)
    end
	f
end