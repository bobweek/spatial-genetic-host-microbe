include("structs.jl")

"""
a neutral microbial mutation event
- `i` = mutating individual
"""
function neutral_mmutation(M₁::mutations,M₂::mutations,i::Int64)
    @unpack_mutations M₁

    # make new locus vector for other mutation table M₂
    νᵥ = sparse(fill(false,size(M₂.ν)[2])')
    M₂.ν = sparse_vcat(M₂.ν,νᵥ)
    dropzeros!(M₂.ν)

    # make locus vector with new mutation
    νᵥ = sparse(fill(false,size(ν)[2])')
    νᵥ[i] = true

    # append to neutral mutation table
    ν = sparse_vcat(ν,νᵥ)
    dropzeros!(ν)

    @pack_mutations! M₁
    return M₁, M₂
end

"""
a causal microbial mutation event
- `i` = mutating individual
"""
function causal_mmutation(M₁::mutations,M₂::mutations,i::Int64)
    @unpack_mutations M₁

    # make new locus vector for other mutation table M₂
    κᵥ = sparse(fill(false,size(M₂.κ)[2])')
    M₂.κ = sparse_vcat(M₂.κ,κᵥ)
    dropzeros!(M₂.κ)

    # make locus vector with new mutation
    κᵥ = sparse(fill(false,size(κ)[2])')
    κᵥ[i] = true

    # append to causal mutation table
    κ = sparse_vcat(κ,κᵥ)
    dropzeros!(κ)

    # append new effect size
    αₚ = rand(Normal(0,σ))
    α = vcat(α,αₚ)
    M₂.α = vcat(M₂.α,αₚ)

    @pack_mutations! M₁
    return M₁, M₂
end

"""
draw neutral microbe mutations
"""
function neutral_mmutations(M₁::mutations,M₂::mutations)
    
    n = size(M₁.ν)[2]
    𝓜 = rand(Binomial(n,M₁.μⁿ))
    who = sample(1:n,𝓜,replace=false)
    for i in who
        M₁, M₂ = neutral_mmutation(M₁,M₂,i)        
    end

    return M₁, M₂

end

"""
draw causal microbe mutations
"""
function causal_mmutations(M₁::mutations,M₂::mutations)

    n = size(M₁.κ)[2]
    𝓜 = rand(Binomial(n,M₁.μᶜ))
    who = sample(1:n,𝓜,replace=false)
    for i in who
        M₁, M₂ = causal_mmutation(M₁,M₂,i)
    end

    return M₁, M₂

end

"""
draw microbial mutations
"""
function mmutate(M₁::mutations,M₂::mutations)

    M₁, M₂ = neutral_mmutations(M₁,M₂)
    M₁, M₂ = causal_mmutations(M₁,M₂)
   
    return M₁, M₂

end

"""
microbial mutations
"""
function mmutates(sys::system)
    @unpack_system sys

    Mₘ, Mₑ = mmutate(Mₘ,Mₑ)
    Mₑ, Mₘ = mmutate(Mₑ,Mₘ)

    @pack_system! sys
    return sys
end

"""
a neutral host mutation event
- `i` = mutating individual
"""
function neutral_hmutation(M::mutations,i::Int64)
    @unpack_mutations M

    # make locus vector with new mutation at i
    νᵥ = sparse(fill(false,size(ν)[2])')
    νᵥ[i] = true

    # append to neutral mutation table
    ν = sparse_vcat(ν,νᵥ)
    dropzeros!(ν)

    @pack_mutations! M
    return M
end

"""
a causal host mutation event
- `i` = mutating individual
"""
function causal_hmutation(M::mutations,i::Int64)
    @unpack_mutations M

    # make locus vector with new mutation at i
    κᵥ = sparse(fill(false,size(κ)[2])')
    κᵥ[i] = true

    # append to causal mutation table
    κ = sparse_vcat(κ,κᵥ)
    dropzeros!(κ)

    # append new effect size
    α = vcat(α,rand(Normal(0,σ)))

    @pack_mutations! M
    return M
end

"""
draw neutral host mutations
"""
function neutral_hmutations(M::mutations)
    
    n = size(M.ν)[2]
    𝓜 = rand(Binomial(n,M.μⁿ))
    who = sample(1:n,𝓜,replace=false)
    for i in who
        neutral_hmutation(M,i)        
    end

    return M

end

"""
draw causal host mutations
"""
function causal_hmutations(M::mutations)

    n = size(M.κ)[2]
    𝓜 = rand(Binomial(n,M.μᶜ))
    who = sample(1:n,𝓜,replace=false)
    for i in who
        causal_hmutation(M,i)
    end

    return M

end

"""
draw host mutations
"""
function hmutate(sys::system)
    @unpack_system sys

    neutral_hmutations(Mₕ)
    causal_hmutations(Mₕ)
   
    @pack_system! sys
    return sys
end

"""
a microbial dispersal event
- `i` = dispersing microbe
"""
function mdisp(M::mutations,i::Int64,Jₑ::Int64)
    @unpack_mutations M
    
    # location index of focal individual
    k = ceil(Int64,i/Jₑ)

    # pick rnd individual in different location
    before = 1:(Jₑ*(k-1))
    after  = (Jₑ*k+1):size(ν)[2]
    j = sample(vcat(before,after))

    # copy focal individuals genome onto rnd ind
    ν[:,j] = ν[:,i]
    κ[:,j] = κ[:,i]

    dropzeros!(ν)
    dropzeros!(κ)
    
    @pack_mutations! M
    return M
end

"""
draw microbial dispersal events
- only environmental microbes
"""
function mdisps(sys::system,π::parameters)
    @unpack_system sys
    @unpack_parameters π

    n = K*Jₑ
    𝓓 = rand(Binomial(n,dₘ))
    who = sample(1:n,𝓓,replace=false)
    for i in who
        mdisp(Mₑ,i,Jₑ)
    end

    @pack_system! sys
    return sys
end

"""
a host dispersal event
- `i` = dispersing host
"""
function hdisp(Mₕ::mutations,Mₘ::mutations,i::Int64,Nₕ::Int64,Jₕ::Int64)
    
    # location index of focal host
    k = ceil(Int64,i/Nₕ)

    # pick rnd host in different location
    before = 1:(Nₕ*(k-1))
    after  = (Nₕ*k+1):size(Mₕ.ν)[2]
    j = sample(vcat(before,after))

    # copy focal hosts genome onto rnd host
    Mₕ.ν[:,j] = Mₕ.ν[:,i]
    Mₕ.κ[:,j] = Mₕ.κ[:,i]
    dropzeros!(Mₕ.ν)
    dropzeros!(Mₕ.κ)
    
    # get focal host microbiome indices
    f = (k-1)*Nₕ + 1 # first host in focal host pop
    𝔦 = (k-1)*Nₕ*Jₕ .+ (i-f)*Jₕ .+ (1:Jₕ)

    # got rnd host microbiome indices
    k = ceil(Int64,j/Nₕ) # population of rnd host
    f = (k-1)*Nₕ + 1     # first host in rnd host pop
    𝔧 = (k-1)*Nₕ*Jₕ .+ (j-f)*Jₕ .+ (1:Jₕ)

    # copy focal host mbiome onto rnd host mbiome
    Mₘ.ν[:,𝔧] = Mₘ.ν[:,𝔦]
    Mₘ.κ[:,𝔧] = Mₘ.κ[:,𝔦]
    dropzeros!(Mₘ.ν)
    dropzeros!(Mₘ.κ)
    
    return Mₕ, Mₘ
end

"""
draw host dispersal events
"""
function hdisps(sys::system,π::parameters)
    @unpack_system sys
    @unpack_parameters π

    n = K*Nₕ
    𝓓 = rand(Binomial(n,dₕ))
    who = sample(1:n,𝓓,replace=false)
    for i in who
        Mₕ, Mₘ = hdisp(Mₕ,Mₘ,i,Nₕ,Jₕ)
    end

    @pack_system! sys
    return sys
end

"""
a shedding event
- `i` = individual being shed
"""
function shed(sys::system,π::parameters,i::Int64)
    @unpack_system sys
    @unpack_parameters π

    # compute location
    k = ceil(Int64,i/(Nₕ*Jₕ))

    # draw microbe from local environment
    j = sample(1:Jₑ) + (k-1)*Jₑ
    # this should be benchmarked against rand(DiscreteUniform(((k-1)*Jₑ+1),(k*Jₑ)))

    # replace environmental microbe mutations with shedded microbe
    Mₑ.ν[:,j] = Mₘ.ν[:,i]
    Mₑ.κ[:,j] = Mₘ.κ[:,i]
    dropzeros!(Mₑ.ν)
    dropzeros!(Mₑ.κ)

    @pack_system! sys
    return sys
end

"""
draw shedding events
"""
function sheds(sys::system,π::parameters)
    @unpack_parameters π

    n = K*Nₕ*Jₕ
    𝓢 = rand(Binomial(n,ψ))
    who = sample(1:n,𝓢,replace=false)
    for i in who
        shed(sys,π,i)
    end

    return sys
end

"""
a environmental acquisition event
"""
function acquire(sys::system,π::parameters,i::Int64)
    @unpack_system sys
    @unpack_parameters π

    # compute location
    k = ceil(Int64,i/(Nₕ*Jₕ))

    # draw microbe from local environment
    j = sample(((k-1)*Jₑ+1):(k*Jₑ)) # this should be benchmarked against rand(DiscreteUniform(((k-1)*Jₑ+1),(k*Jₑ)))

    # replace host microbe mutations with acquired microbe
    Mₘ.ν[:,i] = Mₑ.ν[:,j]
    Mₘ.κ[:,i] = Mₑ.κ[:,j]
    
    @pack_system! sys
    return sys
end

"""
draw environmental acquisition events
"""
function acquires(sys::system,π::parameters)
    @unpack_parameters π

    n = K*Nₕ*Jₕ
    𝓢 = rand(Binomial(n,ψ))
    who = sample(1:n,𝓢,replace=false)
    for i in who
        acquire(sys,π,i)
    end

    return sys
end

"""
a social transmission event
- `i` = individual chosen to be transmitted
"""
function social(Mₘ::mutations,π::parameters,i::Int64)
    @unpack_mutations Mₘ
    @unpack_parameters π

    # compute location
    k = ceil(Int64,i/(Nₕ*Jₕ))

    # first crobe in pop k
    fₘ = (k-1)*Nₕ*Jₕ + 1

    # compute local host index
    h = ceil(Int64,(i-fₘ)/Jₕ)

    # compute index of first microbe in focal host
    ϕ = (k-1)*Nₕ*Jₕ + (h-1)*Jₕ

    # microbes in focal host
    hᵩ = ϕ:(ϕ+Jₕ)

    # draw microbe from any local host other than focal one
    j = sample(setdiff(fₘ:(Nₕ*Jₕ+fₘ-1),hᵩ))

    # transmission
    ν[:,j] = ν[:,i]
    κ[:,j] = κ[:,i]

    @pack_mutations! Mₘ
    return Mₘ
end

"""
draw social transmission events
"""
function socials(sys::system,π::parameters)
    @unpack_system sys
    @unpack_parameters π

    n = K*Nₕ*Jₕ
    𝓣 = rand(Binomial(n,t))
    who = sample(1:n,𝓣,replace=false)
    for i in who
        social(Mₘ,π,i)
    end

    @pack_system! sys
    return sys
end

"""
a distal propagule event (env microbes only!)
- `i` = individual replaced by propagule
"""
function dstl_ppgl(Mₑ::mutations,Mₘ::mutations,i::Int64)
    @unpack_mutations Mₑ
    # a distant microbial propagule lands in the environment
    # and displaces microbe individual i

    # individual i's mutations are wiped clean
    ν[:,i] .= false
    κ[:,i] .= false
    dropzeros!(ν)
    dropzeros!(κ)

    # new propagule mutations are drawn
    𝓜  = rand(Poisson(D))
    𝓜ᶜ = rand(Binomial(𝓜,χ))
    𝓜ⁿ = 𝓜 - 𝓜ᶜ

    # new loci appended to host-associated mutation tables
    n = size(Mₘ.ν)[2]
    νₚ = sparse(fill(false,(𝓜ⁿ,n)))
    κₚ = sparse(fill(false,(𝓜ᶜ,n)))
    Mₘ.ν = sparse_vcat(Mₘ.ν,νₚ)
    Mₘ.κ = sparse_vcat(Mₘ.κ,κₚ)

    # new mutations appended to env mutation tables
    n = size(ν)[2]
    νₚ = sparse(fill(false,(𝓜ⁿ,n)))
    κₚ = sparse(fill(false,(𝓜ᶜ,n)))
    νₚ[:,i] .= true
    κₚ[:,i] .= true
    ν = sparse_vcat(ν,νₚ)
    κ = sparse_vcat(κ,κₚ)
    
    # new causal mutation effects are appended
    αₚ = rand(Normal(0,σ),𝓜ᶜ)
    α = vcat(α,αₚ)
    Mₘ.α = vcat(Mₘ.α,αₚ)

    @pack_mutations! Mₑ
    return Mₑ, Mₘ
end

"""
draw distal propagule events
"""
function dstl_ppgls(sys::system,π::parameters)
    @unpack_system sys
    @unpack_parameters π

    # draw propagules
    n = K*Jₑ
    𝓟 = rand(Binomial(n,p))
    who = sample(1:n,𝓟,replace=false)
    for i in who
        Mₑ, Mₘ = dstl_ppgl(Mₑ,Mₘ,i)        
    end

    @pack_system! sys
    return sys
end

"""
microbiome drift
- `M` = mutation table containing drifting microbiome
- `m` = indices of microbes in drifting microbiome
"""
function mdrift(M::mutations,m::UnitRange{Int64})
    @unpack_mutations M

    # number of microbes
    J = length(m)

    # draw parental indices
    𝒫 = sample(m,J)

    # inheritance of neutral mutations
    𝔐 = SparseMatrixCSC{Bool,Int64}(undef,size(ν)[1],J)
    for i in 1:J
        𝔐[:,i] = ν[:,𝒫[i]]
    end
    ν[:,m] = 𝔐
    dropzeros!(ν)
    
    # inheritance of causal mutations
    𝔐 = SparseMatrixCSC{Bool,Int64}(undef,size(κ)[1],J)
    for i in 1:J
        𝔐[:,i] = κ[:,𝒫[i]]
    end
    κ[:,m] = 𝔐
    dropzeros!(κ)

    @pack_mutations! M
    return M
end

"""
independent drift for each microbiome
"""
function mdrifts(sys::system,π::parameters)
    @unpack_system sys
    @unpack_parameters π

    for k in 1:K

        for h in 1:Nₕ

            # indices of local host microbes
            m = (1:Jₕ) .+ (k-1)*Nₕ*Jₕ .+ (h-1)*Jₕ

            # host-associated microbiome drift
            mdrift(Mₘ,m)

        end

        # indices of local environmental microbes
        m = (1:Jₑ) .+ (k-1)*Jₑ
        
        # environmental microbiome drift
        mdrift(Mₑ,m)

    end
    
    @pack_system! sys
    return sys
end

"""
one microbial generation
"""
function mgen(sys::system,π::parameters)

    ×
    # microbiome drift
    mdrifts(sys,π)
    
    # host-associated microbe mutations
    mmutates(sys)

    # distal propagules land in the environment
    dstl_ppgls(sys,π)

    # microbial transmission via social contact
    socials(sys,π)

    # environmental microbial dispersal
    mdisps(sys,π)

    # host shedding of microbes into the environment
    sheds(sys,π)

    # host acquisition of microbes from the environment
    acquires(sys,π)

    # update microbe generation number
    sys.τₘ += 1

    return sys
end

"""
initialize system and parameters
"""
function sys_init(
    K::Int64 = 3,
    J::Int64 = 10,
    N::Int64 = 10,
    d::Float64 = 0.2,
    ψ::Float64 = 0.1,
    ε::Float64 = 0.1,
    t::Float64 = 0.1,
    p::Float64 = 0.1,
    μ::Float64 = 0.2,
    σ::Float64 = 0.01,
    χ::Float64 = 0.5,
    Γₘ::Int64 = 2,
    Γₕ::Int64 = 10,
    τ::Int64 = 5)

    π = parameters(
        K=K, 
        Jₑ=J, 
        Jₕ=J, 
        Nₕ=N, 
        dₕ=d, 
        dₘ=d, 
        ψ=ψ, 
        ε=ε, 
        t=t, 
        p=p,
        Γₘ=Γₘ,
        Γₕ=Γₕ,
        τ=τ)

    @unpack_parameters π

    Mₕ = mutations()
    Mₘ = mutations()
    Mₑ = mutations()    

    Mₕ.μⁿ = μ
    Mₕ.μᶜ = μ
    Mₘ.μⁿ = μ
    Mₘ.μᶜ = μ
    Mₑ.μⁿ = μ
    Mₑ.μᶜ = μ
    
    Mₕ.σ = σ
    Mₘ.σ = σ
    Mₑ.σ = σ

    Mₕ.χ = χ
    Mₘ.χ = χ
    Mₑ.χ = χ
    
    # microbe mutations

    Mₘ.ν = sparse(fill(false,(1,K*N*J)))
    Mₘ.κ = sparse(fill(false,(1,K*N*J)))
    Mₘ.α = zeros(1)

    Mₑ.ν = sparse(fill(false,(1,K*J)))
    Mₑ.κ = sparse(fill(false,(1,K*J)))
    Mₑ.α = zeros(1)

    Mₘ, Mₑ = mmutate(Mₘ,Mₑ)
    Mₑ, Mₘ = mmutate(Mₑ,Mₘ)

    # host mutations

    Mₕ.ν = sparse(fill(false,(1,K*N)))
    Mₕ.κ = sparse(fill(false,(1,K*N)))
    Mₕ.α = zeros(1)

    neutral_hmutations(Mₕ)
    causal_hmutations(Mₕ)
    
    sys = system(Mₕ=Mₕ,Mₘ=Mₘ,Mₑ=Mₑ)

    mpurge_check(sys)
    hpurge_check(sys)

    return sys, π
end

"""
checks for purged microbial mutations
"""
function mpurge_check(sys::system)
    @unpack_system sys

    # do this for microbes every host generation
    # and at end of sim

    #
    # neutral loci
    #

    # get indices for neutral loci
    𝓛 = eachindex(axes(Mₘ.ν,1))

    # will contain empty loci
    remove = Vector{Int64}(undef,0)

    # iterate over neutral loci
    for ℓ in 𝓛
        # if mutation has been purged
        if !any(Mₘ.ν[ℓ,:]) & !any(Mₑ.ν[ℓ,:])
            # then append locus to remove
            push!(remove,ℓ)
        end
    end

    # drop empty loci/columns from neutral mutation table
    Mₘ.ν = Mₘ.ν[setdiff(1:end,remove),:]
    Mₑ.ν = Mₑ.ν[setdiff(1:end,remove),:]
    dropzeros!(Mₘ.ν)
    dropzeros!(Mₑ.ν)

    #
    # causal loci
    #

    # get indices for causal loci
    𝓛 = eachindex(axes(Mₘ.κ,1))

    # will contain empty loci
    remove = Vector{Int64}(undef,0)

    # iterate over causal loci
    for ℓ in 𝓛
        # if mutation has been purged
        if !any(Mₘ.κ[ℓ,:]) & !any(Mₑ.κ[ℓ,:])
            # then append locus to remove
            push!(remove,ℓ)
        end
    end

    # drop empty loci/columns from neutral mutation table
    Mₘ.κ = Mₘ.κ[setdiff(1:end,remove),:]
    Mₑ.κ = Mₑ.κ[setdiff(1:end,remove),:]
    dropzeros!(Mₘ.κ)
    dropzeros!(Mₑ.κ)

    # remove associated entry from α
    deleteat!(Mₘ.α,remove)
    deleteat!(Mₑ.α,remove)

    @pack_system! sys
    return sys
end

"""
checks for purged host mutations
"""
function hpurge_check(sys::system)
    @unpack_system sys
    @unpack_mutations Mₕ

    # do this for hosts every ten generations
    # and at end of sim

    #
    # neutral loci
    #

    # get indices for neutral loci
    𝓛 = eachindex(axes(ν,1))

    # will contain empty loci
    remove = Vector{Int64}(undef,0)

    # iterate over neutral loci
    for ℓ in 𝓛
        # if mutation has been purged
        if !any(ν'[:,ℓ])
            # then append locus to remove
            push!(remove,ℓ)
        end
    end

    # drop empty loci/columns from neutral mutation table
    ν = ν'[:,setdiff(1:end,remove)]'

    #
    # causal loci
    #

    # get indices for causal loci
    𝓛 = eachindex(axes(κ,1))

    # will contain empty loci
    remove = Vector{Int64}(undef,0)

    # iterate over causal loci
    for ℓ in 𝓛
        # if mutation has been purged
        if !any(κ'[:,ℓ])
            # then append locus to remove
            push!(remove,ℓ)
        end
    end

    # drop empty loci/columns from neutral mutation table
    κ = κ'[:,setdiff(1:end,remove)]'

    # remove associated entry from α
    deleteat!(α,remove)

    @pack_mutations! Mₕ
    @pack_system! sys
    return sys
end

"""
checks for fixed mutations
"""
function fixation_check(M::mutations)
    @unpack_mutations M

    # only do this when
    #   p = 0 because we need to
    #       compare w propagule genomes
    # and/or at the end of the sim
    #   (although it may useful to keep
    #       for reconstructing ancestry)
    # based on these considerations
    #   this fct may never get implemented

    for c in eachcol(ν)
        if all(c)
            # delete column c
        end
    end

    for c in eachcol(κ)
        if all(c)
            # delete column c
            # remove associated entry from α
            # accumulate associated entry in accumulator
        elseif !any(c)
            # delete column c
            # remove associated entry from α
            # don't accumulate
        end
    end
    
    return M
end

"""
computes expected relative fitnesses of host individuals
"""
function fitness(sys::system,π::parameters)
    @unpack_system sys
    @unpack_parameters π

    # compute additive microbial values
    m = vec(sum(reshape((Mₘ.α' * Mₘ.κ)', Jₕ, K*Nₕ), dims=1))

    # compute additive genetic values
    g = (Mₕ.α' * Mₕ.κ)'

    # compute trait values
    z = g₀ .+ g .+ m

    # compute relative fitness matrix
    # an Nₕ×K matrix with columns that sum to 1
    w = Matrix{Float64}(undef,Nₕ,K)
    for k in 1:K
        Z = z[((k-1)*Nₕ+1):(k*Nₕ)]
        w[:,k] = normalize(exp.(-s[k] .* (θ[k] .- Z).^2 ./ 2),1)
    end

    return w
end

"""
local host selection and drift
"""
function hseldift(Mₕ::mutations,Mₘ::mutations,k::Int64,w::Matrix{Float64})
    
    # number of local hosts
    N = length(w[:,k])

    # index of first host in focal pop
    f = 1 + (k-1)*N

    # draw host parent realized fitness
    𝓦 = rand(Multinomial(N,w[:,k]))

    # get host parent indices
    𝒫 = findall(x->x>0,𝓦) .+ (f-1)
    
    # inheritance of neutral host mutations
    𝔬 = 1 # offspring index
    𝔐 = SparseMatrixCSC{Bool,Int64}(undef,size(Mₕ.ν)[1],N)
    for i in eachindex(𝒫)
        for j in 1:𝓦[i]
            𝔐[:,𝔬] = Mₕ.ν[:,𝒫[i]]
            𝔬 += 1
        end
    end
    Mₕ.ν[:,f:(f+N-1)] = 𝔐
    dropzeros!(Mₕ.ν)
    
    # inheritance of causal host mutations
    𝔬 = 1 # offspring index
    𝔐 = SparseMatrixCSC{Bool,Int64}(undef,size(Mₕ.κ)[1],N)
    for i in eachindex(𝒫)
        for j in 1:𝓦[i]
            𝔐[:,𝔬] = Mₕ.κ[:,𝒫[i]]
            𝔬 += 1
        end
    end
    Mₕ.κ[:,f:(f+N-1)] = 𝔐
    dropzeros!(Mₕ.κ)
    
    # do things for host mbiomes
    K = floor(Int64,size(Mₕ.ν)[2]/N)
    J = floor(Int64,size(Mₘ.ν)[2]/(K*N))
    

    # the first stuff below here goes into the loop over host parents
    # i is the (absolute) index of the host parent in, eg, Mₕ.ν

    # inheritance of netural microbe mutations
    𝔬 = 1 # (relative) host offspring index
    𝔐 = SparseMatrixCSC{Bool,Int64}(undef,size(Mₘ.ν)[1],N*J)
    for i in eachindex(𝒫)
        for j in 1:𝓦[i]

            # get (absolute) parent host microbiome indices
            𝔦 = (k-1)*N*J .+ (𝒫[i]-f)*J .+ (1:J)
            
            # get (relative) offspring host microbiome indices
            𝔪 = (1:J) .+ (𝔬-1)*J

            # copy host offspring mbiome into 𝔐
            𝔐[:,𝔪] = Mₘ.ν[:,𝔦]
            
            # next host offspring
            𝔬 += 1

        end
    end
    Mₘ.ν[:,f:(f+N*J-1)] = 𝔐
    dropzeros!(Mₘ.ν)

    # inheritnace of causal microbe mutations
    𝔬 = 1 # (relative) host offspring index
    𝔐 = SparseMatrixCSC{Bool,Int64}(undef,size(Mₘ.κ)[1],N*J)
    for i in eachindex(𝒫)
        for j in 1:𝓦[i]

            # get (absolute) parent host microbiome indices
            𝔦 = (k-1)*N*J .+ (𝒫[i]-f)*J .+ (1:J)
            
            # get (relative) offspring host microbiome indices
            𝔪 = (1:J) .+ (𝔬-1)*J

            # copy host offspring mbiome into 𝔐
            𝔐[:,𝔪] = Mₘ.κ[:,𝔦]
            
            # next host offspring
            𝔬 += 1

        end
    end
    Mₘ.κ[:,f:(f+N*J-1)] = 𝔐
    dropzeros!(Mₘ.κ)

end

"""
independent host drift and selection for each location
"""
function hseldifts(sys::system,π::parameters)
    @unpack_system sys
    @unpack_parameters π

    w = fitness(sys,π)

    for k in 1:K

        # do local host drift
        hseldift(Mₕ,Mₘ,k,w)

    end

    @pack_system! sys
    return sys
end

"""
performs one iteration of simulation model
- equals one host generation
- includes:
    - host dispersal
    - host selection
    - host drift
    - host mutation
    - host development
        - `Γₘ` microbe generations
"""
function hgen(sys::system,π::parameters)
    
    # host dispersal
    hdisps(sys,π)
    
    # host selection and drift
    hseldifts(sys,π)

    # host mutation
    hmutate(sys)

    # host development
    # want this to happen at the very end
    # so that this fct outputs data on adult hosts
    for i in 1:π.Γₘ
        sys = mgen(sys,π)
    end

    # check for purged mutations in microbe genomes
    mpurge_check(sys)

    # if its time to check hosts for purged mutations
    if mod(sys.τₕ,π.τ) == 0
        hpurge_check(sys)
    end

    # update host generation number
    sys.τₕ += 1

    # reset microbe generation number
    sys.τₘ = 0

    return sys
end

"""
run a simulation for `Γₕ` host generations
"""
function sim(sys::system,π::parameters)
    for i in 1:π.Γₕ
        hgen(sys,π)
    end
    return sys    
end
