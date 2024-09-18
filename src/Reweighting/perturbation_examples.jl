#FUNCTIONS FOR PERTURBATIONS

#Perturbing Lennard-Jones interactions
function switching_func(r, switch, cut)
    switching = ((cut^2-r^2)^2*(cut^2-r^2-3*(switch^2-r^2)))/(cut^2-switch^2)^3
    return switching
end

function lennard_jones(r, σ, ϵ) 
    V = 4 * ϵ * ((σ/r)^12 - (σ/r)^6)
    return V
end

function lennard_jones_perturbation(r, σ1, ϵ1, σ2, ϵ2; α = 0.0, β = 0.0, smooth::Bool = false, switch = 8.0, cut = 12.0)
    ϵ_orig = sqrt(ϵ1*ϵ2)
    ϵ_pert = (1 + β) * sqrt(ϵ1*ϵ2)
    σ_orig = (σ1 + σ2)/2
    σ_pert = α * (σ1* + σ2)/2
    V_pert = (lennard_jones(r, σ_pert,  ϵ_pert) - lennard_jones(r,  σ_orig, ϵ_orig))
    if smooth && switch < r <= cut
        V_pert = switching_func(r, switch, cut) * V_pert
    else
        V_pert =  V_pert - (lennard_jones(cut, σ_pert,  ϵ_pert) - lennard_jones(cut,  σ_orig, ϵ_orig))
    end
    return V_pert
end

#Applying a polynomial decay pertubation
poly_decay_perturbation(r, α, cut) = r > cut ? zero(r) : α * ((r/cut)^2 - 1)^2

#Applying a half-gaussian pertubation
gaussian_decay_perturbation(r, alpha, beta) = beta > 0 ? alpha*exp(-beta*r^2) : error("you need to insert positive values for β")
