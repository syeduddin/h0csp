import numpy as np
import emcee
import matplotlib.pyplot as plt
import corner

# Constants
c = 3e5  # Speed of light in km/s

# Effective Hubble function in timescape cosmology
def H_eff(z, H0, f_v0):
    Omega_m_eff = (1 - f_v0)  # Effective matter density
    Omega_Q = 1 - Omega_m_eff  # No Lambda in timescape
    epsilon = 0.6  # Typical value for timescape model
    return H0 * np.sqrt(Omega_m_eff * (1+z)**3 + Omega_Q * (1+z)**(2+epsilon))

# Compute comoving distance with numerical summation
def comoving_distance(z, H0, f_v0):
    z_vals = np.linspace(0, z, 1000)
    dz = z_vals[1] - z_vals[0]
    integrand = 1 / H_eff(z_vals, H0, f_v0)
    return np.sum(integrand) * dz * c  # Trapezoidal approximation

# Luminosity distance and distance modulus
def distance_modulus(z, H0, f_v0):
    dL = (1 + z) * comoving_distance(z, H0, f_v0)
    return 5 * np.log10(dL) + 25

# Define the likelihood function for supernova fitting
def log_likelihood(theta, z_obs, mu_obs, mu_err):
    H0, f_v0, M_B, alpha, beta = theta
    mu_model = np.array([distance_modulus(z, H0, f_v0) + M_B for z in z_obs])
    residuals = (mu_obs - mu_model) / mu_err
    return -0.5 * np.sum(residuals**2)

# Define the prior function
def log_prior(theta):
    H0, f_v0, M_B, alpha, beta = theta
    if 50 < H0 < 80 and 0 < f_v0 < 1 and -20 < M_B < -18 and 0 < alpha < 1 and 1 < beta < 5:
        return 0.0  # Flat priors
    return -np.inf  # Log(0) for values outside range

# Define the posterior probability
def log_probability(theta, z_obs, mu_obs, mu_err):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, z_obs, mu_obs, mu_err)

# Simulated Supernova Data (Replace with real data if available)
z_obs = np.linspace(0.01, 1.5, 40)  # Observed redshifts
true_params = [70, 0.8, -19.3, 0.14, 3.1]  # H0, f_v0, M_B, alpha, beta
mu_obs = np.array([distance_modulus(z, *true_params[:2]) + true_params[2] for z in z_obs])
mu_err = 0.2 * np.ones_like(mu_obs)  # Simulated error

# Initialize MCMC
ndim, nwalkers, steps = 5, 50, 2000
pos = [true_params + 1e-4 * np.random.randn(ndim) for _ in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(z_obs, mu_obs, mu_err))

# Run MCMC
sampler.run_mcmc(pos, steps, progress=True)

# Get the results
samples = sampler.get_chain(discard=500, thin=10, flat=True)
labels = ["H0", "f_v0", "M_B", "alpha", "beta"]

# Plot the corner plot
fig = corner.corner(samples, labels=labels, truths=true_params)
plt.show()

# Print estimated parameters
medians = np.median(samples, axis=0)
print("Estimated Parameters:")
for name, val in zip(labels, medians):
    print(f"{name}: {val:.3f}")
