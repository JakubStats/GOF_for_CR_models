//-----------------------------
#include <math.h> 
#include <TMB.hpp> 

template<class Type>
struct integrand_f {
    int num_occasions; 
    vector<Type> time_effects; 
    Type sigma_individual_effects;
    
    // Evaluate integrand
    Type operator() (Type x) {
        // x is the standard normal variable N(0,1)
        vector<Type> one_minus_p(num_occasions);
        
        for (int j=0; j<num_occasions; j++) {
            // Linear predictor: Time effect + Individual Random Effect (scaled by sigma)
            Type eta = time_effects(j) + sigma_individual_effects * x;
            one_minus_p(j) = 1.0 - invlogit(eta);
        }
        
        Type p_at_least_once = 1.0 - one_minus_p.prod();
        
        // Return [Normal Density] * [Prob(Captured | x)]
        Type ans = dnorm(x, Type(0.0), Type(1.0), false) * p_at_least_once;
        return ans;
    }
};
    
template<class Type>
vector<Type> marginalized_pis(vector<Type> input) { 
    // input order: num_occasions, time_effects, sigma_individual_effects
    int n_occ = CppAD::Integer(input[0]);
    
    // Clean slicing
    vector<Type> time_effects = input.segment(1, n_occ);
    Type sigma_individual_effects = input[1 + n_occ];
    
    integrand_f<Type> f = {n_occ, time_effects, sigma_individual_effects};
    
    // Romberg integration limits
    Type a = -7.0; 
    Type b = 7.0; 
    
    vector<Type> res(1);
    res[0] = romberg::integrate(f, a, b);
    return res;
}

REGISTER_ATOMIC(marginalized_pis)
    
template<class Type>
Type objective_function<Type>::operator() () { 
    using namespace density;

    // Data
    DATA_MATRIX(y); 
    int n_ind = y.rows();
    int n_occ = y.cols();

    // Parameters
    PARAMETER_VECTOR(time_effects);
    PARAMETER_VECTOR(individual_effects);
    PARAMETER(logsigma_individual_effects);
    Type sigma_individual_effects = exp(logsigma_individual_effects);
    
    Type nll = 0.0;

    // 1. Random Effects Likelihood: Individual intercepts ~ N(0, sigma)
    nll -= sum(dnorm(individual_effects, Type(0.0), sigma_individual_effects, true));

    // 2. Conditional Data Likelihood (Bernoulli)
    for (int i=0; i<n_ind; i++) {
        for (int j=0; j<n_occ; j++) {
            Type eta = time_effects(j) + individual_effects(i);
            nll -= dbinom(y(i,j), Type(1.0), invlogit(eta), true);
        }
    }
    
    // 3. Zero-Truncation Adjustment
    vector<Type> prob_cap(n_ind);
    for (int i=0; i<n_ind; i++) {
        // Prepare input vector for atomic function
        vector<Type> input(n_occ + 2); 
	input(0) = Type(n_occ);
	for(int j = 0; j < n_occ; j++) input(j + 1) = time_effects(j);
	input(n_occ + 1) = sigma_individual_effects;
          
        
        Type p_cap = marginalized_pis(input)[0];
  
        
        prob_cap(i) = p_cap;
        
        // Add log(P(Captured)) to NLL to normalize for truncation
        nll += log(p_cap + Type(1e-12));
    }
    
    // Horvitz-Thompson Population Size Estimate
    Type population_size = (1.0 / prob_cap).sum();
     
    // Reporting
    ADREPORT(individual_effects);
    ADREPORT(sigma_individual_effects);
    ADREPORT(time_effects);
    ADREPORT(prob_cap);
    ADREPORT(population_size);
    
    return nll; 
}