#include <math.h> 
#include <TMB.hpp> 

template<class Type>
struct integrand_f {
    int num_occasions; 
    vector<Type> X_i;  
    vector<Type> betas; 
    Type sigma_individual_effects;
    
    Type operator() (Type x) {
        // eta_fixed includes the intercept (beta_0) and covariates
        Type eta_fixed = (X_i * betas).sum();
        
        // In M_h, p is constant over j for a given individual + random effect
        Type eta = eta_fixed + sigma_individual_effects * x;
        Type p = invlogit(eta);
        
        // Probability of being caught at least once: 1 - (1-p)^tau
        Type p_at_least_once = 1.0 - exp(num_occasions * log(1.0 - p));
        
        // Combine standard normal density and capture probability
        return dnorm(x, Type(0.0), Type(1.0), false) * p_at_least_once;
    }
};
    
template<class Type>
vector<Type> marginalized_pis(vector<Type> input) { 
    int num_occasions_int = CppAD::Integer(input[0]);
    int num_X_int = CppAD::Integer(input[1]);
    
    int offset = 2;
    vector<Type> X_i = input.segment(offset, num_X_int);
    offset += num_X_int;
    vector<Type> betas = input.segment(offset, num_X_int);
    offset += num_X_int;
    Type sigma_individual_effects = input[offset];
    
    integrand_f<Type> f = {num_occasions_int, X_i, betas, sigma_individual_effects};
    
    // Romberg integration limits for N(0,1)
    Type a = -5.0;
    Type b = 5.0;
    
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
    DATA_MATRIX(X); // X must contain a column of 1s for the intercept
    int num_individuals = y.rows();
    int num_occasions = y.cols();
    int num_X = X.cols();
    
    // Parameters
    PARAMETER_VECTOR(betas);
    PARAMETER_VECTOR(individual_effects); // Latent variables (random = "individual_effects" in R)
    PARAMETER(logsigma_individual_effects);
    Type sigma_individual_effects = exp(logsigma_individual_effects);
    
    Type nll = 0.0;

    // 1. Random Effects Penalty
    // TMB handles the N(0, sigma) integration, but we use non-centered 
    // parametrization for stability.
    for(int i=0; i<num_individuals; i++){
        nll -= dnorm(individual_effects(i), Type(0.0), Type(0.5), true);
    }

    // 2. Conditional Data Likelihood (Truncated Binomial)
    vector<Type> eta_fixed = X * betas;
    for (int i=0; i<num_individuals; i++) {
        // Scaling the random effect by sigma here
        Type eta = eta_fixed(i) + sigma_individual_effects * individual_effects(i);
        Type p = invlogit(eta);
        
        // Sum captures across row for binomial efficiency
        Type y_i_sum = y.row(i).sum();
        nll -= dbinom(y_i_sum, Type(num_occasions), p, true);
    }
    
    // 3. Zero-Truncation Adjustment (Huggins Denominator)
    vector<Type> probability_captured_once(num_individuals);
    
    for (int i=0; i<num_individuals; i++) {
        vector<Type> X_row = X.row(i);
        // Simplified input: metadata + X_row + betas + sigma
        vector<Type> input(2 + num_X + num_X + 1); 
        input << Type(num_occasions), Type(num_X), X_row, betas, sigma_individual_effects;  
        
        Type p_cap = marginalized_pis(input)[0];
        
        if(p_cap < Type(1e-12)) p_cap = Type(1e-12);
        
        probability_captured_once(i) = p_cap;
        
        // Add log(pi_i) to NLL for the Huggins correction
        nll += log(p_cap);
    }
    
    // Abundance Estimation
    Type population_size = (1.0 / probability_captured_once).sum();
    
    ADREPORT(sigma_individual_effects);
    ADREPORT(betas);
    ADREPORT(probability_captured_once);
    ADREPORT(population_size);
    
    return nll; 
}
