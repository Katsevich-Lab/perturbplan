// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include "var_nb.h"
#include <cmath>
using namespace Rcpp;

//' Compute Test Statistic Distribution for Random Effect Sizes
//' 
//' @description
//' Computes the asymptotic mean and standard deviation of the score test statistic
//' for random effect sizes in perturb-seq experiments. This function handles the
//' case where fold changes vary randomly across perturbations, using the average
//' fold change and its second moment.
//' 
//' @param num_trt_cell Numeric. Number of treatment cells
//' @param num_cntrl_cell Numeric. Number of control cells  
//' @param expression_mean Numeric. Mean baseline expression level
//' @param expression_size Numeric. Size parameter for negative binomial distribution
//' @param avg_fold_change Numeric. Average fold change across perturbations
//' @param avg_fold_change_sq Numeric. Average of squared fold changes (second moment)
//' 
//' @return A list containing:
//' \describe{
//'   \item{mean}{Numeric. Asymptotic mean of the test statistic}
//'   \item{sd}{Numeric. Asymptotic standard deviation of the test statistic}
//' }
//' 
//' @details
//' This function computes the asymptotic distribution of the score test statistic
//' under the assumption of random effect sizes. The key difference from fixed
//' effect sizes is that it accounts for variability in fold changes using the
//' second moment (avg_fold_change_sq).
//' 
//' The computation follows these steps:
//' \enumerate{
//'   \item Compute treatment/control cell proportions
//'   \item Calculate treatment, control, and pooled expression means
//'   \item Compute pooled variance using negative binomial variance formula
//'   \item Calculate denominator squared for the test statistic
//'   \item Compute control group variance
//'   \item Compute treatment group variance (incorporating fold change variability)
//'   \item Calculate final asymptotic mean and standard deviation
//' }
//' 
//' @seealso \code{\link{compute_distribution_teststat_fixed_es_cpp}} for fixed effect sizes
//' 
//' @export
// [[Rcpp::export]]
List compute_distribution_teststat_random_es_cpp(
    double num_trt_cell,
    double num_cntrl_cell,
    double expression_mean,
    double expression_size,
    double avg_fold_change,
    double avg_fold_change_sq
) {
    // Input validation
    if (num_trt_cell <= 0 || num_cntrl_cell <= 0) {
        stop("Cell counts must be positive");
    }
    if (expression_mean <= 0 || expression_size <= 0) {
        stop("Expression parameters must be positive");
    }
    if (avg_fold_change <= 0 || avg_fold_change_sq <= 0) {
        stop("Fold change parameters must be positive");
    }
    if (avg_fold_change_sq < avg_fold_change * avg_fold_change) {
        stop("avg_fold_change_sq must be >= avg_fold_change^2 (second moment constraint)");
    }
    
    // Compute treatment/control cells proportion
    const double num_test_cell = num_trt_cell + num_cntrl_cell;
    const double trt_test_prop = num_trt_cell / num_test_cell;
    const double cntrl_test_prop = 1.0 - trt_test_prop;
    
    // Define treatment/control/pooled mean expression
    const double trt_expression_mean = expression_mean * avg_fold_change;
    const double cntrl_expression_mean = expression_mean;
    const double pooled_expression_mean = trt_expression_mean * trt_test_prop + 
                                         cntrl_expression_mean * cntrl_test_prop;
    
    // Compute the square of the denominator in the score statistic
    const double pooled_var = var_nb(pooled_expression_mean, expression_size);
    const double denominator_sq = pooled_var * (1.0 / num_cntrl_cell + 1.0 / num_trt_cell);
    
    // Guard against division by zero
    if (denominator_sq <= 0) {
        stop("Denominator squared must be positive");
    }
    
    // Compute the asymptotic sd of test stat
    // Control group variance
    const double cntrl_var = var_nb(cntrl_expression_mean, expression_size) / num_cntrl_cell;
    
    // Treatment group variance (incorporating fold change variability)
    const double trt_var = (trt_expression_mean + 
                           (cntrl_expression_mean * cntrl_expression_mean * avg_fold_change_sq / expression_size)) / 
                           num_trt_cell;
    
    // Compute the asymptotic sd
    const double sd = std::sqrt((cntrl_var + trt_var) / denominator_sq);
    
    // Compute the asymptotic mean of test stat
    const double mean = cntrl_expression_mean * (avg_fold_change - 1.0) / std::sqrt(denominator_sq);
    
    // Guard against invalid results
    if (!std::isfinite(mean) || !std::isfinite(sd)) {
        stop("Non-finite result in test statistic computation");
    }
    if (sd <= 0) {
        stop("Standard deviation must be positive");
    }
    
    // Return the mean and sd as a named list
    return List::create(
        Named("mean") = mean,
        Named("sd") = sd
    );
}