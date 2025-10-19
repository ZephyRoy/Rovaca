#ifndef ROVACA_HC_SADDLE_POINT_EXPANSION_H_
#define ROVACA_HC_SADDLE_POINT_EXPANSION_H_

namespace rovaca
{

namespace SaddlePointExpansion
{

double get_stirling_error(double z);

double get_deviance_part(double x, double mu);

double log_binomial_probability(int x, int n, double p, double q);

}  // namespace SaddlePointExpansion

}  // namespace rovaca

#endif  // ROVACA_HC_SADDLE_POINT_EXPANSION_H_
