#include <vector>

namespace water_data_tables
{
	// tables for calculation of constants A, B and C for water properties
	const std::vector<double> tabm = {-0.00751, 0.013624, -0.0781, 0.0, 0.0};
	const std::vector<double> tabn = {0.01193, 0.0851, 1.02766, 0.0, 0.0};
	const std::vector<double> tabo = {1.8316, -7.8119, -3.6231, -0.10733, 1.09192};
	
	// tables for calculation of co2 properties
	const std::vector<double> tabca = {28.9447706, -0.0354581768, -4770.67077, 1.027827E-05,
										33.8126098, 9.04037140E-03, -1.14934031E-03, -0.307405726, 
										-0.0907301486, 9.32713393E-04, 0.0};
	const std::vector<double> tabcb = {-0.411370585, 6.07632013E-04, 97.5347708, 0.0, 0.0, 
										0.0, 0.0, -0.0237622469, 0.0170656236, 0.0, 1.41335834E-05};
	const std::vector<double> tabcc = {3.36389723E-04, -1.9829898E-05, 0.0, 0.0, 0.0, 0.0, 
										0.0, 2.1222083E-03, -5.24873303E-03, 0.0, 0.0};
										
	// tables for water EOS calculation
	const std::vector<double> tabdwt = { -0.127213, 0.645486, 1.03265, -0.070291, 0.639589};
	const std::vector<double> tabewt = {4.221, -3.478, 6.221, 0.5182, -0.4405};
	const std::vector<double> tabfwt = {-11.403, 29.932, 27.952, 0.20684, 0.3768};
	
	// tables for brine temperature coefficient calculation
	const std::vector<double> tdcm2 = {-7.925E-05, -1.93E-06, -3.4254E-04, 0.0, 0.0};
	const std::vector<double> tdcm32 = {1.0998E-03, -2.8755E-03, -3.5819E-03, -0.72877, 1.92016};
	const std::vector<double> tdcm1 = {-7.6402E-03, 3.6963E-02, 4.36083E-02, -0.333661,1.185685};
	const std::vector<double> tdcm12 = {3.746E-4, -3.328E-4, -3.346E-4, 0.0, 0.0};
	
	// tables for brine cp 
	const std::vector<double> tecm2 = {0.0, 0.0, 0.1353, 0.0, 0.0};
	const std::vector<double> tfcm32 = {-1.409, -0.361, -0.2532, 0, 9.216};
	const std::vector<double> tfcm1 = {0, 5.614, 4.6782, -0.307, 2.6069};
	const std::vector<double> tfcm12 = {-0.1127, 0.2047, -0.0452, 0.0, 0.0};
}
namespace air_data_tables
{
	const std::vector<double> tabpsv = {1.2378847E-05,-1.9121316E-2,33.93711047,-6.3431645E3};
	const std::vector<double> tabef = {1.00062,3.14E-08,5.6E-07};
	const std::vector<double> tabz = {1.58123E-06,-2.9331E-08,1.1043E-10,5.707E-06,-2.051E-08,1.9898E-04,-2.376E-06,1.83E-11,-0.765E-08};
	const std::vector<double> tabvisc = {84.986,7.0,113.157,-1,-3.7501E-03,-100.015};
}
namespace utility_functions
{
	double tempfunc(const std::vector<double>& vTab, const double temperature);
}
namespace fluid_properties
{
	void iapws95(const double pressure, const double temperature, const double salinity, double* density, double* compressibility);
	double kestin_viscosity(const double pressure, const double temperature, const double salinity);
	double air_density(const double pressure, const double temperature, const double humidity);
	double air_water_mol_fraction(const double pressure, const double temperature, const double humidity);
	double air_compressibility_factor(const double pressure, const double temperature, const double humidity);
	double air_viscosity(const double pressure, const double temperature, const double humidity);
}
