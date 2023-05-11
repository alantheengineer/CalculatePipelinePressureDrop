#include "enhancedwater.h"
#include <math.h>

double utility_functions::tempfunc(const std::vector<double>& vTab, const double temperature)
{
	double dFuncReturn = 0.0;

	double dTempdiv = temperature * 0.01;
	double dTempsquared = dTempdiv*dTempdiv;

	double dNumerator = vTab[0] * dTempsquared + vTab[1] * dTempdiv + vTab[2];
	double dDemoninator = vTab[3] * dTempsquared + vTab[4] * dTempdiv + 1;

	if (fabs(dDemoninator) > 0.0) dFuncReturn = dNumerator / dDemoninator;

	return dFuncReturn;
}

void fluid_properties::iapws95(const double pressure, const double temperature, const double salinity, double * density, double * compressibility)
{
	// Calculation of water properties using Spivey & McCain method
	// of interpolating the IAPWS - 95 equation of state for brine.
	//
	// JCPT July 2004, Vol 43, No. 7
	// 01 / 05 / 2015 - Alan Tominey

	/* Calculations are performed using metric(�C and MPa). Salinity is ppm for ease of use so convert to mol / kg */

	double molal_salinity = salinity / 54880;

	// calculate reference water density and compressibility
	double ref_water_density = utility_functions::tempfunc(water_data_tables::tabdwt, temperature); // water density at 70 MPa and given temperature
	double pressure_coeff = utility_functions::tempfunc(water_data_tables::tabewt, temperature); // pressure coefficient at given temperature
	double temperature_coeff = utility_functions::tempfunc(water_data_tables::tabfwt, temperature); // temperature coefficient at given temperature

	double inv_ewt = 1.0;
	if (fabs(pressure_coeff) > 0.0) inv_ewt = 1.0 / pressure_coeff;

	// calculate compressibility of water at given t and p
	double iwpt = inv_ewt*log(pressure_coeff*(pressure*0.0142857)+temperature_coeff);
	double iwrt = inv_ewt*log(pressure_coeff + temperature_coeff);

	// calculate density at given t and p
	double water_density = ref_water_density * exp(iwpt - iwrt);

	if (salinity > 0.0)
	{
		double sqrt_molal=sqrt(molal_salinity);

		// coefficients for brine at 70 MPa
		double dcm2 = utility_functions::tempfunc(water_data_tables::tdcm2, temperature);
		double dcm32 = utility_functions::tempfunc(water_data_tables::tdcm32, temperature);
		double dcm1 = utility_functions::tempfunc(water_data_tables::tdcm1, temperature);
		double dcm12 = utility_functions::tempfunc(water_data_tables::tdcm12, temperature);

		double brine_poly = dcm2*(molal_salinity*molal_salinity) + dcm32*pow(molal_salinity, 1.5) + dcm1*molal_salinity + dcm12*sqrt_molal;
		double ref_brine_density = ref_water_density + brine_poly;

		// coefficients for brine compressibility at t and p and s
		double ecmt = utility_functions::tempfunc(water_data_tables::tecm2, temperature);
		double ebtcm = pressure_coeff + ecmt*molal_salinity;
		double inv_etcm = 1.0 / ebtcm;

		double fcm32 = utility_functions::tempfunc(water_data_tables::tfcm32, temperature);
		double fcm1 = utility_functions::tempfunc(water_data_tables::tfcm1, temperature);
		double fcm12 = utility_functions::tempfunc(water_data_tables::tfcm12, temperature);

		brine_poly = fcm32*pow(molal_salinity, 1.5) + fcm1 * molal_salinity + fcm12 * sqrt_molal;
		double fbtcm = temperature_coeff + brine_poly;
		double ibptcm = inv_etcm*log(ebtcm*(pressure*0.0142857) + fbtcm);
		double ibrtcm = inv_etcm*log(ebtcm + fbtcm);

		double brine_density = ref_brine_density * exp(ibptcm - ibrtcm);

		water_density = brine_density;
		pressure_coeff = ebtcm;
		temperature_coeff = fbtcm;
	}

	*density = water_density*1000.0; // returned in kg/m^3
	*compressibility = 1.42857e-2 * (1 / (pressure_coeff*(pressure*1.42857E-2) + temperature_coeff));
}

double fluid_properties::kestin_viscosity(const double pressure, const double temperature, const double salinity)
{
	// Calculation of brine viscosity using modified Kestin.et.al
	//
	// Modified to use Chen et.al correlation for theoretical brine
	// viscosity at STP tuned to IAPWS - 2008, considered valid + / -1.5%
	// between 0 - 180 �C
	// Pressure correction applied from Kestin, pressure up to 500 MPa
	//
	// Journal of the Acoustical Society of America, 62, 1129 - 1135, 1977
	// J.Phys.Chem.Ref.Data Vol. 10. No. 1, 1981
	// 23 / 02 / 2016 - Alan Tominey

	/* Calculations are performed using metric(�C and MPa). Salinity is ppm for ease of use so convert to mol / kg */

	double molal_salinity = salinity / 54880;
	double ww_salinity = molal_salinity*0.05544; // converting to kg/kg for other calculations
	double coefficient1 = (0.0031*molal_salinity - 0.0514)*molal_salinity + 0.9991;

	// Calculate McCain & Spivey correction factor
	double mcSLNT = log(temperature*8.0e-03);
	double mcSPVA = 0.068*(pressure*0.01) + 0.0173;
	double mcSPVB = 0.0273*(pressure*0.01) - 1.0531;
	double mcSPVF = (mcSLNT*mcSPVA + mcSPVB)*mcSLNT;

	if (mcSPVF > 10.0) mcSPVF = 10.0;

	mcSPVF = exp(mcSPVF);

	// Calculate theoretical viscosity at T and ref. pressure (0)
	double tmetpl = temperature + 64.993;
	tmetpl = (0.157 * tmetpl * tmetpl) - 91.296;
	double iawvis = ((4.42844E-5 + (1 / tmetpl))) * 1.0E06; //- uPa.s

	// calculate pressure correction parameters
	double betaW = (((-1.05e-8*temperature + 4.47e-6)*temperature - 6.97e-4)*temperature + 5.74e-2)*temperature - 1.297;
	double expc = 0.545 + (2.80E-3*temperature) - betaW;
	double satcon = (3.60E-5*temperature + 2.8E-3)*temperature + 6.044;

	double consm = (salinity / 54880) / satcon;
	double redprc = ((0.5*consm - 2)*consm + 2.5)*consm;

	double pcoeff = expc*redprc + betaW;

	double calcmu = iawvis*(1 + pcoeff*(pressure*0.001));

	return (calcmu*0.001)*coefficient1;
}

double fluid_properties::air_density(const double pressure /*MPa*/, const double temperature /*deg C*/, const double humidity /*%*/)
{
	double temperature_abs = temperature + 273.15;
	double pressure_pascals = pressure * 1E06;

	double water_mol_fraction = air_water_mol_fraction(pressure, temperature, humidity);

	double Z = air_compressibility_factor(pressure, temperature,humidity);

	double xc = 0.0004;
	double density = ((3.48349 + 1.44*(xc - 0.0004))*0.001)*(pressure_pascals / (Z*temperature_abs))*(1.0 - 0.378*water_mol_fraction);

	return density;
}

double fluid_properties::air_water_mol_fraction(const double pressure, const double temperature, const double humidity)
{
	double temperature_abs = temperature + 273.15;
	double pressure_pascals = pressure * 1E06;

	double saturated_vapour_pres = exp(air_data_tables::tabpsv[0] * pow(temperature_abs, 2)
		+ air_data_tables::tabpsv[1] * temperature_abs
		+ air_data_tables::tabpsv[2]
		+ air_data_tables::tabpsv[3] / temperature_abs);

	double enhancement_fact = air_data_tables::tabef[0] + air_data_tables::tabef[1] * pressure_pascals + air_data_tables::tabef[2] * pow(temperature, 2);

	double water_mol_fraction = (humidity / 100.0)*(saturated_vapour_pres / pressure_pascals)*enhancement_fact;

	return water_mol_fraction;
}

double fluid_properties::air_compressibility_factor(const double pressure, const double temperature, const double humidity)
{
	double temperature_abs = temperature + 273.15;
	double pressure_pascals = pressure * 1E06;

	double water_mol_fraction = air_water_mol_fraction(pressure, temperature, humidity);

	double Z = 1.0 - (pressure_pascals / temperature_abs)*(air_data_tables::tabz[0]
		+ air_data_tables::tabz[1] * temperature
		+ air_data_tables::tabz[2] * pow(temperature, 2)
		+ (air_data_tables::tabz[3] + air_data_tables::tabz[4] * temperature)*water_mol_fraction
		+ (air_data_tables::tabz[5] + air_data_tables::tabz[6] * temperature)*pow(water_mol_fraction, 2))
		+ pow(pressure_pascals / temperature_abs, 2)*(air_data_tables::tabz[7] + air_data_tables::tabz[8] * pow(water_mol_fraction, 2));

	return Z;
}

double fluid_properties::air_viscosity(const double pressure, const double temperature, const double humidity)
{
	double temperature_abs = temperature + 273.15;
	double pressure_pascals = pressure * 1E06;

	double water_mol_fraction = air_water_mol_fraction(pressure, temperature, humidity);

	double air_viscosity = 1e-02*(air_data_tables::tabvisc[0]
		+ air_data_tables::tabvisc[1] * temperature_abs
		+ (air_data_tables::tabvisc[2] + air_data_tables::tabvisc[3] * temperature_abs)*water_mol_fraction
		+ air_data_tables::tabvisc[4] * pow(temperature_abs, 2)
		+ air_data_tables::tabvisc[5] * pow(water_mol_fraction, 2));

	return air_viscosity;
}
