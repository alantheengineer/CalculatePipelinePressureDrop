#include <math.h>
#include <iostream>
#include "pressureloss_calcs.h"
#include "enhancedwater.h"

using namespace std;
using namespace fluid_properties;

double pressure_loss::friction_factor(const double * relativeRoughness, const double * reynoldsNumber)
{
	double ff = 0.002;
	double rn = *reynoldsNumber;
	if (rn < 3000.0)
	{
		if (rn < 700.0)
		{
			ff = 64.0/rn;
		}
		else
		{
			ff = 64.0 / rn + 1.008E-05*(rn - 700);
		}
	}
	else
	{
		for (int i = 0; i < 5; i++)
		{
			double rr = *relativeRoughness;
			double x = log10(0.27*rr + 2.52 / rn / sqrt(ff));
			ff = 0.25 / x / x;
		}
	}

	return ff;
}

double pressure_loss::calc_pressure_loss(const pipe_segment* p_segment,
																		const double* p_pressureIn,
																		const double* p_temperature,
																		const double* p_rate,
																		const double* p_property,
																	  const int* p_fluid_type,
																	  const bool* b_print_out)
{
	double dDensity = 0.0;
	double dStdDensity = 0.0;
	double dPStd = 0.101325348872;
	double dTStd = 21.0;
	double dCompressibility = 0.0;
	double dViscosity = 0.0;
	double dViscosity_metric = 0.0;

	double d_current_pressure = *p_pressureIn;

	double d_calculation_step_length = p_segment->d_segment_length * 0.1;
	int i_num_steps = 10;
	// always taking at least 10 steps along a pipe_segment
	if (d_calculation_step_length < 0.1) {
		d_calculation_step_length = p_segment->d_segment_length; // unless it is just a metre in length
		i_num_steps = 1;
	}

	if (p_segment->i_segment_type == 1) i_num_steps = 1;

	for (int i_calc = 0; i_calc < i_num_steps; ++i_calc) {

		if (*p_fluid_type == 0)
		{
			// std
			iapws95(dPStd, dTStd, *p_property, &dStdDensity, &dCompressibility);
			// conditions
			iapws95(d_current_pressure, *p_temperature, *p_property, &dDensity, &dCompressibility);
			dViscosity = kestin_viscosity(d_current_pressure, *p_temperature, *p_property);
			dViscosity_metric = dViscosity * 0.001;
		}

		if (*p_fluid_type == 1)
		{
			// std
			dStdDensity = air_density(dPStd, dTStd, *p_property);
			// conditions
			dDensity = air_density(d_current_pressure, *p_temperature, *p_property);
			dCompressibility = air_compressibility_factor(d_current_pressure, *p_temperature, *p_property);
			dViscosity = air_viscosity(d_current_pressure, *p_temperature, *p_property);

			dViscosity_metric = dViscosity * 0.000001;
		}

		double dMassRate = dStdDensity * (*p_rate);
		double dVelocity = (dMassRate / dDensity) / (3.14159265359 * pow(0.5*p_segment->d_segment_diameter,2));

		dVelocity *= p_segment->d_rate_multiplier;

		double dTotalGradient;
		double dFF;

		if (p_segment->i_segment_type == 0) {

			double dRelativeRoughness = (0.001*p_segment->d_segment_roughness) / p_segment->d_segment_diameter;
			double dReynoldsNumber = (dDensity*dVelocity*p_segment->d_segment_diameter) / dViscosity_metric;

			dFF = friction_factor(&dRelativeRoughness, &dReynoldsNumber);

			double dFrictionGradient = dFF*dDensity*pow(dVelocity, 2) / (2 * p_segment->d_segment_diameter);
			double dGravityHead = dDensity*sin(p_segment->d_segment_angle) * 9.81;
			double dAccelerationHead = dDensity*dVelocity*dVelocity / (d_current_pressure*9.81*1.E06);

			if (dAccelerationHead > 0.99) dAccelerationHead = 0.99;

			double dAccelerationGradient = (dFrictionGradient + dGravityHead)*dAccelerationHead;

			dTotalGradient = dGravityHead + dFrictionGradient + dAccelerationGradient;

		} else {
			double d_head_lost = p_segment->d_segment_k_value * (dVelocity*dVelocity)/(2*9.81);

			d_current_pressure -= (d_head_lost * 9.81 * dDensity)*1.0E-06;

			if (*b_print_out) {
				cout << "fitting: " << d_current_pressure << " MPa " << dVelocity << " m/s " << d_head_lost << " m " << dDensity << " kg/m3" << endl;
			}
			return d_current_pressure;
		}

		d_current_pressure -= (dTotalGradient*d_calculation_step_length) * 1.0E-06;

		if (*b_print_out) {
			cout << "pipe: " << d_current_pressure << " MPa " << dVelocity << " m/s " << dTotalGradient << " Pa/m " << dDensity << " kg/m3 " << dFF << " - friction factor " << endl;
		}
	}

	return d_current_pressure;
}
