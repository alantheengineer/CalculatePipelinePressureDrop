#pragma once
#include <vector>

namespace pressure_loss
{
	struct pipe_segment {
	  int i_segment_type;
	  double d_segment_length;
	  double d_segment_roughness;
	  double d_segment_diameter;
	  double d_segment_angle;
	  double d_segment_k_value;
		double d_rate_multiplier;
	};

	double friction_factor(const double* relativeRoughness, const double* reynoldsNumber);
	double calc_pressure_loss(const pipe_segment* p_segment, const double* p_pressureIn, const double* p_temperature, const double* p_rate, const double* p_property, const int* p_fluid_type, const bool* b_print_out);
}
