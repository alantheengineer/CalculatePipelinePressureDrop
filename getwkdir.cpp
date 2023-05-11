#include<iostream>
#include<fstream>
#include<filesystem>
#include<vector>

#include "pressureloss_calcs.h"

using namespace std;
using namespace pressure_loss;

namespace fs = filesystem;

// for string delimiter
vector<string> split (string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}

int main(){
  string this_dir = fs::current_path();

  vector<pipe_segment> v_pipeline;

  cout << "reading files in: " << this_dir << endl;

  for (auto& p: fs::directory_iterator(this_dir)) {
    string this_file = p.path();
    // setup inlet conditions
    double d_inlet_pressure, d_pressure; // text file enters kPa, turns into MPa internally
    double d_inlet_temperature, d_temperature; // deg C
    double d_inlet_rate, d_rate; // m3/s
    double d_salinity = 0.0; // ppm
    double d_humidity = 0.0; // %
    double d_fluid_property = 0.0;
    int i_fluid_type = -1;
    bool b_output_gradient = false;

    // looking for .pipe file to read
    if(this_file.substr(this_file.find_last_of(".") + 1) == "pipe") {
      string line;
      ifstream myfile(this_file);
      if (myfile.is_open()){
        // print data for file...
        cout << "loading data and calculating: " << endl;
        cout << this_file << endl;
        while (getline(myfile, line)){
          // skip comment or empty lines
          if (line.length() == 0 || line.at(0) == '#') continue;

          // split the string using the space delimiter
          string delimiter = " ";
          vector<string> v = split(line, delimiter);

          // iterate through the vector to extract data
          for (auto i = 0; i < v.size(); ++i) {
            if (v[i] == "fitting") {
              pipe_segment* seg = new pipe_segment();
              seg->i_segment_type = 1;
              // set default values
              seg->d_segment_diameter = 0.1;
              for (auto j = i; j < v.size(); ++j) {
                if (v[j] == "k") seg->d_segment_k_value = atof(v[++j].c_str());
                if (v[j] == "diameter") seg->d_segment_diameter = atof(v[++j].c_str());
                if (v[j] == "multiplier") seg->d_rate_multiplier = atof(v[++j].c_str());
              }
              v_pipeline.push_back(*seg);
            }
            if (v[i] == "pipe") {
              pipe_segment* seg = new pipe_segment();
              seg->i_segment_type = 0;
              // set default values
              seg->d_segment_length = 1.0;
              seg->d_segment_diameter = 0.1;
              seg->d_segment_angle = 0.0;
              seg->d_segment_roughness = 0.01524;
              seg->d_rate_multiplier = 1.0;
              for (auto j = i; j < v.size(); ++j) {
                if (v[j] == "length") seg->d_segment_length = atof(v[++j].c_str());
                if (v[j] == "diameter") seg->d_segment_diameter = atof(v[++j].c_str());
                if (v[j] == "roughness") seg->d_segment_roughness = atof(v[++j].c_str());
                if (v[j] == "angle") seg->d_segment_angle = atof(v[++j].c_str());
                if (v[j] == "multiplier") seg->d_rate_multiplier = atof(v[++j].c_str());
              }
              v_pipeline.push_back(*seg);
            }
            if (v[i] == "print") b_output_gradient = true;
            if (v[i] == "pressure") d_inlet_pressure = d_pressure = atof(v[++i].c_str());
            if (v[i] == "temperature") d_inlet_temperature = d_temperature = atof(v[++i].c_str());
            if (v[i] == "rate") d_inlet_rate = d_rate = atof(v[++i].c_str());
            if (v[i] == "humidity") d_fluid_property = atof(v[++i].c_str());
            if (v[i] == "salinity") d_fluid_property = atof(v[++i].c_str());
            if (v[i] == "fluid") {
              if (v[++i] == "air"){
                i_fluid_type = 1;
              } else {
                i_fluid_type = 0;
              }
            }
          }
        }

        string s_fluid = (i_fluid_type == 1) ? "air" : "water";
        cout << "Inlet pressure: " << d_pressure << " MPa" << endl;
        cout << "Inlet temperature: " << d_temperature << " deg C" << endl;
        cout << "Fluid rate: " << d_rate << " m3/s of " << s_fluid << endl;

        if (b_output_gradient) {
          cout << "" << endl;
        }

        for (vector<pipe_segment>::iterator it = v_pipeline.begin(); it != v_pipeline.end(); ++it) {
          d_pressure =calc_pressure_loss(&(*it),&d_pressure,&d_temperature,&d_rate,&d_fluid_property,&i_fluid_type,&b_output_gradient);
        }
        double d_total_pressure_loss = (d_inlet_pressure - d_pressure)*1.0E06;
        cout << "outlet pressure of pipe: " << d_pressure << " MPa " << endl;
        cout << "total pressure loss: " << d_total_pressure_loss << " Pa" << endl;
        myfile.close();
      } else {
        cout << "Unable to open file";
      }
    }
  }

  return 1;
}
