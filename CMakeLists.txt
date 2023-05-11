cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)

project(calc_pressure_loss)

link_libraries(stdc++fs)

add_executable(
  calc_pressure_loss.exe
  getwkdir.cpp
  pressureloss_calcs.cpp
  enhancedwater.cpp
  )
