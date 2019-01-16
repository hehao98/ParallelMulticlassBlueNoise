CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/RainbowNoiseLibrary)

ADD_LIBRARY(RainbowNoiseSpectrum

DeformableLFT.cpp
DeformableLFT.hpp
Deformer.cpp
Deformer.hpp
IsotropicDeformer.cpp
IsotropicDeformer.hpp
LFT.cpp
LFT.hpp
LinearLFT.cpp
LinearLFT.hpp
PFT.cpp
PFT.hpp
RandomLFT.cpp
RandomLFT.hpp
SFT.cpp
SFT.hpp
)
