#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include <pybind11/complex.h>
#include <pybind11/operators.h>
#include <pybind11/pytypes.h>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j) {
    return i + j;
}

#include <vector>
#include <cmath>

std::vector<double> generate_sine_wave(double freq, double sample_rate, int samples) {
    std::vector<double> result(samples);
    for (int i = 0; i < samples; ++i) {
        result[i] = sin(2 * 3.14159265358979323846 * freq * i / sample_rate);
    }
    return result;
}

std::vector<double> generate_cosine_wave(double freq, double sample_rate, int samples) {
    std::vector<double> result(samples);
    for (int i = 0; i < samples; ++i) {
        result[i] = cos(2 * 3.14159265358979323846 * freq * i / sample_rate);
    }
    return result;
}

std::vector<double> generate_square_wave(double freq, double sample_rate, int samples) {
    std::vector<double> result(samples);
    for (int i = 0; i < samples; ++i) {
        double t = static_cast<double>(i) / sample_rate;
        result[i] = (sin(2 * M_PI * freq * t) >= 0) ? 1.0 : -1.0;
    }
    return result;
}

std::vector<double> generate_sawtooth_wave(double freq, double sample_rate, int samples) {
    std::vector<double> result(samples);
    for (int i = 0; i < samples; ++i) {
        double t = static_cast<double>(i) / sample_rate;
        result[i] = 2.0 * (t * freq - floor(t * freq + 0.5));
    }
    return result;
}

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
		   generate_sine_wave
    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

	m.def("generate_sine_wave", &generate_sine_wave,
		R"pbdoc(
		Generate a sine wave

		Args:
			freq (float): Frequency of the sine wave.
			sample_rate (float): Sample rate.
			samples (int): Number of samples to generate.

		Returns:
			list: Generated sine wave samples.
		)pbdoc");
        
    m.def("generate_cosine_wave", &generate_cosine_wave,
        R"pbdoc(
        Generate a cosine wave

        Args:
            freq (float): Frequency of the cosine wave.
            sample_rate (float): Sample rate.
            samples (int): Number of samples to generate.

        Returns:
            list: Generated cosine wave samples.
    )pbdoc");
        
    m.def("generate_square_wave", &generate_square_wave,
        R"pbdoc(
        Generate a square wave

        Args:
            freq (float): Frequency of the wave.
            sample_rate (float): Sample rate.
            samples (int): Number of samples.

        Returns:
            list: Square wave samples.
    )pbdoc");
        

    m.def("generate_sawtooth_wave", &generate_sawtooth_wave,
        R"pbdoc(
        Generate a sawtooth wave

        Args:
            freq (float): Frequency of the wave.
            sample_rate (float): Sample rate.
            samples (int): Number of samples.

        Returns:
            list: Sawtooth wave samples.
    )pbdoc");
        



    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
