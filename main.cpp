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

#include <complex>

std::vector<std::complex<double>> dft(const std::vector<double>& input) {
    int N = input.size();
    std::vector<std::complex<double>> output(N);
    for (int k = 0; k < N; ++k) {
        std::complex<double> sum = 0.0;
        for (int n = 0; n < N; ++n) {
            double angle = -2 * M_PI * k * n / N;
            sum += std::polar(input[n], angle);
        }
        output[k] = sum;
    }
    return output;
}

std::vector<double> idft(const std::vector<std::complex<double>>& input) {
    int N = input.size();
    std::vector<double> output(N);
    for (int n = 0; n < N; ++n) {
        std::complex<double> sum = 0.0;
        for (int k = 0; k < N; ++k) {
            double angle = 2 * M_PI * k * n / N;
            sum += input[k] * std::polar(1.0, angle);
        }
        output[n] = sum.real() / N;
    }
    return output;
}

std::vector<double> filter_1d(const std::vector<double>& signal, const std::vector<double>& kernel) {
    int n = signal.size();
    int k = kernel.size();
    std::vector<double> result(n, 0.0);
    int offset = k / 2;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            int idx = i + j - offset;
            if (idx >= 0 && idx < n) {
                result[i] += signal[idx] * kernel[j];
            }
        }
    }
    return result;
}

std::vector<std::vector<double>> filter_2d(const std::vector<std::vector<double>>& input, const std::vector<std::vector<double>>& kernel) {
    int rows = input.size();
    int cols = input[0].size();
    int k_rows = kernel.size();
    int k_cols = kernel[0].size();
    int row_offset = k_rows / 2;
    int col_offset = k_cols / 2;

    std::vector<std::vector<double>> output(rows, std::vector<double>(cols, 0.0));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double sum = 0.0;
            for (int ki = 0; ki < k_rows; ++ki) {
                for (int kj = 0; kj < k_cols; ++kj) {
                    int ni = i + ki - row_offset;
                    int nj = j + kj - col_offset;
                    if (ni >= 0 && ni < rows && nj >= 0 && nj < cols) {
                        sum += input[ni][nj] * kernel[ki][kj];
                    }
                }
            }
            output[i][j] = sum;
        }
    }

    return output;
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

    m.def("dft", &dft, "Compute DFT of real signal");
    m.def("idft", &idft, "Compute inverse DFT");

    m.def("filter_1d", &filter_1d, "Apply 1D filter (convolution)");

    m.def("filter_2d", &filter_2d, "Apply 2D filter to 2D signal");



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
