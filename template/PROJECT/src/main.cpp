//command to rebuild this project: pip install --no-cache-dir --force-reinstall .
#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include <cmath>
#include <vector>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

constexpr double PI_CONST = 3.14159265358979;

int add(int i, int j)
{
	return i + j;
}

std::vector<double> linspace(double start, double end, size_t num)
{
	std::vector<double> result(num);
	double step = (end - start) / (num - 1);
	for (size_t i = 0; i < num; i++)
	{
		result[i] = start + i * step;
	}
	return result;
}

void sinPattern(double frequency, double start, double end, size_t samples)
{
	std::vector<double> timePoints = linspace(start, end, samples);
	std::vector<double> y;
	for (double timePoint : timePoints)
	{
		y.push_back(std::sin(2 * PI_CONST * frequency * timePoint));
	}
	matplot::plot(timePoints, y);
	matplot::title("Sine Wave");
	matplot::xlabel("Time");
	matplot::ylabel("Value");
	matplot::show();
}

void cosPattern(double frequency, double start, double end, size_t samples)
{
	auto timePoints = linspace(start, end, samples);
	std::vector<double> y;
	for (double timePoint : timePoints)
	{
		y.push_back(std::cos(2 * PI_CONST * frequency * timePoint));
	}
	matplot::plot(timePoints, y);
	matplot::title("Cosine Wave");
	matplot::xlabel("Time");
	matplot::ylabel("Value");
	matplot::show();
}

void squarePattern(double frequency, double start, double end, size_t samples)
{
	auto timePoints = linspace(start, end, samples);
	std::vector<double> y;
	for (double timePoint : timePoints)
	{
		y.push_back(std::sin(2 * PI_CONST * frequency * timePoint) >= 0 ? 1.0 : -1.0);
	}
	matplot::plot(timePoints, y);
	matplot::title("Square Wave");
	matplot::xlabel("Time");
	matplot::ylabel("Value");
	matplot::show();
}

void sawtoothPattern(double frequency, double start, double end, size_t samples)
{
	auto timePoints = linspace(start, end, samples);
	std::vector<double> y;
	for (double timePoint : timePoints)
	{
		double value = 2.0 * (timePoint * frequency - std::floor(0.5 + timePoint * frequency));
		y.push_back(value);
	}
	matplot::plot(timePoints, y);
	matplot::title("Sawtooth Wave");
	matplot::xlabel("Time");
	matplot::ylabel("Value");
	matplot::show();
}

namespace py = pybind11;

PYBIND11_MODULE(_core, m)
{
	m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

	m.def("add", &add,
		R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
		)pbdoc"
	);

	m.def("sinPattern", &sinPattern,
		R"pbdoc(
        Draw a sin wave with given arguments
		)pbdoc"
	);

	m.def("subtract", [](int i, int j) { return i - j; },
		R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
		)pbdoc"
	);
}
