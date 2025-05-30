//command to rebuild this project: pip install --no-cache-dir --force-reinstall .
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <matplot/matplot.h>
#include <cmath>
#include <vector>
#include <utility>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

constexpr double PI_CONST = 3.14159265358979;

std::vector<std::pair<double, double>> sinPattern(double frequency, double start, double end, size_t samples)
{
	std::vector<std::pair<double, double>> result;
	result.reserve(samples);
	double step = (end - start) / (samples - 1);
	for (size_t i = 0; i < samples; ++i)
	{
		double timePoint = start + i * step;
		double y = std::sin(2 * PI_CONST * frequency * timePoint);
		result.emplace_back(timePoint, y);
	}
	return result;
}

std::vector<std::pair<double, double>> cosPattern(double frequency, double start, double end, size_t samples)
{
	std::vector<std::pair<double, double>> result;
	result.reserve(samples);
	double step = (end - start) / (samples - 1);
	for (size_t i = 0; i < samples; ++i)
	{
		double timePoint = start + i * step;
		double y = std::cos(2 * PI_CONST * frequency * timePoint);
		result.emplace_back(timePoint, y);
	}
	return result;
}

std::vector<std::pair<double, double>> squarePattern(double frequency, double start, double end, size_t samples)
{
	std::vector<std::pair<double, double>> result;
	result.reserve(samples);
	double step = (end - start) / (samples - 1);
	for (size_t i = 0; i < samples; ++i)
	{
		double timePoint = start + i * step;
		double y = std::sin(2 * PI_CONST * frequency * timePoint) >= 0 ? 1.0 : -1.0;
		result.emplace_back(timePoint, y);
	}
	return result;
}

std::vector<std::pair<double, double>> sawtoothPattern(double frequency, double start, double end, size_t samples)
{
	std::vector<std::pair<double, double>> result;
	result.reserve(samples);
	double step = (end - start) / (samples - 1);
	for (size_t i = 0; i < samples; ++i)
	{
		double timePoint = start + i * step;
		double value = 2.0 * (timePoint * frequency - std::floor(0.5 + timePoint * frequency));
		result.emplace_back(timePoint, value);
	}
	return result;
}

std::vector<std::pair<double, std::pair<double, double>>> discreteFourierTransform(const std::vector<std::pair<double, double>>& patternSamples)
{
	std::vector<std::pair<double, std::pair<double, double>>> dftResult;
	size_t N = patternSamples.size();
	dftResult.reserve(N);
	for (size_t k = 0; k < N; k++)
	{
		double real = 0.0;
		double imag = 0.0;
		for (size_t n = 0; n < N; n++)
		{
			double angle = 2 * PI_CONST * k * n / N;
			real += patternSamples[n].second * std::cos(angle);
			imag += patternSamples[n].second * -std::sin(angle);
		}
		dftResult.emplace_back(std::make_pair(static_cast<double>(k), std::make_pair(real, imag)));
	}
	return dftResult;
}

void showDFTPattern(const std::vector<std::pair<double, std::pair<double, double>>>& dftPattern)
{
	std::vector<double> x(dftPattern.size()), yReal(dftPattern.size()), yImag(dftPattern.size());
	for (size_t i = 0; i < dftPattern.size(); ++i)
	{
		x[i] = dftPattern[i].first;
		yReal[i] = dftPattern[i].second.first;
		yImag[i] = dftPattern[i].second.second;
	}
	matplot::figure();
	matplot::plot(x, yReal);
	matplot::title("Real Part of DFT");
	matplot::xlabel("Frequency Index");
	matplot::ylabel("Amplitude");
	matplot::show();
}

void inverseDiscreteFourierTransform(const std::vector<std::pair<double, std::pair<double, double>>>& realImagPattern)
{
	size_t N = realImagPattern.size();
	std::vector<double> functionSamples(N);
	std::vector<double> x(N);
	for (size_t n = 0; n < N; n++)
	{
		x[n] = realImagPattern[n].first;
		functionSamples[n] = 0.0;
		for (size_t k = 0; k < N; k++)
		{
			double angle = 2 * PI_CONST * k * n / N;
			functionSamples[n] += realImagPattern[k].second.first * std::cos(angle) - realImagPattern[k].second.second * std::sin(angle);
		}
		functionSamples[n] /= N;
	}
	matplot::figure();
	matplot::plot(x, functionSamples);
	matplot::title("Inverse DFT");
	matplot::xlabel("Sample Index");
	matplot::ylabel("Amplitude");
	matplot::show();
}

void showPattern(const std::vector<std::pair<double, double>>& patternSamples)
{
	std::vector<double> x(patternSamples.size()), y(patternSamples.size());
	for (size_t i = 0; i < patternSamples.size(); ++i)
	{
		x[i] = patternSamples[i].first;
		y[i] = patternSamples[i].second;
	}
	matplot::figure();
	matplot::plot(x, y);
	matplot::title("Pattern");
	matplot::xlabel("Time");
	matplot::ylabel("Amplitude");
	matplot::show();
}

namespace py = pybind11;

// Register std::pair<double, double> as a Python tuple
PYBIND11_MODULE(_core, m)
{
	py::class_<std::pair<double, double>>(m, "PairDoubleDouble")
		.def(py::init<>())
		.def(py::init<double, double>())
		.def_readwrite("first", &std::pair<double, double>::first)
		.def_readwrite("second", &std::pair<double, double>::second)
		.def("__repr__",
			[](const std::pair<double, double>& p) {
				return "(" + std::to_string(p.first) + ", " + std::to_string(p.second) + ")";
			}
		);

	py::bind_vector<std::vector<std::pair<double, double>>>(m, "VectorPairDoubleDouble");

	m.doc() = R"pbdoc(
        Pybind11 math visualisation plugin
        ----------------------------------
    )pbdoc";

	m.def("sinPattern", &sinPattern, py::arg("frequency"), py::arg("start"), py::arg("end"), py::arg("samples"),
		R"pbdoc(
        Draw a sinus wave with given arguments
    )pbdoc"
	);

	m.def("cosPattern", &cosPattern, py::arg("frequency"), py::arg("start"), py::arg("end"), py::arg("samples"),
		R"pbdoc(
        Draw a cosinus wave with given arguments
    )pbdoc"
	);

	m.def("squarePattern", &squarePattern, py::arg("frequency"), py::arg("start"), py::arg("end"), py::arg("samples"),
		R"pbdoc(
        Draw a square wave with given arguments
    )pbdoc"
	);

	m.def("sawtoothPattern", &sawtoothPattern, py::arg("frequency"), py::arg("start"), py::arg("end"), py::arg("samples"),
		R"pbdoc(
        Draw a sawtooth wave with given arguments
    )pbdoc"
	);

	m.def("discreteFourierTransform", &discreteFourierTransform, py::arg("patternSamples"),
		R"pbdoc(
    Perform discrete Fourier transform on given samples
    )pbdoc"
	);

	m.def("showDFTPattern", &showDFTPattern, py::arg("dftPattern"),
		R"pbdoc(
	Show the discrete Fourier transform pattern in a plot
	)pbdoc"
	);

	m.def("inverseDiscreteFourierTransform", &inverseDiscreteFourierTransform, py::arg("realImagPattern"),
		R"pbdoc(
    Perform inverse discrete Fourier transform on given real and imaginary parts
    )pbdoc"
	);
	m.def("showPattern", &showPattern, py::arg("patternSamples"),
		R"pbdoc(
		Show the pattern samples in a plot
	)pbdoc"
	);

}
