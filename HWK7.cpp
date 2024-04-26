#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>
#include <random>

// trapezoid 
double trapezoid_integration(int n, std::function<double(double)> f, double a, double b) {
    double h = (b - a) / n;
    double integral = 0.5 * (f(a) + f(b));

    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        integral += f(x);
    }

    integral *= h;
    return integral;
}

// monte-carlo
double monte_carlo_integration(int n, std::function<double(double)> f, double a, double b) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(a, b);

    double integral = 0.0;

    for (int i = 0; i < n; ++i) {
        double x = dis(gen);
        integral += f(x);
    }

    integral *= (b - a) / n;
    return integral;
}

// Integral 1: e^(-x^2)
double integrand1(double x) {
    return exp(-x * x);
}

// Integral 2: sin^2(1 / (x(2-x) + epsilon))
double integrand2(double x) {
    double epsilon = 1e-12;
    return pow(sin(1 / (x * (2 - x) + epsilon)), 2);
}

int main() {
    std::cout << "\nIntegral 1:\n";
    std::cout << std::setw(8) << "N" << std::setw(20) << "Trapezoid" << std::setw(20) << "Monte Carlo" << std::endl;

    for (int N = 8; N <= 1024; N *= 2) {
        double trapezoid_result1 = trapezoid_integration(N, integrand1, -5, 5);
        double monte_carlo_result1 = monte_carlo_integration(N, integrand1, -5, 5);

        std::cout << std::setw(8) << N << std::setw(20) << trapezoid_result1 << std::setw(20) << monte_carlo_result1 << std::endl;
    }

    std::cout << "\nIntegral 2:\n";
    std::cout << std::setw(8) << "N" << std::setw(20) << "Trapezoid" << std::setw(20) << "Monte Carlo" << std::endl;

    for (int N = 8; N <= 1024; N *= 2) {
        double trapezoid_result2 = trapezoid_integration(N, integrand2, 0, 2);
        double monte_carlo_result2 = monte_carlo_integration(N, integrand2, 0, 2);

        std::cout << std::setw(8) << N << std::setw(20) << trapezoid_result2 << std::setw(20) << monte_carlo_result2 << std::endl;
    }

    return 0;
}
