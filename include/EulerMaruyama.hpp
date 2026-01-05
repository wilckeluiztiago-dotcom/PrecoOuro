/*
 * EulerMaruyama.hpp - Header Euler-Maruyama - C++23
 * Autor: Luiz Tiago Wilcke
 */

#ifndef EULER_MARUYAMA_HPP
#define EULER_MARUYAMA_HPP

#include <vector>
#include <random>
#include <functional>

namespace analise_ouro {

using FuncaoDrift = std::function<double(double x, double t)>;
using FuncaoDifusao = std::function<double(double x, double t)>;

class EulerMaruyama {
private:
    double dt;
    int numPassos;
    mutable std::mt19937_64 gerador;

public:
    EulerMaruyama();
    EulerMaruyama(double deltaT, int nPassos);
    ~EulerMaruyama();

    void definirParametros(double deltaT, int nPassos);
    [[nodiscard]] auto resolver(double x0, FuncaoDrift drift, FuncaoDifusao difusao) const -> std::vector<double>;
    [[nodiscard]] auto resolverGBM(double s0, double mu, double sigma) const -> std::vector<double>;
    [[nodiscard]] auto resolverOU(double x0, double theta, double mu, double sigma) const -> std::vector<double>;
    [[nodiscard]] auto resolverCIR(double x0, double a, double b, double sigma) const -> std::vector<double>;
    [[nodiscard]] auto resolverMultiplasTrajetorias(double x0, FuncaoDrift drift, FuncaoDifusao difusao, int numTrajetorias) const -> std::vector<std::vector<double>>;
    [[nodiscard]] auto calcularErroForte(const std::vector<double>& numerica, const std::vector<double>& exata) const -> double;
    [[nodiscard]] auto calcularErroFraco(const std::vector<std::vector<double>>& numericas, const std::vector<double>& mediasExatas) const -> double;
    [[nodiscard]] auto obterPassoTempo() const noexcept -> double;
    [[nodiscard]] auto obterNumPassos() const noexcept -> int;
};

} // namespace analise_ouro

#endif // EULER_MARUYAMA_HPP
