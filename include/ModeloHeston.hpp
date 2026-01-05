/*
 * ModeloHeston.hpp - Header Heston - C++23
 * Autor: Luiz Tiago Wilcke
 */

#ifndef MODELO_HESTON_HPP
#define MODELO_HESTON_HPP

#include <vector>
#include <random>
#include <utility>

namespace analise_ouro {

class ModeloHeston {
private:
    double precoInicial, volInicial, taxaLivreRisco;
    double kappa, theta, xi, rho;
    double dt;
    mutable std::mt19937_64 gerador;

public:
    ModeloHeston();
    ModeloHeston(double s0, double v0, double r, double kappa_, double theta_, double xi_, double rho_);
    ~ModeloHeston();

    void definirParametros(double s0, double v0, double r, double kappa_, double theta_, double xi_, double rho_);

    [[nodiscard]] auto simularTrajetoria(int numPassos) const -> std::pair<std::vector<double>, std::vector<double>>;
    [[nodiscard]] auto simularTrajetoriaQE(int numPassos) const -> std::pair<std::vector<double>, std::vector<double>>;
    [[nodiscard]] auto simularMultiplasTrajetorias(int numTrajetorias, int numPassos) const -> std::vector<std::pair<std::vector<double>, std::vector<double>>>;
    [[nodiscard]] auto precoCallMonteCarlo(double strike, double tempoExpiracao, int numSimulacoes) const -> double;
    [[nodiscard]] auto verificarCondicaoFeller() const -> bool;

    [[nodiscard]] auto obterKappa() const noexcept -> double;
    [[nodiscard]] auto obterTheta() const noexcept -> double;
    [[nodiscard]] auto obterXi() const noexcept -> double;
    [[nodiscard]] auto obterRho() const noexcept -> double;
};

} // namespace analise_ouro

#endif // MODELO_HESTON_HPP
