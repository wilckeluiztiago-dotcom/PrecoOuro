/*
 * SimulacaoMonteCarlo.hpp - Header Monte Carlo - C++23
 * Autor: Luiz Tiago Wilcke
 */

#ifndef SIMULACAO_MONTE_CARLO_HPP
#define SIMULACAO_MONTE_CARLO_HPP

#include <vector>
#include <random>
#include <utility>

namespace analise_ouro {

struct ResultadoMonteCarlo {
    double media;
    double variancia;
    double erroPadrao;
    double precoDescontado;
    std::pair<double, double> intervaloConfianca95;
};

class SimulacaoMonteCarlo {
private:
    int numSimulacoes, numPassos;
    double dt;
    mutable std::mt19937_64 gerador;

    [[nodiscard]] auto calcularEstatisticas(const std::vector<double>& payoffs, double r, double T) const -> ResultadoMonteCarlo;

public:
    SimulacaoMonteCarlo();
    SimulacaoMonteCarlo(int nSim, int nPassos);
    ~SimulacaoMonteCarlo();

    void definirParametros(int nSim, int nPassos);
    [[nodiscard]] auto simularGBM(double s0, double mu, double sigma) const -> std::vector<std::vector<double>>;
    [[nodiscard]] auto calcularPrecoOpcaoCall(double s0, double K, double r, double sigma, double T) const -> ResultadoMonteCarlo;
    [[nodiscard]] auto calcularPrecoOpcaoPut(double s0, double K, double r, double sigma, double T) const -> ResultadoMonteCarlo;
    [[nodiscard]] auto simularVaR(double s0, double mu, double sigma, double nivelConfianca = 0.95) const -> double;
    [[nodiscard]] auto simularCVaR(double s0, double mu, double sigma, double nivelConfianca = 0.95) const -> double;
    [[nodiscard]] auto calcularPercentis(const std::vector<std::vector<double>>& trajetorias, int passo) const -> std::vector<double>;
    [[nodiscard]] auto obterNumSimulacoes() const noexcept -> int;
    [[nodiscard]] auto obterNumPassos() const noexcept -> int;
};

} // namespace analise_ouro

#endif // SIMULACAO_MONTE_CARLO_HPP
