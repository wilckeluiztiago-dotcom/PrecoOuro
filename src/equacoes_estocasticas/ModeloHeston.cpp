/*
 * ModeloHeston.cpp - Modelo de Volatilidade Estocástica de Heston - C++23
 * Autor: Luiz Tiago Wilcke
 * Projeto: Análise do Preço do Ouro
 */

#include "ModeloHeston.hpp"
#include <cmath>
#include <random>
#include <algorithm>
#include <numbers>

namespace analise_ouro {

ModeloHeston::ModeloHeston() 
    : precoInicial(1800.0), volInicial(0.04), taxaLivreRisco(0.05),
      kappa(2.0), theta(0.04), xi(0.3), rho(-0.7), dt(1.0/252.0) {
    auto semente = std::random_device{}();
    gerador.seed(semente);
}

ModeloHeston::ModeloHeston(double s0, double v0, double r, double kappa_, double theta_, double xi_, double rho_)
    : precoInicial(s0), volInicial(v0), taxaLivreRisco(r),
      kappa(kappa_), theta(theta_), xi(xi_), rho(rho_), dt(1.0/252.0) {
    auto semente = std::random_device{}();
    gerador.seed(semente);
}

ModeloHeston::~ModeloHeston() = default;

void ModeloHeston::definirParametros(double s0, double v0, double r, double kappa_, double theta_, double xi_, double rho_) {
    precoInicial = s0;
    volInicial = v0;
    taxaLivreRisco = r;
    kappa = kappa_;
    theta = theta_;
    xi = xi_;
    rho = rho_;
}

[[nodiscard]] auto ModeloHeston::simularTrajetoria(int numPassos) const -> std::pair<std::vector<double>, std::vector<double>> {
    std::vector<double> precos(numPassos + 1);
    std::vector<double> variancia(numPassos + 1);
    
    precos[0] = precoInicial;
    variancia[0] = volInicial;
    
    std::mt19937_64 gen = gerador;
    std::normal_distribution<double> Z(0.0, 1.0);
    double sqrtDt = std::sqrt(dt);
    
    // dS = r*S*dt + sqrt(v)*S*dW1
    // dv = kappa*(theta-v)*dt + xi*sqrt(v)*dW2
    // corr(dW1, dW2) = rho
    
    for (int t = 1; t <= numPassos; ++t) {
        double Z1 = Z(gen);
        double Z2 = rho * Z1 + std::sqrt(1.0 - rho*rho) * Z(gen);
        
        double v = std::max(variancia[t-1], 0.0);
        double sqrtV = std::sqrt(v);
        
        // Euler para preço
        double dS = taxaLivreRisco * precos[t-1] * dt + sqrtV * precos[t-1] * Z1 * sqrtDt;
        precos[t] = precos[t-1] + dS;
        precos[t] = std::max(precos[t], 0.0);
        
        // Euler para variância (com reflexão para manter positivo)
        double dv = kappa * (theta - v) * dt + xi * sqrtV * Z2 * sqrtDt;
        variancia[t] = std::max(variancia[t-1] + dv, 0.0);
    }
    
    return {precos, variancia};
}

[[nodiscard]] auto ModeloHeston::simularTrajetoriaQE(int numPassos) const -> std::pair<std::vector<double>, std::vector<double>> {
    // Esquema QE (Quadratic Exponential) para melhor precisão
    std::vector<double> precos(numPassos + 1);
    std::vector<double> variancia(numPassos + 1);
    
    precos[0] = precoInicial;
    variancia[0] = volInicial;
    
    std::mt19937_64 gen = gerador;
    std::normal_distribution<double> Z(0.0, 1.0);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    
    for (int t = 1; t <= numPassos; ++t) {
        double v = variancia[t-1];
        
        double m = theta + (v - theta) * std::exp(-kappa * dt);
        double s2 = v * xi * xi * std::exp(-kappa * dt) * (1.0 - std::exp(-kappa * dt)) / kappa +
                    theta * xi * xi * std::pow(1.0 - std::exp(-kappa * dt), 2) / (2.0 * kappa);
        
        double psi = s2 / (m * m);
        double psiC = 1.5;
        
        if (psi <= psiC) {
            double b2 = 2.0 / psi - 1.0 + std::sqrt(2.0 / psi) * std::sqrt(2.0 / psi - 1.0);
            double a = m / (1.0 + b2);
            double b = std::sqrt(b2);
            double Zv = Z(gen);
            variancia[t] = a * std::pow(b + Zv, 2);
        } else {
            double p = (psi - 1.0) / (psi + 1.0);
            double beta = (1.0 - p) / m;
            double u = U(gen);
            variancia[t] = (u <= p) ? 0.0 : std::log((1.0 - p) / (1.0 - u)) / beta;
        }
        
        double k0 = -rho * kappa * theta * dt / xi;
        double k1 = 0.5 * dt * (kappa * rho / xi - 0.5) - rho / xi;
        double k2 = 0.5 * dt * (kappa * rho / xi - 0.5) + rho / xi;
        double k3 = 0.5 * dt * (1.0 - rho * rho);
        double k4 = 0.5 * dt * (1.0 - rho * rho);
        
        double Zs = Z(gen);
        double logS = std::log(precos[t-1]) + taxaLivreRisco * dt + k0 + k1 * v + k2 * variancia[t]
                      + std::sqrt(k3 * v + k4 * variancia[t]) * Zs;
        precos[t] = std::exp(logS);
    }
    
    return {precos, variancia};
}

[[nodiscard]] auto ModeloHeston::simularMultiplasTrajetorias(int numTrajetorias, int numPassos) const -> std::vector<std::pair<std::vector<double>, std::vector<double>>> {
    std::vector<std::pair<std::vector<double>, std::vector<double>>> trajetorias(numTrajetorias);
    for (int i = 0; i < numTrajetorias; ++i) {
        trajetorias[i] = simularTrajetoria(numPassos);
    }
    return trajetorias;
}

[[nodiscard]] auto ModeloHeston::precoCallMonteCarlo(double strike, double tempoExpiracao, int numSimulacoes) const -> double {
    int numPassos = static_cast<int>(tempoExpiracao / dt);
    double soma = 0.0;
    
    for (int i = 0; i < numSimulacoes; ++i) {
        auto [precos, vol] = simularTrajetoria(numPassos);
        double ST = precos.back();
        soma += std::max(ST - strike, 0.0);
    }
    
    return std::exp(-taxaLivreRisco * tempoExpiracao) * soma / numSimulacoes;
}

[[nodiscard]] auto ModeloHeston::verificarCondicaoFeller() const -> bool {
    return 2.0 * kappa * theta > xi * xi;
}

[[nodiscard]] auto ModeloHeston::obterKappa() const noexcept -> double { return kappa; }
[[nodiscard]] auto ModeloHeston::obterTheta() const noexcept -> double { return theta; }
[[nodiscard]] auto ModeloHeston::obterXi() const noexcept -> double { return xi; }
[[nodiscard]] auto ModeloHeston::obterRho() const noexcept -> double { return rho; }

} // namespace analise_ouro
