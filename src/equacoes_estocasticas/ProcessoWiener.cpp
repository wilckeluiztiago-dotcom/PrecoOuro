/*
 * ProcessoWiener.cpp - Processo de Wiener e Cálculo de Itô - C++23
 * Autor: Luiz Tiago Wilcke
 * Projeto: Análise do Preço do Ouro
 */

#include "ProcessoWiener.hpp"
#include <cmath>
#include <random>
#include <stdexcept>
#include <algorithm>

namespace analise_ouro {

ProcessoWiener::ProcessoWiener() : dt(1.0/252.0) {
    auto semente = std::random_device{}();
    gerador.seed(semente);
}

ProcessoWiener::ProcessoWiener(double deltaT) : dt(deltaT) {
    auto semente = std::random_device{}();
    gerador.seed(semente);
}

ProcessoWiener::~ProcessoWiener() = default;

void ProcessoWiener::definirPassoTempo(double deltaT) {
    if (deltaT <= 0) throw std::invalid_argument("dt deve ser positivo");
    dt = deltaT;
}

void ProcessoWiener::definirSemente(unsigned long semente) {
    gerador.seed(semente);
}

[[nodiscard]] auto ProcessoWiener::simularTrajetoria(int numPassos) const -> std::vector<double> {
    std::vector<double> W(numPassos + 1, 0.0);
    
    std::mt19937_64 gen = gerador;
    std::normal_distribution<double> Z(0.0, 1.0);
    double sqrtDt = std::sqrt(dt);
    
    for (int t = 1; t <= numPassos; ++t) {
        W[t] = W[t-1] + Z(gen) * sqrtDt;
    }
    
    return W;
}

[[nodiscard]] auto ProcessoWiener::simularIncrementos(int numPassos) const -> std::vector<double> {
    std::vector<double> dW(numPassos);
    
    std::mt19937_64 gen = gerador;
    std::normal_distribution<double> Z(0.0, 1.0);
    double sqrtDt = std::sqrt(dt);
    
    for (int t = 0; t < numPassos; ++t) {
        dW[t] = Z(gen) * sqrtDt;
    }
    
    return dW;
}

[[nodiscard]] auto ProcessoWiener::simularMultiplasTrajetorias(int numTrajetorias, int numPassos) const -> std::vector<std::vector<double>> {
    std::vector<std::vector<double>> trajetorias(numTrajetorias);
    for (int i = 0; i < numTrajetorias; ++i) {
        trajetorias[i] = simularTrajetoria(numPassos);
    }
    return trajetorias;
}

[[nodiscard]] auto ProcessoWiener::simularWienerCorrelacionado(int numPassos, double correlacao) const -> std::pair<std::vector<double>, std::vector<double>> {
    auto W1 = simularTrajetoria(numPassos);
    auto Z = simularIncrementos(numPassos);
    
    std::vector<double> W2(numPassos + 1, 0.0);
    double sqrtDt = std::sqrt(dt);
    
    for (int t = 1; t <= numPassos; ++t) {
        double dW1 = W1[t] - W1[t-1];
        double dZ = Z[t-1];
        double dW2 = correlacao * dW1 + std::sqrt(1.0 - correlacao*correlacao) * dZ;
        W2[t] = W2[t-1] + dW2;
    }
    
    return {W1, W2};
}

[[nodiscard]] auto ProcessoWiener::calcularVarianciaQuadratica(const std::vector<double>& W) const -> double {
    double soma = 0.0;
    for (size_t i = 1; i < W.size(); ++i) {
        double dW = W[i] - W[i-1];
        soma += dW * dW;
    }
    return soma;
}

[[nodiscard]] auto ProcessoWiener::integralIto(const std::vector<double>& f, const std::vector<double>& W) const -> double {
    if (f.size() != W.size()) throw std::invalid_argument("Tamanhos incompatíveis");
    
    double integral = 0.0;
    for (size_t i = 0; i < W.size() - 1; ++i) {
        double dW = W[i+1] - W[i];
        integral += f[i] * dW;
    }
    return integral;
}

[[nodiscard]] auto ProcessoWiener::integralStratonovich(const std::vector<double>& f, const std::vector<double>& W) const -> double {
    if (f.size() != W.size()) throw std::invalid_argument("Tamanhos incompatíveis");
    
    double integral = 0.0;
    for (size_t i = 0; i < W.size() - 1; ++i) {
        double dW = W[i+1] - W[i];
        double fMedio = (f[i] + f[i+1]) / 2.0;
        integral += fMedio * dW;
    }
    return integral;
}

[[nodiscard]] auto ProcessoWiener::aplicarLemaIto(const std::vector<double>& S, FuncaoDerivada df, FuncaoDerivada d2f) const -> std::vector<double> {
    std::vector<double> dF(S.size() - 1);
    
    for (size_t i = 0; i < S.size() - 1; ++i) {
        double dfS = df(S[i]);
        double d2fS = d2f(S[i]);
        double dS = S[i+1] - S[i];
        
        // Lema de Itô: dF = df/dS * dS + 0.5 * d²f/dS² * (dS)² (termo de volatilidade)
        dF[i] = dfS * dS + 0.5 * d2fS * dS * dS;
    }
    
    return dF;
}

[[nodiscard]] auto ProcessoWiener::calcularCovarianciaEmpírica(const std::vector<double>& W1, const std::vector<double>& W2, int passo) const -> double {
    if (static_cast<size_t>(passo) >= W1.size() || static_cast<size_t>(passo) >= W2.size()) {
        throw std::invalid_argument("Passo fora do range");
    }
    return W1[passo] * W2[passo];
}

[[nodiscard]] auto ProcessoWiener::propriedadeMartingale(const std::vector<std::vector<double>>& trajetorias, int passo) const -> double {
    // E[W(t+s) | W(t)] = W(t) - verifica propriedade martingale
    double soma = 0.0;
    for (const auto& traj : trajetorias) {
        if (static_cast<size_t>(passo) < traj.size()) {
            soma += traj[passo];
        }
    }
    return soma / trajetorias.size();
}

[[nodiscard]] auto ProcessoWiener::obterPassoTempo() const noexcept -> double { 
    return dt; 
}

} // namespace analise_ouro
