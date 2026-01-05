/*
 * ModeloCIR.cpp - Modelo Cox-Ingersoll-Ross - C++23
 * Autor: Luiz Tiago Wilcke
 * Projeto: Análise do Preço do Ouro
 */

#include "ModeloCIR.hpp"
#include <cmath>
#include <random>
#include <algorithm>

namespace analise_ouro {

ModeloCIR::ModeloCIR() 
    : taxaInicial(0.05), velocidadeReversao(0.5), mediaLongoPrazo(0.05), volatilidade(0.1), dt(1.0/252.0) {
    auto semente = std::random_device{}();
    gerador.seed(semente);
}

ModeloCIR::ModeloCIR(double r0, double a, double b, double sigma)
    : taxaInicial(r0), velocidadeReversao(a), mediaLongoPrazo(b), volatilidade(sigma), dt(1.0/252.0) {
    auto semente = std::random_device{}();
    gerador.seed(semente);
}

ModeloCIR::~ModeloCIR() = default;

void ModeloCIR::definirParametros(double r0, double a, double b, double sigma) {
    taxaInicial = r0;
    velocidadeReversao = a;
    mediaLongoPrazo = b;
    volatilidade = sigma;
}

[[nodiscard]] auto ModeloCIR::simularTrajetoria(int numPassos) const -> std::vector<double> {
    std::vector<double> taxas(numPassos + 1);
    taxas[0] = taxaInicial;
    
    std::mt19937_64 gen = gerador;
    std::normal_distribution<double> Z(0.0, 1.0);
    double sqrtDt = std::sqrt(dt);
    
    // dr = a(b - r)dt + sigma*sqrt(r)*dW
    for (int t = 1; t <= numPassos; ++t) {
        double dW = Z(gen) * sqrtDt;
        double r = std::max(taxas[t-1], 0.0);
        double dr = velocidadeReversao * (mediaLongoPrazo - r) * dt + volatilidade * std::sqrt(r) * dW;
        taxas[t] = std::max(r + dr, 0.0);
    }
    
    return taxas;
}

[[nodiscard]] auto ModeloCIR::simularTrajetoriaExata(int numPassos) const -> std::vector<double> {
    std::vector<double> taxas(numPassos + 1);
    taxas[0] = taxaInicial;
    
    std::mt19937_64 gen = gerador;
    
    double expA = std::exp(-velocidadeReversao * dt);
    double c = volatilidade * volatilidade * (1.0 - expA) / (4.0 * velocidadeReversao);
    double d = 4.0 * velocidadeReversao * mediaLongoPrazo / (volatilidade * volatilidade);
    
    for (int t = 1; t <= numPassos; ++t) {
        double lambda = 4.0 * velocidadeReversao * expA * taxas[t-1] / (volatilidade * volatilidade * (1.0 - expA));
        
        // Chi-quadrado não-central
        std::chi_squared_distribution<double> chiSq(d);
        std::poisson_distribution<int> poisson(lambda / 2.0);
        
        int N = poisson(gen);
        double chiVal = chiSq(gen);
        
        std::chi_squared_distribution<double> chiSq2(2.0 * N);
        if (N > 0) chiVal += chiSq2(gen);
        
        taxas[t] = c * chiVal;
    }
    
    return taxas;
}

[[nodiscard]] auto ModeloCIR::verificarCondicaoFeller() const -> bool {
    return 2.0 * velocidadeReversao * mediaLongoPrazo >= volatilidade * volatilidade;
}

[[nodiscard]] auto ModeloCIR::mediaCondicional(double r0, double t) const -> double {
    return r0 * std::exp(-velocidadeReversao * t) + mediaLongoPrazo * (1.0 - std::exp(-velocidadeReversao * t));
}

[[nodiscard]] auto ModeloCIR::varianciaCondicional(double r0, double t) const -> double {
    double expA = std::exp(-velocidadeReversao * t);
    return r0 * volatilidade * volatilidade * expA * (1.0 - expA) / velocidadeReversao +
           mediaLongoPrazo * volatilidade * volatilidade * (1.0 - expA) * (1.0 - expA) / (2.0 * velocidadeReversao);
}

[[nodiscard]] auto ModeloCIR::precoBondZeroCoupon(double taxaAtual, double maturidade) const -> double {
    double h = std::sqrt(velocidadeReversao * velocidadeReversao + 2.0 * volatilidade * volatilidade);
    double expHT = std::exp(h * maturidade);
    
    double A = std::pow(2.0 * h * std::exp((velocidadeReversao + h) * maturidade / 2.0) /
                        ((velocidadeReversao + h) * (expHT - 1.0) + 2.0 * h),
                        2.0 * velocidadeReversao * mediaLongoPrazo / (volatilidade * volatilidade));
    double B = 2.0 * (expHT - 1.0) / ((velocidadeReversao + h) * (expHT - 1.0) + 2.0 * h);
    
    return A * std::exp(-B * taxaAtual);
}

[[nodiscard]] auto ModeloCIR::obterVelocidadeReversao() const noexcept -> double { return velocidadeReversao; }
[[nodiscard]] auto ModeloCIR::obterMediaLongoPrazo() const noexcept -> double { return mediaLongoPrazo; }
[[nodiscard]] auto ModeloCIR::obterVolatilidade() const noexcept -> double { return volatilidade; }

} // namespace analise_ouro
