/*
 * ModeloBlackScholes.cpp - Modelo Black-Scholes para opções - C++23
 * Autor: Luiz Tiago Wilcke
 * Projeto: Análise do Preço do Ouro
 */

#include "ModeloBlackScholes.hpp"
#include <cmath>
#include <random>
#include <stdexcept>
#include <algorithm>
#include <numbers>

namespace analise_ouro {

ModeloBlackScholes::ModeloBlackScholes() 
    : precoAtivo(1800.0), taxaLivreRisco(0.05), volatilidade(0.2), dt(1.0/252.0) {
    auto semente = std::random_device{}();
    gerador.seed(semente);
}

ModeloBlackScholes::ModeloBlackScholes(double s0, double r, double sigma)
    : precoAtivo(s0), taxaLivreRisco(r), volatilidade(sigma), dt(1.0/252.0) {
    auto semente = std::random_device{}();
    gerador.seed(semente);
}

ModeloBlackScholes::~ModeloBlackScholes() = default;

void ModeloBlackScholes::definirParametros(double s0, double r, double sigma) {
    precoAtivo = s0;
    taxaLivreRisco = r;
    volatilidade = sigma;
}

[[nodiscard]] auto ModeloBlackScholes::cdfNormal(double x) const -> double {
    return 0.5 * std::erfc(-x / std::numbers::sqrt2);
}

[[nodiscard]] auto ModeloBlackScholes::pdfNormal(double x) const -> double {
    return std::exp(-0.5 * x * x) / std::sqrt(2.0 * std::numbers::pi);
}

[[nodiscard]] auto ModeloBlackScholes::calcularD1(double K, double T) const -> double {
    return (std::log(precoAtivo / K) + (taxaLivreRisco + 0.5 * volatilidade * volatilidade) * T) / (volatilidade * std::sqrt(T));
}

[[nodiscard]] auto ModeloBlackScholes::calcularD2(double K, double T) const -> double {
    return calcularD1(K, T) - volatilidade * std::sqrt(T);
}

[[nodiscard]] auto ModeloBlackScholes::precoCall(double strike, double tempoExpiracao) const -> double {
    double d1 = calcularD1(strike, tempoExpiracao);
    double d2 = calcularD2(strike, tempoExpiracao);
    
    return precoAtivo * cdfNormal(d1) - strike * std::exp(-taxaLivreRisco * tempoExpiracao) * cdfNormal(d2);
}

[[nodiscard]] auto ModeloBlackScholes::precoPut(double strike, double tempoExpiracao) const -> double {
    double d1 = calcularD1(strike, tempoExpiracao);
    double d2 = calcularD2(strike, tempoExpiracao);
    
    return strike * std::exp(-taxaLivreRisco * tempoExpiracao) * cdfNormal(-d2) - precoAtivo * cdfNormal(-d1);
}

[[nodiscard]] auto ModeloBlackScholes::delta(double strike, double tempoExpiracao, bool ehCall) const -> double {
    double d1 = calcularD1(strike, tempoExpiracao);
    return ehCall ? cdfNormal(d1) : cdfNormal(d1) - 1.0;
}

[[nodiscard]] auto ModeloBlackScholes::gamma(double strike, double tempoExpiracao) const -> double {
    double d1 = calcularD1(strike, tempoExpiracao);
    return pdfNormal(d1) / (precoAtivo * volatilidade * std::sqrt(tempoExpiracao));
}

[[nodiscard]] auto ModeloBlackScholes::vega(double strike, double tempoExpiracao) const -> double {
    double d1 = calcularD1(strike, tempoExpiracao);
    return precoAtivo * pdfNormal(d1) * std::sqrt(tempoExpiracao);
}

[[nodiscard]] auto ModeloBlackScholes::theta(double strike, double tempoExpiracao, bool ehCall) const -> double {
    double d1 = calcularD1(strike, tempoExpiracao);
    double d2 = calcularD2(strike, tempoExpiracao);
    
    double termo1 = -precoAtivo * pdfNormal(d1) * volatilidade / (2.0 * std::sqrt(tempoExpiracao));
    double desconto = std::exp(-taxaLivreRisco * tempoExpiracao);
    
    if (ehCall) {
        return termo1 - taxaLivreRisco * strike * desconto * cdfNormal(d2);
    } else {
        return termo1 + taxaLivreRisco * strike * desconto * cdfNormal(-d2);
    }
}

[[nodiscard]] auto ModeloBlackScholes::rho(double strike, double tempoExpiracao, bool ehCall) const -> double {
    double d2 = calcularD2(strike, tempoExpiracao);
    double desconto = std::exp(-taxaLivreRisco * tempoExpiracao);
    
    if (ehCall) {
        return strike * tempoExpiracao * desconto * cdfNormal(d2);
    } else {
        return -strike * tempoExpiracao * desconto * cdfNormal(-d2);
    }
}

[[nodiscard]] auto ModeloBlackScholes::volatidadeImplicita(double precoOpcao, double strike, double tempoExpiracao, bool ehCall) const -> double {
    // Newton-Raphson
    double sigma = 0.2;
    const int maxIter = 100;
    const double tol = 1e-6;
    
    ModeloBlackScholes modelo(precoAtivo, taxaLivreRisco, sigma);
    
    for (int i = 0; i < maxIter; ++i) {
        modelo.definirParametros(precoAtivo, taxaLivreRisco, sigma);
        double preco = ehCall ? modelo.precoCall(strike, tempoExpiracao) : modelo.precoPut(strike, tempoExpiracao);
        double vegaVal = modelo.vega(strike, tempoExpiracao);
        
        if (std::abs(vegaVal) < 1e-10) break;
        
        double diff = preco - precoOpcao;
        if (std::abs(diff) < tol) return sigma;
        
        sigma -= diff / vegaVal;
        sigma = std::clamp(sigma, 0.001, 5.0);
    }
    
    return sigma;
}

[[nodiscard]] auto ModeloBlackScholes::simularTrajetoriaPreco(int numPassos) const -> std::vector<double> {
    std::vector<double> precos(numPassos + 1);
    precos[0] = precoAtivo;
    
    std::mt19937_64 gen = gerador;
    std::normal_distribution<double> Z(0.0, 1.0);
    double sqrtDt = std::sqrt(dt);
    
    for (int t = 1; t <= numPassos; ++t) {
        double dW = Z(gen) * sqrtDt;
        double expoente = (taxaLivreRisco - 0.5 * volatilidade * volatilidade) * dt + volatilidade * dW;
        precos[t] = precos[t-1] * std::exp(expoente);
    }
    
    return precos;
}

[[nodiscard]] auto ModeloBlackScholes::precoCallMonteCarlo(double strike, double tempoExpiracao, int numSimulacoes) const -> double {
    int numPassos = static_cast<int>(tempoExpiracao / dt);
    double soma = 0.0;
    
    for (int i = 0; i < numSimulacoes; ++i) {
        auto traj = simularTrajetoriaPreco(numPassos);
        double ST = traj.back();
        soma += std::max(ST - strike, 0.0);
    }
    
    return std::exp(-taxaLivreRisco * tempoExpiracao) * soma / numSimulacoes;
}

[[nodiscard]] auto ModeloBlackScholes::precoPutMonteCarlo(double strike, double tempoExpiracao, int numSimulacoes) const -> double {
    int numPassos = static_cast<int>(tempoExpiracao / dt);
    double soma = 0.0;
    
    for (int i = 0; i < numSimulacoes; ++i) {
        auto traj = simularTrajetoriaPreco(numPassos);
        double ST = traj.back();
        soma += std::max(strike - ST, 0.0);
    }
    
    return std::exp(-taxaLivreRisco * tempoExpiracao) * soma / numSimulacoes;
}

} // namespace analise_ouro
