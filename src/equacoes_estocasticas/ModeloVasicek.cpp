/*
 * ModeloVasicek.cpp - Modelo de Taxas de Juros de Vasicek - C++23
 * Autor: Luiz Tiago Wilcke
 * Projeto: Análise do Preço do Ouro
 */

#include "ModeloVasicek.hpp"
#include <cmath>
#include <random>
#include <algorithm>
#include <numbers>

namespace analise_ouro {

ModeloVasicek::ModeloVasicek() 
    : taxaInicial(0.05), velocidadeReversao(0.5), mediaLongoPrazo(0.05), volatilidade(0.01), dt(1.0/252.0) {
    auto semente = std::random_device{}();
    gerador.seed(semente);
}

ModeloVasicek::ModeloVasicek(double r0, double a, double b, double sigma)
    : taxaInicial(r0), velocidadeReversao(a), mediaLongoPrazo(b), volatilidade(sigma), dt(1.0/252.0) {
    auto semente = std::random_device{}();
    gerador.seed(semente);
}

ModeloVasicek::~ModeloVasicek() = default;

void ModeloVasicek::definirParametros(double r0, double a, double b, double sigma) {
    taxaInicial = r0;
    velocidadeReversao = a;
    mediaLongoPrazo = b;
    volatilidade = sigma;
}

[[nodiscard]] auto ModeloVasicek::simularTrajetoria(int numPassos) const -> std::vector<double> {
    std::vector<double> taxas(numPassos + 1);
    taxas[0] = taxaInicial;
    
    std::mt19937_64 gen = gerador;
    std::normal_distribution<double> Z(0.0, 1.0);
    double sqrtDt = std::sqrt(dt);
    
    // dr = a(b - r)dt + sigma*dW
    for (int t = 1; t <= numPassos; ++t) {
        double dW = Z(gen) * sqrtDt;
        double dr = velocidadeReversao * (mediaLongoPrazo - taxas[t-1]) * dt + volatilidade * dW;
        taxas[t] = taxas[t-1] + dr;
    }
    
    return taxas;
}

[[nodiscard]] auto ModeloVasicek::simularTrajetoriaExata(int numPassos) const -> std::vector<double> {
    std::vector<double> taxas(numPassos + 1);
    taxas[0] = taxaInicial;
    
    std::mt19937_64 gen = gerador;
    std::normal_distribution<double> Z(0.0, 1.0);
    
    double expA = std::exp(-velocidadeReversao * dt);
    double varCondicional = (volatilidade * volatilidade / (2.0 * velocidadeReversao)) * (1.0 - expA * expA);
    double dpCondicional = std::sqrt(varCondicional);
    
    for (int t = 1; t <= numPassos; ++t) {
        double mediaCondicional = taxas[t-1] * expA + mediaLongoPrazo * (1.0 - expA);
        taxas[t] = mediaCondicional + dpCondicional * Z(gen);
    }
    
    return taxas;
}

[[nodiscard]] auto ModeloVasicek::mediaCondicional(double r0, double t) const -> double {
    return r0 * std::exp(-velocidadeReversao * t) + mediaLongoPrazo * (1.0 - std::exp(-velocidadeReversao * t));
}

[[nodiscard]] auto ModeloVasicek::varianciaCondicional(double t) const -> double {
    return (volatilidade * volatilidade / (2.0 * velocidadeReversao)) * (1.0 - std::exp(-2.0 * velocidadeReversao * t));
}

[[nodiscard]] auto ModeloVasicek::precoBondZeroCoupon(double taxaAtual, double maturidade) const -> double {
    double B = (1.0 - std::exp(-velocidadeReversao * maturidade)) / velocidadeReversao;
    double A = std::exp((mediaLongoPrazo - volatilidade * volatilidade / (2.0 * velocidadeReversao * velocidadeReversao)) 
                        * (B - maturidade) - volatilidade * volatilidade * B * B / (4.0 * velocidadeReversao));
    return A * std::exp(-B * taxaAtual);
}

[[nodiscard]] auto ModeloVasicek::taxaForward(double taxaAtual, double t, double T) const -> double {
    double Pt = precoBondZeroCoupon(taxaAtual, t);
    double PT = precoBondZeroCoupon(taxaAtual, T);
    return -std::log(PT / Pt) / (T - t);
}

[[nodiscard]] auto ModeloVasicek::estimarParametros(const std::vector<double>& taxasHistoricas) const -> std::tuple<double, double, double> {
    size_t n = taxasHistoricas.size();
    if (n < 3) return {0.0, 0.0, 0.0};
    
    double somaR = 0.0, somaR2 = 0.0, somaRLag = 0.0, somaRRLag = 0.0;
    for (size_t i = 1; i < n; ++i) {
        somaR += taxasHistoricas[i];
        somaR2 += taxasHistoricas[i] * taxasHistoricas[i];
        somaRLag += taxasHistoricas[i-1];
        somaRRLag += taxasHistoricas[i] * taxasHistoricas[i-1];
    }
    
    double nf = static_cast<double>(n - 1);
    double beta = (somaRRLag - somaR * somaRLag / nf) / (somaRLag * somaRLag / nf - somaR2 + somaR * somaR / nf);
    double aEst = -std::log(beta) / dt;
    double bEst = (somaR / nf - beta * somaRLag / nf) / (1.0 - beta);
    
    double somaRes2 = 0.0;
    for (size_t i = 1; i < n; ++i) {
        double pred = taxasHistoricas[i-1] * beta + bEst * (1.0 - beta);
        double res = taxasHistoricas[i] - pred;
        somaRes2 += res * res;
    }
    double sigmaEst = std::sqrt(somaRes2 / (nf - 2) * 2.0 * aEst / (1.0 - beta * beta));
    
    return {aEst, bEst, sigmaEst};
}

[[nodiscard]] auto ModeloVasicek::obterVelocidadeReversao() const noexcept -> double { return velocidadeReversao; }
[[nodiscard]] auto ModeloVasicek::obterMediaLongoPrazo() const noexcept -> double { return mediaLongoPrazo; }
[[nodiscard]] auto ModeloVasicek::obterVolatilidade() const noexcept -> double { return volatilidade; }

} // namespace analise_ouro
