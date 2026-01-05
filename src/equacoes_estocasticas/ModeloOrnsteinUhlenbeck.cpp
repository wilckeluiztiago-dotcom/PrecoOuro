/*
 * ModeloOrnsteinUhlenbeck.cpp - Processo de Ornstein-Uhlenbeck - C++23
 * Autor: Luiz Tiago Wilcke
 * Projeto: Análise do Preço do Ouro
 */

#include "ModeloOrnsteinUhlenbeck.hpp"
#include <cmath>
#include <random>
#include <algorithm>
#include <numeric>

namespace analise_ouro {

ModeloOrnsteinUhlenbeck::ModeloOrnsteinUhlenbeck() 
    : valorInicial(0.0), velocidadeReversao(1.0), mediaLongoPrazo(0.0), volatilidade(0.1), dt(1.0/252.0) {
    auto semente = std::random_device{}();
    gerador.seed(semente);
}

ModeloOrnsteinUhlenbeck::ModeloOrnsteinUhlenbeck(double x0, double theta, double mu, double sigma)
    : valorInicial(x0), velocidadeReversao(theta), mediaLongoPrazo(mu), volatilidade(sigma), dt(1.0/252.0) {
    auto semente = std::random_device{}();
    gerador.seed(semente);
}

ModeloOrnsteinUhlenbeck::~ModeloOrnsteinUhlenbeck() = default;

void ModeloOrnsteinUhlenbeck::definirParametros(double x0, double theta, double mu, double sigma) {
    valorInicial = x0;
    velocidadeReversao = theta;
    mediaLongoPrazo = mu;
    volatilidade = sigma;
}

[[nodiscard]] auto ModeloOrnsteinUhlenbeck::simularTrajetoria(int numPassos) const -> std::vector<double> {
    std::vector<double> trajetoria(numPassos + 1);
    trajetoria[0] = valorInicial;
    
    std::mt19937_64 gen = gerador;
    std::normal_distribution<double> Z(0.0, 1.0);
    double sqrtDt = std::sqrt(dt);
    
    // dX = theta*(mu - X)dt + sigma*dW
    for (int t = 1; t <= numPassos; ++t) {
        double dW = Z(gen) * sqrtDt;
        double dX = velocidadeReversao * (mediaLongoPrazo - trajetoria[t-1]) * dt + volatilidade * dW;
        trajetoria[t] = trajetoria[t-1] + dX;
    }
    
    return trajetoria;
}

[[nodiscard]] auto ModeloOrnsteinUhlenbeck::simularTrajetoriaExata(int numPassos) const -> std::vector<double> {
    std::vector<double> trajetoria(numPassos + 1);
    trajetoria[0] = valorInicial;
    
    std::mt19937_64 gen = gerador;
    std::normal_distribution<double> Z(0.0, 1.0);
    
    double expTheta = std::exp(-velocidadeReversao * dt);
    double varCondicional = (volatilidade * volatilidade / (2.0 * velocidadeReversao)) * (1.0 - expTheta * expTheta);
    double dpCondicional = std::sqrt(varCondicional);
    
    for (int t = 1; t <= numPassos; ++t) {
        double mediaCondicionalT = trajetoria[t-1] * expTheta + mediaLongoPrazo * (1.0 - expTheta);
        trajetoria[t] = mediaCondicionalT + dpCondicional * Z(gen);
    }
    
    return trajetoria;
}

[[nodiscard]] auto ModeloOrnsteinUhlenbeck::simularMultiplasTrajetorias(int numTrajetorias, int numPassos) const -> std::vector<std::vector<double>> {
    std::vector<std::vector<double>> trajetorias(numTrajetorias);
    for (int i = 0; i < numTrajetorias; ++i) {
        trajetorias[i] = simularTrajetoriaExata(numPassos);
    }
    return trajetorias;
}

[[nodiscard]] auto ModeloOrnsteinUhlenbeck::mediaCondicional(double x0, double t) const -> double {
    return x0 * std::exp(-velocidadeReversao * t) + mediaLongoPrazo * (1.0 - std::exp(-velocidadeReversao * t));
}

[[nodiscard]] auto ModeloOrnsteinUhlenbeck::varianciaCondicional(double t) const -> double {
    return (volatilidade * volatilidade / (2.0 * velocidadeReversao)) * (1.0 - std::exp(-2.0 * velocidadeReversao * t));
}

[[nodiscard]] auto ModeloOrnsteinUhlenbeck::varianciaEstacionaria() const -> double {
    return volatilidade * volatilidade / (2.0 * velocidadeReversao);
}

[[nodiscard]] auto ModeloOrnsteinUhlenbeck::autoCorrelacao(double lag) const -> double {
    return std::exp(-velocidadeReversao * lag);
}

[[nodiscard]] auto ModeloOrnsteinUhlenbeck::meiaVida() const -> double {
    return std::log(2.0) / velocidadeReversao;
}

[[nodiscard]] auto ModeloOrnsteinUhlenbeck::estimarParametros(const std::vector<double>& dados) const -> std::tuple<double, double, double> {
    size_t n = dados.size();
    if (n < 3) return {0.0, 0.0, 0.0};
    
    double somaX = 0.0, somaX2 = 0.0, somaXLag = 0.0, somaXXLag = 0.0;
    for (size_t i = 1; i < n; ++i) {
        somaX += dados[i];
        somaX2 += dados[i] * dados[i];
        somaXLag += dados[i-1];
        somaXXLag += dados[i] * dados[i-1];
    }
    
    double nf = static_cast<double>(n - 1);
    double mediaX = somaX / nf;
    double mediaXLag = somaXLag / nf;
    
    double rho = (somaXXLag - nf * mediaX * mediaXLag) / (somaX2 - nf * mediaX * mediaX);
    rho = std::clamp(rho, -0.999, 0.999);
    
    double thetaEst = -std::log(rho) / dt;
    double muEst = mediaX * (1.0 - rho) / (1.0 - rho);
    
    double somaRes2 = 0.0;
    for (size_t i = 1; i < n; ++i) {
        double pred = dados[i-1] * rho + muEst * (1.0 - rho);
        double res = dados[i] - pred;
        somaRes2 += res * res;
    }
    double sigmaEst = std::sqrt(somaRes2 / (nf - 2) * 2.0 * thetaEst / (1.0 - rho * rho));
    
    return {thetaEst, muEst, sigmaEst};
}

[[nodiscard]] auto ModeloOrnsteinUhlenbeck::obterVelocidadeReversao() const noexcept -> double { return velocidadeReversao; }
[[nodiscard]] auto ModeloOrnsteinUhlenbeck::obterMediaLongoPrazo() const noexcept -> double { return mediaLongoPrazo; }
[[nodiscard]] auto ModeloOrnsteinUhlenbeck::obterVolatilidade() const noexcept -> double { return volatilidade; }

} // namespace analise_ouro
