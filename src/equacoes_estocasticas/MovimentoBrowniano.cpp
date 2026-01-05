/*
 * MovimentoBrowniano.cpp - Movimento Browniano Geométrico - C++23
 * Autor: Luiz Tiago Wilcke
 * Projeto: Análise do Preço do Ouro
 */

#include "MovimentoBrowniano.hpp"
#include <cmath>
#include <random>
#include <stdexcept>
#include <numbers>

namespace analise_ouro {

MovimentoBrowniano::MovimentoBrowniano() : precoInicial(100.0), drift(0.05), volatilidade(0.2), dt(1.0/252.0) {
    auto semente = std::random_device{}();
    gerador.seed(semente);
}

MovimentoBrowniano::MovimentoBrowniano(double s0, double mu, double sigma, double deltaT)
    : precoInicial(s0), drift(mu), volatilidade(sigma), dt(deltaT) {
    auto semente = std::random_device{}();
    gerador.seed(semente);
}

MovimentoBrowniano::~MovimentoBrowniano() = default;

void MovimentoBrowniano::definirParametros(double s0, double mu, double sigma) {
    precoInicial = s0;
    drift = mu;
    volatilidade = sigma;
}

void MovimentoBrowniano::definirPassoTempo(double deltaT) {
    if (deltaT <= 0) throw std::invalid_argument("dt deve ser positivo");
    dt = deltaT;
}

void MovimentoBrowniano::definirSemente(unsigned long semente) {
    gerador.seed(semente);
}

[[nodiscard]] auto MovimentoBrowniano::gerarNormalPadrao() -> double {
    std::normal_distribution<double> dist(0.0, 1.0);
    return dist(gerador);
}

[[nodiscard]] auto MovimentoBrowniano::simularTrajetoria(int numPassos) const -> std::vector<double> {
    std::vector<double> trajetoria(numPassos + 1);
    trajetoria[0] = precoInicial;
    
    std::mt19937_64 gen = gerador;
    std::normal_distribution<double> Z(0.0, 1.0);
    
    double sqrtDt = std::sqrt(dt);
    
    // dS = mu*S*dt + sigma*S*dW (GBM)
    for (int t = 1; t <= numPassos; ++t) {
        double dW = Z(gen) * sqrtDt;
        double dS = drift * trajetoria[t-1] * dt + volatilidade * trajetoria[t-1] * dW;
        trajetoria[t] = trajetoria[t-1] + dS;
        trajetoria[t] = std::max(trajetoria[t], 0.0);
    }
    
    return trajetoria;
}

[[nodiscard]] auto MovimentoBrowniano::simularTrajetoriaExata(int numPassos) const -> std::vector<double> {
    std::vector<double> trajetoria(numPassos + 1);
    trajetoria[0] = precoInicial;
    
    std::mt19937_64 gen = gerador;
    std::normal_distribution<double> Z(0.0, 1.0);
    
    double sqrtDt = std::sqrt(dt);
    
    // Solução exata: S(t) = S(0) * exp((mu - sigma²/2)*t + sigma*W(t))
    for (int t = 1; t <= numPassos; ++t) {
        double dW = Z(gen) * sqrtDt;
        double expoente = (drift - 0.5 * volatilidade * volatilidade) * dt + volatilidade * dW;
        trajetoria[t] = trajetoria[t-1] * std::exp(expoente);
    }
    
    return trajetoria;
}

[[nodiscard]] auto MovimentoBrowniano::simularMultiplasTrajetorias(int numTrajetorias, int numPassos) const -> std::vector<std::vector<double>> {
    std::vector<std::vector<double>> trajetorias(numTrajetorias);
    
    for (int i = 0; i < numTrajetorias; ++i) {
        trajetorias[i] = simularTrajetoriaExata(numPassos);
    }
    
    return trajetorias;
}

[[nodiscard]] auto MovimentoBrowniano::calcularMediaEstatistica(const std::vector<std::vector<double>>& trajetorias, int passo) const -> double {
    double soma = 0.0;
    for (const auto& traj : trajetorias) {
        if (static_cast<size_t>(passo) < traj.size()) {
            soma += traj[passo];
        }
    }
    return soma / trajetorias.size();
}

[[nodiscard]] auto MovimentoBrowniano::calcularVarianciaEstatistica(const std::vector<std::vector<double>>& trajetorias, int passo) const -> double {
    double media = calcularMediaEstatistica(trajetorias, passo);
    double soma = 0.0;
    for (const auto& traj : trajetorias) {
        if (static_cast<size_t>(passo) < traj.size()) {
            double diff = traj[passo] - media;
            soma += diff * diff;
        }
    }
    return soma / (trajetorias.size() - 1);
}

[[nodiscard]] auto MovimentoBrowniano::valorEsperadoTeorico(double tempo) const -> double {
    return precoInicial * std::exp(drift * tempo);
}

[[nodiscard]] auto MovimentoBrowniano::varianciaTeorica(double tempo) const -> double {
    double E_S2 = precoInicial * precoInicial * std::exp((2*drift + volatilidade*volatilidade) * tempo);
    double E_S = valorEsperadoTeorico(tempo);
    return E_S2 - E_S * E_S;
}

[[nodiscard]] auto MovimentoBrowniano::calcularRetornosLogaritmicos(const std::vector<double>& trajetoria) const -> std::vector<double> {
    std::vector<double> retornos(trajetoria.size() - 1);
    for (size_t i = 1; i < trajetoria.size(); ++i) {
        if (trajetoria[i-1] > 0 && trajetoria[i] > 0) {
            retornos[i-1] = std::log(trajetoria[i] / trajetoria[i-1]);
        } else {
            retornos[i-1] = 0.0;
        }
    }
    return retornos;
}

[[nodiscard]] auto MovimentoBrowniano::estimarParametros(const std::vector<double>& precos) const -> std::pair<double, double> {
    auto retornos = calcularRetornosLogaritmicos(precos);
    if (retornos.empty()) return {0.0, 0.0};
    
    double soma = 0.0;
    for (const auto& r : retornos) soma += r;
    double mediaRetorno = soma / retornos.size();
    
    double varRetorno = 0.0;
    for (const auto& r : retornos) {
        double d = r - mediaRetorno;
        varRetorno += d * d;
    }
    varRetorno /= (retornos.size() - 1);
    
    double sigmaEst = std::sqrt(varRetorno / dt);
    double muEst = mediaRetorno / dt + 0.5 * sigmaEst * sigmaEst;
    
    return {muEst, sigmaEst};
}

[[nodiscard]] auto MovimentoBrowniano::obterDrift() const noexcept -> double { return drift; }
[[nodiscard]] auto MovimentoBrowniano::obterVolatilidade() const noexcept -> double { return volatilidade; }
[[nodiscard]] auto MovimentoBrowniano::obterPrecoInicial() const noexcept -> double { return precoInicial; }

} // namespace analise_ouro
