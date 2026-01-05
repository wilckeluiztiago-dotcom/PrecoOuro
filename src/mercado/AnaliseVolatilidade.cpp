/*
 * AnaliseVolatilidade.cpp - Análise de Volatilidade do Ouro - C++23
 * Autor: Luiz Tiago Wilcke
 * Projeto: Análise do Preço do Ouro
 */

#include "AnaliseVolatilidade.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace analise_ouro {

AnaliseVolatilidade::AnaliseVolatilidade() : serieCarregada(false), tamanhoSerie(0) {}

AnaliseVolatilidade::AnaliseVolatilidade(const std::vector<double>& precos) { carregarPrecos(precos); }

AnaliseVolatilidade::~AnaliseVolatilidade() = default;

void AnaliseVolatilidade::carregarPrecos(const std::vector<double>& precos) {
    if (precos.size() < 2) throw std::invalid_argument("Precisa de pelo menos 2 preços");
    precosOriginais = precos;
    tamanhoSerie = precos.size();
    calcularRetornos();
    serieCarregada = true;
}

void AnaliseVolatilidade::calcularRetornos() {
    retornosLog.resize(tamanhoSerie - 1);
    for (size_t i = 1; i < tamanhoSerie; ++i) {
        retornosLog[i - 1] = std::log(precosOriginais[i] / precosOriginais[i - 1]);
    }
}

[[nodiscard]] auto AnaliseVolatilidade::volatilidadeHistorica(int janela) const -> std::vector<double> {
    verificarDados();
    if (janela < 2) janela = 20;
    
    std::vector<double> volatilidades(retornosLog.size() - janela + 1);
    double fatorAnualizacao = std::sqrt(252.0);
    
    for (size_t i = 0; i <= retornosLog.size() - janela; ++i) {
        double soma = 0.0, soma2 = 0.0;
        for (int j = 0; j < janela; ++j) {
            soma += retornosLog[i + j];
            soma2 += retornosLog[i + j] * retornosLog[i + j];
        }
        double media = soma / janela;
        double var = (soma2 - janela * media * media) / (janela - 1);
        volatilidades[i] = std::sqrt(var) * fatorAnualizacao;
    }
    return volatilidades;
}

[[nodiscard]] auto AnaliseVolatilidade::volatilidadeRealizada(int janela) const -> double {
    verificarDados();
    
    double soma2 = 0.0;
    int n = std::min(janela, static_cast<int>(retornosLog.size()));
    int inicio = retornosLog.size() - n;
    
    for (size_t i = inicio; i < retornosLog.size(); ++i) {
        soma2 += retornosLog[i] * retornosLog[i];
    }
    
    return std::sqrt(soma2 / n * 252.0);
}

[[nodiscard]] auto AnaliseVolatilidade::volatilidadeEWMA(double lambda) const -> std::vector<double> {
    verificarDados();
    
    std::vector<double> volatilidades(retornosLog.size());
    double var = retornosLog[0] * retornosLog[0];
    volatilidades[0] = std::sqrt(var * 252.0);
    
    for (size_t i = 1; i < retornosLog.size(); ++i) {
        var = lambda * var + (1.0 - lambda) * retornosLog[i - 1] * retornosLog[i - 1];
        volatilidades[i] = std::sqrt(var * 252.0);
    }
    return volatilidades;
}

[[nodiscard]] auto AnaliseVolatilidade::volatilidadeParkinson(int janela) const -> double {
    verificarDados();
    if (tamanhoSerie < 4) return 0.0;
    
    // Aproximação usando preços de fechamento
    double soma = 0.0;
    int n = std::min(janela, static_cast<int>(tamanhoSerie) - 1);
    
    for (size_t i = tamanhoSerie - n; i < tamanhoSerie; ++i) {
        double logRatio = std::log(precosOriginais[i] / precosOriginais[i - 1]);
        soma += logRatio * logRatio;
    }
    
    return std::sqrt(soma / (4.0 * n * std::log(2.0)) * 252.0);
}

[[nodiscard]] auto AnaliseVolatilidade::volatilidadeGarmanKlass() const -> double {
    verificarDados();
    
    // Usando preços de fechamento como aproximação
    double soma = 0.0;
    for (size_t i = 1; i < tamanhoSerie; ++i) {
        double u = std::log(precosOriginais[i] / precosOriginais[i - 1]);
        soma += u * u;
    }
    
    return std::sqrt(soma / (tamanhoSerie - 1) * 252.0);
}

[[nodiscard]] auto AnaliseVolatilidade::detectarClusteringVolatilidade() const -> std::vector<int> {
    verificarDados();
    
    auto vols = volatilidadeHistorica(20);
    if (vols.size() < 10) return {};
    
    double mediaVol = std::accumulate(vols.begin(), vols.end(), 0.0) / vols.size();
    double dpVol = 0.0;
    for (const auto& v : vols) {
        dpVol += (v - mediaVol) * (v - mediaVol);
    }
    dpVol = std::sqrt(dpVol / (vols.size() - 1));
    
    std::vector<int> clusters;
    double limiar = mediaVol + 1.5 * dpVol;
    
    for (size_t i = 0; i < vols.size(); ++i) {
        if (vols[i] > limiar) {
            clusters.push_back(static_cast<int>(i));
        }
    }
    return clusters;
}

[[nodiscard]] auto AnaliseVolatilidade::calcularConeVolatilidade(const std::vector<int>& janelas) const -> std::vector<std::tuple<int, double, double, double>> {
    verificarDados();
    
    std::vector<std::tuple<int, double, double, double>> cone;
    
    for (int janela : janelas) {
        auto vols = volatilidadeHistorica(janela);
        if (vols.empty()) continue;
        
        std::sort(vols.begin(), vols.end());
        double minimo = vols.front();
        double maximo = vols.back();
        double mediana = vols[vols.size() / 2];
        
        cone.emplace_back(janela, minimo, mediana, maximo);
    }
    return cone;
}

[[nodiscard]] auto AnaliseVolatilidade::calcularVolatilidadeImplicita(double precoOpcao, double precoAtivo, double strike, double r, double T, bool ehCall) const -> double {
    double sigma = 0.2;
    const int maxIter = 100;
    const double tol = 1e-6;
    
    for (int i = 0; i < maxIter; ++i) {
        double d1 = (std::log(precoAtivo / strike) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
        double d2 = d1 - sigma * std::sqrt(T);
        
        double N_d1 = 0.5 * std::erfc(-d1 / std::sqrt(2.0));
        double N_d2 = 0.5 * std::erfc(-d2 / std::sqrt(2.0));
        double n_d1 = std::exp(-0.5 * d1 * d1) / std::sqrt(2.0 * M_PI);
        
        double preco = ehCall ? precoAtivo * N_d1 - strike * std::exp(-r * T) * N_d2
                              : strike * std::exp(-r * T) * (1 - N_d2) - precoAtivo * (1 - N_d1);
        double vega = precoAtivo * n_d1 * std::sqrt(T);
        
        if (std::abs(vega) < 1e-10) break;
        double diff = preco - precoOpcao;
        if (std::abs(diff) < tol) return sigma;
        sigma -= diff / vega;
        sigma = std::clamp(sigma, 0.01, 5.0);
    }
    return sigma;
}

void AnaliseVolatilidade::verificarDados() const {
    if (!serieCarregada) throw std::runtime_error("Dados não carregados");
}

[[nodiscard]] auto AnaliseVolatilidade::obterRetornos() const -> const std::vector<double>& { return retornosLog; }

} // namespace analise_ouro
