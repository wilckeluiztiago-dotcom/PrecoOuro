/*
 * IndicadoresTecnicos.cpp - Indicadores Técnicos para Ouro - C++23
 * Autor: Luiz Tiago Wilcke
 * Projeto: Análise do Preço do Ouro
 */

#include "IndicadoresTecnicos.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace analise_ouro {

IndicadoresTecnicos::IndicadoresTecnicos() : dadosCarregados(false), tamanhoSerie(0) {}
IndicadoresTecnicos::IndicadoresTecnicos(const std::vector<double>& precos) { carregarPrecos(precos); }
IndicadoresTecnicos::~IndicadoresTecnicos() = default;

void IndicadoresTecnicos::carregarPrecos(const std::vector<double>& precos) {
    if (precos.size() < 2) throw std::invalid_argument("Precisa de pelo menos 2 preços");
    precosOriginais = precos;
    tamanhoSerie = precos.size();
    dadosCarregados = true;
}

[[nodiscard]] auto IndicadoresTecnicos::rsi(int periodo) const -> std::vector<double> {
    verificarDados();
    if (periodo < 1) periodo = 14;
    
    std::vector<double> resultado(tamanhoSerie, 0.0);
    std::vector<double> ganhos(tamanhoSerie - 1, 0.0);
    std::vector<double> perdas(tamanhoSerie - 1, 0.0);
    
    for (size_t i = 1; i < tamanhoSerie; ++i) {
        double delta = precosOriginais[i] - precosOriginais[i - 1];
        if (delta > 0) ganhos[i - 1] = delta;
        else perdas[i - 1] = -delta;
    }
    
    double ganhoMedio = 0.0, perdaMedia = 0.0;
    for (int i = 0; i < periodo && i < static_cast<int>(ganhos.size()); ++i) {
        ganhoMedio += ganhos[i];
        perdaMedia += perdas[i];
    }
    ganhoMedio /= periodo;
    perdaMedia /= periodo;
    
    for (size_t i = periodo; i < tamanhoSerie; ++i) {
        ganhoMedio = (ganhoMedio * (periodo - 1) + ganhos[i - 1]) / periodo;
        perdaMedia = (perdaMedia * (periodo - 1) + perdas[i - 1]) / periodo;
        
        double rs = (perdaMedia > 1e-10) ? ganhoMedio / perdaMedia : 100.0;
        resultado[i] = 100.0 - 100.0 / (1.0 + rs);
    }
    return resultado;
}

[[nodiscard]] auto IndicadoresTecnicos::macd(int periodoRapido, int periodoLento, int periodoSinal) const -> std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> {
    verificarDados();
    
    auto emaRapida = ema(periodoRapido);
    auto emaLenta = ema(periodoLento);
    
    std::vector<double> linhaMACD(tamanhoSerie);
    for (size_t i = 0; i < tamanhoSerie; ++i) {
        linhaMACD[i] = emaRapida[i] - emaLenta[i];
    }
    
    // EMA do MACD
    std::vector<double> sinal(tamanhoSerie);
    double alfa = 2.0 / (periodoSinal + 1.0);
    sinal[0] = linhaMACD[0];
    for (size_t i = 1; i < tamanhoSerie; ++i) {
        sinal[i] = alfa * linhaMACD[i] + (1.0 - alfa) * sinal[i - 1];
    }
    
    std::vector<double> histograma(tamanhoSerie);
    for (size_t i = 0; i < tamanhoSerie; ++i) {
        histograma[i] = linhaMACD[i] - sinal[i];
    }
    
    return {linhaMACD, sinal, histograma};
}

[[nodiscard]] auto IndicadoresTecnicos::bollingerBands(int periodo, double numDesvios) const -> std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> {
    verificarDados();
    
    std::vector<double> bandaSuperior(tamanhoSerie);
    std::vector<double> mediaMovel(tamanhoSerie);
    std::vector<double> bandaInferior(tamanhoSerie);
    
    for (size_t i = 0; i < tamanhoSerie; ++i) {
        int inicio = std::max(0, static_cast<int>(i) - periodo + 1);
        int n = static_cast<int>(i) - inicio + 1;
        
        double soma = 0.0;
        for (int j = inicio; j <= static_cast<int>(i); ++j) {
            soma += precosOriginais[j];
        }
        double media = soma / n;
        
        double somaQuad = 0.0;
        for (int j = inicio; j <= static_cast<int>(i); ++j) {
            double diff = precosOriginais[j] - media;
            somaQuad += diff * diff;
        }
        double dp = std::sqrt(somaQuad / n);
        
        mediaMovel[i] = media;
        bandaSuperior[i] = media + numDesvios * dp;
        bandaInferior[i] = media - numDesvios * dp;
    }
    
    return {bandaSuperior, mediaMovel, bandaInferior};
}

[[nodiscard]] auto IndicadoresTecnicos::stochastic(int periodoPrincipal, int periodoSuavizacao) const -> std::pair<std::vector<double>, std::vector<double>> {
    verificarDados();
    
    std::vector<double> percentK(tamanhoSerie, 0.0);
    std::vector<double> percentD(tamanhoSerie, 0.0);
    
    for (size_t i = periodoPrincipal - 1; i < tamanhoSerie; ++i) {
        double minimo = precosOriginais[i];
        double maximo = precosOriginais[i];
        
        for (int j = 0; j < periodoPrincipal; ++j) {
            size_t idx = i - j;
            minimo = std::min(minimo, precosOriginais[idx]);
            maximo = std::max(maximo, precosOriginais[idx]);
        }
        
        if (maximo - minimo > 1e-10) {
            percentK[i] = 100.0 * (precosOriginais[i] - minimo) / (maximo - minimo);
        }
    }
    
    // %D é a média móvel de %K
    for (size_t i = periodoSuavizacao - 1; i < tamanhoSerie; ++i) {
        double soma = 0.0;
        for (int j = 0; j < periodoSuavizacao; ++j) {
            soma += percentK[i - j];
        }
        percentD[i] = soma / periodoSuavizacao;
    }
    
    return {percentK, percentD};
}

[[nodiscard]] auto IndicadoresTecnicos::atr(int periodo) const -> std::vector<double> {
    verificarDados();
    
    std::vector<double> tr(tamanhoSerie - 1);
    for (size_t i = 1; i < tamanhoSerie; ++i) {
        double highLow = std::abs(precosOriginais[i] - precosOriginais[i - 1]);
        tr[i - 1] = highLow;
    }
    
    std::vector<double> resultado(tamanhoSerie, 0.0);
    double soma = 0.0;
    for (int i = 0; i < periodo && i < static_cast<int>(tr.size()); ++i) {
        soma += tr[i];
    }
    resultado[periodo] = soma / periodo;
    
    for (size_t i = periodo + 1; i < tamanhoSerie; ++i) {
        resultado[i] = (resultado[i - 1] * (periodo - 1) + tr[i - 1]) / periodo;
    }
    return resultado;
}

[[nodiscard]] auto IndicadoresTecnicos::ema(int periodo) const -> std::vector<double> {
    verificarDados();
    
    std::vector<double> resultado(tamanhoSerie);
    double alfa = 2.0 / (periodo + 1.0);
    resultado[0] = precosOriginais[0];
    
    for (size_t i = 1; i < tamanhoSerie; ++i) {
        resultado[i] = alfa * precosOriginais[i] + (1.0 - alfa) * resultado[i - 1];
    }
    return resultado;
}

[[nodiscard]] auto IndicadoresTecnicos::adx(int periodo) const -> std::vector<double> {
    verificarDados();
    
    std::vector<double> resultado(tamanhoSerie, 0.0);
    // Simplificado - usa variações de preço
    for (size_t i = periodo; i < tamanhoSerie; ++i) {
        double soma = 0.0;
        for (int j = 0; j < periodo; ++j) {
            soma += std::abs(precosOriginais[i - j] - precosOriginais[i - j - 1]);
        }
        resultado[i] = soma / periodo;
    }
    return resultado;
}

void IndicadoresTecnicos::verificarDados() const {
    if (!dadosCarregados) throw std::runtime_error("Dados não carregados");
}

} // namespace analise_ouro
