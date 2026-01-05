/*
 * AnaliseRisco.cpp - Análise de Risco para Investimento em Ouro - C++23
 * Autor: Luiz Tiago Wilcke
 * Projeto: Análise do Preço do Ouro
 */

#include "AnaliseRisco.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace analise_ouro {

AnaliseRisco::AnaliseRisco() : dadosCarregados(false), tamanhoSerie(0) {}
AnaliseRisco::AnaliseRisco(const std::vector<double>& precos) { carregarPrecos(precos); }
AnaliseRisco::~AnaliseRisco() = default;

void AnaliseRisco::carregarPrecos(const std::vector<double>& precos) {
    if (precos.size() < 2) throw std::invalid_argument("Série muito curta");
    precosOriginais = precos;
    tamanhoSerie = precos.size();
    calcularRetornos();
    dadosCarregados = true;
}

void AnaliseRisco::calcularRetornos() {
    retornos.resize(tamanhoSerie - 1);
    for (size_t i = 1; i < tamanhoSerie; ++i) {
        retornos[i - 1] = (precosOriginais[i] - precosOriginais[i - 1]) / precosOriginais[i - 1];
    }
}

[[nodiscard]] auto AnaliseRisco::valorEmRisco(double nivelConfianca, int janela) const -> double {
    verificarDados();
    
    std::vector<double> retornosJanela;
    int inicio = std::max(0, static_cast<int>(retornos.size()) - janela);
    for (size_t i = inicio; i < retornos.size(); ++i) {
        retornosJanela.push_back(retornos[i]);
    }
    
    std::sort(retornosJanela.begin(), retornosJanela.end());
    int indice = static_cast<int>((1.0 - nivelConfianca) * retornosJanela.size());
    return -retornosJanela[indice];
}

[[nodiscard]] auto AnaliseRisco::valorEmRiscoCondicional(double nivelConfianca, int janela) const -> double {
    verificarDados();
    
    std::vector<double> retornosJanela;
    int inicio = std::max(0, static_cast<int>(retornos.size()) - janela);
    for (size_t i = inicio; i < retornos.size(); ++i) {
        retornosJanela.push_back(retornos[i]);
    }
    
    std::sort(retornosJanela.begin(), retornosJanela.end());
    int numCauda = static_cast<int>((1.0 - nivelConfianca) * retornosJanela.size());
    
    double soma = 0.0;
    for (int i = 0; i < numCauda; ++i) {
        soma += retornosJanela[i];
    }
    return -soma / numCauda;
}

[[nodiscard]] auto AnaliseRisco::valorEmRiscoParametrico(double nivelConfianca) const -> double {
    verificarDados();
    
    double media = std::accumulate(retornos.begin(), retornos.end(), 0.0) / retornos.size();
    double var = 0.0;
    for (const auto& r : retornos) {
        var += (r - media) * (r - media);
    }
    var /= (retornos.size() - 1);
    double dp = std::sqrt(var);
    
    // Quantil da normal padrão
    double z = 0.0;
    if (nivelConfianca >= 0.99) z = 2.326;
    else if (nivelConfianca >= 0.95) z = 1.645;
    else if (nivelConfianca >= 0.90) z = 1.282;
    else z = 1.0;
    
    return -(media - z * dp);
}

[[nodiscard]] auto AnaliseRisco::indiceSharte(double taxaLivreRisco) const -> double {
    verificarDados();
    
    double retornoMedio = std::accumulate(retornos.begin(), retornos.end(), 0.0) / retornos.size();
    double var = 0.0;
    for (const auto& r : retornos) {
        var += (r - retornoMedio) * (r - retornoMedio);
    }
    double dp = std::sqrt(var / (retornos.size() - 1));
    
    double excessoRetorno = retornoMedio * 252 - taxaLivreRisco;
    double volAnualizada = dp * std::sqrt(252.0);
    
    return (std::abs(volAnualizada) > 1e-10) ? excessoRetorno / volAnualizada : 0.0;
}

[[nodiscard]] auto AnaliseRisco::indiceSortino(double taxaLivreRisco) const -> double {
    verificarDados();
    
    double retornoMedio = std::accumulate(retornos.begin(), retornos.end(), 0.0) / retornos.size();
    double taxaDiariaLivreRisco = taxaLivreRisco / 252.0;
    
    double somaQuadNeg = 0.0;
    int contNeg = 0;
    for (const auto& r : retornos) {
        if (r < taxaDiariaLivreRisco) {
            double diff = r - taxaDiariaLivreRisco;
            somaQuadNeg += diff * diff;
            contNeg++;
        }
    }
    
    double downside = (contNeg > 0) ? std::sqrt(somaQuadNeg / contNeg) : 1e-10;
    double excessoRetorno = retornoMedio * 252 - taxaLivreRisco;
    double downsideAnualizado = downside * std::sqrt(252.0);
    
    return (std::abs(downsideAnualizado) > 1e-10) ? excessoRetorno / downsideAnualizado : 0.0;
}

[[nodiscard]] auto AnaliseRisco::maxDrawdown() const -> double {
    verificarDados();
    
    double pico = precosOriginais[0];
    double maxDD = 0.0;
    
    for (const auto& preco : precosOriginais) {
        if (preco > pico) pico = preco;
        double dd = (pico - preco) / pico;
        if (dd > maxDD) maxDD = dd;
    }
    return maxDD;
}

[[nodiscard]] auto AnaliseRisco::calcularDrawdowns() const -> std::vector<double> {
    verificarDados();
    
    std::vector<double> drawdowns(tamanhoSerie);
    double pico = precosOriginais[0];
    
    for (size_t i = 0; i < tamanhoSerie; ++i) {
        if (precosOriginais[i] > pico) pico = precosOriginais[i];
        drawdowns[i] = (pico - precosOriginais[i]) / pico;
    }
    return drawdowns;
}

[[nodiscard]] auto AnaliseRisco::indiceCalmar(double taxaLivreRisco) const -> double {
    verificarDados();
    
    double retornoAnual = std::accumulate(retornos.begin(), retornos.end(), 0.0) / retornos.size() * 252;
    double mdd = maxDrawdown();
    
    return (std::abs(mdd) > 1e-10) ? (retornoAnual - taxaLivreRisco) / mdd : 0.0;
}

[[nodiscard]] auto AnaliseRisco::betaModificado(const std::vector<double>& retornosMercado) const -> double {
    verificarDados();
    
    size_t n = std::min(retornos.size(), retornosMercado.size());
    if (n < 2) return 0.0;
    
    double mediaAtivo = std::accumulate(retornos.begin(), retornos.begin() + n, 0.0) / n;
    double mediaMercado = std::accumulate(retornosMercado.begin(), retornosMercado.begin() + n, 0.0) / n;
    
    double cov = 0.0, varMercado = 0.0;
    for (size_t i = 0; i < n; ++i) {
        cov += (retornos[i] - mediaAtivo) * (retornosMercado[i] - mediaMercado);
        varMercado += (retornosMercado[i] - mediaMercado) * (retornosMercado[i] - mediaMercado);
    }
    
    return (std::abs(varMercado) > 1e-10) ? cov / varMercado : 0.0;
}

[[nodiscard]] auto AnaliseRisco::informationRatio(const std::vector<double>& retornosBenchmark) const -> double {
    verificarDados();
    
    size_t n = std::min(retornos.size(), retornosBenchmark.size());
    if (n < 2) return 0.0;
    
    std::vector<double> excessos(n);
    for (size_t i = 0; i < n; ++i) {
        excessos[i] = retornos[i] - retornosBenchmark[i];
    }
    
    double mediaExcesso = std::accumulate(excessos.begin(), excessos.end(), 0.0) / n;
    double var = 0.0;
    for (const auto& e : excessos) {
        var += (e - mediaExcesso) * (e - mediaExcesso);
    }
    double trackingError = std::sqrt(var / (n - 1)) * std::sqrt(252.0);
    
    return (std::abs(trackingError) > 1e-10) ? mediaExcesso * 252 / trackingError : 0.0;
}

void AnaliseRisco::verificarDados() const {
    if (!dadosCarregados) throw std::runtime_error("Dados não carregados");
}

[[nodiscard]] auto AnaliseRisco::obterRetornos() const -> const std::vector<double>& { return retornos; }

} // namespace analise_ouro
