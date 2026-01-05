/*
 * PrevisaoPrecos.cpp - Previsão de Preços do Ouro - C++23
 * Autor: Luiz Tiago Wilcke
 * Projeto: Análise do Preço do Ouro
 */

#include "PrevisaoPrecos.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <stdexcept>

namespace analise_ouro {

PrevisaoPrecos::PrevisaoPrecos() : dadosCarregados(false), tamanhoSerie(0) {
    auto semente = std::random_device{}();
    gerador.seed(semente);
}

PrevisaoPrecos::PrevisaoPrecos(const std::vector<double>& precos) {
    auto semente = std::random_device{}();
    gerador.seed(semente);
    carregarPrecos(precos);
}

PrevisaoPrecos::~PrevisaoPrecos() = default;

void PrevisaoPrecos::carregarPrecos(const std::vector<double>& precos) {
    if (precos.size() < 10) throw std::invalid_argument("Série muito curta");
    precosOriginais = precos;
    tamanhoSerie = precos.size();
    dadosCarregados = true;
}

[[nodiscard]] auto PrevisaoPrecos::preverMediaMovel(int horizonte, int janela) const -> std::vector<double> {
    verificarDados();
    
    std::vector<double> previsoes(horizonte);
    std::vector<double> historia(precosOriginais);
    
    for (int h = 0; h < horizonte; ++h) {
        double soma = 0.0;
        int inicio = std::max(0, static_cast<int>(historia.size()) - janela);
        for (size_t i = inicio; i < historia.size(); ++i) {
            soma += historia[i];
        }
        double previsao = soma / std::min(janela, static_cast<int>(historia.size()));
        previsoes[h] = previsao;
        historia.push_back(previsao);
    }
    return previsoes;
}

[[nodiscard]] auto PrevisaoPrecos::preverSuavizacaoExponencial(int horizonte, double alfa) const -> std::vector<double> {
    verificarDados();
    
    double nivel = precosOriginais[0];
    for (size_t i = 1; i < tamanhoSerie; ++i) {
        nivel = alfa * precosOriginais[i] + (1.0 - alfa) * nivel;
    }
    
    std::vector<double> previsoes(horizonte, nivel);
    return previsoes;
}

[[nodiscard]] auto PrevisaoPrecos::preverHoltLinear(int horizonte, double alfa, double beta) const -> std::vector<double> {
    verificarDados();
    
    double nivel = precosOriginais[0];
    double tendencia = precosOriginais[1] - precosOriginais[0];
    
    for (size_t t = 1; t < tamanhoSerie; ++t) {
        double nivelAnterior = nivel;
        nivel = alfa * precosOriginais[t] + (1.0 - alfa) * (nivel + tendencia);
        tendencia = beta * (nivel - nivelAnterior) + (1.0 - beta) * tendencia;
    }
    
    std::vector<double> previsoes(horizonte);
    for (int h = 0; h < horizonte; ++h) {
        previsoes[h] = nivel + (h + 1) * tendencia;
    }
    return previsoes;
}

[[nodiscard]] auto PrevisaoPrecos::preverMonteCarlo(int horizonte, int numSimulacoes) const -> ResultadoPrevisao {
    verificarDados();
    
    // Calcula retornos e parâmetros
    std::vector<double> retornos(tamanhoSerie - 1);
    for (size_t i = 1; i < tamanhoSerie; ++i) {
        retornos[i - 1] = std::log(precosOriginais[i] / precosOriginais[i - 1]);
    }
    
    double mediaRetorno = std::accumulate(retornos.begin(), retornos.end(), 0.0) / retornos.size();
    double var = 0.0;
    for (const auto& r : retornos) {
        var += (r - mediaRetorno) * (r - mediaRetorno);
    }
    double sigma = std::sqrt(var / (retornos.size() - 1));
    
    std::vector<std::vector<double>> simulacoes(numSimulacoes);
    std::mt19937_64 gen = gerador;
    std::normal_distribution<double> Z(0.0, 1.0);
    
    for (int sim = 0; sim < numSimulacoes; ++sim) {
        simulacoes[sim].resize(horizonte);
        double preco = precosOriginais.back();
        
        for (int h = 0; h < horizonte; ++h) {
            double retorno = (mediaRetorno - 0.5 * sigma * sigma) + sigma * Z(gen);
            preco *= std::exp(retorno);
            simulacoes[sim][h] = preco;
        }
    }
    
    ResultadoPrevisao resultado;
    resultado.previsaoMedia.resize(horizonte);
    resultado.intervaloInferior.resize(horizonte);
    resultado.intervaloSuperior.resize(horizonte);
    
    for (int h = 0; h < horizonte; ++h) {
        std::vector<double> valoresH(numSimulacoes);
        for (int sim = 0; sim < numSimulacoes; ++sim) {
            valoresH[sim] = simulacoes[sim][h];
        }
        std::sort(valoresH.begin(), valoresH.end());
        
        resultado.previsaoMedia[h] = std::accumulate(valoresH.begin(), valoresH.end(), 0.0) / numSimulacoes;
        resultado.intervaloInferior[h] = valoresH[static_cast<int>(0.025 * numSimulacoes)];
        resultado.intervaloSuperior[h] = valoresH[static_cast<int>(0.975 * numSimulacoes)];
    }
    
    return resultado;
}

[[nodiscard]] auto PrevisaoPrecos::preverArima(int horizonte, int p, int d, int q) const -> std::vector<double> {
    verificarDados();
    
    // Diferenciação
    std::vector<double> serie = precosOriginais;
    for (int i = 0; i < d; ++i) {
        std::vector<double> diff(serie.size() - 1);
        for (size_t j = 1; j < serie.size(); ++j) {
            diff[j - 1] = serie[j] - serie[j - 1];
        }
        serie = diff;
    }
    
    double media = std::accumulate(serie.begin(), serie.end(), 0.0) / serie.size();
    
    // Coeficientes AR (simplificado)
    std::vector<double> phi(p, 0.1);
    
    // Previsão
    std::vector<double> previsoesDiff(horizonte);
    for (int h = 0; h < horizonte; ++h) {
        double pred = media;
        for (int i = 0; i < p && static_cast<int>(serie.size()) - 1 - i >= 0; ++i) {
            pred += phi[i] * (serie[serie.size() - 1 - i] - media);
        }
        previsoesDiff[h] = pred;
        serie.push_back(pred);
    }
    
    // Integração
    std::vector<double> previsoes(horizonte);
    double ultimoValor = precosOriginais.back();
    for (int h = 0; h < horizonte; ++h) {
        ultimoValor += previsoesDiff[h];
        previsoes[h] = ultimoValor;
    }
    
    return previsoes;
}

[[nodiscard]] auto PrevisaoPrecos::calcularErroPrevisao(const std::vector<double>& reais, const std::vector<double>& previstos) const -> MetricasErro {
    if (reais.size() != previstos.size() || reais.empty()) {
        throw std::invalid_argument("Tamanhos incompatíveis");
    }
    
    MetricasErro metricas;
    double somaE = 0.0, somaE2 = 0.0, somaAE = 0.0, somaAPE = 0.0;
    size_t n = reais.size();
    
    for (size_t i = 0; i < n; ++i) {
        double erro = reais[i] - previstos[i];
        somaE += erro;
        somaE2 += erro * erro;
        somaAE += std::abs(erro);
        if (std::abs(reais[i]) > 1e-10) {
            somaAPE += std::abs(erro / reais[i]);
        }
    }
    
    metricas.me = somaE / n;
    metricas.mse = somaE2 / n;
    metricas.rmse = std::sqrt(metricas.mse);
    metricas.mae = somaAE / n;
    metricas.mape = somaAPE / n * 100.0;
    
    return metricas;
}

void PrevisaoPrecos::verificarDados() const {
    if (!dadosCarregados) throw std::runtime_error("Dados não carregados");
}

} // namespace analise_ouro
