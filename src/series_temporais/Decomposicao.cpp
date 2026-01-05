/*
 * Decomposicao.cpp
 * Módulo de decomposição de séries temporais
 * Usa recursos avançados de C++23
 * 
 * Autor: Luiz Tiago Wilcke
 * Projeto: Análise do Preço do Ouro
 */

#include "Decomposicao.hpp"
#include "MediaMovel.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <ranges>

namespace analise_ouro {

// Construtor padrão
Decomposicao::Decomposicao() 
    : serieCarregada(false), tamanhoSerie(0), periodoSazonal(12) {}

// Construtor com dados e período
Decomposicao::Decomposicao(const std::vector<double>& dados, size_t periodo) 
    : periodoSazonal(periodo) {
    carregarDados(dados);
}

// Destrutor
Decomposicao::~Decomposicao() = default;

// Carrega dados
void Decomposicao::carregarDados(const std::vector<double>& novosDados) {
    if (novosDados.empty()) {
        throw std::invalid_argument("Dados não podem estar vazios");
    }
    
    serieOriginal = novosDados;
    tamanhoSerie = novosDados.size();
    serieCarregada = true;
}

// Define período sazonal
void Decomposicao::definirPeriodo(size_t novoPeriodo) {
    if (novoPeriodo == 0) {
        throw std::invalid_argument("Período deve ser maior que zero");
    }
    periodoSazonal = novoPeriodo;
}

// Decomposição aditiva: Y = T + S + R
[[nodiscard]] auto Decomposicao::decomposicaoAditiva() const -> ComponentesDecomposicao {
    verificarDados();
    
    ComponentesDecomposicao resultado;
    resultado.original = serieOriginal;
    
    // 1. Extrai tendência usando média móvel centrada
    resultado.tendencia = extrairTendencia();
    
    // 2. Remove tendência para obter série dessazonalizada
    std::vector<double> semTendencia(tamanhoSerie);
    for (size_t i = 0; i < tamanhoSerie; ++i) {
        if (std::isnan(resultado.tendencia[i])) {
            semTendencia[i] = std::nan("");
        } else {
            semTendencia[i] = serieOriginal[i] - resultado.tendencia[i];
        }
    }
    
    // 3. Extrai componente sazonal
    resultado.sazonal = extrairSazonalidade(semTendencia);
    
    // 4. Resíduo = Original - Tendência - Sazonal
    resultado.residuo.resize(tamanhoSerie);
    for (size_t i = 0; i < tamanhoSerie; ++i) {
        if (std::isnan(resultado.tendencia[i])) {
            resultado.residuo[i] = std::nan("");
        } else {
            resultado.residuo[i] = serieOriginal[i] - resultado.tendencia[i] - resultado.sazonal[i];
        }
    }
    
    return resultado;
}

// Decomposição multiplicativa: Y = T * S * R
[[nodiscard]] auto Decomposicao::decomposicaoMultiplicativa() const -> ComponentesDecomposicao {
    verificarDados();
    
    // Verifica se há valores não-positivos (não pode usar log)
    for (const auto& val : serieOriginal) {
        if (val <= 0) {
            throw std::runtime_error("Decomposição multiplicativa requer valores positivos");
        }
    }
    
    ComponentesDecomposicao resultado;
    resultado.original = serieOriginal;
    
    // 1. Extrai tendência
    resultado.tendencia = extrairTendencia();
    
    // 2. Série sem tendência (divisão)
    std::vector<double> semTendencia(tamanhoSerie);
    for (size_t i = 0; i < tamanhoSerie; ++i) {
        if (std::isnan(resultado.tendencia[i]) || resultado.tendencia[i] <= 0) {
            semTendencia[i] = std::nan("");
        } else {
            semTendencia[i] = serieOriginal[i] / resultado.tendencia[i];
        }
    }
    
    // 3. Extrai sazonalidade
    resultado.sazonal = extrairSazonalidadeMultiplicativa(semTendencia);
    
    // 4. Resíduo = Original / (Tendência * Sazonal)
    resultado.residuo.resize(tamanhoSerie);
    for (size_t i = 0; i < tamanhoSerie; ++i) {
        double ts = resultado.tendencia[i] * resultado.sazonal[i];
        if (std::isnan(ts) || std::abs(ts) < 1e-15) {
            resultado.residuo[i] = std::nan("");
        } else {
            resultado.residuo[i] = serieOriginal[i] / ts;
        }
    }
    
    return resultado;
}

// Extrai tendência usando média móvel centrada
[[nodiscard]] auto Decomposicao::extrairTendencia() const -> std::vector<double> {
    verificarDados();
    
    std::vector<double> tendencia(tamanhoSerie, std::nan(""));
    
    if (tamanhoSerie < periodoSazonal) {
        return tendencia;
    }
    
    MediaMovel mm(serieOriginal);
    
    // Média móvel centrada
    if (periodoSazonal % 2 == 0) {
        // Para período par: usa média móvel 2x
        auto sma1 = mm.sma(periodoSazonal);
        if (sma1.size() < 2) return tendencia;
        
        MediaMovel mm2(sma1);
        auto sma2 = mm2.sma(2);
        
        // Posiciona corretamente
        size_t offset = periodoSazonal / 2;
        for (size_t i = 0; i < sma2.size(); ++i) {
            if (offset + i < tamanhoSerie) {
                tendencia[offset + i] = sma2[i];
            }
        }
    } else {
        // Para período ímpar: média móvel simples centrada
        auto sma = mm.sma(periodoSazonal);
        size_t offset = periodoSazonal / 2;
        for (size_t i = 0; i < sma.size(); ++i) {
            if (offset + i < tamanhoSerie) {
                tendencia[offset + i] = sma[i];
            }
        }
    }
    
    return tendencia;
}

// Extrai sazonalidade (aditiva)
[[nodiscard]] auto Decomposicao::extrairSazonalidade(const std::vector<double>& semTendencia) const -> std::vector<double> {
    std::vector<double> sazonal(tamanhoSerie);
    
    // Calcula índices sazonais médios
    std::vector<double> indicesSazonais(periodoSazonal, 0.0);
    std::vector<int> contagem(periodoSazonal, 0);
    
    for (size_t i = 0; i < tamanhoSerie; ++i) {
        if (!std::isnan(semTendencia[i])) {
            size_t posicaoSazonal = i % periodoSazonal;
            indicesSazonais[posicaoSazonal] += semTendencia[i];
            contagem[posicaoSazonal]++;
        }
    }
    
    // Calcula médias
    for (size_t i = 0; i < periodoSazonal; ++i) {
        if (contagem[i] > 0) {
            indicesSazonais[i] /= contagem[i];
        }
    }
    
    // Normaliza para soma zero
    double somaIndices = std::accumulate(indicesSazonais.begin(), indicesSazonais.end(), 0.0);
    double ajuste = somaIndices / static_cast<double>(periodoSazonal);
    for (auto& idx : indicesSazonais) {
        idx -= ajuste;
    }
    
    // Aplica índices sazonais
    for (size_t i = 0; i < tamanhoSerie; ++i) {
        sazonal[i] = indicesSazonais[i % periodoSazonal];
    }
    
    return sazonal;
}

// Extrai sazonalidade (multiplicativa)
[[nodiscard]] auto Decomposicao::extrairSazonalidadeMultiplicativa(const std::vector<double>& semTendencia) const -> std::vector<double> {
    std::vector<double> sazonal(tamanhoSerie);
    
    // Calcula índices sazonais médios
    std::vector<double> indicesSazonais(periodoSazonal, 0.0);
    std::vector<int> contagem(periodoSazonal, 0);
    
    for (size_t i = 0; i < tamanhoSerie; ++i) {
        if (!std::isnan(semTendencia[i])) {
            size_t posicaoSazonal = i % periodoSazonal;
            indicesSazonais[posicaoSazonal] += semTendencia[i];
            contagem[posicaoSazonal]++;
        }
    }
    
    // Calcula médias
    for (size_t i = 0; i < periodoSazonal; ++i) {
        if (contagem[i] > 0) {
            indicesSazonais[i] /= contagem[i];
        } else {
            indicesSazonais[i] = 1.0;
        }
    }
    
    // Normaliza para produto = 1 (média geométrica = 1)
    double somaIndices = std::accumulate(indicesSazonais.begin(), indicesSazonais.end(), 0.0);
    double fator = static_cast<double>(periodoSazonal) / somaIndices;
    for (auto& idx : indicesSazonais) {
        idx *= fator;
    }
    
    // Aplica índices
    for (size_t i = 0; i < tamanhoSerie; ++i) {
        sazonal[i] = indicesSazonais[i % periodoSazonal];
    }
    
    return sazonal;
}

// Decomposição STL (Seasonal-Trend decomposition using LOESS) - simplificada
[[nodiscard]] auto Decomposicao::decomposicaoSTL(int numIteracoes) const -> ComponentesDecomposicao {
    verificarDados();
    
    ComponentesDecomposicao resultado;
    resultado.original = serieOriginal;
    resultado.tendencia.resize(tamanhoSerie, 0.0);
    resultado.sazonal.resize(tamanhoSerie, 0.0);
    resultado.residuo.resize(tamanhoSerie, 0.0);
    
    // Inicializa tendência
    MediaMovel mm(serieOriginal);
    resultado.tendencia = mm.smaComPadding(periodoSazonal);
    
    // Iterações STL
    for (int iter = 0; iter < numIteracoes; ++iter) {
        // 1. Remove tendência
        std::vector<double> detrendado(tamanhoSerie);
        for (size_t i = 0; i < tamanhoSerie; ++i) {
            detrendado[i] = serieOriginal[i] - resultado.tendencia[i];
        }
        
        // 2. Extrai sazonalidade
        resultado.sazonal = extrairSazonalidade(detrendado);
        
        // 3. Remove sazonalidade
        std::vector<double> dessazonalizado(tamanhoSerie);
        for (size_t i = 0; i < tamanhoSerie; ++i) {
            dessazonalizado[i] = serieOriginal[i] - resultado.sazonal[i];
        }
        
        // 4. Atualiza tendência com LOESS simplificado (usa média móvel ponderada)
        MediaMovel mmDess(dessazonalizado);
        resultado.tendencia = mmDess.smaComPadding(periodoSazonal);
    }
    
    // Resíduo final
    for (size_t i = 0; i < tamanhoSerie; ++i) {
        resultado.residuo[i] = serieOriginal[i] - resultado.tendencia[i] - resultado.sazonal[i];
    }
    
    return resultado;
}

// Suavização exponencial de Holt-Winters
[[nodiscard]] auto Decomposicao::holtWinters(double alfa, double beta, double gamma, bool multiplicativo) const -> std::vector<double> {
    verificarDados();
    
    if (tamanhoSerie < 2 * periodoSazonal) {
        throw std::runtime_error("Série muito curta para Holt-Winters");
    }
    
    std::vector<double> nivel(tamanhoSerie);
    std::vector<double> tendencia(tamanhoSerie);
    std::vector<double> sazonal(tamanhoSerie + periodoSazonal);
    std::vector<double> previsao(tamanhoSerie);
    
    // Inicialização
    double somaInicial = 0.0;
    for (size_t i = 0; i < periodoSazonal; ++i) {
        somaInicial += serieOriginal[i];
    }
    nivel[0] = somaInicial / static_cast<double>(periodoSazonal);
    
    // Tendência inicial
    double somaTend = 0.0;
    for (size_t i = 0; i < periodoSazonal; ++i) {
        somaTend += (serieOriginal[periodoSazonal + i] - serieOriginal[i]);
    }
    tendencia[0] = somaTend / (static_cast<double>(periodoSazonal) * static_cast<double>(periodoSazonal));
    
    // Sazonalidade inicial
    for (size_t i = 0; i < periodoSazonal; ++i) {
        if (multiplicativo) {
            sazonal[i] = serieOriginal[i] / nivel[0];
        } else {
            sazonal[i] = serieOriginal[i] - nivel[0];
        }
    }
    
    // Aplicação recursiva
    for (size_t t = 1; t < tamanhoSerie; ++t) {
        size_t sIdx = t % periodoSazonal;
        
        if (multiplicativo) {
            nivel[t] = alfa * (serieOriginal[t] / sazonal[sIdx]) + 
                       (1.0 - alfa) * (nivel[t-1] + tendencia[t-1]);
            tendencia[t] = beta * (nivel[t] - nivel[t-1]) + 
                          (1.0 - beta) * tendencia[t-1];
            sazonal[t + periodoSazonal] = gamma * (serieOriginal[t] / nivel[t]) + 
                                          (1.0 - gamma) * sazonal[sIdx];
            previsao[t] = (nivel[t] + tendencia[t]) * sazonal[sIdx];
        } else {
            nivel[t] = alfa * (serieOriginal[t] - sazonal[sIdx]) + 
                       (1.0 - alfa) * (nivel[t-1] + tendencia[t-1]);
            tendencia[t] = beta * (nivel[t] - nivel[t-1]) + 
                          (1.0 - beta) * tendencia[t-1];
            sazonal[t + periodoSazonal] = gamma * (serieOriginal[t] - nivel[t]) + 
                                          (1.0 - gamma) * sazonal[sIdx];
            previsao[t] = nivel[t] + tendencia[t] + sazonal[sIdx];
        }
    }
    
    return previsao;
}

// Verifica dados
void Decomposicao::verificarDados() const {
    if (!serieCarregada || tamanhoSerie == 0) {
        throw std::runtime_error("Dados não carregados");
    }
}

// Getters
[[nodiscard]] auto Decomposicao::obterTamanho() const noexcept -> size_t {
    return tamanhoSerie;
}

[[nodiscard]] auto Decomposicao::obterPeriodo() const noexcept -> size_t {
    return periodoSazonal;
}

} // namespace analise_ouro
