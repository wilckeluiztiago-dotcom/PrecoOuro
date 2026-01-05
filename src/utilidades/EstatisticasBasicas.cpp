/*
 * EstatisticasBasicas.cpp
 * Módulo de estatísticas descritivas para análise financeira
 * 
 * Autor: Luiz Tiago Wilcke
 * Projeto: Análise do Preço do Ouro
 */

#include "EstatisticasBasicas.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <limits>

namespace analise_ouro {

// Construtor padrão
EstatisticasBasicas::EstatisticasBasicas() : dadosCarregados(false), tamanhoDados(0) {}

// Construtor com dados
EstatisticasBasicas::EstatisticasBasicas(const std::vector<double>& dados) {
    carregarDados(dados);
}

// Destrutor
EstatisticasBasicas::~EstatisticasBasicas() {}

// Carrega dados para análise
void EstatisticasBasicas::carregarDados(const std::vector<double>& novosDados) {
    if (novosDados.empty()) {
        throw std::invalid_argument("Dados não podem estar vazios");
    }
    
    serieOriginal = novosDados;
    dadosOrdenados = novosDados;
    std::sort(dadosOrdenados.begin(), dadosOrdenados.end());
    
    tamanhoDados = novosDados.size();
    dadosCarregados = true;
    
    // Pré-calcula estatísticas básicas
    calcularEstatisticasBasicas();
}

// Calcula estatísticas básicas
void EstatisticasBasicas::calcularEstatisticasBasicas() {
    mediaCache = 0.0;
    for (const auto& valor : serieOriginal) {
        mediaCache += valor;
    }
    mediaCache /= tamanhoDados;
    
    varianciaCache = 0.0;
    for (const auto& valor : serieOriginal) {
        double diff = valor - mediaCache;
        varianciaCache += diff * diff;
    }
    varianciaCache /= (tamanhoDados - 1); // Variância amostral
}

// Média aritmética
double EstatisticasBasicas::media() const {
    verificarDados();
    return mediaCache;
}

// Média ponderada
double EstatisticasBasicas::mediaPonderada(const std::vector<double>& pesos) const {
    verificarDados();
    
    if (pesos.size() != tamanhoDados) {
        throw std::invalid_argument("Número de pesos deve ser igual ao número de dados");
    }
    
    double somaPesos = 0.0;
    double somaPonderada = 0.0;
    
    for (size_t i = 0; i < tamanhoDados; i++) {
        somaPonderada += serieOriginal[i] * pesos[i];
        somaPesos += pesos[i];
    }
    
    return somaPonderada / somaPesos;
}

// Mediana
double EstatisticasBasicas::mediana() const {
    verificarDados();
    
    size_t meio = tamanhoDados / 2;
    if (tamanhoDados % 2 == 0) {
        return (dadosOrdenados[meio - 1] + dadosOrdenados[meio]) / 2.0;
    }
    return dadosOrdenados[meio];
}

// Moda (valor mais frequente - aproximado para dados contínuos)
double EstatisticasBasicas::moda() const {
    verificarDados();
    
    // Usa histograma para encontrar moda em dados contínuos
    int numBins = static_cast<int>(std::sqrt(tamanhoDados));
    if (numBins < 5) numBins = 5;
    
    double minVal = dadosOrdenados.front();
    double maxVal = dadosOrdenados.back();
    double larguraBin = (maxVal - minVal) / numBins;
    
    std::vector<int> frequencias(numBins, 0);
    
    for (const auto& valor : serieOriginal) {
        int bin = static_cast<int>((valor - minVal) / larguraBin);
        if (bin >= numBins) bin = numBins - 1;
        frequencias[bin]++;
    }
    
    auto maxIt = std::max_element(frequencias.begin(), frequencias.end());
    int binModa = std::distance(frequencias.begin(), maxIt);
    
    return minVal + (binModa + 0.5) * larguraBin;
}

// Variância amostral
double EstatisticasBasicas::variancia() const {
    verificarDados();
    return varianciaCache;
}

// Variância populacional
double EstatisticasBasicas::varianciaPopulacional() const {
    verificarDados();
    return varianciaCache * (tamanhoDados - 1) / tamanhoDados;
}

// Desvio padrão amostral
double EstatisticasBasicas::desvioPadrao() const {
    verificarDados();
    return std::sqrt(varianciaCache);
}

// Coeficiente de variação
double EstatisticasBasicas::coeficienteVariacao() const {
    verificarDados();
    if (std::abs(mediaCache) < 1e-10) {
        throw std::runtime_error("Média próxima de zero, CV indefinido");
    }
    return desvioPadrao() / mediaCache;
}

// Erro padrão da média
double EstatisticasBasicas::erroPadrao() const {
    verificarDados();
    return desvioPadrao() / std::sqrt(static_cast<double>(tamanhoDados));
}

// Assimetria (Skewness)
double EstatisticasBasicas::assimetria() const {
    verificarDados();
    
    double dp = desvioPadrao();
    if (dp < 1e-10) return 0.0;
    
    double soma = 0.0;
    for (const auto& valor : serieOriginal) {
        double z = (valor - mediaCache) / dp;
        soma += z * z * z;
    }
    
    double n = static_cast<double>(tamanhoDados);
    return (n / ((n - 1) * (n - 2))) * soma;
}

// Curtose
double EstatisticasBasicas::curtose() const {
    verificarDados();
    
    double dp = desvioPadrao();
    if (dp < 1e-10) return 0.0;
    
    double soma = 0.0;
    for (const auto& valor : serieOriginal) {
        double z = (valor - mediaCache) / dp;
        soma += z * z * z * z;
    }
    
    double n = static_cast<double>(tamanhoDados);
    double termo1 = (n * (n + 1)) / ((n - 1) * (n - 2) * (n - 3));
    double termo2 = (3 * (n - 1) * (n - 1)) / ((n - 2) * (n - 3));
    
    return termo1 * soma - termo2;
}

// Excesso de curtose
double EstatisticasBasicas::excessoCurtose() const {
    return curtose(); // Já calcula excesso de curtose
}

// Percentil genérico
double EstatisticasBasicas::percentil(double p) const {
    verificarDados();
    
    if (p < 0 || p > 100) {
        throw std::invalid_argument("Percentil deve estar entre 0 e 100");
    }
    
    if (p == 0) return dadosOrdenados.front();
    if (p == 100) return dadosOrdenados.back();
    
    double indice = (p / 100.0) * (tamanhoDados - 1);
    size_t inferior = static_cast<size_t>(std::floor(indice));
    size_t superior = static_cast<size_t>(std::ceil(indice));
    
    if (inferior == superior) {
        return dadosOrdenados[inferior];
    }
    
    double fracao = indice - inferior;
    return dadosOrdenados[inferior] * (1 - fracao) + dadosOrdenados[superior] * fracao;
}

// Quartis
double EstatisticasBasicas::quartil1() const { return percentil(25); }
double EstatisticasBasicas::quartil2() const { return percentil(50); }
double EstatisticasBasicas::quartil3() const { return percentil(75); }

// Intervalo interquartil
double EstatisticasBasicas::intervaloInterquartil() const {
    return quartil3() - quartil1();
}

// Valores extremos
double EstatisticasBasicas::minimo() const {
    verificarDados();
    return dadosOrdenados.front();
}

double EstatisticasBasicas::maximo() const {
    verificarDados();
    return dadosOrdenados.back();
}

double EstatisticasBasicas::amplitude() const {
    return maximo() - minimo();
}

// Soma
double EstatisticasBasicas::soma() const {
    verificarDados();
    return std::accumulate(serieOriginal.begin(), serieOriginal.end(), 0.0);
}

// Covariância entre duas séries
double EstatisticasBasicas::covariancia(const std::vector<double>& outraSerie) const {
    verificarDados();
    
    if (outraSerie.size() != tamanhoDados) {
        throw std::invalid_argument("Séries devem ter o mesmo tamanho");
    }
    
    double mediaOutra = 0.0;
    for (const auto& val : outraSerie) {
        mediaOutra += val;
    }
    mediaOutra /= tamanhoDados;
    
    double cov = 0.0;
    for (size_t i = 0; i < tamanhoDados; i++) {
        cov += (serieOriginal[i] - mediaCache) * (outraSerie[i] - mediaOutra);
    }
    
    return cov / (tamanhoDados - 1);
}

// Correlação de Pearson
double EstatisticasBasicas::correlacaoPearson(const std::vector<double>& outraSerie) const {
    verificarDados();
    
    if (outraSerie.size() != tamanhoDados) {
        throw std::invalid_argument("Séries devem ter o mesmo tamanho");
    }
    
    EstatisticasBasicas estatOutra(outraSerie);
    double cov = covariancia(outraSerie);
    double dpProduto = desvioPadrao() * estatOutra.desvioPadrao();
    
    if (dpProduto < 1e-10) {
        return 0.0;
    }
    
    return cov / dpProduto;
}

// Correlação de Spearman (ranks)
double EstatisticasBasicas::correlacaoSpearman(const std::vector<double>& outraSerie) const {
    verificarDados();
    
    if (outraSerie.size() != tamanhoDados) {
        throw std::invalid_argument("Séries devem ter o mesmo tamanho");
    }
    
    // Calcula ranks
    std::vector<double> ranks1 = calcularRanks(serieOriginal);
    std::vector<double> ranks2 = calcularRanks(outraSerie);
    
    EstatisticasBasicas estatRanks1(ranks1);
    return estatRanks1.correlacaoPearson(ranks2);
}

// Calcula ranks para correlação de Spearman
std::vector<double> EstatisticasBasicas::calcularRanks(const std::vector<double>& serie) const {
    size_t n = serie.size();
    std::vector<std::pair<double, size_t>> pareado(n);
    
    for (size_t i = 0; i < n; i++) {
        pareado[i] = {serie[i], i};
    }
    
    std::sort(pareado.begin(), pareado.end());
    
    std::vector<double> ranks(n);
    for (size_t i = 0; i < n; i++) {
        ranks[pareado[i].second] = static_cast<double>(i + 1);
    }
    
    return ranks;
}

// Teste de Jarque-Bera para normalidade
double EstatisticasBasicas::testeJarqueBera() const {
    verificarDados();
    
    double n = static_cast<double>(tamanhoDados);
    double s = assimetria();
    double k = excessoCurtose();
    
    return (n / 6.0) * (s * s + (k * k) / 4.0);
}

// Retornos simples
std::vector<double> EstatisticasBasicas::retornosSimples() const {
    verificarDados();
    
    std::vector<double> retornos(tamanhoDados - 1);
    for (size_t i = 1; i < tamanhoDados; i++) {
        if (std::abs(serieOriginal[i-1]) < 1e-10) {
            retornos[i-1] = 0.0;
        } else {
            retornos[i-1] = (serieOriginal[i] - serieOriginal[i-1]) / serieOriginal[i-1];
        }
    }
    
    return retornos;
}

// Retornos logarítmicos
std::vector<double> EstatisticasBasicas::retornosLogaritmicos() const {
    verificarDados();
    
    std::vector<double> retornos(tamanhoDados - 1);
    for (size_t i = 1; i < tamanhoDados; i++) {
        if (serieOriginal[i] <= 0 || serieOriginal[i-1] <= 0) {
            retornos[i-1] = 0.0;
        } else {
            retornos[i-1] = std::log(serieOriginal[i] / serieOriginal[i-1]);
        }
    }
    
    return retornos;
}

// Verifica se dados foram carregados
void EstatisticasBasicas::verificarDados() const {
    if (!dadosCarregados || tamanhoDados == 0) {
        throw std::runtime_error("Dados não carregados");
    }
}

// Getters
size_t EstatisticasBasicas::tamanho() const { return tamanhoDados; }
const std::vector<double>& EstatisticasBasicas::obterDados() const { return serieOriginal; }

} // namespace analise_ouro
