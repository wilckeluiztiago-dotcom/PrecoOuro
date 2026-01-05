/*
 * CorrelacaoAtivos.cpp - Correlação do Ouro com Outros Ativos - C++23
 * Autor: Luiz Tiago Wilcke
 * Projeto: Análise do Preço do Ouro
 */

#include "CorrelacaoAtivos.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace analise_ouro {

CorrelacaoAtivos::CorrelacaoAtivos() : dadosCarregados(false) {}
CorrelacaoAtivos::~CorrelacaoAtivos() = default;

void CorrelacaoAtivos::adicionarSerie(const std::string& nome, const std::vector<double>& precos) {
    if (precos.size() < 2) throw std::invalid_argument("Série muito curta");
    
    std::vector<double> retornos(precos.size() - 1);
    for (size_t i = 1; i < precos.size(); ++i) {
        retornos[i - 1] = std::log(precos[i] / precos[i - 1]);
    }
    
    seriesRetornos[nome] = retornos;
    nomesAtivos.push_back(nome);
    dadosCarregados = true;
}

[[nodiscard]] auto CorrelacaoAtivos::correlacaoPearson(const std::string& ativo1, const std::string& ativo2) const -> double {
    verificarDados();
    
    auto it1 = seriesRetornos.find(ativo1);
    auto it2 = seriesRetornos.find(ativo2);
    
    if (it1 == seriesRetornos.end() || it2 == seriesRetornos.end()) {
        throw std::invalid_argument("Ativo não encontrado");
    }
    
    const auto& r1 = it1->second;
    const auto& r2 = it2->second;
    size_t n = std::min(r1.size(), r2.size());
    
    double soma1 = 0.0, soma2 = 0.0, soma12 = 0.0, soma1_2 = 0.0, soma2_2 = 0.0;
    for (size_t i = 0; i < n; ++i) {
        soma1 += r1[i];
        soma2 += r2[i];
        soma12 += r1[i] * r2[i];
        soma1_2 += r1[i] * r1[i];
        soma2_2 += r2[i] * r2[i];
    }
    
    double nf = static_cast<double>(n);
    double num = nf * soma12 - soma1 * soma2;
    double den = std::sqrt((nf * soma1_2 - soma1 * soma1) * (nf * soma2_2 - soma2 * soma2));
    
    return (std::abs(den) > 1e-10) ? num / den : 0.0;
}

[[nodiscard]] auto CorrelacaoAtivos::matrizCorrelacao() const -> std::vector<std::vector<double>> {
    verificarDados();
    
    size_t n = nomesAtivos.size();
    std::vector<std::vector<double>> matriz(n, std::vector<double>(n, 0.0));
    
    for (size_t i = 0; i < n; ++i) {
        matriz[i][i] = 1.0;
        for (size_t j = i + 1; j < n; ++j) {
            double corr = correlacaoPearson(nomesAtivos[i], nomesAtivos[j]);
            matriz[i][j] = corr;
            matriz[j][i] = corr;
        }
    }
    return matriz;
}

[[nodiscard]] auto CorrelacaoAtivos::correlacaoRolante(const std::string& ativo1, const std::string& ativo2, int janela) const -> std::vector<double> {
    verificarDados();
    
    auto it1 = seriesRetornos.find(ativo1);
    auto it2 = seriesRetornos.find(ativo2);
    if (it1 == seriesRetornos.end() || it2 == seriesRetornos.end()) {
        throw std::invalid_argument("Ativo não encontrado");
    }
    
    const auto& r1 = it1->second;
    const auto& r2 = it2->second;
    size_t n = std::min(r1.size(), r2.size());
    
    std::vector<double> correlacoes(n - janela + 1);
    
    for (size_t i = 0; i <= n - janela; ++i) {
        double s1 = 0.0, s2 = 0.0, s12 = 0.0, s1_2 = 0.0, s2_2 = 0.0;
        for (int j = 0; j < janela; ++j) {
            s1 += r1[i + j];
            s2 += r2[i + j];
            s12 += r1[i + j] * r2[i + j];
            s1_2 += r1[i + j] * r1[i + j];
            s2_2 += r2[i + j] * r2[i + j];
        }
        double jf = static_cast<double>(janela);
        double num = jf * s12 - s1 * s2;
        double den = std::sqrt((jf * s1_2 - s1 * s1) * (jf * s2_2 - s2 * s2));
        correlacoes[i] = (std::abs(den) > 1e-10) ? num / den : 0.0;
    }
    return correlacoes;
}

[[nodiscard]] auto CorrelacaoAtivos::beta(const std::string& ativo, const std::string& mercado) const -> double {
    verificarDados();
    
    auto itAtivo = seriesRetornos.find(ativo);
    auto itMercado = seriesRetornos.find(mercado);
    if (itAtivo == seriesRetornos.end() || itMercado == seriesRetornos.end()) {
        throw std::invalid_argument("Ativo não encontrado");
    }
    
    const auto& ra = itAtivo->second;
    const auto& rm = itMercado->second;
    size_t n = std::min(ra.size(), rm.size());
    
    double mediaA = std::accumulate(ra.begin(), ra.begin() + n, 0.0) / n;
    double mediaM = std::accumulate(rm.begin(), rm.begin() + n, 0.0) / n;
    
    double cov = 0.0, varM = 0.0;
    for (size_t i = 0; i < n; ++i) {
        cov += (ra[i] - mediaA) * (rm[i] - mediaM);
        varM += (rm[i] - mediaM) * (rm[i] - mediaM);
    }
    
    return (std::abs(varM) > 1e-10) ? cov / varM : 0.0;
}

[[nodiscard]] auto CorrelacaoAtivos::covariancia(const std::string& ativo1, const std::string& ativo2) const -> double {
    verificarDados();
    
    auto it1 = seriesRetornos.find(ativo1);
    auto it2 = seriesRetornos.find(ativo2);
    if (it1 == seriesRetornos.end() || it2 == seriesRetornos.end()) {
        throw std::invalid_argument("Ativo não encontrado");
    }
    
    const auto& r1 = it1->second;
    const auto& r2 = it2->second;
    size_t n = std::min(r1.size(), r2.size());
    
    double media1 = std::accumulate(r1.begin(), r1.begin() + n, 0.0) / n;
    double media2 = std::accumulate(r2.begin(), r2.begin() + n, 0.0) / n;
    
    double cov = 0.0;
    for (size_t i = 0; i < n; ++i) {
        cov += (r1[i] - media1) * (r2[i] - media2);
    }
    return cov / (n - 1);
}

[[nodiscard]] auto CorrelacaoAtivos::matrizCovariancia() const -> std::vector<std::vector<double>> {
    verificarDados();
    
    size_t n = nomesAtivos.size();
    std::vector<std::vector<double>> matriz(n, std::vector<double>(n, 0.0));
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            double cov = covariancia(nomesAtivos[i], nomesAtivos[j]);
            matriz[i][j] = cov;
            matriz[j][i] = cov;
        }
    }
    return matriz;
}

void CorrelacaoAtivos::verificarDados() const {
    if (!dadosCarregados || nomesAtivos.empty()) throw std::runtime_error("Dados não carregados");
}

[[nodiscard]] auto CorrelacaoAtivos::obterNomesAtivos() const -> const std::vector<std::string>& { return nomesAtivos; }

} // namespace analise_ouro
