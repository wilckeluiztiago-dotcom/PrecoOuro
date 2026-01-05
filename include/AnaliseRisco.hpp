/*
 * AnaliseRisco.hpp - Header Risco - C++23
 * Autor: Luiz Tiago Wilcke
 */

#ifndef ANALISE_RISCO_HPP
#define ANALISE_RISCO_HPP

#include <vector>

namespace analise_ouro {

class AnaliseRisco {
private:
    std::vector<double> precosOriginais;
    std::vector<double> retornos;
    bool dadosCarregados;
    size_t tamanhoSerie;
    
    void calcularRetornos();
    void verificarDados() const;

public:
    AnaliseRisco();
    explicit AnaliseRisco(const std::vector<double>& precos);
    ~AnaliseRisco();

    void carregarPrecos(const std::vector<double>& precos);
    [[nodiscard]] auto valorEmRisco(double nivelConfianca = 0.95, int janela = 252) const -> double;
    [[nodiscard]] auto valorEmRiscoCondicional(double nivelConfianca = 0.95, int janela = 252) const -> double;
    [[nodiscard]] auto valorEmRiscoParametrico(double nivelConfianca = 0.95) const -> double;
    [[nodiscard]] auto indiceSharte(double taxaLivreRisco = 0.02) const -> double;
    [[nodiscard]] auto indiceSortino(double taxaLivreRisco = 0.02) const -> double;
    [[nodiscard]] auto maxDrawdown() const -> double;
    [[nodiscard]] auto calcularDrawdowns() const -> std::vector<double>;
    [[nodiscard]] auto indiceCalmar(double taxaLivreRisco = 0.02) const -> double;
    [[nodiscard]] auto betaModificado(const std::vector<double>& retornosMercado) const -> double;
    [[nodiscard]] auto informationRatio(const std::vector<double>& retornosBenchmark) const -> double;
    [[nodiscard]] auto obterRetornos() const -> const std::vector<double>&;
};

} // namespace analise_ouro

#endif // ANALISE_RISCO_HPP
