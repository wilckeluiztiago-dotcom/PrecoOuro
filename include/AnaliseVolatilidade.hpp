/*
 * AnaliseVolatilidade.hpp - Header Volatilidade - C++23
 * Autor: Luiz Tiago Wilcke
 */

#ifndef ANALISE_VOLATILIDADE_HPP
#define ANALISE_VOLATILIDADE_HPP

#include <vector>
#include <tuple>

namespace analise_ouro {

class AnaliseVolatilidade {
private:
    std::vector<double> precosOriginais;
    std::vector<double> retornosLog;
    bool serieCarregada;
    size_t tamanhoSerie;
    
    void calcularRetornos();
    void verificarDados() const;

public:
    AnaliseVolatilidade();
    explicit AnaliseVolatilidade(const std::vector<double>& precos);
    ~AnaliseVolatilidade();

    void carregarPrecos(const std::vector<double>& precos);
    [[nodiscard]] auto volatilidadeHistorica(int janela = 20) const -> std::vector<double>;
    [[nodiscard]] auto volatilidadeRealizada(int janela = 20) const -> double;
    [[nodiscard]] auto volatilidadeEWMA(double lambda = 0.94) const -> std::vector<double>;
    [[nodiscard]] auto volatilidadeParkinson(int janela = 20) const -> double;
    [[nodiscard]] auto volatilidadeGarmanKlass() const -> double;
    [[nodiscard]] auto detectarClusteringVolatilidade() const -> std::vector<int>;
    [[nodiscard]] auto calcularConeVolatilidade(const std::vector<int>& janelas) const -> std::vector<std::tuple<int, double, double, double>>;
    [[nodiscard]] auto calcularVolatilidadeImplicita(double precoOpcao, double precoAtivo, double strike, double r, double T, bool ehCall) const -> double;
    [[nodiscard]] auto obterRetornos() const -> const std::vector<double>&;
};

} // namespace analise_ouro

#endif // ANALISE_VOLATILIDADE_HPP
