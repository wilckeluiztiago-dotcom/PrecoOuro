/*
 * PrevisaoPrecos.hpp - Header Previsão - C++23
 * Autor: Luiz Tiago Wilcke
 */

#ifndef PREVISAO_PRECOS_HPP
#define PREVISAO_PRECOS_HPP

#include <vector>
#include <random>

namespace analise_ouro {

struct ResultadoPrevisao {
    std::vector<double> previsaoMedia;
    std::vector<double> intervaloInferior;
    std::vector<double> intervaloSuperior;
};

struct MetricasErro {
    double me;      // Mean Error
    double mse;     // Mean Squared Error
    double rmse;    // Root Mean Squared Error
    double mae;     // Mean Absolute Error
    double mape;    // Mean Absolute Percentage Error
};

class PrevisaoPrecos {
private:
    std::vector<double> precosOriginais;
    bool dadosCarregados;
    size_t tamanhoSerie;
    mutable std::mt19937_64 gerador;
    
    void verificarDados() const;

public:
    PrevisaoPrecos();
    explicit PrevisaoPrecos(const std::vector<double>& precos);
    ~PrevisaoPrecos();

    void carregarPrecos(const std::vector<double>& precos);
    [[nodiscard]] auto preverMediaMovel(int horizonte, int janela = 20) const -> std::vector<double>;
    [[nodiscard]] auto preverSuavizacaoExponencial(int horizonte, double alfa = 0.3) const -> std::vector<double>;
    [[nodiscard]] auto preverHoltLinear(int horizonte, double alfa = 0.3, double beta = 0.1) const -> std::vector<double>;
    [[nodiscard]] auto preverMonteCarlo(int horizonte, int numSimulacoes = 10000) const -> ResultadoPrevisao;
    [[nodiscard]] auto preverArima(int horizonte, int p = 1, int d = 1, int q = 0) const -> std::vector<double>;
    [[nodiscard]] auto calcularErroPrevisao(const std::vector<double>& reais, const std::vector<double>& previstos) const -> MetricasErro;
};

} // namespace analise_ouro

#endif // PREVISAO_PRECOS_HPP
