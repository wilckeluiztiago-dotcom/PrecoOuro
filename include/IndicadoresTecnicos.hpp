/*
 * IndicadoresTecnicos.hpp - Header Indicadores - C++23
 * Autor: Luiz Tiago Wilcke
 */

#ifndef INDICADORES_TECNICOS_HPP
#define INDICADORES_TECNICOS_HPP

#include <vector>
#include <tuple>
#include <utility>

namespace analise_ouro {

class IndicadoresTecnicos {
private:
    std::vector<double> precosOriginais;
    bool dadosCarregados;
    size_t tamanhoSerie;
    
    void verificarDados() const;

public:
    IndicadoresTecnicos();
    explicit IndicadoresTecnicos(const std::vector<double>& precos);
    ~IndicadoresTecnicos();

    void carregarPrecos(const std::vector<double>& precos);
    [[nodiscard]] auto rsi(int periodo = 14) const -> std::vector<double>;
    [[nodiscard]] auto macd(int periodoRapido = 12, int periodoLento = 26, int periodoSinal = 9) const -> std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>;
    [[nodiscard]] auto bollingerBands(int periodo = 20, double numDesvios = 2.0) const -> std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>;
    [[nodiscard]] auto stochastic(int periodoPrincipal = 14, int periodoSuavizacao = 3) const -> std::pair<std::vector<double>, std::vector<double>>;
    [[nodiscard]] auto atr(int periodo = 14) const -> std::vector<double>;
    [[nodiscard]] auto ema(int periodo) const -> std::vector<double>;
    [[nodiscard]] auto adx(int periodo = 14) const -> std::vector<double>;
};

} // namespace analise_ouro

#endif // INDICADORES_TECNICOS_HPP
