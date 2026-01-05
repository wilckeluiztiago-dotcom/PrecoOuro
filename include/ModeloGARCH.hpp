/*
 * ModeloGARCH.hpp - Header GARCH - C++23
 * Autor: Luiz Tiago Wilcke
 */

#ifndef MODELO_GARCH_HPP
#define MODELO_GARCH_HPP

#include <vector>
#include <utility>

namespace analise_ouro {

class ModeloGARCH {
private:
    std::vector<double> serieOriginal;
    std::vector<double> retornos;
    std::vector<double> volatilidades;
    std::vector<double> alfa;
    std::vector<double> beta;
    double omega;
    double mediaRetornos;
    bool serieCarregada;
    size_t tamanhoSerie;
    int ordemP, ordemQ;
    bool estimado;
    
    [[nodiscard]] auto calcularRetornos() const -> std::vector<double>;
    [[nodiscard]] auto calcularLogVerossimilhanca() const -> std::pair<double, std::vector<double>>;
    void verificarDados() const;

public:
    ModeloGARCH();
    ModeloGARCH(const std::vector<double>& dados, int p = 1, int q = 1);
    ~ModeloGARCH();

    void carregarDados(const std::vector<double>& novosDados);
    void definirOrdens(int p, int q);
    [[nodiscard]] auto estimar() -> bool;
    [[nodiscard]] auto preverVolatilidade(int horizontePrevisao) const -> std::vector<double>;
    [[nodiscard]] auto obterVolatilidades() const -> const std::vector<double>&;
    [[nodiscard]] auto obterOmega() const noexcept -> double;
    [[nodiscard]] auto obterAlfa() const -> const std::vector<double>&;
    [[nodiscard]] auto obterBeta() const -> const std::vector<double>&;
};

} // namespace analise_ouro

#endif // MODELO_GARCH_HPP
