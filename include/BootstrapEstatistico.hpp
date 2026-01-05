/*
 * BootstrapEstatistico.hpp - Header Bootstrap - C++23
 * Autor: Luiz Tiago Wilcke
 */

#ifndef BOOTSTRAP_ESTATISTICO_HPP
#define BOOTSTRAP_ESTATISTICO_HPP

#include <vector>
#include <random>
#include <functional>
#include <utility>

namespace analise_ouro {

using FuncaoEstatistica = std::function<double(const std::vector<double>&)>;

class BootstrapEstatistico {
private:
    int numReamostras;
    mutable std::mt19937_64 gerador;

    [[nodiscard]] auto gerarAmostra(const std::vector<double>& dados) const -> std::vector<double>;

public:
    BootstrapEstatistico();
    explicit BootstrapEstatistico(int nReamostras);
    ~BootstrapEstatistico();

    void definirNumReamostras(int n);
    [[nodiscard]] auto intervaloConfiancaMedia(const std::vector<double>& dados, double nivelConfianca = 0.95) const -> std::pair<double, double>;
    [[nodiscard]] auto intervaloConfiancaVariancia(const std::vector<double>& dados, double nivelConfianca = 0.95) const -> std::pair<double, double>;
    [[nodiscard]] auto intervaloConfiancaMediana(const std::vector<double>& dados, double nivelConfianca = 0.95) const -> std::pair<double, double>;
    [[nodiscard]] auto intervaloConfianca(const std::vector<double>& dados, FuncaoEstatistica estatistica, double nivelConfianca = 0.95) const -> std::pair<double, double>;
    [[nodiscard]] auto erroPadraoBootstrap(const std::vector<double>& dados, FuncaoEstatistica estatistica) const -> double;
    [[nodiscard]] auto testeHipotese(const std::vector<double>& dados, double valorHipotese, FuncaoEstatistica estatistica) const -> double;
    [[nodiscard]] auto obterNumReamostras() const noexcept -> int;
};

} // namespace analise_ouro

#endif // BOOTSTRAP_ESTATISTICO_HPP
