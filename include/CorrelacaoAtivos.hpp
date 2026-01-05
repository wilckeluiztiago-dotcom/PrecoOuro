/*
 * CorrelacaoAtivos.hpp - Header Correlação - C++23
 * Autor: Luiz Tiago Wilcke
 */

#ifndef CORRELACAO_ATIVOS_HPP
#define CORRELACAO_ATIVOS_HPP

#include <vector>
#include <string>
#include <map>

namespace analise_ouro {

class CorrelacaoAtivos {
private:
    std::map<std::string, std::vector<double>> seriesRetornos;
    std::vector<std::string> nomesAtivos;
    bool dadosCarregados;
    
    void verificarDados() const;

public:
    CorrelacaoAtivos();
    ~CorrelacaoAtivos();

    void adicionarSerie(const std::string& nome, const std::vector<double>& precos);
    [[nodiscard]] auto correlacaoPearson(const std::string& ativo1, const std::string& ativo2) const -> double;
    [[nodiscard]] auto matrizCorrelacao() const -> std::vector<std::vector<double>>;
    [[nodiscard]] auto correlacaoRolante(const std::string& ativo1, const std::string& ativo2, int janela) const -> std::vector<double>;
    [[nodiscard]] auto beta(const std::string& ativo, const std::string& mercado) const -> double;
    [[nodiscard]] auto covariancia(const std::string& ativo1, const std::string& ativo2) const -> double;
    [[nodiscard]] auto matrizCovariancia() const -> std::vector<std::vector<double>>;
    [[nodiscard]] auto obterNomesAtivos() const -> const std::vector<std::string>&;
};

} // namespace analise_ouro

#endif // CORRELACAO_ATIVOS_HPP
