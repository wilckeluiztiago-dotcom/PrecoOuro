/*
 * RelatorioNumerico.hpp - Header Relatórios - C++23
 * Autor: Luiz Tiago Wilcke
 */

#ifndef RELATORIO_NUMERICO_HPP
#define RELATORIO_NUMERICO_HPP

#include <vector>
#include <string>

namespace analise_ouro {

class RelatorioNumerico {
private:
    std::string diretorioSaida;
    [[nodiscard]] auto obterDataHora() const -> std::string;

public:
    RelatorioNumerico();
    explicit RelatorioNumerico(const std::string& dir);
    ~RelatorioNumerico();

    void definirDiretorioSaida(const std::string& dir);
    void gerarRelatorioEstatisticas(const std::vector<double>& dados, const std::string& nomeArquivo) const;
    void gerarRelatorioRisco(double var95, double var99, double cvar95, double sharpe, double sortino, double maxDD, const std::string& nomeArquivo) const;
    void gerarRelatorioModelo(const std::string& nomeModelo, const std::vector<double>& parametros, const std::vector<std::string>& nomesParam, double aic, double bic, const std::string& nomeArquivo) const;
    void gerarRelatorioPrevisao(const std::vector<double>& previsoes, const std::vector<double>& icInferior, const std::vector<double>& icSuperior, const std::string& nomeArquivo) const;
};

} // namespace analise_ouro

#endif // RELATORIO_NUMERICO_HPP
