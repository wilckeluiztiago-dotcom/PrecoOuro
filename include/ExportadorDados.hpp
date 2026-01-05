/*
 * ExportadorDados.hpp - Header Exportação - C++23
 * Autor: Luiz Tiago Wilcke
 */

#ifndef EXPORTADOR_DADOS_HPP
#define EXPORTADOR_DADOS_HPP

#include <vector>
#include <string>
#include <map>

namespace analise_ouro {

class ExportadorDados {
private:
    std::string diretorioSaida;
    std::string separador;

public:
    ExportadorDados();
    explicit ExportadorDados(const std::string& dir);
    ~ExportadorDados();

    void definirDiretorioSaida(const std::string& dir);
    void definirSeparador(const std::string& sep);
    
    void exportarCSV(const std::vector<double>& dados, const std::string& cabecalho, const std::string& nomeArquivo) const;
    void exportarCSVMultiplo(const std::vector<std::vector<double>>& dados, const std::vector<std::string>& cabecalhos, const std::string& nomeArquivo) const;
    void exportarMatrizCSV(const std::vector<std::vector<double>>& matriz, const std::vector<std::string>& cabecalhosLinha, const std::vector<std::string>& cabecalhosColuna, const std::string& nomeArquivo) const;
    void exportarJSON(const std::vector<double>& dados, const std::string& nomeVariavel, const std::string& nomeArquivo) const;
    void exportarJSONCompleto(const std::map<std::string, std::vector<double>>& dados, const std::string& nomeArquivo) const;
    void exportarGnuplotDat(const std::vector<double>& x, const std::vector<double>& y, const std::string& nomeArquivo) const;
    void exportarResumoTexto(const std::map<std::string, double>& metricas, const std::string& titulo, const std::string& nomeArquivo) const;
};

} // namespace analise_ouro

#endif // EXPORTADOR_DADOS_HPP
