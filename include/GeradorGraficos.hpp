/*
 * GeradorGraficos.hpp - Header Gráficos - C++23
 * Autor: Luiz Tiago Wilcke
 */

#ifndef GERADOR_GRAFICOS_HPP
#define GERADOR_GRAFICOS_HPP

#include <vector>
#include <string>
#include <fstream>

namespace analise_ouro {

class GeradorGraficos {
private:
    int largura, altura;
    std::string diretorioSaida;
    
    void gerarScriptBase(std::ofstream& script) const;

public:
    GeradorGraficos();
    explicit GeradorGraficos(const std::string& dir);
    ~GeradorGraficos();

    void definirDiretorioSaida(const std::string& dir);
    void definirTamanho(int w, int h);
    
    void plotarSerie(const std::vector<double>& dados, const std::string& titulo, const std::string& nomeArquivo) const;
    void plotarMultiplasSeries(const std::vector<std::vector<double>>& dados, const std::vector<std::string>& legendas, const std::string& titulo, const std::string& nomeArquivo) const;
    void plotarHistograma(const std::vector<double>& dados, int numBins, const std::string& titulo, const std::string& nomeArquivo) const;
    void plotarAcf(const std::vector<double>& acf, const std::string& titulo, const std::string& nomeArquivo) const;
    void plotarCandlestick(const std::vector<double>& abertura, const std::vector<double>& maxima, const std::vector<double>& minima, const std::vector<double>& fechamento, const std::string& titulo, const std::string& nomeArquivo) const;
};

} // namespace analise_ouro

#endif // GERADOR_GRAFICOS_HPP
