/*
 * MatrizOperacoes.hpp
 * Header do módulo de operações matriciais
 * 
 * Autor: Luiz Tiago Wilcke
 */

#ifndef MATRIZ_OPERACOES_HPP
#define MATRIZ_OPERACOES_HPP

#include <vector>
#include <string>

namespace analise_ouro {

class MatrizOperacoes {
private:
    std::vector<std::vector<double>> dados;
    size_t linhas;
    size_t colunas;

public:
    // Construtores e destrutor
    MatrizOperacoes();
    MatrizOperacoes(size_t numLinhas, size_t numColunas, double valorInicial = 0.0);
    explicit MatrizOperacoes(const std::vector<std::vector<double>>& matrizDados);
    ~MatrizOperacoes();

    // Métodos estáticos para criação
    static MatrizOperacoes identidade(size_t n);
    static MatrizOperacoes zeros(size_t numLinhas, size_t numColunas);
    static MatrizOperacoes uns(size_t numLinhas, size_t numColunas);

    // Acesso aos elementos
    double& operator()(size_t i, size_t j);
    const double& operator()(size_t i, size_t j) const;

    // Operações básicas
    MatrizOperacoes operator+(const MatrizOperacoes& outra) const;
    MatrizOperacoes operator-(const MatrizOperacoes& outra) const;
    MatrizOperacoes operator*(const MatrizOperacoes& outra) const;
    MatrizOperacoes operator*(double escalar) const;
    MatrizOperacoes transposta() const;

    // Decomposições
    void decomposicaoLU(MatrizOperacoes& L, MatrizOperacoes& U) const;
    MatrizOperacoes decomposicaoCholesky() const;

    // Propriedades
    double determinante() const;
    MatrizOperacoes inversa() const;
    double traco() const;
    double normaFrobenius() const;

    // Sistemas lineares
    std::vector<double> resolverSistema(const std::vector<double>& b) const;

    // Getters
    size_t obterLinhas() const;
    size_t obterColunas() const;

    // Utilidades
    std::string paraString(int precisao = 4) const;
};

} // namespace analise_ouro

#endif // MATRIZ_OPERACOES_HPP
