/*
 * MatrizOperacoes.cpp
 * Módulo de operações matriciais para cálculos financeiros
 * 
 * Autor: Luiz Tiago Wilcke
 * Projeto: Análise do Preço do Ouro
 */

#include "MatrizOperacoes.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <sstream>

namespace analise_ouro {

// Construtor padrão - matriz vazia
MatrizOperacoes::MatrizOperacoes() : linhas(0), colunas(0) {}

// Construtor com dimensões
MatrizOperacoes::MatrizOperacoes(size_t numLinhas, size_t numColunas, double valorInicial) 
    : linhas(numLinhas), colunas(numColunas) {
    dados.resize(linhas, std::vector<double>(colunas, valorInicial));
}

// Construtor com dados
MatrizOperacoes::MatrizOperacoes(const std::vector<std::vector<double>>& matrizDados) {
    if (matrizDados.empty()) {
        linhas = 0;
        colunas = 0;
        return;
    }
    
    linhas = matrizDados.size();
    colunas = matrizDados[0].size();
    
    // Verifica consistência das dimensões
    for (const auto& linha : matrizDados) {
        if (linha.size() != colunas) {
            throw std::invalid_argument("Matriz com linhas de tamanhos diferentes");
        }
    }
    
    dados = matrizDados;
}

// Destrutor
MatrizOperacoes::~MatrizOperacoes() {}

// Cria matriz identidade
MatrizOperacoes MatrizOperacoes::identidade(size_t n) {
    MatrizOperacoes resultado(n, n, 0.0);
    for (size_t i = 0; i < n; i++) {
        resultado.dados[i][i] = 1.0;
    }
    return resultado;
}

// Cria matriz de zeros
MatrizOperacoes MatrizOperacoes::zeros(size_t numLinhas, size_t numColunas) {
    return MatrizOperacoes(numLinhas, numColunas, 0.0);
}

// Cria matriz de uns
MatrizOperacoes MatrizOperacoes::uns(size_t numLinhas, size_t numColunas) {
    return MatrizOperacoes(numLinhas, numColunas, 1.0);
}

// Acesso aos elementos
double& MatrizOperacoes::operator()(size_t i, size_t j) {
    if (i >= linhas || j >= colunas) {
        throw std::out_of_range("Índice fora dos limites da matriz");
    }
    return dados[i][j];
}

const double& MatrizOperacoes::operator()(size_t i, size_t j) const {
    if (i >= linhas || j >= colunas) {
        throw std::out_of_range("Índice fora dos limites da matriz");
    }
    return dados[i][j];
}

// Soma de matrizes
MatrizOperacoes MatrizOperacoes::operator+(const MatrizOperacoes& outra) const {
    if (linhas != outra.linhas || colunas != outra.colunas) {
        throw std::invalid_argument("Dimensões incompatíveis para soma");
    }
    
    MatrizOperacoes resultado(linhas, colunas);
    for (size_t i = 0; i < linhas; i++) {
        for (size_t j = 0; j < colunas; j++) {
            resultado.dados[i][j] = dados[i][j] + outra.dados[i][j];
        }
    }
    return resultado;
}

// Subtração de matrizes
MatrizOperacoes MatrizOperacoes::operator-(const MatrizOperacoes& outra) const {
    if (linhas != outra.linhas || colunas != outra.colunas) {
        throw std::invalid_argument("Dimensões incompatíveis para subtração");
    }
    
    MatrizOperacoes resultado(linhas, colunas);
    for (size_t i = 0; i < linhas; i++) {
        for (size_t j = 0; j < colunas; j++) {
            resultado.dados[i][j] = dados[i][j] - outra.dados[i][j];
        }
    }
    return resultado;
}

// Multiplicação de matrizes
MatrizOperacoes MatrizOperacoes::operator*(const MatrizOperacoes& outra) const {
    if (colunas != outra.linhas) {
        throw std::invalid_argument("Dimensões incompatíveis para multiplicação");
    }
    
    MatrizOperacoes resultado(linhas, outra.colunas, 0.0);
    for (size_t i = 0; i < linhas; i++) {
        for (size_t j = 0; j < outra.colunas; j++) {
            for (size_t k = 0; k < colunas; k++) {
                resultado.dados[i][j] += dados[i][k] * outra.dados[k][j];
            }
        }
    }
    return resultado;
}

// Multiplicação por escalar
MatrizOperacoes MatrizOperacoes::operator*(double escalar) const {
    MatrizOperacoes resultado(linhas, colunas);
    for (size_t i = 0; i < linhas; i++) {
        for (size_t j = 0; j < colunas; j++) {
            resultado.dados[i][j] = dados[i][j] * escalar;
        }
    }
    return resultado;
}

// Transposta
MatrizOperacoes MatrizOperacoes::transposta() const {
    MatrizOperacoes resultado(colunas, linhas);
    for (size_t i = 0; i < linhas; i++) {
        for (size_t j = 0; j < colunas; j++) {
            resultado.dados[j][i] = dados[i][j];
        }
    }
    return resultado;
}

// Decomposição LU
void MatrizOperacoes::decomposicaoLU(MatrizOperacoes& L, MatrizOperacoes& U) const {
    if (linhas != colunas) {
        throw std::invalid_argument("Matriz deve ser quadrada para decomposição LU");
    }
    
    size_t n = linhas;
    L = zeros(n, n);
    U = zeros(n, n);
    
    for (size_t i = 0; i < n; i++) {
        // Matriz U superior
        for (size_t k = i; k < n; k++) {
            double soma = 0.0;
            for (size_t j = 0; j < i; j++) {
                soma += L.dados[i][j] * U.dados[j][k];
            }
            U.dados[i][k] = dados[i][k] - soma;
        }
        
        // Matriz L inferior
        for (size_t k = i; k < n; k++) {
            if (i == k) {
                L.dados[i][i] = 1.0;
            } else {
                double soma = 0.0;
                for (size_t j = 0; j < i; j++) {
                    soma += L.dados[k][j] * U.dados[j][i];
                }
                if (std::abs(U.dados[i][i]) < 1e-10) {
                    throw std::runtime_error("Divisão por zero na decomposição LU");
                }
                L.dados[k][i] = (dados[k][i] - soma) / U.dados[i][i];
            }
        }
    }
}

// Decomposição de Cholesky
MatrizOperacoes MatrizOperacoes::decomposicaoCholesky() const {
    if (linhas != colunas) {
        throw std::invalid_argument("Matriz deve ser quadrada para Cholesky");
    }
    
    size_t n = linhas;
    MatrizOperacoes L = zeros(n, n);
    
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j <= i; j++) {
            double soma = 0.0;
            
            if (j == i) {
                for (size_t k = 0; k < j; k++) {
                    soma += L.dados[j][k] * L.dados[j][k];
                }
                double valor = dados[j][j] - soma;
                if (valor < 0) {
                    throw std::runtime_error("Matriz não é positiva definida");
                }
                L.dados[j][j] = std::sqrt(valor);
            } else {
                for (size_t k = 0; k < j; k++) {
                    soma += L.dados[i][k] * L.dados[j][k];
                }
                if (std::abs(L.dados[j][j]) < 1e-10) {
                    throw std::runtime_error("Divisão por zero em Cholesky");
                }
                L.dados[i][j] = (dados[i][j] - soma) / L.dados[j][j];
            }
        }
    }
    
    return L;
}

// Determinante
double MatrizOperacoes::determinante() const {
    if (linhas != colunas) {
        throw std::invalid_argument("Matriz deve ser quadrada para calcular determinante");
    }
    
    size_t n = linhas;
    if (n == 1) return dados[0][0];
    if (n == 2) return dados[0][0] * dados[1][1] - dados[0][1] * dados[1][0];
    
    // Usando decomposição LU
    MatrizOperacoes L, U;
    decomposicaoLU(L, U);
    
    double det = 1.0;
    for (size_t i = 0; i < n; i++) {
        det *= U.dados[i][i];
    }
    
    return det;
}

// Inversa usando eliminação de Gauss-Jordan
MatrizOperacoes MatrizOperacoes::inversa() const {
    if (linhas != colunas) {
        throw std::invalid_argument("Matriz deve ser quadrada para inversão");
    }
    
    size_t n = linhas;
    MatrizOperacoes aumentada(n, 2 * n);
    
    // Copia matriz original e identidade
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            aumentada.dados[i][j] = dados[i][j];
            aumentada.dados[i][j + n] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    // Eliminação de Gauss-Jordan
    for (size_t i = 0; i < n; i++) {
        // Pivoteamento parcial
        size_t maxLinha = i;
        for (size_t k = i + 1; k < n; k++) {
            if (std::abs(aumentada.dados[k][i]) > std::abs(aumentada.dados[maxLinha][i])) {
                maxLinha = k;
            }
        }
        std::swap(aumentada.dados[i], aumentada.dados[maxLinha]);
        
        if (std::abs(aumentada.dados[i][i]) < 1e-10) {
            throw std::runtime_error("Matriz singular, não pode ser invertida");
        }
        
        // Normaliza linha do pivô
        double pivo = aumentada.dados[i][i];
        for (size_t j = 0; j < 2 * n; j++) {
            aumentada.dados[i][j] /= pivo;
        }
        
        // Elimina outras linhas
        for (size_t k = 0; k < n; k++) {
            if (k != i) {
                double fator = aumentada.dados[k][i];
                for (size_t j = 0; j < 2 * n; j++) {
                    aumentada.dados[k][j] -= fator * aumentada.dados[i][j];
                }
            }
        }
    }
    
    // Extrai inversa
    MatrizOperacoes resultado(n, n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            resultado.dados[i][j] = aumentada.dados[i][j + n];
        }
    }
    
    return resultado;
}

// Traço da matriz
double MatrizOperacoes::traco() const {
    if (linhas != colunas) {
        throw std::invalid_argument("Matriz deve ser quadrada para calcular traço");
    }
    
    double soma = 0.0;
    for (size_t i = 0; i < linhas; i++) {
        soma += dados[i][i];
    }
    return soma;
}

// Norma de Frobenius
double MatrizOperacoes::normaFrobenius() const {
    double soma = 0.0;
    for (size_t i = 0; i < linhas; i++) {
        for (size_t j = 0; j < colunas; j++) {
            soma += dados[i][j] * dados[i][j];
        }
    }
    return std::sqrt(soma);
}

// Resolve sistema linear Ax = b
std::vector<double> MatrizOperacoes::resolverSistema(const std::vector<double>& b) const {
    if (linhas != colunas || linhas != b.size()) {
        throw std::invalid_argument("Dimensões incompatíveis para sistema linear");
    }
    
    MatrizOperacoes inv = inversa();
    std::vector<double> resultado(linhas, 0.0);
    
    for (size_t i = 0; i < linhas; i++) {
        for (size_t j = 0; j < colunas; j++) {
            resultado[i] += inv.dados[i][j] * b[j];
        }
    }
    
    return resultado;
}

// Getters
size_t MatrizOperacoes::obterLinhas() const { return linhas; }
size_t MatrizOperacoes::obterColunas() const { return colunas; }

// Representação em string
std::string MatrizOperacoes::paraString(int precisao) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precisao);
    
    for (size_t i = 0; i < linhas; i++) {
        oss << "[ ";
        for (size_t j = 0; j < colunas; j++) {
            oss << std::setw(10) << dados[i][j] << " ";
        }
        oss << "]\n";
    }
    
    return oss.str();
}

} // namespace analise_ouro
