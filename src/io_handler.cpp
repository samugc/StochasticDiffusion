#include "io_handler.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>

using namespace std;
using namespace Eigen;

// LEER DATOS

vector<vector<double>> readCSV(const string& filename) {
    vector<vector<double>> data;  // Almacenará el contenido del CSV
    ifstream file(filename);  // Abrir archivo CSV
    string line;  // Línea que leemos del archivo

    // Verificamos si el archivo está abierto correctamente
    if (!file.is_open()) {
        cerr << "No se puede abrir el archivo: " << filename << endl;
        return data;
    }

    bool isFirstLine = true;  // Variable para identificar la primera línea

    // Leer cada línea del archivo
    while (getline(file, line)) {
        // Saltar la primera línea (encabezado)
        if (isFirstLine) {
            isFirstLine = false;
            continue;
        }

        stringstream ss(line);  // Usamos un stringstream para dividir la línea
        string value;  // Valor de cada campo en la línea
        vector<double> row;  // Una fila del CSV

        // Dividimos la línea por comas y almacenamos los valores
        while (getline(ss, value, ';')) {
            size_t pos = value.find(','); // Encuentra la posición de la coma
            if (pos != string::npos) { // Si la coma se encuentra
                value.replace(pos, 1, "."); // Reemplaza la coma por un punto
            }
            row.push_back(stod(value));
        }

        // Almacenamos la fila en el vector de datos
        data.push_back(row);
    }

    // Cerramos el archivo
    file.close();
    return data;
}

// MANIPULAR VECTORES

VectorXd createEvaluatedVector(int tam, double value) {
    VectorXd v(tam + 1); 
    v(0) = 0;            
    for (int i = 1; i <= tam; ++i) {
        v(i) = value;    
    }
    return v;
}

VectorXd removeFirst(const VectorXd& vec) {
    return vec.tail(vec.size() - 1);
}

VectorXd removeLast(const Eigen::VectorXd& vec) {
    return vec.head(vec.size() - 1);
}

// MANIPULAR MATRIZ

MatrixXd convertToMatrix(const vector<vector<double>>& data) {
    // Creamos una matriz Eigen de tamaño adecuado
    MatrixXd matrix(data.size(), data[0].size());

    // Llenamos la matriz con los datos
    for (size_t i = 0; i < data.size(); ++i) {
        for (size_t j = 0; j < data[i].size(); ++j) {
            matrix(i, j) = data[i][j];
        }
    }
    
    return matrix;
}

// IMPRIMIR

void printMatrix(const MatrixXd& matrix) {
    cout <<  setprecision(10) << matrix << endl;
}

void printVector(const Eigen::VectorXd& v) {
    cout << "[ ";
    copy(v.data(), v.data() + v.size(), ostream_iterator<double>(cout, ", "));
    cout << "\b\b ]" << endl; 
}

void printRows(const MatrixXd& matrix, int numRows) {
    for (int i = 0; i < numRows && i < matrix.rows(); ++i) {
        cout <<  setprecision(10) << matrix.row(i) << endl;  // Imprime la i-ésima fila
    }
}

void printVectorRange(const vector<vector<double>>& data, size_t endRow, size_t endCol) {
    for (size_t i = 0; i <endRow && i < data.size(); ++i) {
        for (size_t j = 0; j < endCol && j < data[i].size(); ++j) {
            cout <<  setprecision(10) << data[i][j] << " ";  // Imprimir el valor en la fila i y columna j
        }
        cout << endl; 
    }
}

void printPoly(const VectorXd& poly) {
    cout << "f(x) = ";
    for (int i = poly.size() - 1; i >= 0; --i) {
        double coef = poly[i];

        // Mostrar el signo
        if (i != poly.size() - 1) {
            cout << (coef >= 0 ? " + " : " - ");
        } else if (coef < 0) {
            cout << "-";
        }

        // Valor absoluto del coeficiente
        double absCoef = abs(coef);

        // No mostrar el 1 en x^n
        if (!(absCoef == 1.0 && i > 0)) {
            cout << absCoef;
        }

        // Mostrar la x^n
        if (i > 1) {
            cout << "*x^" << i;
        } else if (i == 1) {
            cout << "*x";
        }
    }
    cout << endl;
}
