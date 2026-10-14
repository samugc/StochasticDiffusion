#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <filesystem>
#include "io_handler.hpp"
#include "model_data.hpp"
#include "model_functions.hpp"
#include "visualization.hpp"
#include "regression.hpp"
#include "estimations.hpp"
#include "equations.hpp"
#include "optimization.hpp"
#include "evolutive.hpp"

using namespace std;
using namespace Eigen;

/**************************************************************************************************/
//                             OBTENEMOS EL TIEMPO Y LA VARIABLE X                                //
/**************************************************************************************************/

// Necesitamos añadir una fila y columna de 0
MatrixXd add0(const MatrixXd& matrix){
    
    int rows = matrix.rows();
    int cols = matrix.cols();

    // Crear una nueva matriz con una fila y una columna extra (inicializadas en cero)
    MatrixXd extended = MatrixXd::Zero(rows + 1, cols + 1);

    // Copiar la matriz original dentro de la nueva (desplazada una fila y una columna)
    extended.block(1, 1, rows, cols) = matrix;

    return extended;
}

// Obtenemos la matriz del tiempo
MatrixXd obtainTime(const MatrixXd& transposedMatrix){

    // Obtener dimensiones de la matriz transpuesta
    int filas = transposedMatrix.rows();

    // Extraemos los valores del tiempo
    RowVectorXd t_values = transposedMatrix.row(0);

    // Extendemos la fila de tiempos el número de trayectorias que haya
    MatrixXd t_matrix = t_values.replicate(filas, 1);
    
    // We add a row and column of 0
    MatrixXd t_extended = add0(t_matrix);

    return t_extended;
}

// Obtenemos la matriz de la variable aleatoria X
MatrixXd obtainX(const MatrixXd& transposedMatrix){

    // Eliminamos la primera fila de la matriz transpuesta
    MatrixXd x_matrix = transposedMatrix.bottomRows(transposedMatrix.rows() - 1);

    // We add a row and column of 0
    MatrixXd x_extended = add0(x_matrix);

    return x_extended;
}

/**************************************************************************************************/
//                                      VARIABLE ALEATORIA V                                      //
/**************************************************************************************************/

// Obtenemos la matriz de la variable aleatoria V
MatrixXd obtainV(){
    // Crear la matriz `v` de tamaño (d, n)
    MatrixXd v = MatrixXd::Zero(d + 1, n(1));

    // Llenar la primera fila de `v`
    for (int i = 1; i <= d; i++) {
        v(0, i) = x(i, 1);
    }

    // Llenar el resto de `v`
    for (int i = 1; i <= d; i++) {
        for (int j = 1; j <= (n(i) - 1); j++) {  
            v(i, j) = (1.0/pow(delta<double>(i, j + 1, j), 0.5)) * log(x(i, j + 1) / x(i, j));
        }
    }

    return v;
}

/**************************************************************************************************/
//                                           LÍMITES                                              //
/**************************************************************************************************/

// Calculamos los límites
void computeBoundsFromData(
    int degree,
    MatrixXd& beta_bounds,   
    Vector2d& eta_bounds,    
    Vector2d& sigma2_bounds, 
    const VectorXd& estimationEtaBeta
) {
    VectorXd m_arith = arithmeticMean();

    double m1 = m_arith(1);
    double mN = m_arith(N);
    double eta_hat = 1.0 / (mN / m1 - 1.0);

    // Construir los vectores para regresión polinómica
    VectorXd y(N-2);
    y(0)=0;
    for (int j = 1; j <= N-3; ++j) {
        double ratio = (mN / m_arith(j) - 1.0) * eta_hat;
        y(j) = -log(ratio);
    }

    VectorXd newY = removeFirst(y);

    VectorXd newT2 = removeLast(t_values);
    VectorXd newT1 = removeLast(newT2);
    VectorXd newT = removeLast(newT1);

    // Regresión polinómica
    VectorXd beta_est(degree);
    beta_est = polynomialRegression(newT, newY, degree);

    beta_bounds = MatrixXd(2, p);
    for (int j = 0; j < p; ++j) {
        double val = estimationEtaBeta[j + 1]; // salto η
        beta_bounds(0, j) = val * 0.9;
        beta_bounds(1, j) = val * 1.1;
    }

    double eta_min = std::numeric_limits<double>::infinity();
    double eta_max = 0.0;
    for (int i = 1; i <= d; ++i) {
        double ratio = x(i, N) / x(i, 1);
        if (ratio > 1.0) {
            double inv = 1.0 / (ratio - 1.0);
            eta_min = std::min(eta_min, inv);
            eta_max = std::max(eta_max, inv);
        }
    }
    eta_bounds << eta_min, eta_max;

    sigma2_bounds << 0.0, 0.01;
}

/**************************************************************************************************/
//                                           ERRORES                                              //
/**************************************************************************************************/

double computeRAE(const VectorXd& sample, const VectorXd& estimated) {
    double sum = 0.0;
    for (int i = 0; i < N; ++i) {
        sum += abs((sample(i) - estimated(i)) / sample(i));
    }
    return sum / N;
}

/**************************************************************************************************/
//                                            MAIN                                                //
/**************************************************************************************************/

int main(int argc, char* argv[]){

    bool savePlot = false;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];

        if (arg == "--plot") {
            savePlot = true;
        } else if (arg == "--help" || arg == "-h") {
            std::cout
                << "Usage: stochastic_diffusion [--plot]\n"
                << "\n"
                << "Runs a fast parameter-estimation demo on the supplied\n"
                << "stochastic diffusion trajectories.\n"
                << "\n"
                << "Options:\n"
                << "  --plot   Save figures/model_fit.png\n";
            return 0;
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            return 2;
        }
    }

    // Leer CSV
    vector<vector<double>> data = readCSV("data/TrayectoriasSimuladas_Multisigmoidal_Articulo.csv");
    
    // Convertir a Eigen matriz
    MatrixXd eigenMatrix = convertToMatrix(data);
    MatrixXd transposedMatrix = eigenMatrix.transpose();

    // Definir d y n
    d = transposedMatrix.rows()-1;
    int col = transposedMatrix.cols();
    n = createEvaluatedVector(d,col);
    N = n[1];
    
    // Calcular la matriz tiempo y de la variable aleatoria X
    t_values = transposedMatrix.row(0).transpose();
    t = obtainTime(transposedMatrix);
    x = obtainX(transposedMatrix);

    // Obtenemos V 
    v = obtainV();

    /// Coeficientes de beta
    p = 3;

    /// Vector de valores iniciales
    VectorXd estimationEtaBeta(p);
    estimationEtaBeta = estimateEtaBeta(p);

    double val_mu1 = obtainMu0();
    double val_sigma12 = obtainSigma0();

    VectorXd alpha(2);
    alpha << val_mu1, val_sigma12; 
    
    int tam_var = p+2;

    double val_eta = exp(estimationEtaBeta[0]);
    double val_beta1 = estimationEtaBeta[1];
    double val_beta2 = estimationEtaBeta[2];
    double val_beta3 = estimationEtaBeta[3];
    double val_sigma = 0.01;

    VectorXd epsilon(tam_var);
    epsilon << val_eta, val_beta1, val_beta2, val_beta3, val_sigma;  

    // Ejecutamos Newton-Raphson
    //VectorXd solution = newton_raphson(epsilon);
    //cout << solution << endl;

    // Calculamos los límites
    MatrixXd beta_bounds;
    Vector2d eta_bounds, sigma2_bounds;

    computeBoundsFromData(p, beta_bounds, eta_bounds, sigma2_bounds, estimationEtaBeta);

    MatrixXd bounds(2, tam_var); // 2 filas: min y max

    // Insertar límites: orden = [η, β1, ..., βp, σ²]
    bounds(0, 0) = eta_bounds(0); // mínimo η
    bounds(1, 0) = eta_bounds(1); // máximo η

    for (int j = 0; j < beta_bounds.cols(); ++j) {
        bounds(0, j + 1) = beta_bounds(0, j); // mínimos β
        bounds(1, j + 1) = beta_bounds(1, j); // máximos β
    }

    bounds(0, 1) = beta_bounds(0, 0); // mínimos β
    bounds(1, 1) = beta_bounds(1, 0); // máximos β

    bounds(0, 2) = beta_bounds(1, 1); // mínimos β
    bounds(1, 2) = beta_bounds(0, 1); // máximos β

    bounds(0, 3) = beta_bounds(0, 2); // mínimos β
    bounds(1, 3) = beta_bounds(1, 2); // máximos β

    bounds(0, tam_var - 1) = sigma2_bounds(0); // mínimo σ²
    bounds(1, tam_var - 1) = sigma2_bounds(1); // máximo σ²

    // Valores reales
    double real_eta = exp(-1);
    double real_beta1 = 0.1;
    double real_beta2 = -0.009;
    double real_beta3 = 0.0002;
    double real_sigma = 0.01*0.01;

    VectorXd epsilon_real(tam_var);
    epsilon_real << real_eta, real_beta1, real_beta2, real_beta3, real_sigma; 


    // -------------------------------------------------------------------------
    // Fast portfolio demo
    // -------------------------------------------------------------------------

    const VectorXd sampleMean = removeFirst(arithmeticMean());

    const VectorXd estimatedMean =
        media(val_eta, val_beta1, val_beta2, val_beta3);

    const VectorXd referenceMean =
        media(real_eta, real_beta1, real_beta2, real_beta3);

    auto rmse = [](const VectorXd& observed, const VectorXd& predicted) {
        if (observed.size() != predicted.size() || observed.size() == 0) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        return std::sqrt(
            (observed - predicted).squaredNorm()
            / static_cast<double>(observed.size())
        );
    };

    const double estimateRmse = rmse(sampleMean, estimatedMean);
    const double referenceRmse = rmse(sampleMean, referenceMean);

    std::cout << "\n";
    std::cout << "StochasticDiffusion\n";
    std::cout << "============================================================\n";
    std::cout << "Parameter estimation for a stochastic multi-sigmoidal\n";
    std::cout << "diffusion model from simulated trajectories.\n\n";

    std::cout << "Dataset\n";
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Trajectories      : " << d << "\n";
    std::cout << "Time observations : " << N << "\n\n";

    std::cout << std::fixed << std::setprecision(6);

    std::cout << "Regression-based parameter estimate\n";
    std::cout << "------------------------------------------------------------\n";
    std::cout
        << std::left  << std::setw(12) << "Parameter"
        << std::right << std::setw(14) << "Estimate"
        << std::setw(14) << "Reference"
        << std::setw(14) << "Abs. error"
        << "\n";

    std::cout << std::string(54, '-') << "\n";

    auto printParameter = [](const std::string& name,
                             double estimate,
                             double reference) {
        std::cout
            << std::left  << std::setw(12) << name
            << std::right << std::setw(14) << estimate
            << std::setw(14) << reference
            << std::setw(14) << std::abs(estimate - reference)
            << "\n";
    };

    printParameter("eta",   val_eta,   real_eta);
    printParameter("beta1", val_beta1, real_beta1);
    printParameter("beta2", val_beta2, real_beta2);
    printParameter("beta3", val_beta3, real_beta3);

    const double maxRelativeParameterError = std::max({
        std::abs((val_eta   - real_eta)   / real_eta),
        std::abs((val_beta1 - real_beta1) / real_beta1),
        std::abs((val_beta2 - real_beta2) / real_beta2),
        std::abs((val_beta3 - real_beta3) / real_beta3)
    }) * 100.0;

    std::cout
        << "\nMaximum relative shape-parameter error: "
        << std::setprecision(2)
        << maxRelativeParameterError
        << "%\n";

    std::cout << std::setprecision(6);

    std::cout << "\nMean-trajectory fit\n";
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Estimated model RMSE : " << estimateRmse << "\n";
    std::cout << "Reference model RMSE  : " << referenceRmse << "\n";

    if (savePlot) {
        std::filesystem::create_directories("figures");

        const bool modelPlotOk = saveDemoFitPlot(
            t_values,
            sampleMean,
            estimatedMean,
            referenceMean,
            "figures/model_fit.png"
        );

        const bool trajectoryPlotOk = saveTrajectoryFitPlot(
            x,
            t_values,
            sampleMean,
            estimatedMean,
            referenceMean,
            "figures/trajectory_fit.png"
        );

        VectorXd estimatedParameters(4);
        estimatedParameters
            << val_eta, val_beta1, val_beta2, val_beta3;

        VectorXd referenceParameters(4);
        referenceParameters
            << real_eta, real_beta1, real_beta2, real_beta3;

        const bool parameterPlotOk = saveParameterRecoveryPlot(
            estimatedParameters,
            referenceParameters,
            "figures/parameter_recovery.png"
        );

        if (!modelPlotOk || !trajectoryPlotOk || !parameterPlotOk) {
            std::cerr
                << "\nCould not generate the figures. "
                << "Check that gnuplot is installed.\n";
            return 3;
        }

        std::cout << "\nFigures saved:\n";
        std::cout << "  figures/model_fit.png\n";
        std::cout << "  figures/trajectory_fit.png\n";
        std::cout << "  figures/parameter_recovery.png\n";
    }

    std::cout << "\nDemo completed successfully.\n\n";

    // Historical curve estimates were used only for thesis visualisation.
    // The plotting block is intentionally omitted from the portfolio executable
    // because its precomputed reference vectors are not part of this repository.



    // Ejecución del algoritmo
    
    /*
    cout<<"NEWTON-RAPHSON"<<endl;
    
    ofstream outfile("resultados_newtonRaphson.txt");
    if (!outfile.is_open()) {
        cerr << "Error al abrir el archivo de salida." << endl;
        return 1;
    }

    for (int i = 0; i < 30; ++i) {
        Result result = newton_raphsonR(epsilon);
        outfile << "Ejecución " << i + 1 << ":\n"
                << "Solución: " << result.solucion.transpose() << "\n"
                << "Fitness: " << result.fitness << "\n"
                << "Tiempo: " << result.time << " segundos\n\n"
                << std::flush;
    }
    */
    return 0;
}
