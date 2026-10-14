// visualization.hpp
#ifndef VISUALIZATION_HPP
#define VISUALIZATION_HPP

#include <iostream>
#include <Eigen/Dense>
#include <string>
#include <vector>
#include "model_data.hpp"
#include "model_functions.hpp"

// Función que grafica los puntos (x, y) con un título específico.
void graph(const Eigen::VectorXd& x, const Eigen::VectorXd& y, const std::string& titulo);

// Función que grafica una función definida por una cadena de texto entre dos puntos.
void graphFunction(const std::string& function, const double x_min, const double x_max, const std::string& titulo);

// Función que grafica una función definida por un puntero a función entre dos puntos.
void graphFunction(const std::function<double(double)>& func, double x_min, double x_max, const std::string& titulo);

// Función que grafica puntos (x, y) y una función definida por una cadena de texto entre dos puntos.
void graphPointsAndFunction(const Eigen::VectorXd& x, const Eigen::VectorXd& y, const std::function<double(double)>& func, double x_min, double x_max, const std::string& titulo);

// Función que grafica múltiples trayectorias.
void graphMultiplePaths(const Eigen::VectorXd& x, const std::vector<Eigen::VectorXd>& allPaths, const std::string& titulo);

// Función que grafica todas las trayectorias.
void graphAllPaths();

// Función que grafica la media de los datos
Eigen::VectorXd media(double eta, double beta1, double beta2, double beta3);

// Función que grafica todos los datos en gris junto a su media.
void graphMultiplePathsWithMean(const Eigen::MatrixXd& x, const Eigen::VectorXd& t_values, double val_eta, double val_beta1, double val_beta2, double val_beta3, const std::string& titulo);

// Función que gradica tres curvas: la teórica, la estimada y la de la muestra.
void graphThreeCurves(const Eigen::VectorXd& x, const Eigen::VectorXd& m_teorica, const Eigen::VectorXd& m_estimada, const Eigen::VectorXd& m_sample, const std::string& titulo);

// Función que grafica múltiples curvas, cada una con su propio nombre y datos.
void graphMultipleCurves(const Eigen::VectorXd& x, const std::vector<std::pair<std::string, Eigen::VectorXd>>& curves,const std::string& titulo);

// Función que grafica múltiples curvas, cada una con su propio nombre y datos, que incluye NR.
void graphMultipleCurvesNR(const Eigen::VectorXd& x, const std::vector<std::pair<std::string, Eigen::VectorXd>>& curves,const std::string& titulo);


// Save the fast portfolio demo as a PNG figure.
bool saveDemoFitPlot(
    const Eigen::VectorXd& time,
    const Eigen::VectorXd& sampleMean,
    const Eigen::VectorXd& estimatedMean,
    const Eigen::VectorXd& referenceMean,
    const std::string& outputPath);


// Save a portfolio figure with representative stochastic trajectories
// together with the observed mean and fitted diffusion model.
bool saveTrajectoryFitPlot(
    const Eigen::MatrixXd& trajectories,
    const Eigen::VectorXd& time,
    const Eigen::VectorXd& sampleMean,
    const Eigen::VectorXd& estimatedMean,
    const Eigen::VectorXd& referenceMean,
    const std::string& outputPath);


// Save relative recovery errors for the four shape parameters.
bool saveParameterRecoveryPlot(
    const Eigen::VectorXd& estimated,
    const Eigen::VectorXd& reference,
    const std::string& outputPath);

#endif
