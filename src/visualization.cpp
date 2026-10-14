#include "visualization.hpp"

using namespace std;
using namespace Eigen;

void graph(const VectorXd& x, const VectorXd& y, const std::string& titulo) {
    if (x.size() != y.size()) {
        std::cerr << "Error: vectores x e y deben tener el mismo tamaño.\n";
        return;
    }

    FILE* pipe = popen("gnuplot -persistent", "w");
    if (!pipe) {
        std::cerr << "Error al abrir gnuplot.\n";
        return;
    }

    fprintf(pipe, "set title '%s'\n", titulo.c_str());
    fprintf(pipe, "set xlabel 'x'\n");
    fprintf(pipe, "set ylabel 'y'\n");
    fprintf(pipe, "set yrange [0:20]\n");  // Límite del eje y
    fprintf(pipe, "plot '-' with points pointtype 7 pointsize 0.5 title '%s'\n", titulo.c_str());

    for (int i = 0; i < x.size(); ++i) {
        fprintf(pipe, "%f %f\n", x(i), y(i));
    }

    fprintf(pipe, "e\n");
    pclose(pipe);
}

void graphFunction(const std::string& function, const double x_min, const double x_max, const std::string& titulo) {
    // Open a pipe to gnuplot
    FILE* pipe = popen("gnuplot -persistent", "w");
    if (!pipe) {
        std::cerr << "Error opening gnuplot.\n";
        return;
    }

    // Set up plot labels and title
    fprintf(pipe, "set title '%s'\n", titulo.c_str());
    fprintf(pipe, "set xlabel 'x'\n");
    fprintf(pipe, "set ylabel 'y'\n");

    // Plot the function with a range for x values
    fprintf(pipe, "plot [%f:%f] %s title '%s' with lines\n", x_min, x_max, function.c_str(), titulo.c_str());

    // Close the pipe
    pclose(pipe);
}

void graphFunction(const std::function<double(double)>& func, double x_min, double x_max, const std::string& titulo) {
    FILE* pipe = popen("gnuplot -persistent", "w");
    if (!pipe) {
        std::cerr << "Error opening gnuplot.\n";
        return;
    }

    // Configuración del gráfico
    fprintf(pipe, "set title '%s'\n", titulo.c_str());
    fprintf(pipe, "set xlabel 'x'\n");
    fprintf(pipe, "set ylabel 'y'\n");
    fprintf(pipe, "set grid\n");

    // Solo graficar la función
    fprintf(pipe, "plot '-' with lines lw 2 title 'Funcion'\n");

    // Evaluar y enviar puntos de la función
    for (double x_val = x_min; x_val <= x_max; x_val += 0.05) {
        double y_val = func(x_val);
        fprintf(pipe, "%f %f\n", x_val, y_val);
    }
    fprintf(pipe, "e\n");

    pclose(pipe);
}

void graphPointsAndFunction(const VectorXd& x, const VectorXd& y, const std::function<double(double)>& func, const double x_min, const double x_max, const std::string& titulo) {
    FILE* pipe = popen("gnuplot -persistent", "w");
    if (!pipe) {
        std::cerr << "Error opening gnuplot.\n";
        return;
    }

    // Set up plot labels and title
    fprintf(pipe, "set title '%s'\n", titulo.c_str());
    fprintf(pipe, "set xlabel 't_j'\n");
    fprintf(pipe, "set ylabel 'σ_j^2-σ_0^2'\n");

    fprintf(pipe, "set key top left\n");

    // Plot points and the function in one command
    fprintf(pipe, "plot '-' with points pointtype 7 pointsize 0.1 title 'Puntos', ");
    fprintf(pipe, "'-' with lines title 'Regresión polinomial' \n");

    // Data points
    for (int i = 0; i < x.size(); ++i) {
        fprintf(pipe, "%f %f\n", x(i), y(i));
    }
    fprintf(pipe, "e\n"); // End the first plot section

    // Function plot: Plot the function y = func(x) over a range of x
    for (double x_val = x_min; x_val <= x_max; x_val += 0.1) {
        double y_val = func(x_val);  // Evaluate the function
        fprintf(pipe, "%f %f\n", x_val, y_val);
    }
    fprintf(pipe, "e\n"); // End the second plot section

    // Close the pipe
    pclose(pipe);
}

void graphMultiplePaths(const Eigen::VectorXd& x, const std::vector<Eigen::VectorXd>& allPaths, const std::string& titulo) {
    FILE* pipe = popen("gnuplot -persistent", "w");
    if (!pipe) {
        std::cerr << "Error al abrir gnuplot.\n";
        return;
    }

    fprintf(pipe, "set title '%s'\n", titulo.c_str());
    fprintf(pipe, "set xlabel 'Tiempo'\n");
    fprintf(pipe, "set ylabel 'Valor'\n");

    // Desactivar leyenda
    fprintf(pipe, "unset key\n");

    // Comando plot para múltiples trayectorias sin leyenda
    fprintf(pipe, "plot ");
    for (size_t i = 0; i < allPaths.size(); ++i) {
        if (i > 0) fprintf(pipe, ", ");
        fprintf(pipe, "'-' with linespoints pointtype 7 pointsize 0.25");  // Sin título
    }
    fprintf(pipe, "\n");

    // Datos para cada trayectoria
    for (const auto& path : allPaths) {
        for (int j = 0; j < path.size(); ++j) {
            fprintf(pipe, "%f %f\n", x(j), path(j));
        }
        fprintf(pipe, "e\n");
    }

    pclose(pipe);
}

void graphAllPaths() {
    vector<VectorXd> allPaths;  

    for (int i = 1; i <= d; ++i) {
        VectorXd fullPath = v.row(i).transpose();
        VectorXd path = fullPath.segment(1, fullPath.size() - 1);  // Exclude the first element

        allPaths.push_back(path);  
    }

    graphMultiplePaths(t_values, allPaths, "All Paths");
}

VectorXd media(double eta,double beta1,double beta2, double beta3) {
    // Inicializamos el vector de salida
    int tam = static_cast<int>(t_values.size());
    VectorXd m(tam);  // Asumiendo que queremos almacenar 502 valores
    double theta[] = {0.0, beta1, beta2, beta3};
    
    // Asumiendo que t es un objeto que tiene la dimensión adecuada para pasar a evaluate_polynomial
    // Aquí debes asegurarte de que t esté correctamente definido

    for (int i = 0; i < tam; ++i) {
        double Q0 = evaluate_polynomial<double>(theta, t_values(0,0));  // Ajusta el valor de t según sea necesario
        double Q = evaluate_polynomial<double>(theta, t_values(i,0));     // Ajusta el valor de t según sea necesario

        // Almacenar el valor en la posición i-1 de m
        m(i) = 5 * (eta + exp(-Q0)) / (eta + exp(-Q));
    }

    return m;
}

void graphMultiplePathsWithMean(const Eigen::MatrixXd& x,
                                 const Eigen::VectorXd& t_values,
                                 double val_eta,
                                 double val_beta1,
                                 double val_beta2,
                                 double val_beta3,
                                 const std::string& titulo) {
    FILE* pipe = popen("gnuplot -persistent", "w");
    if (!pipe) {
        std::cerr << "Error al abrir gnuplot.\n";
        return;
    }

    int d = x.rows() - 1;  // Asumes que la fila 0 es especial y se ignora

    // Construir las trayectorias
    std::vector<Eigen::VectorXd> allPaths;
    for (int i = 1; i <= d; ++i) {
        Eigen::VectorXd fullPath = x.row(i).transpose();
        Eigen::VectorXd path = fullPath.segment(1, fullPath.size() - 1);  // Excluir el primer valor
        allPaths.push_back(path);
    }

    // Calcular la media
    Eigen::VectorXd med = media(val_eta, val_beta1, val_beta2, val_beta3);

    // Configurar Gnuplot
    fprintf(pipe, "set title '%s'\n", titulo.c_str());
    fprintf(pipe, "set xlabel 'Tiempo'\n");
    fprintf(pipe, "set ylabel 'Valor'\n");
    fprintf(pipe, "unset key\n");

    // Comando plot para todas las trayectorias + una media
    fprintf(pipe, "plot ");
    for (size_t i = 0; i < allPaths.size(); ++i) {
        if (i > 0) fprintf(pipe, ", ");
        fprintf(pipe, "'-' with lines lc rgb '#aaaaaa' lw 1");  // gris suave
    }
    fprintf(pipe, ", '-' with lines lc rgb '#d63044' lw 2.5\n"); // media destacada en rojo

    // Datos de trayectorias
    for (const auto& path : allPaths) {
        for (int j = 0; j < path.size(); ++j) {
            fprintf(pipe, "%f %f\n", t_values(j), path(j));
        }
        fprintf(pipe, "e\n");
    }

    // Datos de la trayectoria media
    for (int j = 0; j < med.size(); ++j) {
        fprintf(pipe, "%f %f\n", t_values(j), med(j));
    }
    fprintf(pipe, "e\n");

    pclose(pipe);
}

void graphThreeCurves(
    const VectorXd& x,
    const VectorXd& m_teorica,
    const VectorXd& m_estimada,
    const VectorXd& m_sample,
    const std::string& titulo)
{
    if (x.size() != m_teorica.size() || x.size() != m_estimada.size() || x.size() != m_sample.size()) {
        std::cerr << "Error: todos los vectores deben tener el mismo tamaño.\n";
        return;
    }

    FILE* pipe = popen("gnuplot -persistent", "w");
    if (!pipe) {
        std::cerr << "Error al abrir gnuplot.\n";
        return;
    }

    fprintf(pipe, "set title '%s'\n", titulo.c_str());
    fprintf(pipe, "set xlabel 'Tiempo'\n");
    fprintf(pipe, "set ylabel 'Media'\n");
    fprintf(pipe, "set grid\n");
    fprintf(pipe, "set key left top \n"); 

    // Ajusta el rango Y si lo deseas, o comenta la línea siguiente
    // fprintf(pipe, "set yrange [0:1]\n"); 

    // Definimos las tres curvas con estilos diferentes
    fprintf(pipe, "plot '-' with lines lw 2 linecolor rgb 'black' title 'Teórica', "
                  "'-' with lines lw 2 dt 2 linecolor rgb 'red' title 'Estimada', "
                  "'-' with lines lw 2 dt 3 linecolor rgb 'blue' title 'Muestral'\n");

    // Datos de la curva teórica
    for (int i = 0; i < x.size(); ++i) {
        fprintf(pipe, "%f %f\n", x(i), m_teorica(i));
    }
    fprintf(pipe, "e\n");

    // Datos de la curva estimada
    for (int i = 0; i < x.size(); ++i) {
        fprintf(pipe, "%f %f\n", x(i), m_estimada(i));
    }
    fprintf(pipe, "e\n");

    // Datos de la curva muestral
    for (int i = 0; i < x.size(); ++i) {
        fprintf(pipe, "%f %f\n", x(i), m_sample(i));
    }
    fprintf(pipe, "e\n");

    pclose(pipe);
}

void graphMultipleCurves(const VectorXd& x,
                         const std::vector<std::pair<std::string, VectorXd>>& curvas,
                         const std::string& titulo) {
    FILE* pipe = popen("gnuplot -persistent", "w");
    if (!pipe) {
        std::cerr << "Error al abrir gnuplot.\n";
        return;
    }

    std::vector<std::string> colores = {
        "black", "gray", "red", "blue", "green", "orange", "violet", "cyan", "magenta"
    };

    fprintf(pipe, "set title '%s'\n", titulo.c_str());
    fprintf(pipe, "set xlabel 'Tiempo'\n");
    fprintf(pipe, "set ylabel 'Media'\n");
    fprintf(pipe, "set key left top\n");
    fprintf(pipe, "set grid\n");

    // Comando plot: línea sólida para la teórica, luego discontinuas para el resto
    fprintf(pipe, "plot ");
    for (size_t i = 0; i < curvas.size(); ++i) {
        if (i > 0) fprintf(pipe, ", ");

        std::string estilo = (i == 0) 
            ? "with lines lw 3"  // Teórica: línea sólida
            : "with lines dt 2 lw 2";  // Resto: línea discontinua

        std::string color = colores[i % colores.size()];
        fprintf(pipe, "'-' %s linecolor rgb '%s' title '%s'",
                estilo.c_str(), color.c_str(), curvas[i].first.c_str());
    }
    fprintf(pipe, "\n");

    // Datos
    for (const auto& [nombre, valores] : curvas) {
        for (int i = 0; i < x.size(); ++i)
            fprintf(pipe, "%f %f\n", x(i), valores(i));
        fprintf(pipe, "e\n");
    }

    pclose(pipe);
}

void graphMultipleCurvesNR(const VectorXd& x,
                         const std::vector<std::pair<std::string, VectorXd>>& curvas,
                         const std::string& titulo) {
    FILE* pipe = popen("gnuplot -persistent", "w");
    if (!pipe) {
        std::cerr << "Error al abrir gnuplot.\n";
        return;
    }

    std::vector<std::string> colores = {
        "black",         // Teórica
        "grey",          // NR: gris oscuro
        "brown",         // BL
        "red",           // SA
        "blue",          // ILS
        "green",         // EA
        "orange",        // EA+E
        "violet",        // MA+E
        "cyan",          // DE
        "magenta"        // por si acaso
    };

    fprintf(pipe, "set title '%s'\n", titulo.c_str());
    fprintf(pipe, "set xlabel 'Tiempo'\n");
    fprintf(pipe, "set ylabel 'Media'\n");
    fprintf(pipe, "set key left top\n");
    fprintf(pipe, "set grid\n");

    // Comando plot
    fprintf(pipe, "plot ");
    for (size_t i = 0; i < curvas.size(); ++i) {
        if (i > 0) fprintf(pipe, ", ");
        std::string estilo = (i == 0)
            ? "with lines lw 3"
            : "with lines dt 2 lw 2";
        std::string color = colores[i % colores.size()];
        fprintf(pipe, "'-' %s linecolor rgb '%s' title '%s'",
                estilo.c_str(), color.c_str(), curvas[i].first.c_str());
    }
    fprintf(pipe, "\n");

    // Datos
    for (const auto& [nombre, valores] : curvas) {
        for (int i = 0; i < x.size(); ++i)
            fprintf(pipe, "%f %f\n", x(i), valores(i));
        fprintf(pipe, "e\n");
    }

    pclose(pipe);
}


bool saveDemoFitPlot(
    const Eigen::VectorXd& time,
    const Eigen::VectorXd& sampleMean,
    const Eigen::VectorXd& estimatedMean,
    const Eigen::VectorXd& referenceMean,
    const std::string& outputPath)
{
    if (time.size() != sampleMean.size() ||
        time.size() != estimatedMean.size() ||
        time.size() != referenceMean.size()) {
        std::cerr << "Cannot generate plot: vector sizes differ.\n";
        return false;
    }

    FILE* pipe = popen("gnuplot", "w");

    if (!pipe) {
        std::cerr << "Could not start gnuplot.\n";
        return false;
    }

    fprintf(
        pipe,
        "set terminal pngcairo size 1400,800 enhanced font 'Arial,13'\n"
    );

    fprintf(pipe, "set output '%s'\n", outputPath.c_str());

    fprintf(
        pipe,
        "set title 'Stochastic diffusion model fit' font ',18'\n"
    );

    fprintf(pipe, "set xlabel 'Time'\n");
    fprintf(pipe, "set ylabel 'Mean trajectory'\n");
    fprintf(pipe, "set grid back\n");
    fprintf(pipe, "set border linewidth 1.2\n");
    fprintf(pipe, "set tics out\n");
    fprintf(pipe, "set key left top box opaque\n");

    fprintf(
        pipe,
        "plot "
        "'-' with points pt 7 ps 0.55 lc rgb '#4C78A8' "
        "title 'Observed sample mean', "
        "'-' with lines lw 3 lc rgb '#E45756' "
        "title 'Estimated model', "
        "'-' with lines lw 2 dt 2 lc rgb '#333333' "
        "title 'Reference model'\n"
    );

    for (int i = 0; i < time.size(); ++i) {
        fprintf(pipe, "%.10f %.10f\n", time(i), sampleMean(i));
    }
    fprintf(pipe, "e\n");

    for (int i = 0; i < time.size(); ++i) {
        fprintf(pipe, "%.10f %.10f\n", time(i), estimatedMean(i));
    }
    fprintf(pipe, "e\n");

    for (int i = 0; i < time.size(); ++i) {
        fprintf(pipe, "%.10f %.10f\n", time(i), referenceMean(i));
    }
    fprintf(pipe, "e\n");

    fflush(pipe);
    return pclose(pipe) == 0;
}


bool saveTrajectoryFitPlot(
    const Eigen::MatrixXd& trajectories,
    const Eigen::VectorXd& time,
    const Eigen::VectorXd& sampleMean,
    const Eigen::VectorXd& estimatedMean,
    const Eigen::VectorXd& referenceMean,
    const std::string& outputPath)
{
    if (time.size() != sampleMean.size() ||
        time.size() != estimatedMean.size() ||
        time.size() != referenceMean.size()) {
        std::cerr << "Cannot generate trajectory plot: vector sizes differ.\n";
        return false;
    }

    if (trajectories.rows() < 2 ||
        trajectories.cols() < time.size() + 1) {
        std::cerr << "Cannot generate trajectory plot: invalid trajectory matrix.\n";
        return false;
    }

    FILE* pipe = popen("gnuplot", "w");

    if (!pipe) {
        std::cerr << "Could not start gnuplot.\n";
        return false;
    }

    fprintf(
        pipe,
        "set terminal pngcairo size 1400,800 enhanced font 'Arial,13'\n"
    );

    fprintf(pipe, "set output '%s'\n", outputPath.c_str());

    fprintf(
        pipe,
        "set title 'Stochastic trajectories and fitted diffusion model' font ',18'\n"
    );

    fprintf(pipe, "set xlabel 'Time'\n");
    fprintf(pipe, "set ylabel 'Process value'\n");
    fprintf(pipe, "set grid back lc rgb '#dddddd'\n");
    fprintf(pipe, "set border linewidth 1.2\n");
    fprintf(pipe, "set tics out\n");
    fprintf(pipe, "set key left top box opaque\n");

    const int available = static_cast<int>(trajectories.rows()) - 1;
    const int maxDisplayed = 20;

    int step = available / maxDisplayed;
    if (step < 1) {
        step = 1;
    }

    int displayed = 0;
    for (int i = 1; i <= available && displayed < maxDisplayed; i += step) {
        ++displayed;
    }

    fprintf(pipe, "plot ");

    for (int k = 0; k < displayed; ++k) {
        if (k > 0) {
            fprintf(pipe, ", ");
        }

        if (k == 0) {
            fprintf(
                pipe,
                "'-' with lines lw 0.7 lc rgb '#c7c7c7' "
                "title 'Individual trajectories'"
            );
        } else {
            fprintf(
                pipe,
                "'-' with lines lw 0.7 lc rgb '#c7c7c7' notitle"
            );
        }
    }

    fprintf(
        pipe,
        ", '-' with lines lw 2.5 lc rgb '#4C78A8' "
        "title 'Observed sample mean'"
    );

    fprintf(
        pipe,
        ", '-' with lines lw 3 lc rgb '#E45756' "
        "title 'Estimated model'"
    );

    fprintf(
        pipe,
        ", '-' with lines lw 2 dt 2 lc rgb '#333333' "
        "title 'Reference model'\n"
    );

    int written = 0;
    for (int i = 1; i <= available && written < displayed; i += step) {
        for (int j = 0; j < time.size(); ++j) {
            fprintf(
                pipe,
                "%.10f %.10f\n",
                time(j),
                trajectories(i, j + 1)
            );
        }
        fprintf(pipe, "e\n");
        ++written;
    }

    for (int i = 0; i < time.size(); ++i) {
        fprintf(pipe, "%.10f %.10f\n", time(i), sampleMean(i));
    }
    fprintf(pipe, "e\n");

    for (int i = 0; i < time.size(); ++i) {
        fprintf(pipe, "%.10f %.10f\n", time(i), estimatedMean(i));
    }
    fprintf(pipe, "e\n");

    for (int i = 0; i < time.size(); ++i) {
        fprintf(pipe, "%.10f %.10f\n", time(i), referenceMean(i));
    }
    fprintf(pipe, "e\n");

    fflush(pipe);
    return pclose(pipe) == 0;
}


bool saveParameterRecoveryPlot(
    const Eigen::VectorXd& estimated,
    const Eigen::VectorXd& reference,
    const std::string& outputPath)
{
    if (estimated.size() != 4 || reference.size() != 4) {
        std::cerr << "Parameter recovery plot expects four shape parameters.\n";
        return false;
    }

    double errors[4];
    double maxError = 0.0;

    for (int i = 0; i < 4; ++i) {
        if (reference(i) == 0.0) {
            std::cerr << "Cannot compute relative error for zero reference.\n";
            return false;
        }

        errors[i] =
            std::abs((estimated(i) - reference(i)) / reference(i)) * 100.0;

        if (errors[i] > maxError) {
            maxError = errors[i];
        }
    }

    FILE* pipe = popen("gnuplot", "w");

    if (!pipe) {
        std::cerr << "Could not start gnuplot.\n";
        return false;
    }

    const double yMax = std::max(5.0, maxError * 1.35);

    fprintf(
        pipe,
        "set terminal pngcairo size 1100,700 enhanced font 'Arial,13'\n"
    );

    fprintf(pipe, "set output '%s'\n", outputPath.c_str());

    fprintf(
        pipe,
        "set title 'Recovery of model parameters from stochastic data' "
        "font ',18'\n"
    );

    fprintf(pipe, "set ylabel 'Relative error (%%)'\n");
    fprintf(pipe, "set xrange [-0.6:3.6]\n");
    fprintf(pipe, "set yrange [0:%.6f]\n", yMax);

    fprintf(
        pipe,
        "set xtics ('eta' 0, 'beta1' 1, 'beta2' 2, 'beta3' 3)\n"
    );

    fprintf(pipe, "set grid ytics back lc rgb '#dddddd'\n");
    fprintf(pipe, "set border linewidth 1.2\n");
    fprintf(pipe, "set tics out\n");
    fprintf(pipe, "set boxwidth 0.58\n");
    fprintf(pipe, "set style fill solid 0.85 border -1\n");
    fprintf(pipe, "unset key\n");

    fprintf(
        pipe,
        "plot "
        "'-' using 1:2 with boxes lc rgb '#4C78A8', "
        "'-' using 1:2:3 with labels center offset 0,0.8\n"
    );

    for (int i = 0; i < 4; ++i) {
        fprintf(pipe, "%d %.10f\n", i, errors[i]);
    }
    fprintf(pipe, "e\n");

    for (int i = 0; i < 4; ++i) {
        fprintf(
            pipe,
            "%d %.10f \"%.2f%%\"\n",
            i,
            errors[i],
            errors[i]
        );
    }
    fprintf(pipe, "e\n");

    fflush(pipe);
    return pclose(pipe) == 0;
}
