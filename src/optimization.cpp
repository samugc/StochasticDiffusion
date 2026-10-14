#include "optimization.hpp"

using namespace std;
using namespace Eigen;

double logLikelihood(const VectorXd&, const VectorXd& epsilon) {

    const double sigma2 = epsilon(4);

    double lGorro = -(nTotal<double>()/2.0)*log(sigma2) - (Z1<double>() + phi(epsilon) - 2*varGamma(epsilon)) / (2*sigma2);
    
    if (!std::isfinite(lGorro) || std::isnan(lGorro) || std::abs(lGorro) > 1e10) {
        return -1e10;  // Penaliza valores inválidos
    }

    return lGorro;
}

bool accept(double currentVal, double candidateVal, double T) {
    if (candidateVal > currentVal) return true;  // Mejor: mayor valor (maximizamos)
    double p = std::exp((candidateVal - currentVal) / T); // Probabilidad para peor solución
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dis(0.0, 1.0);
    return dis(gen) < p;
}

VectorXd generateNeighbor(const VectorXd& current, const MatrixXd& bounds, double stepScale) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    VectorXd neighbor = current;

    for (int i = 0; i < current.size(); ++i) {
        double range = bounds(1, i) - bounds(0, i);
        std::uniform_real_distribution<> dist(-stepScale * range, stepScale * range);
        double step = dist(gen);
        neighbor(i) = std::min(bounds(1, i), std::max(bounds(0, i), current(i) + step));
    }

    return neighbor;
}

// Generador aleatorio
std::mt19937 genRandom(std::random_device{}());

double randDouble(double a, double b) {
    std::uniform_real_distribution<> dis(a, b);
    return dis(genRandom);
}

Result localSearch(
    const VectorXd& alpha,
    const VectorXd& initial,
    const MatrixXd& bounds,
    int maxIter,
    double baseStep,
    double minStep)
{
    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    VectorXd current = initial;
    double currentVal = logLikelihood(alpha, current);
    int dim = current.size();
    std::mt19937 gen(std::random_device{}());

    for (int iter = 0; iter < maxIter; ++iter) {
        bool improved = false;
        std::vector<int> indices(dim);
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), gen);

        for (int i : indices) {
            double range = bounds(1, i) - bounds(0, i);
            double step = baseStep * range;

            VectorXd temp = current;
            double bestVal = currentVal;
            double bestCoord = current(i);

            for (double delta : {-step, step}) {
                double candidate = current(i) + delta;
                candidate = std::min(bounds(1, i), std::max(bounds(0, i), candidate));
                temp(i) = candidate;

                double val = logLikelihood(alpha, temp);
                if (val > bestVal) {
                    bestVal = val;
                    bestCoord = candidate;
                    improved = true;
                }
            }

            current(i) = bestCoord;
            currentVal = bestVal;
        }

        std::cout << "Iteración " << iter << ", LogL = " << currentVal << std::endl;

        if (!improved && baseStep > minStep) {
            baseStep *= 0.5;
            std::cout << "↓ Refinando stepSize → " << baseStep << std::endl;
        }

        if (!improved && baseStep <= minStep) {
            std::cout << "Sin mejoras. Finaliza.\n";
            break;
        }
    }

    auto end = high_resolution_clock::now();
    double elapsed = duration<double>(end - start).count();

    std::cout << "Tiempo ejecución coordinate descent: " << elapsed << " s\n";

    return Result{current, currentVal, elapsed};
}

Result simulatedAnnealing(
    const VectorXd& alpha,
    const VectorXd& epsilon_init,
    const MatrixXd& bounds,
    double T0,
    double gamma,
    int L,
    int maxIter,
    double T_min,
    double stepSize,
    int restarts)
{
    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    VectorXd bestGlobal = epsilon_init;
    double bestGlobalVal = logLikelihood(alpha, bestGlobal);

    for (int restart = 0; restart < restarts; ++restart) {
        VectorXd current = bestGlobal;
        double currentVal = bestGlobalVal;
        VectorXd bestLocal = current;
        double bestLocalVal = currentVal;

        double T = T0;
        int iter = 0, unchanged = 0;

        while (iter < maxIter && T > T_min && unchanged < L) {
            bool changed = false;

            for (int l = 0; l < L; ++l) {
                VectorXd candidate = generateNeighbor(current, bounds, stepSize);
                double candidateVal = logLikelihood(alpha, candidate);

                if (!std::isfinite(candidateVal)) continue;

                if (candidateVal > currentVal || randDouble(0.0, 1.0) < exp((candidateVal - currentVal) / T)) {
                    current = candidate;
                    currentVal = candidateVal;
                    changed = true;

                    if (currentVal > bestLocalVal) {
                        bestLocal = current;
                        bestLocalVal = currentVal;
                        if (bestLocalVal > bestGlobalVal) {
                            bestGlobal = bestLocal;
                            bestGlobalVal = bestLocalVal;
                        }
                    }
                }
            }

            if (!changed) ++unchanged;
            else unchanged = 0;

            if (iter % 10 == 0) {
                std::cout << "Iteración " << iter << ", T = " << T << ", LogL actual = " << currentVal
                          << ", LogL mejor = " << bestGlobalVal << ", ∆ = " << bestGlobalVal - currentVal << std::endl;
            }

            T *= gamma;
            ++iter;
        }

        std::cout << "\n[Restart " << restart+1 << "] Mejor LogL local = " << bestLocalVal << "\n";
        stepSize *= 0.5;
    }

    auto end = high_resolution_clock::now();
    double elapsed = duration<double>(end - start).count();

    std::cout << "\nTiempo de ejecución SA: " << elapsed << " segundos.\n";
    std::cout << "Mejor solución encontrada SA:\n" << bestGlobal.transpose() << std::endl;

    Result resultado;
    resultado.solucion = bestGlobal;
    resultado.fitness = bestGlobalVal;
    resultado.time = elapsed;

    return resultado;
}

Result iteratedLocalSearch(
    const VectorXd& alpha,
    const VectorXd& initial,
    const MatrixXd& bounds,
    int ilsMaxIter,
    double perturbationSize,
    int lsMaxIter,
    double lsStepSize)
{
    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    Result resultLS = localSearch(alpha, initial, bounds, lsMaxIter, lsStepSize);
    VectorXd best = resultLS.solucion;
    double bestVal = resultLS.fitness;

    for (int iter = 0; iter < ilsMaxIter; ++iter) {
        VectorXd perturbed = best;
        for (int i = 0; i < perturbed.size(); ++i) {
            double noise = perturbationSize * (2.0 * randDouble(-1.0, 1.0));
            perturbed[i] += noise;
            perturbed[i] = std::max(bounds(0, i), std::min(bounds(1, i), perturbed[i]));
        }

        Result resultCand = localSearch(alpha, perturbed, bounds, lsMaxIter, lsStepSize);
        VectorXd candidate = resultCand.solucion;
        double candidateVal = resultCand.fitness;

        if (candidateVal > bestVal) {
            best = candidate;
            bestVal = candidateVal;
        }
    }

    auto end = high_resolution_clock::now();
    double elapsedSec = duration<double>(end - start).count();

    std::cout << "Tiempo de ejecución ILS: " << elapsedSec << " segundos.\n";
    std::cout << "Mejor solución encontrada ILS:\n" << best.transpose() << std::endl;

    return Result{best, bestVal, elapsedSec};
}
