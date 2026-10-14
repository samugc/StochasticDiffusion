#include "evolutive.hpp"

using namespace std;
using namespace Eigen;
using namespace std::chrono;

// Generador aleatorio
std::mt19937 gen(std::random_device{}());

Individual randomIndividual(const MatrixXd& bounds, const VectorXd& alpha) {
    int dim = bounds.cols();
    VectorXd genes(dim);
    for (int i = 0; i < dim; ++i)
        genes[i] = randDouble(bounds(0, i), bounds(1, i));

    double fitness = logLikelihood(alpha, genes);
    return {genes, fitness};
}

Individual crossover(const Individual& p1, const Individual& p2,
                     const MatrixXd& bounds, const VectorXd& alpha, double alpha_blx) {
    int dim = p1.genes.size();
    VectorXd childGenes(dim);
    for (int i = 0; i < dim; ++i) {
        double cMin = std::min(p1.genes[i], p2.genes[i]);
        double cMax = std::max(p1.genes[i], p2.genes[i]);
        double I = cMax - cMin;
        double low = std::max(bounds(0, i), cMin - alpha_blx * I);
        double high = std::min(bounds(1, i), cMax + alpha_blx * I);
        childGenes[i] = randDouble(low, high);
    }
    double fitness = logLikelihood(alpha, childGenes);
    return {childGenes, fitness};
}

void mutate(Individual& ind, const MatrixXd& bounds, double mutationRate, double mutationStep) {
    std::normal_distribution<> noise(0.0, mutationStep);
    for (int i = 0; i < ind.genes.size(); ++i) {
        if (randDouble(0.0, 1.0) < mutationRate) {
            ind.genes[i] += noise(gen);
            // Respetar los límites
            ind.genes[i] = std::min(std::max(ind.genes[i], bounds(0, i)), bounds(1, i));
        }
    }
}

Individual tournament(const std::vector<Individual>& pop) {
    int i = rand() % pop.size();
    int j = rand() % pop.size();
    return pop[i].fitness > pop[j].fitness ? pop[i] : pop[j];
}

Result evolutionaryAlgorithm(
    const VectorXd& alpha,
    const MatrixXd& bounds,
    int popSize,
    int generations,
    double crossoverRate,
    double mutationRate,
    double mutationStep)
{
    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    std::vector<Individual> population;

    // Inicializar población
    for (int i = 0; i < popSize; ++i)
        population.push_back(randomIndividual(bounds, alpha));

    // Evolución
    for (int gen = 0; gen < generations; ++gen) {
        std::vector<Individual> newPop;

        while (static_cast<int>(newPop.size()) < popSize) {
            Individual parent1 = tournament(population);
            Individual parent2 = tournament(population);

            Individual child;
            if (randDouble(0.0, 1.0) < crossoverRate)
                child = crossover(parent1, parent2, bounds, alpha);
            else
                child = parent1;

            mutate(child, bounds, mutationRate, mutationStep);
            child.fitness = logLikelihood(alpha, child.genes);

            newPop.push_back(child);
        }

        population = std::move(newPop);
        std::sort(population.begin(), population.end()); // menor a mayor fitness
    }

    auto end = high_resolution_clock::now(); 
    double elapsed = duration<double>(end - start).count();

    const Individual& best = population.back(); // mayor fitness (asumiendo orden ascendente)

    std::cout << "Tiempo de ejecución EA: " << elapsed << " segundos.\n";
    std::cout << "Mejor solución encontrada EA:\n" << best.genes.transpose() << std::endl;

    return Result{best.genes, best.fitness, elapsed};
}

Result evolutionaryAlgorithmWithElitism(
    const VectorXd& alpha,
    const MatrixXd& bounds,
    int popSize,
    int generations,
    double crossoverRate,
    double mutationRate,
    double mutationStep)
{
    using namespace std::chrono;
    auto start = high_resolution_clock::now(); 

    std::vector<Individual> population;
    for (int i = 0; i < popSize; ++i)
        population.push_back(randomIndividual(bounds, alpha));

    for (int gen = 0; gen < generations; ++gen) {
        std::sort(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
            return a.fitness > b.fitness; // orden descendente: mejor primero
        });

        // ELITISMO: conservar los mejores
        int eliteSize = 2;
        std::vector<Individual> newPopulation(population.begin(), population.begin() + eliteSize);

        while ((int)newPopulation.size() < popSize) {
            Individual parent1 = tournament(population);
            Individual parent2 = tournament(population);

            Individual child;
            if (randDouble(0.0, 1.0) < crossoverRate)
                child = crossover(parent1, parent2, bounds, alpha);
            else
                child = parent1;

            mutate(child, bounds, mutationRate, mutationStep);
            child.fitness = logLikelihood(alpha, child.genes);

            newPopulation.push_back(child);
        }

        population = std::move(newPopulation);
    }

    std::sort(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
        return a.fitness > b.fitness; // asegurar que el mejor está al frente
    });

    auto end = high_resolution_clock::now(); 
    double elapsed = duration<double>(end - start).count();

    const Individual& best = population.front();

    std::cout << "Tiempo de ejecución EA con elitismo: " << elapsed << " segundos.\n";
    std::cout << "Mejor solución encontrada:\n" << best.genes.transpose() << std::endl;

    return Result{best.genes, best.fitness, elapsed};
}

Result memeticAlgorithmWithElitism(
    const VectorXd& alpha,
    const MatrixXd& bounds,
    int popSize,
    int generations,
    double crossoverRate,
    double mutationRate,
    double mutationStep,
    int localSearchInterval,
    int localSearchSteps)
{
    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    std::vector<Individual> population;
    for (int i = 0; i < popSize; ++i)
        population.push_back(randomIndividual(bounds, alpha));

    for (int gen = 0; gen < generations; ++gen) {
        std::sort(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
            return a.fitness > b.fitness;
        });

        int eliteSize = 2;
        std::vector<Individual> newPopulation(population.begin(), population.begin() + eliteSize);

        while ((int)newPopulation.size() < popSize) {
            Individual parent1 = tournament(population);
            Individual parent2 = tournament(population);

            Individual child;
            if (randDouble(0.0, 1.0) < crossoverRate)
                child = crossover(parent1, parent2, bounds, alpha);
            else
                child = parent1;

            mutate(child, bounds, mutationRate, mutationStep);
            child.fitness = logLikelihood(alpha, child.genes);

            newPopulation.push_back(child);
        }

        population = std::move(newPopulation);

        // Búsqueda local periódica
        if (gen % localSearchInterval == 0 && gen > 0) {
            Individual& best = population[0];
            Result local = localSearch(alpha, best.genes, bounds, localSearchSteps, 0.01);
            best.genes = local.solucion;
            best.fitness = local.fitness;
        }
    }

    std::sort(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
        return a.fitness > b.fitness;
    });

    auto end = high_resolution_clock::now();
    double elapsed = duration<double>(end - start).count();

    const Individual& best = population.front();

    std::cout << "Tiempo de ejecución memético: " << elapsed << " segundos.\n";
    std::cout << "Mejor solución encontrada:\n" << best.genes.transpose() << std::endl;

    return Result{best.genes, best.fitness, elapsed};
}

Result differentialEvolution(
    const VectorXd& alpha,
    const MatrixXd& bounds,
    int popSize,
    int generations,
    double F,
    double CR)
{
    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    int D = bounds.cols(); // Dimensión
    std::vector<Individual> population(popSize);

    // Inicialización aleatoria
    for (int i = 0; i < popSize; ++i) {
        population[i].genes = VectorXd(D);
        for (int j = 0; j < D; ++j) {
            population[i].genes(j) = randDouble(bounds(0, j), bounds(1, j));
        }
        population[i].fitness = logLikelihood(alpha, population[i].genes);
    }

    for (int gen = 0; gen < generations; ++gen) {
        std::vector<Individual> newPopulation = population;

        for (int i = 0; i < popSize; ++i) {
            int a, b, c;
            do { a = rand() % popSize; } while (a == i);
            do { b = rand() % popSize; } while (b == i || b == a);
            do { c = rand() % popSize; } while (c == i || c == a || c == b);

            VectorXd mutant = population[a].genes + F * (population[b].genes - population[c].genes);

            // Crossover
            VectorXd trial = population[i].genes;
            for (int j = 0; j < D; ++j) {
                if (randDouble(0.0, 1.0) < CR) {
                    trial(j) = std::clamp(mutant(j), bounds(0, j), bounds(1, j));
                }
            }

            double trialFitness = logLikelihood(alpha, trial);

            // Selección
            if (trialFitness > population[i].fitness) {
                newPopulation[i].genes = trial;
                newPopulation[i].fitness = trialFitness;
            }
        }

        // Elitismo opcional
        auto bestOld = *std::max_element(population.begin(), population.end(),
                          [](const Individual& a, const Individual& b) { return a.fitness < b.fitness; });
        auto bestNew = *std::max_element(newPopulation.begin(), newPopulation.end(),
                          [](const Individual& a, const Individual& b) { return a.fitness < b.fitness; });
        if (bestOld.fitness > bestNew.fitness) {
            int worstIdx = std::min_element(newPopulation.begin(), newPopulation.end(),
                              [](const Individual& a, const Individual& b) { return a.fitness < b.fitness; }) - newPopulation.begin();
            newPopulation[worstIdx] = bestOld;
        }

        population = std::move(newPopulation);
    }

    auto end = high_resolution_clock::now();
    double elapsed = duration<double>(end - start).count();

    auto best = *std::max_element(population.begin(), population.end(),
                  [](const Individual& a, const Individual& b) { return a.fitness < b.fitness; });

    return {best.genes, best.fitness, elapsed};
}

VectorXd adaptiveMemeticWithSmartRestarts(
    const VectorXd& alpha,
    const MatrixXd& bounds,
    int popSize,
    int generations,
    double crossoverRate,
    double mutationRate,
    double mutationStep,
    int localSearchInterval,
    int localSearchSteps,
    int maxNoImprovement)
{
    auto start = high_resolution_clock::now();

    std::vector<Individual> population;
    for (int i = 0; i < popSize; ++i)
        population.push_back(randomIndividual(bounds, alpha));

    int noImprovementCount = 0;
    double bestVal = -1e10;
    VectorXd bestGenes;

    for (int gen = 0; gen < generations; ++gen) {
        std::sort(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
            return a.fitness > b.fitness;
        });

        if (population[0].fitness > bestVal) {
            bestVal = population[0].fitness;
            bestGenes = population[0].genes;
            noImprovementCount = 0;
        } else {
            ++noImprovementCount;
        }

        // Reinicio inteligente si no mejora
        if (noImprovementCount >= maxNoImprovement) {
            for (int i = 1; i < popSize; ++i)
                population[i] = randomIndividual(bounds, alpha);
            noImprovementCount = 0;
        }

        int eliteSize = 2;
        std::vector<Individual> newPop(population.begin(), population.begin() + eliteSize);

        while ((int)newPop.size() < popSize) {
            Individual p1 = tournament(population);
            Individual p2 = tournament(population);

            Individual child = (randDouble(0.0, 1.0) < crossoverRate)
                ? crossover(p1, p2, bounds, alpha)
                : p1;

            mutate(child, bounds, mutationRate, mutationStep);
            child.fitness = logLikelihood(alpha, child.genes);
            newPop.push_back(child);
        }

        if (gen % localSearchInterval == 0 && gen > 0) {
            newPop[0].genes = localSearch(alpha, newPop[0].genes, bounds, localSearchSteps, 0.01).solucion;
            newPop[0].fitness = logLikelihood(alpha, newPop[0].genes);
        }

        population = std::move(newPop);
    }

    auto end = high_resolution_clock::now();
    std::cout << "Tiempo de ejecución Memético Adaptativo: " << duration<double>(end - start).count() << " segundos.\n";
    return bestGenes;
}
