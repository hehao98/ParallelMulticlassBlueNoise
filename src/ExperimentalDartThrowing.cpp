
/**
 * Dart throwing algorithm for multi-class samples
 *
 * It is based on codes from Li-Yi Wei[SIGGRAPH Multiclass Blue Noise]
 *
 * @author He, Hao
 */

#include <algorithm>
#include <fstream>
#include <iostream>

using namespace std;

#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "Exception.hpp"
#include "Grid.hpp"
#include "Timer.hpp"

#include "UniformDistanceField.hpp"

#include "UniformConflictChecker.hpp"

#include "Math.hpp"
#include "RMatrix.hpp"
#include "Random.hpp"
#include "Utility.hpp"

#include "SampleRecord.hpp"
#include "SequentialCounter.hpp"
#include "RandomCounter.hpp"

//#define DEBUG
//#define _RECORD_SAMPLE_HISTORY

int RandomClass(const vector<float> &cdf)
{
    const float value = Random::UniformRandom() * cdf[cdf.size() - 1];

    int selection = 0;

    for (unsigned int i = 0; i < cdf.size(); i++)
    {
        if (value > cdf[i])
        {
            selection = i + 1;
        }
        else
        {
            break;
        }
    }

    return selection % cdf.size();
}

struct MulticlassParameters
{
    // Values calculated from input paramters
    int dimension;
    int num_classes;
    int class_probability_all_int;
    vector<float> class_probability;
    vector<int> priority_values;
    vector<float> r_values;
    Array<float> r_matrix;
    float patience_factor;
    int k_number;
    vector<float> domain_size;
    UniformConflictChecker *conflict_checker;

    // Sample domain related information
    Grid::DomainSpec grid_domain_spec;
    Grid *grid;
    float min_r_value;
    int num_trials;
};

struct SamplingContext
{
    int total_num_samples;
    vector<int> num_samples_per_class;
    vector<int> current_failures_per_class;
};

int parseInput(MulticlassParameters &params, int argc, char **argv);
int setDomainSize(MulticlassParameters &params, vector<float> domain_size_spec);
vector<const Sample *> play(MulticlassParameters &params,
                            SamplingContext &context);

int main(int argc, char **argv)
{
    try
    {
        MulticlassParameters *params = new MulticlassParameters();
        SamplingContext *context = new SamplingContext();

        // Get parameters from input
        if (parseInput(*params, argc, argv) < 0)
        {
            cerr << "Something wrong with input, quiting..." << std::endl;
            return 1;
        }

        Timer timer;
        timer.Start();

        // Run the algorithm
        vector<const Sample *> samples = play(*params, *context);

        timer.Stop();
        const double total_time = timer.ElapsedTime();
        cerr << "total time " << total_time << endl;
        cerr << "# samples per second " << samples.size() / total_time << endl;

        // Output result
        for (unsigned int i = 0; i < samples.size(); i++)
        {
            cout << samples[i]->id;
            for (int j = 0; j < samples[i]->coordinate.Dimension(); j++)
            {
                cout << " " << samples[i]->coordinate[j];
            }
            cout << endl;
        }

        // Clean up
        for (auto &sample : samples) {
            delete sample;
        }
        samples.clear();
        delete params;
        delete context;
        return 0;
    }
    catch (Exception &e)
    {
        cerr << "Error: " << e.Message() << endl;
        return 1;
    }
}

int parseInput(MulticlassParameters &params, int argc, char **argv)
{
    // input arguments
    int min_argc = 7, argCtr = 0;

    if (argc < min_argc)
    {
        cerr << "Usage: " << argv[0]
             << " dimension num_classes (positive for optimal rmatrix computation, "
                "negative for uniform off-diagonal entries) priority (c "
                "integers with low values indicating higher priority) r_values "
                "(c*(c+1)/2 numbers in row major order of the upper matrix, or "
                "only c diagonal entries) k_number domain_size "
                "(dimension numbers)";
        cerr << endl;

        return -1;
    }

    params.dimension = atoi(argv[++argCtr]);
    if (params.dimension <= 0)
    {
        cerr << "dimension must be > 0" << endl;
        return -1;
    }

    char *num_classes_spec = argv[++argCtr];
    if (strlen(num_classes_spec) <= 0)
    {
        cerr << "empty specification for num_classes" << endl;
        return -1;
    }

    int num_classes_i = atoi(num_classes_spec);
    RMatrix::Method rmatrix_method =
        (num_classes_spec[0] == '0'
             ? RMatrix::GEOMETRIC
             : (num_classes_i > 0 ? RMatrix::OPTIMAL : RMatrix::UNIFORM));

    params.num_classes = num_classes_i > 0 ? num_classes_i : -num_classes_i;
    if (params.num_classes <= 0)
    {
        cerr << "num_classes must be > 0" << endl;
        return -1;
    }

    int class_probability_all_int = 1;
    vector<float> class_probability(params.num_classes);
    for (unsigned int i = 0; i < class_probability.size(); i++)
    {
        class_probability[i] = atof(argv[++argCtr]);
        if (floor(class_probability[i]) != class_probability[i])
        {
            class_probability_all_int = 0;
        }
    }

    vector<int> priority_values(params.num_classes);
    for (unsigned int i = 0; i < priority_values.size(); i++)
    {
        priority_values[i] = class_probability_all_int
                                 ? static_cast<int>(floor(class_probability[i]))
                                 : 0;

        if ((priority_values[i] < 0) ||
            (priority_values[i] >= params.num_classes))
        {
            cerr << "priority number must be within [0 num_classes-1]" << endl;
            return -1;
        }
    }

    if (class_probability_all_int)
    {
        for (unsigned int i = 0; i < class_probability.size(); i++)
        {
            class_probability[i] = (i + 1.0) / params.num_classes;
        }
    }
    else
    {
        for (unsigned int i = 1; i < class_probability.size(); i++)
        {
            class_probability[i] += class_probability[i - 1];
        }
    }

    cerr << "class_probablility cdf: ";
    for (unsigned int i = 0; i < class_probability.size(); i++)
    {
        cerr << class_probability[i] << " ";
    }
    cerr << endl;

    min_argc +=
        (params.dimension - 1) +
        (params.num_classes - 1); // we now know the # of domain_size arguments
    // and # of priority arguments

    int num_r_values = argc - min_argc + 1;

    if ((num_r_values != params.num_classes) &&
        (num_r_values != params.num_classes * (params.num_classes + 1) / 2))
    {
        cerr << "wrong number of r_values";
        return -1;
    }

    vector<float> input_r_values;
    for (int i = 0; i < num_r_values; i++)
    {
        input_r_values.push_back(atof(argv[++argCtr]));
    }

    params.r_matrix = RMatrix::BuildRMatrix(params.dimension, params.num_classes,
                                            rmatrix_method, input_r_values);

    vector<float> r_values;
    params.r_matrix.Get(r_values);

    for (unsigned int i = 0; i < r_values.size(); i++)
    {
        if (r_values[i] <= 0)
        {
            cerr << "bad r_matrix" << endl;
            return -1;
        }
    }

    int k_number = atof(argv[++argCtr]);

    const float patience_factor = max(
        static_cast<double>(1.0), abs(static_cast<double>(floor(k_number))));

    if (patience_factor != 0)
    {
        cerr << "patience_factor: " << patience_factor << endl;
    }

    // distance field
    vector<float> class_r_values(params.num_classes);
    {
        for (unsigned int i = 0; i < class_r_values.size(); i++)
        {
            if (!params.r_matrix.Get(vector<int>(2, i), class_r_values[i]))
            {
                throw Exception("cannot get value for class_r_values");
            }
        }
    }

    params.r_values = r_values;
    params.k_number = k_number;
    params.class_probability_all_int = class_probability_all_int;
    params.class_probability = class_probability;
    params.priority_values = priority_values;
    params.patience_factor = patience_factor;
    params.conflict_checker = new UniformConflictChecker(params.r_matrix);

    float min_r_value = params.r_values[0];
    for (unsigned int i = 0; i < params.r_values.size(); i++)
    {
        if (params.r_values[i] < min_r_value)
            min_r_value = params.r_values[i];
    }

    if (min_r_value <= 0)
    {
        cerr << "min_r_value <= 0" << endl;
        return -1;
    }
    params.min_r_value = min_r_value;

    vector<float> domain_size_spec;
    while (((argCtr + 1) < argc) &&
           (domain_size_spec.size() < params.dimension))
    {
        domain_size_spec.push_back(atof(argv[++argCtr]));
    }

    if (domain_size_spec.size() != params.dimension)
    {
        // should probably use assert...
        throw Exception("domain_size_spec.size() != dimension");
    }

    for (unsigned int i = 0; i < domain_size_spec.size(); i++)
    {
        if (domain_size_spec[i] <= 0)
        {
            cerr << "domain_size_spec[" << i << "] <= 0 - " << domain_size_spec[i]
                 << endl;
            return -1;
        }
    }
    params.domain_size = domain_size_spec;

    setDomainSize(params, domain_size_spec);

    return 0;
}

int setDomainSize(MulticlassParameters &params, vector<float> domain_size_spec)
{
    params.grid_domain_spec = Grid::BuildDomainSpec(domain_size_spec);
    params.grid = new Grid(params.grid_domain_spec, params.min_r_value);

    int num_trials = params.k_number;
    int total_num_grid_cells = 1;
    {
        Grid grid_local(params.grid_domain_spec, params.min_r_value);
        const vector<int> num_grid_cells = grid_local.NumCells();
        for (unsigned int i = 0; i < num_grid_cells.size(); i++)
        {
            num_trials *= num_grid_cells[i];
            total_num_grid_cells *= num_grid_cells[i];
        }
    }

    cerr << num_trials << " trials" << endl;

    params.num_trials = num_trials;
}

void generateSample(MulticlassParameters &params, SamplingContext &context,
                     vector<int> gridIndex)
{
    Sample *sample = new Sample(params.dimension);

    auto &total_num_samples = context.total_num_samples;
    auto &num_samples_per_class = context.num_samples_per_class;
    auto &current_failures_per_class = context.current_failures_per_class;


    vector<int> sample_ids;

    sample->id = RandomClass(params.class_probability);
    sample_ids.push_back(sample->id);

    if (sample->coordinate.Dimension() !=
        params.grid_domain_spec.domain_size.size())
    {
        throw Exception(
            "sample->coordinate.size() != grid_domain_spec.domain_size.size()");
    }

    vector<float> sample_domain_min_corner(params.dimension);
    vector<float> sample_domain_max_corner(params.dimension);

    for (int i = 0; i < params.dimension; i++)
    {
        sample_domain_min_corner[i] = gridIndex[i] * params.grid->CellSize();
        sample_domain_max_corner[i] = (gridIndex[i] + 1) * params.grid->CellSize();
    }
    for (int i = 0; i < sample->coordinate.Dimension(); i++)
    {
        sample->coordinate[i] =
            Random::UniformRandom() *
                (sample_domain_max_corner[i] - sample_domain_min_corner[i]) +
            sample_domain_min_corner[i];
    }

    if (sample_ids.size() <= 0)
    {
        throw Exception("sample_ids.size() <= 0");
    }

    SampleRecord::Status sample_fate = SampleRecord::UNKNOWN;

    for (unsigned int which_id = 0; (which_id < sample_ids.size()) &&
                                    (sample_fate != SampleRecord::ACCEPTED);
         which_id++)
    {
        sample->id = sample_ids[which_id];

        if (!params.grid->Inside(*sample))
        {
            sample_fate = SampleRecord::OUTSIDE;
        }
        else if (!params.grid->Conflict(*sample, *params.conflict_checker))
        {
            sample_fate = SampleRecord::ACCEPTED;
        }
        else // has conflict
        {
            sample_fate = SampleRecord::REJECTED;
        }

        if (sample_fate == SampleRecord::ACCEPTED)
        {
            if (!params.grid->Add(*sample))
            {
                throw Exception("cannot add sample");
            }
            else
            {
                total_num_samples++;
                num_samples_per_class[sample->id] += 1;
                sample = 0;
            }
            break;
        }
    }

    if (sample != 0)
    {
        delete sample;
        sample = 0;
    }
}

vector<const Sample *> play(MulticlassParameters &params,
                            SamplingContext &context)
{
    // Initialize random number engine
    const unsigned long random_seed = time(0);
    Random::InitRandomNumberGenerator(random_seed);

    // Initialize algorithm context
    context.total_num_samples = 0;
    context.num_samples_per_class = vector<int>(params.num_classes, 0);
    context.current_failures_per_class = vector<int>(params.num_classes, 0);

    RandomCounter counter(params.dimension, 0, int(params.domain_size[0]/params.grid->CellSize()));

    do {
        vector<int> gridIndex;
        counter.Get(gridIndex);
        generateSample(params, context, gridIndex);
    } while (counter.Next());
    // output
    vector<const Sample *> samples;
    params.grid->GetSamples(samples);

    if (context.total_num_samples != samples.size())
    {
        throw Exception("total_num_samples != samples.size()");
    }

    cerr << samples.size() << " samples" << endl;

    // done
    return samples;
}
