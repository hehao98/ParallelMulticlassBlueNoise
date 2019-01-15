
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

//#define DEBUG

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
    int dimension;
    int num_classes;
    int class_probability_all_int;
    vector<float> class_probability;
    vector<int> priority_values;
    vector<float> r_values;
    Array<float> r_matrix;
    float patience_factor;
    int k_number;

    UniformConflictChecker *conflict_checker;
    Grid::DomainSpec grid_domain_spec;
    Grid *grid;
    float min_r_value;

    int num_trials;
    int target_num_samples;
    vector<int> target_num_samples_per_class;
};

struct SamplingContext
{
    vector<PriorityGroup> priority_groups;
    int total_num_samples;
    vector<int> num_samples_per_class;
    vector<int> current_failures_per_class;
};

int parseInput(MulticlassParameters &params, int argc, char **argv);
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

        // Run the algorithm
        vector<const Sample *> samples = play(*params, *context);

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
                "negative for uniform off-diagonal entries) priority (either c "
                "integers with low values indicating higher priority, or c "
                "floating points indicating class selection probability) r_values "
                "(c*(c+1)/2 numbers in row major order of the upper matrix, or "
                "only c diagonal entries) k_number (positive integer for the usual "
                "k number, negative integer for target number of samples, [0 1] "
                "float for rho-number, or positive float for specifying both the "
                "k-number/patience-factor and the rho-number) domain_size "
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

#ifdef DEBUG
    cerr << "input_r_values: ";
    for (unsigned int i = 0; i < input_r_values.size(); i++)
    {
        cerr << input_r_values[i] << " ";
    }
    cerr << endl;
#endif

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

#ifdef DEBUG
    cerr << "r_values: ";
    for (unsigned int i = 0; i < r_values.size(); i++)
    {
        cerr << r_values[i] << " ";
    }
    cerr << endl;
#endif

    float k_rho_number = atof(argv[++argCtr]);
    float rho_number = k_rho_number - floor(k_rho_number);
    int k_number =
        (k_rho_number == floor(k_rho_number) ? static_cast<int>(k_rho_number)
                                             : 0);

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

    const float patience_factor = max(
        static_cast<double>(1.0), abs(static_cast<double>(floor(k_rho_number))));

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

    // init grid
    params.grid_domain_spec = Grid::BuildDomainSpec(domain_size_spec);

    params.r_values = r_values;
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

    params.k_number = k_number;
    params.class_probability_all_int = class_probability_all_int;
    params.class_probability = class_probability;
    params.priority_values = priority_values;
    params.patience_factor = patience_factor;
    params.min_r_value = min_r_value;
    params.conflict_checker = new UniformConflictChecker(params.r_matrix);
    params.grid = new Grid(params.grid_domain_spec, min_r_value);

    int num_trials = k_number;
    int total_num_grid_cells = 1;
    {
        Grid grid_local(params.grid_domain_spec, params.min_r_value);
        const vector<int> num_grid_cells = grid_local.NumCells();
        for (unsigned int i = 0; i < num_grid_cells.size(); i++)
        {
            num_trials *= num_grid_cells[i];
            total_num_grid_cells *= num_grid_cells[i];
        }

#ifdef DEBUG
        cerr << "num_grid_cells: ";
        for (unsigned int i = 0; i < num_grid_cells.size(); i++)
        {
            cerr << num_grid_cells[i] << " ";
        }
        cerr << endl;
#endif
    }

    int target_num_samples = -k_number;

    vector<int> target_num_samples_per_class(params.num_classes, 0);
    if ((rho_number > 0) && (rho_number < 1))
    {
        cerr << "rho_number " << rho_number << endl;

        for (int i = 0; i < params.num_classes; i++)
        {
            float value = 0;
            params.r_matrix.Get(vector<int>(2, i), value);

            target_num_samples_per_class[i] =
                Math::ComputeMaxNumSamples(params.dimension, value / rho_number);

            for (unsigned int j = 0; j < domain_size_spec.size(); j++)
            {
                target_num_samples_per_class[i] *= domain_size_spec[j];
            }
        }

        cerr << "target number of samples per class :";
        for (unsigned int i = 0; i < target_num_samples_per_class.size(); i++)
        {
            cerr << " " << target_num_samples_per_class[i];
        }
        cerr << endl;

        target_num_samples = 0;
        for (unsigned int i = 0; i < target_num_samples_per_class.size(); i++)
        {
            target_num_samples += target_num_samples_per_class[i];
        }

        cerr << target_num_samples << " target samples" << endl;
    }
    else
    {
        if (params.k_number >= 0)
        {
            cerr << num_trials << " trials" << endl;
        }
        else
        {
            cerr << target_num_samples << " target samples" << endl;
        }
    }

    params.num_trials = num_trials;
    params.target_num_samples = target_num_samples;
    params.target_num_samples_per_class = target_num_samples_per_class;

    return 0;
}

void generateSample(MulticlassParameters &params, SamplingContext &context,
                    int num_samples, int which_priority_group)
{
    Sample *sample = new Sample(params.dimension);

    const PriorityGroup &priority_group =
        context.priority_groups[which_priority_group];
    auto &total_num_samples = context.total_num_samples;
    auto &num_samples_per_class = context.num_samples_per_class;
    auto &current_failures_per_class = context.current_failures_per_class;

    if (priority_group.classes.size() <= 0)
    {
        throw Exception("priority_group.classes.size() <= 0");
    }

    vector<int> sample_ids;

    if (params.k_number == 0)
    {
        // pick the most under filled class
        float min_fill_ratio = 2.0;
        int min_fill_id = -1;

        for (unsigned int i = 0; i < priority_group.classes.size(); i++)
        {
            const int current_id = priority_group.classes[i];
            const float fill_ratio = num_samples_per_class[current_id] * 1.0 /
                                     params.target_num_samples_per_class[current_id];
            if (min_fill_ratio >= fill_ratio)
            {
                min_fill_ratio = fill_ratio;
                min_fill_id = current_id;
            }
        }

        if (min_fill_id < 0)
        {
            throw Exception("min_fill_id < 0");
        }
        else
        {
            sample->id = min_fill_id;
            sample_ids.push_back(sample->id);
        }
    }
    else if (params.class_probability_all_int)
    {
        sample->id = static_cast<int>(floor(Random::UniformRandom() *
                                            priority_group.classes.size())) %
                     priority_group.classes.size();
        sample->id = priority_group.classes[sample->id];
        sample_ids.push_back(sample->id);
    }
    else
    {
        if (context.priority_groups.size() != 1)
        {
            throw Exception("priority_groups.size() != 1");
        }

        sample->id = RandomClass(params.class_probability);
        sample_ids.push_back(sample->id);
    }

#ifdef DEBUG // useful debug
    cerr << "num_samples_per_class :";
    for (unsigned int i = 0; i < num_samples_per_class.size(); i++)
    {
        cerr << " " << num_samples_per_class[i];
    }
    cerr << " "; // << endl;
    cerr << "throw sample in class " << sample->id << endl;
#endif

    if (sample->coordinate.Dimension() !=
        params.grid_domain_spec.domain_size.size())
    {
        throw Exception(
            "sample->coordinate.size() != grid_domain_spec.domain_size.size()");
    }

    vector<float> sample_domain_min_corner(params.dimension);
    vector<float> sample_domain_max_corner(params.dimension);

    const int has_tried_hard_enough =
        (current_failures_per_class[sample->id] >
         params.patience_factor *
             params.target_num_samples_per_class[sample->id]);

    for (int i = 0; i < params.dimension; i++)
    {
        sample_domain_min_corner[i] = 0;
        sample_domain_max_corner[i] = params.grid_domain_spec.domain_size[i] *
                                      params.grid_domain_spec.cell_size;
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

#ifdef DEBUG
            if (current_failures_per_class[sample->id] > 0) // debug
            {
                cerr << "current_failures_per_class[" << sample->id
                     << "]: " << current_failures_per_class[sample->id] << endl;
            }
#endif

            if ((params.k_number == 0) && has_tried_hard_enough)
            {
                // cerr << "try killing" << endl; // debug
                vector<const Sample *> neighbors;
                if (!params.grid->GetConflicts(*sample, *params.conflict_checker,
                                               neighbors))
                {
                    throw Exception("cannot get conflicts");
                }

                int neighbors_all_removable = 1;

                const float sample_fill_ratio =
                    num_samples_per_class[sample->id] * 1.0 /
                    params.target_num_samples_per_class[sample->id];
                float sample_r_value = 0;
                params.r_matrix.Get(vector<int>(2, sample->id), sample_r_value);

                for (unsigned int j = 0; j < neighbors.size(); j++)
                {
                    const Sample &current_neighbor = *neighbors[j];

                    float neighbor_r_value = 0;
                    params.r_matrix.Get(vector<int>(2, current_neighbor.id),
                                        neighbor_r_value);

                    if (neighbor_r_value >= sample_r_value)
                    {
                        neighbors_all_removable = 0;
                        break;
                    }

                    float current_neighbor_fill_ratio = 0;

                    for (unsigned int i = 0; i < priority_group.classes.size(); i++)
                    {
                        const int current_id = priority_group.classes[i];
                        if (current_neighbor.id == current_id)
                        {
                            current_neighbor_fill_ratio =
                                num_samples_per_class[current_id] * 1.0 /
                                params.target_num_samples_per_class[current_id];
                        }
                    }

                    if ((current_neighbor_fill_ratio <= 0) ||
                        (current_neighbor_fill_ratio < sample_fill_ratio))
                    {
                        neighbors_all_removable = 0;
                        break;
                    }
                }

                if (neighbors_all_removable)
                {
                    for (unsigned int j = 0; j < neighbors.size(); j++)
                    {
                        if (!params.grid->Remove(*neighbors[j]))
                        {
                            throw Exception("cannot remove neighbor");
                        }
                        num_samples--;
                        total_num_samples--;
                        num_samples_per_class[neighbors[j]->id] -= 1;

#ifdef _RECORD_SAMPLE_HISTORY
                        const SampleRecord new_record(neighbors[j]->id,
                                                      SampleRecord::KILLED);
                        sample_history.push_back(new_record);
#endif

                        delete neighbors[j];
                    }

                    neighbors.clear();
                    if (!params.grid->Conflict(*sample, *params.conflict_checker))
                    {
                        sample_fate = SampleRecord::ACCEPTED;
                    }
                    else
                    {
                        throw Exception("still conflict after removing neighbors");
                    }
                }
                else
                {
                    // cerr << "cannot remove neighbors" << endl; // debug
                }
            }
        }

        if (sample_fate == SampleRecord::ACCEPTED)
        {
            current_failures_per_class[sample->id] = 0;
        }
        else
        {
            current_failures_per_class[sample->id] += 1;
        }

        switch (sample_fate)
        {
        case SampleRecord::ACCEPTED:
        {
#ifdef _RECORD_SAMPLE_HISTORY
            const SampleRecord new_record(sample->id, SampleRecord::ACCEPTED);
            sample_history.push_back(new_record);
#endif
            if (!params.grid->Add(*sample))
            {
                throw Exception("cannot add sample");
                delete sample;
                sample = 0;
            }
            else
            {
                num_samples++;
                total_num_samples++;
                num_samples_per_class[sample->id] += 1;
                sample = 0;
            }
            break;
        }

        case SampleRecord::REJECTED:
        {
#ifdef _RECORD_SAMPLE_HISTORY
            const SampleRecord new_record(sample->id, SampleRecord::REJECTED);
            sample_history.push_back(new_record);
#endif
            break;
        }

        case SampleRecord::OUTSIDE:
            // do nothing
            break;

        default:
            throw Exception("unknown sample fate");
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
    context.priority_groups = BuildPriorityGroups(
        params.dimension, params.priority_values, params.r_matrix,
        params.num_trials, params.target_num_samples);
    context.total_num_samples = 0;
    context.num_samples_per_class = vector<int>(params.num_classes, 0);
    context.current_failures_per_class = vector<int>(params.num_classes, 0);

#ifdef _RECORD_SAMPLE_HISTORY
    // record class numbers for each sample drawn
    vector<SampleRecord> sample_history;
#endif

    Timer timer;

    timer.Start();

    for (unsigned int which_priority_group = 0;
         which_priority_group < context.priority_groups.size();
         which_priority_group++)
    {
        const PriorityGroup &priority_group =
            context.priority_groups[which_priority_group];

#ifdef DEBUG
        cerr << "classes: (";
        for (unsigned int i = 0; i < priority_group.classes.size(); i++)
        {
            cerr << priority_group.classes[i] << " ";
        }
        cerr << "), k_number: " << params.k_number
             << ", num_trials: " << priority_group.num_trials
             << ", num_samples: " << priority_group.target_num_samples << endl;
#endif

        int num_samples = 0;
        for (unsigned int i = 0; i < priority_group.classes.size(); i++)
        {
            num_samples += context.num_samples_per_class[priority_group.classes[i]];
        }

        for (int trial = 0;
             ((params.k_number > 0) && (trial < priority_group.num_trials)) ||
             ((params.k_number <= 0) &&
              (num_samples < priority_group.target_num_samples));
             trial++)
        {
            generateSample(params, context, num_samples, which_priority_group);
        }
    }

    timer.Stop();
    const double total_time = timer.ElapsedTime();
    cerr << "total time " << total_time << endl;

    // output
    vector<const Sample *> samples;
    params.grid->GetSamples(samples);

    if (context.total_num_samples != samples.size())
    {
        throw Exception("total_num_samples != samples.size()");
    }

    cerr << samples.size() << " samples" << endl;

#ifdef _RECORD_SAMPLE_HISTORY
    if (history_file_name)
    {
        if (!WriteSampleHistory(history_file_name, sample_history))
        {
            cerr << "cannot write to " << history_file_name << endl;
        }
    }
    else
    {
        cerr << "null history file name" << endl;
    }
#endif

    cerr << "# samples per second " << samples.size() / total_time << endl;

    // done
    return samples;
}
