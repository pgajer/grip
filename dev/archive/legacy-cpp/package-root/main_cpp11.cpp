#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <memory>
#include "DrawGraph.hpp"
#include "MesaPlot.hpp"

enum engine_t {
    MISH_v5 = 12,
    MISH_v6 = 13
};

struct configuration_t {
    size_tt numOfVert = 5;
    size_tt rounds = 30;
    size_tt finalRounds = 30;
    size_tt tinitFactor = 6;
    size_tt p1 = 0;
    size_tt p2 = 0;
    size_tt numOfItr = 1;
    size_tt dim = 3;
    engine_t engf = MISH_v5;
    size_tt slow = 0;
    float roundsPar = 1.0f;
    coord_t edge = 32;
    bool readFromFile = false;
    bool printToFile = false;
    std::string fileName;
    std::string posFile;
    std::string graphStr = "complete";
    size_tt par = 3;
    size_tt displayPar = 1;
    int width = 700;
    int height = 700;
    size_tt numOfInitVert = 15;
    size_tt numOfNbrs = 10;
    float r = 0.15f;
    float s = 3.0f;
    size_tt T = 4;
    bool color = false;
};

void print_dg_config_pars(const configuration_t& config) { // const DrawGraph& dg
    std::cout << "DrawGraph configuration_t Parameters:\n";
    std::cout << "------------------------------------\n";
    std::cout << "Number of vertices: " << config.numOfVert << "\n";
    std::cout << "Dimension: " << config.dim << "\n";
    std::cout << "Rounds: " << config.rounds << "\n";
    std::cout << "Final rounds: " << config.finalRounds << "\n";
    std::cout << "Initial temperature factor: " << config.tinitFactor << "\n";
    std::cout << "Engine function: " << config.engf << "\n";
    std::cout << "Number of initial vertices: " << config.numOfInitVert << "\n";
    std::cout << "Number of neighbors: " << config.numOfNbrs << "\n";
    std::cout << "r parameter: " << config.r << "\n";
    std::cout << "s parameter: " << config.s << "\n";
    std::cout << "Display parameter: " << config.displayPar << "\n";
    std::cout << "------------------------------------\n";

    // Add any additional parameters from the DrawGraph object
    // that are not part of the configuration_t struct.
    // std::cout << "Some DrawGraph specific parameter: " << dg.getSomeParameter() << "\n";
}

void process_command_line_args(int argc, char** argv, configuration_t& config);
std::unique_ptr<Graph> create_graph(const configuration_t& config);
void run_engine(DrawGraph& dg, engine_t engf);

int main(int argc, char** argv) {
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

    configuration_t config;
    process_command_line_args(argc, argv, config);

    std::unique_ptr<Graph> graph = create_graph(config);

    // print_dg_config_pars(config);

    DrawGraph dg(*graph,
                 config.dim,
                 config.rounds,
                 config.finalRounds,
                 config.tinitFactor,
                 static_cast<size_tt>(config.engf),
                 config.numOfInitVert,
                 config.numOfNbrs,
                 config.r,
                 config.s,
                 config.displayPar);

    if (config.displayPar == 1) {
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        std::cout << "Elapsed time: " << elapsed.count() << " s\n";

        MesaPlot mp(graph.get(), &dg, static_cast<size_tt>(config.engf), config.slow, config.width, config.height, config.color);
        mp.draw_Graph(argc, argv);

    } else if (config.displayPar == 0) {

        run_engine(dg, config.engf);

        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        std::cout << "Elapsed time: " << elapsed.count() << " s\n";

        MesaPlot mp(graph.get(), &dg, 0, config.slow, config.width, config.height, config.color);
        mp.draw_Graph(argc, argv);
    }

    return 0;
}

void process_command_line_args(int argc, char** argv, configuration_t& config) {
    for (int arg = 1; arg < argc; ++arg) {
        std::string option = argv[arg];
        if (option[0] == '-') {
            switch (option[1]) {
                case '#': config.numOfVert = std::stoul(argv[++arg]); break;
                case 'r': config.rounds = std::stoul(argv[++arg]); break;
                case 'R': config.finalRounds = std::stoul(argv[++arg]); break;
                case 'n': config.numOfItr = std::stoul(argv[++arg]); break;
                case 'i': config.numOfInitVert = std::stoul(argv[++arg]); break;
                case 'b': config.numOfNbrs = std::stoul(argv[++arg]); break;
                case 'c': config.color = static_cast<bool>(std::stoi(argv[++arg])); break;
                case 'S': config.s = std::stof(argv[++arg]); break;
                case 'T': config.T = std::stoul(argv[++arg]); break;
                case 'd': config.dim = std::stoul(argv[++arg]); break;
                case 'k':
                    if (std::string(argv[arg + 1]) == "mish_v5") config.engf = MISH_v5;
                    else if (std::string(argv[arg + 1]) == "mish_v6") config.engf = MISH_v6;
                    else config.engf = static_cast<engine_t>(std::stoul(argv[++arg]));
                    break;
                case 'e': config.edge = std::stoul(argv[++arg]); break;
                case 'w': config.width = std::stoi(argv[++arg]); break;
                case 't': config.tinitFactor = std::stoul(argv[++arg]); break;
                case 'g': config.graphStr = argv[++arg]; break;
                case 's': config.slow = std::stoul(argv[++arg]); break;
                case 'P': config.par = std::stoul(argv[++arg]); break;
                case 'p':
                    config.printToFile = true;
                    if (arg + 1 < argc && argv[arg + 1][0] != '-') config.posFile = argv[++arg];
                    break;
                case 'D': config.displayPar = std::stoul(argv[++arg]); break;
                default:
                    throw std::runtime_error("Illegal option: " + std::string(1, option[1]));
            }
        }
    }
}

std::unique_ptr<Graph> create_graph(const configuration_t& config) {
    std::unique_ptr<Graph> graph(new Graph());

    if (config.graphStr == "complete") {
        graph->complete_Graph(config.numOfVert);
    } else if (config.graphStr == "random") {
        graph->rand_Graph(config.numOfVert, config.T);
    } else if (config.graphStr == "path") {
        graph->path_Graph(config.numOfVert);
    } else if (config.graphStr == "cycle") {
        graph->cycle_Graph(config.numOfVert);
    } else if (config.graphStr == "hypercube") {
        graph->hyper_Cube(config.numOfVert);
    } else if (config.graphStr == "cylinder") {
        graph->square_Cylinder(config.numOfVert, config.T);
    } else if (config.graphStr == "torus") {
        graph->torus(config.numOfVert, config.T);
    } else if (config.graphStr == "mesh") {
        graph->mesh(config.numOfVert);
    } else if (config.graphStr == "twistedtorus") {
        graph->twistedTorus(config.numOfVert, config.T, config.p1, config.p2);
    } else if (config.graphStr == "tree") {
        graph->tree(config.numOfVert, config.T);
    } else if (config.graphStr == "meshX") {
        graph->meshX(config.numOfVert);
    } else if (config.graphStr == "meshT") {
        graph->meshT(config.numOfVert);
    } else if (config.graphStr == "meshTX") {
        graph->meshTX(config.numOfVert);
    } else if (config.graphStr == "sierpinski") {
        graph->sierpinski(config.numOfVert, config.T);
    } else if (config.graphStr == "read_graph") {
        graph->read_Graph_From_IWB_File(config.graphStr.c_str());
    } else {
        // Print the value of config.graphStr
        std::cout << "Unrecognized graph type: " << config.graphStr << std::endl;

        // State that the string is not found
        std::cout << "The specified graph type was not found." << std::endl;

        // List allowed strings
        std::cout << "Allowed graph types are: " << std::endl;
        std::cout << "complete, random, path, cycle, hypercube, cylinder, torus, mesh, "
                  << "twistedtorus, tree, meshX, meshT, meshTX, sierpinski, read_graph" << std::endl;

        // Optionally, you might want to throw an exception or return a null pointer here
        // to indicate that graph creation failed
        return nullptr;
    }

    return graph;
}

void run_engine(DrawGraph& dg, engine_t engf) {
    switch (engf) {
        case MISH_v5:
            dg.mish_engine_v5();
            break;
        case MISH_v6:
            dg.mish_engine_v6();
            break;
        default:
            throw std::runtime_error("Unknown engine type");
    }
}
