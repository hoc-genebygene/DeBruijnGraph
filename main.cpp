#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <memory>
#include <random>
#include <stack>
#include <string>
#include <vector>

namespace {
std::vector<std::string> GenerateKmers(int K, const std::string & input_string) {
    std::vector<std::string> kmers;

    for (int k = 0; k <= input_string.size() - K; ++k) {
        kmers.push_back(input_string.substr(k, K));
    }

    return kmers;
}

std::string GetKMerPrefix(const std::string & kmer) {
    return kmer.substr(0, kmer.size() - 1);
}

std::string GetKMerSuffix(const std::string & kmer) {
    return kmer.substr(1, kmer.size() - 1);
}
}

class DeBruijnGraph {
private:
    struct Edge;

    struct Node {
        std::vector<Edge *> in_edges;
        std::vector<std::unique_ptr<Edge>> out_edges;

        std::string data;
    };

    struct Edge {
        Node * to_node;
        uint64_t weight = 0;
        uint64_t times_visited = 0;
    };

public:
    DeBruijnGraph() {
        //kmer_node_map_.max_load_factor(2.0f);
    }

    Node * FindNode(const std::string & kmer) {
        auto it = kmer_node_map_.find(kmer);

        if (it != kmer_node_map_.end()) {
            return it->second.get();
        } else {
            return nullptr;
        }
    }

    Node * AddNode(const std::string & kmer) {
        std::unique_ptr<Node> node = std::make_unique<Node>();
        node->data = kmer;
        Node * node_raw = node.get();
        kmer_node_map_[kmer] = std::move(node);

        return node_raw;
    }

    // Find Edge between nodes
    Edge * FindEdge(Node * prefix_node, Node * suffix_node) {
        for (auto & edge : prefix_node->out_edges) {
            if (edge->to_node == suffix_node) {
                return edge.get();
            }
        }

        return nullptr;
    }

    void ConnectNodes(Node * prefix_node, Node * suffix_node) {
        auto edge = FindEdge(prefix_node, suffix_node);

        if (edge == nullptr) {
            auto edge = std::make_unique<Edge>();
            edge->to_node = suffix_node;
            edge->weight = 1;

            suffix_node->in_edges.push_back(edge.get());
            prefix_node->out_edges.push_back(std::move(edge));
        } else {
            ++(edge->weight);
        }

    }

    std::vector<Node *> FindSemiBalancedNodes() {
        std::vector<Node *> semi_balanced_nodes;

        for (auto & kmer_node_pair : kmer_node_map_) {
            int64_t num_in_edges = 0;
            for (auto in_edge : kmer_node_pair.second->in_edges) {
                num_in_edges += in_edge->weight;
            }

            int64_t num_out_edges = 0;
            for (auto & out_edge : kmer_node_pair.second->out_edges) {
                num_out_edges += out_edge->weight;
            }

            if (num_in_edges != num_out_edges) {
                assert(std::abs(num_in_edges - num_out_edges) == 1);

                semi_balanced_nodes.push_back(kmer_node_pair.second.get());
            }
        }

        return semi_balanced_nodes;
    }

    void BuildGraph(const std::vector<std::string> & kmers) {
        for (const auto & kmer : kmers) {

            auto prefix = GetKMerPrefix(kmer);
            auto suffix = GetKMerSuffix(kmer);

            // Find the nodes
            auto prefix_node = FindNode(prefix);

            // Node doesn't exist so generate it
            if (prefix_node == nullptr) {
                prefix_node = AddNode(prefix);
            }

            auto suffix_node = FindNode(suffix);

            if (suffix_node == nullptr) {
                suffix_node = AddNode(suffix);
            }

            ConnectNodes(prefix_node, suffix_node);
        }
    }

    std::vector<Node *> FindEulerianPath() {
        auto semi_balanced_nodes = FindSemiBalancedNodes();
        assert(semi_balanced_nodes.size() == 2 || semi_balanced_nodes.size() == 0); // 2 means Eulerian path, 0 means Eulerian cycle

        Node * start_node = nullptr;
        Node * stop_node = nullptr;
        // Pick the start node with more out_edges than in_edges, it's pretty obvious why this must be the start node...

        if (semi_balanced_nodes.size() == 2) {
            if (semi_balanced_nodes[0]->out_edges.size() > semi_balanced_nodes[0]->in_edges.size()) {
                start_node = semi_balanced_nodes[0];
                stop_node = semi_balanced_nodes[1];
            } else {
                start_node = semi_balanced_nodes[1];
                stop_node = semi_balanced_nodes[0];
            }
        } else {
            start_node = kmer_node_map_.begin()->second.get();
        }

        // Now begin Hierholzer's algorithm
        std::stack<Node *> temp_path;
        std::stack<Node *> final_path;

        auto AllOutEdgesVisited = [](Node * node) -> bool {
            for (auto & out_edge : node->out_edges) {
                if (out_edge->times_visited < out_edge->weight) {
                    return false;
                }
            }

            return true;
        };

        Node * current_node = start_node;
        while (current_node != start_node || !AllOutEdgesVisited(current_node)) {
            // Find a cycle (or we hit the stop_node)
            do {
                for (auto & out_edge : current_node->out_edges) {
                    if (out_edge->times_visited < out_edge->weight) {
                        ++(out_edge->times_visited);
                        temp_path.push(current_node);
                        current_node = out_edge->to_node;
                        break;
                    }
                }
            }
            while (!AllOutEdgesVisited(current_node));

            // Now backtrack
            while (current_node != start_node || !AllOutEdgesVisited(current_node)) {
                final_path.push(temp_path.top());
                temp_path.pop();
                current_node = temp_path.top();
            }
        }

        // Push everything remaining on temp_path to final_path
        while (temp_path.size() != 0) {
            final_path.push(temp_path.top());
            temp_path.pop();
        }

        std::vector<Node *> eulerian_path;
        while (final_path.size() != 0) {
            eulerian_path.push_back(final_path.top());
            final_path.pop();
        }
        if (semi_balanced_nodes.size() == 2) {
            eulerian_path.push_back(stop_node);
        }

        return eulerian_path;
    }

    void PrintNodes() {
        for (auto & kmer_node_pair : kmer_node_map_) {
            for (auto & out_edge : kmer_node_pair.second->out_edges) {
                std::cout << kmer_node_pair.second->data << " -> " << out_edge->to_node->data << " " << out_edge->weight << std::endl;
            }
        }
    }

private:
    std::unordered_map<std::string, std::unique_ptr<Node>> kmer_node_map_;
};

std::string GenerateRandomNucleotideString(size_t length) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 3);

    std::vector<char> vec(length);

    std::generate_n(vec.begin(), length, [&dis, &gen]() {
        switch(dis(gen)) {
            case 0:
                return 'A';
            case 1:
                return 'C';
            case 2:
                return 'G';
            case 3:
                return 'T';
            default:
                throw std::logic_error("Random number generator too large");
        }
    });

    return std::string(vec.begin(), vec.end());
}

int main() {
std::string input_string = GenerateRandomNucleotideString(250 * 1000000);
    //std::string input_string = "a_long_long_long_time";
    //std::cout << input_string << std::endl;
    int K = 25;
    std::vector<std::string> kmers = GenerateKmers(K, input_string);

    DeBruijnGraph graph;
    {
        auto start = std::chrono::steady_clock::now();
        graph.BuildGraph(kmers);
        auto stop = std::chrono::steady_clock::now();

        std::cout << "Time to build DeBruijnGraph: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " ms" << std::endl;
    }
    //graph.PrintNodes();

    {
        auto start = std::chrono::steady_clock::now();
        auto eulerian_path = graph.FindEulerianPath();
        auto stop = std::chrono::steady_clock::now();
    //    std::cout << "Eulerian path: " << std::endl;
    //    for (auto node : eulerian_path) {
    //        std::cout << node->data << " ";
    //    }
    //    std::cout << std::endl;
        std::cout << "Time to find Eulerian path: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " ms" << std::endl;
    }

    return 0;
}
