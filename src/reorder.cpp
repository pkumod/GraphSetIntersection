#include "porder.hpp"
#include "util.hpp"

POrder porder;

int main(int argc, char* argv[])
{
    if (argc == 1) {
        printf("Please provide graph file path.\n");
        quit();
    }

    std::string graph_file_path(argv[1]);
    std::string graph_name = extract_filename(graph_file_path);
    std::string order_option = "gro";
    std::string reordered_graph_file_path = graph_name + "_GRO.txt";
    std::string reordered_newid_file_path = graph_name + "_GRO_newID.txt";
    int i;
    if ((i = arg_pos((char *)"-order", argc, argv)) > 0)
        order_option = std::string(argv[i + 1]);
    else if ((i = arg_pos((char *)"-ratio", argc, argv)) > 0) {
        EdgeVector edge_vec = load_graph(graph_file_path);
        porder.load_org_graph(edge_vec);
        porder.comp_ratio();
        return 0;
    }
        
    struct timeval time_start;
    struct timeval time_end;

    gettimeofday(&time_start, NULL);
    EdgeVector edge_vec = load_graph(graph_file_path);
    gettimeofday(&time_end, NULL);
    double read_time = (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;
    printf("read_time=%.3fms\n", read_time);

    gettimeofday(&time_start, NULL);
    porder.load_org_graph(edge_vec);
    gettimeofday(&time_end, NULL);
    double build_time = (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;
    printf("build_time=%.3fms\n", build_time);

    // porder.leaf_node_count();
    // porder.select_bignode(0.8);
    EdgeVector new_edge_vec;
    gettimeofday(&time_start, NULL);

    if (order_option == "hybrid") {
        reordered_graph_file_path = graph_name + "_Horder.txt";
        reordered_newid_file_path = graph_name + "_Horder_newID.txt";
        new_edge_vec = porder.hybrid_bfsdeg();
    } else if (order_option == "mloggapa") {
        reordered_graph_file_path = graph_name + "_MLOGGAPAorder.txt";
        reordered_newid_file_path = graph_name + "_MLOGGAPAorder_newID.txt";
        new_edge_vec = porder.mloggapa_order();       
    } else if (order_option == "metis") {
        reordered_graph_file_path = graph_name + "_METISorder.txt";
        reordered_newid_file_path = graph_name + "_METIS_newID.txt";
        new_edge_vec = porder.metis_order();
    } else if (order_option == "slashburn") {
        reordered_graph_file_path = graph_name + "_SBorder.txt";
        reordered_newid_file_path = graph_name + "_SB_newID.txt";
        new_edge_vec = porder.slashburn_order();        
    } else if (order_option == "bfsr") {
        reordered_graph_file_path = graph_name + "_BFSRorder.txt";
        reordered_newid_file_path = graph_name + "_BFSR_newID.txt";
        new_edge_vec = porder.bfsr_order();        
    } else if (order_option == "dfs") {
        reordered_graph_file_path = graph_name + "_DFSorder.txt";
        reordered_newid_file_path = graph_name + "_DFS_newID.txt";
        new_edge_vec = porder.dfs_order();        
    } else {
        order_option = "gro";
        reordered_graph_file_path = graph_name + "_GRO.txt";
        reordered_newid_file_path = graph_name + "_GRO_newID.txt";
        new_edge_vec = porder.greedy_mheap(); 
    }

    printf("order_algo=%s graph_file=%s\n", order_option.c_str(), graph_name.c_str());
    gettimeofday(&time_end, NULL);
    double reorder_time = (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;    
    if (reorder_time > 10000.0) printf("reorder_time=%.3fs\n", reorder_time / 1000.0);
    else printf("reorder_time=%.3fms\n", reorder_time);

    porder.comp_ratio();

    gettimeofday(&time_start, NULL);
    save_graph(reordered_graph_file_path, new_edge_vec);
    save_newid(reordered_newid_file_path, porder.org2newid);
    gettimeofday(&time_end, NULL);
    double write_time = (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;
    printf("write_time=%.3fms\n", write_time);

    return 0;
}


