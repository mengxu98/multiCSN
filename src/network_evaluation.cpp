#include <Rcpp.h>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <vector>
#include <cmath>

using namespace Rcpp;

// Length-prefix the first endpoint so arbitrary punctuation cannot merge edges.
static std::string metric_edge_key(const std::string &reg, const std::string &tar)
{
    return std::to_string(reg.size()) + ":" + reg + tar;
}

// [[Rcpp::export]]
DataFrame prepare_calculate_metrics(DataFrame network_table, DataFrame ground_truth)
{
    CharacterVector network_reg = as<CharacterVector>(network_table["regulator"]);
    CharacterVector network_tar = as<CharacterVector>(network_table["target"]);
    NumericVector network_weight = abs(as<NumericVector>(network_table["weight"]));

    CharacterVector truth_reg = as<CharacterVector>(ground_truth["regulator"]);
    CharacterVector truth_tar = as<CharacterVector>(ground_truth["target"]);

    std::unordered_set<std::string> gt_gene_set;
    std::vector<std::string> gt_genes;
    std::unordered_set<std::string> truth_map;
    for (R_xlen_t i = 0; i < truth_reg.length(); i++)
    {
        if (truth_reg[i] == NA_STRING || truth_tar[i] == NA_STRING)
            continue;
        if (truth_reg[i] == truth_tar[i])
            continue;

        std::string reg = as<std::string>(truth_reg[i]);
        std::string tar = as<std::string>(truth_tar[i]);
        if (!gt_gene_set.count(reg))
        {
            gt_gene_set.insert(reg);
            gt_genes.push_back(reg);
        }
        if (!gt_gene_set.count(tar))
        {
            gt_gene_set.insert(tar);
            gt_genes.push_back(tar);
        }

        std::string key = metric_edge_key(reg, tar);
        truth_map.insert(key);
    }

    if (gt_genes.size() < 2)
    {
        return DataFrame::create(
            _["regulator"] = CharacterVector(0),
            _["target"] = CharacterVector(0),
            _["weight"] = NumericVector(0),
            _["label"] = IntegerVector(0));
    }

    std::sort(gt_genes.begin(), gt_genes.end());

    std::unordered_map<std::string, double> pred_weights;
    for (R_xlen_t i = 0; i < network_reg.length(); i++)
    {
        if (network_reg[i] == NA_STRING || network_tar[i] == NA_STRING)
            continue;

        std::string reg = as<std::string>(network_reg[i]);
        std::string tar = as<std::string>(network_tar[i]);
        if (reg == tar)
            continue;
        if (!gt_gene_set.count(reg) || !gt_gene_set.count(tar))
            continue;
        if (!std::isfinite(network_weight[i]))
            continue;

        std::string key = metric_edge_key(reg, tar);
        double weight = network_weight[i];
        auto it = pred_weights.find(key);
        if (it == pred_weights.end() || weight > it->second)
        {
            pred_weights[key] = weight;
        }
    }

    struct EdgeRow
    {
        std::string regulator;
        std::string target;
        double weight;
        int label;
    };
    std::vector<EdgeRow> rows;
    rows.reserve(gt_genes.size() * (gt_genes.size() - 1));

    for (size_t i = 0; i < gt_genes.size(); i++)
    {
        for (size_t j = 0; j < gt_genes.size(); j++)
        {
            if (i == j)
                continue;

            std::string key = metric_edge_key(gt_genes[i], gt_genes[j]);
            double weight = 0.0;
            auto pred_it = pred_weights.find(key);
            if (pred_it != pred_weights.end())
            {
                weight = pred_it->second;
            }
            rows.push_back(
                EdgeRow{
                    gt_genes[i],
                    gt_genes[j],
                    weight,
                    truth_map.count(key) > 0 ? 1 : 0});
        }
    }

    std::stable_sort(
        rows.begin(),
        rows.end(),
        [](const EdgeRow &a, const EdgeRow &b)
        {
            return a.weight > b.weight;
        });

    R_xlen_t n = static_cast<R_xlen_t>(rows.size());
    CharacterVector out_reg(n);
    CharacterVector out_tar(n);
    NumericVector out_weight(n);
    IntegerVector out_label(n);

    for (R_xlen_t i = 0; i < n; i++)
    {
        out_reg[i] = rows[i].regulator;
        out_tar[i] = rows[i].target;
        out_weight[i] = rows[i].weight;
        out_label[i] = rows[i].label;
    }

    DataFrame result = DataFrame::create(
        _["regulator"] = out_reg,
        _["target"] = out_tar,
        _["weight"] = out_weight,
        _["label"] = out_label);
    result.attr("row.names") = IntegerVector::create(NA_INTEGER, -n);

    return result;
}
