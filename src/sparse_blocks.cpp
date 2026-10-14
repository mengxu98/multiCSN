#include <Rcpp.h>
#include <algorithm>
#include <limits>

// [[Rcpp::export(rng = false)]]
Rcpp::S4 gather_sparse_block(
    Rcpp::IntegerVector pointers,
    Rcpp::IntegerVector indices,
    Rcpp::NumericVector values,
    Rcpp::IntegerVector selected,
    int n_rows
) {
  const int columns = selected.size();
  size_t total = 0;
  for (int column = 0; column < columns; ++column) {
    const int source = selected[column] - 1;
    total += pointers[source + 1] - pointers[source];
  }
  if (total > static_cast<size_t>(std::numeric_limits<int>::max())) {
    Rcpp::stop("Sparse block exceeds the compressed matrix index limit.");
  }
  Rcpp::IntegerVector output_indices(total);
  Rcpp::NumericVector output_values(total);
  Rcpp::IntegerVector output_pointers(columns + 1);
  int offset = 0;
  for (int column = 0; column < columns; ++column) {
    const int source = selected[column] - 1;
    const int start = pointers[source];
    const int length = pointers[source + 1] - start;
    std::copy_n(indices.begin() + start, length, output_indices.begin() + offset);
    std::copy_n(values.begin() + start, length, output_values.begin() + offset);
    offset += length;
    output_pointers[column + 1] = offset;
  }
  Rcpp::S4 block("dgCMatrix");
  block.slot("i") = output_indices;
  block.slot("p") = output_pointers;
  block.slot("x") = output_values;
  block.slot("Dim") = Rcpp::IntegerVector::create(n_rows, columns);
  return block;
}

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector positive_sparse_row(
    Rcpp::IntegerVector pointers, Rcpp::IntegerVector indices,
    Rcpp::NumericVector values, int row
) {
  std::vector<int> output;
  for (int pos = pointers[row - 1]; pos < pointers[row]; ++pos) {
    if (Rcpp::NumericVector::is_na(values[pos])) {
      output.push_back(NA_INTEGER);
    } else if (values[pos] > 0) {
      output.push_back(indices[pos] + 1);
    }
  }
  return Rcpp::wrap(output);
}

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector sparse_candidate_union(
    Rcpp::IntegerVector domain_p, Rcpp::IntegerVector domain_j,
    Rcpp::NumericVector domain_x, Rcpp::IntegerVector gate_p,
    Rcpp::IntegerVector gate_j, Rcpp::NumericVector gate_x,
    Rcpp::IntegerVector gate_row, int row, int n_regulators
) {
  std::vector<unsigned char> selected(n_regulators, 0);
  for (int pos = domain_p[row - 1]; pos < domain_p[row]; ++pos) {
    if (Rcpp::NumericVector::is_na(domain_x[pos]) || domain_x[pos] == 0) continue;
    const int gate = gate_row[domain_j[pos]];
    if (gate == NA_INTEGER) continue;
    for (int inner = gate_p[gate - 1]; inner < gate_p[gate]; ++inner) {
      if (gate_x[inner] > 0) selected[gate_j[inner]] = 1;
    }
  }
  std::vector<int> output;
  for (int regulator = 0; regulator < n_regulators; ++regulator) {
    if (selected[regulator]) output.push_back(regulator + 1);
  }
  return Rcpp::wrap(output);
}
