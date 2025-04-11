#include "query_all.hpp"

template <class Mat>
QueryAll<Mat>::QueryAll() {
    n_intervals = 0;
    n_cpgs = 0;
}

template <class Mat>
QueryAll<Mat>::QueryAll(
    std::vector<std::string>& bedfile_vec,
    std::vector<std::string>& regions,
    const BSType type,
    const int valInd,
    const bool merged,
    const bool sparse,
    const int prealloc,
    const int nthreads
) {

    n_cpgs = 0;
    chr_id = 0;
    resize_count = 0;
    n_samples = bedfile_vec.size();
    n_intervals = regions.size();
    sample_names = Rcpp::CharacterVector(bedfile_vec.size());
    cpg_map = khmap_init();
    is_merged = merged;
    stop_invalid_argument = false;

    stop_invalid_argument = false;
    stop_out_of_range = false;

    bitmat.resize(prealloc, bedfile_vec.size());

    setup_logger("iscream::query_all");

    spdlog::info("Querying {0} regions from {1} bedfiles\n", regions.size(), bedfile_vec.size());
    Progress bar(bedfile_vec.size(), true); 

    spdlog::stopwatch sw;
    #if defined(_OPENMP)
        #pragma omp parallel for num_threads(nthreads)
    #endif
    for (int bedfile_n = 0; bedfile_n < bedfile_vec.size(); bedfile_n++) {
        if (stop_invalid_argument | stop_out_of_range) continue;
        if ( !Progress::check_abort() ) {
            spdlog::debug("Querying {}", bedfile_vec[bedfile_n]);
            MultiRegionQuery cpgs_in_file = query_intervals(bedfile_vec[bedfile_n].c_str(), regions);
            for (RegionQuery cpgs_in_interval : cpgs_in_file) {
                if (stop_invalid_argument | stop_out_of_range) continue;
                #pragma omp critical
                {
                    try {
                        populate_matrix(cpgs_in_interval, bedfile_n, type, valInd);
                    } catch (std::invalid_argument const& ex) {
                        stop_invalid_argument = true;
                        Rprintf("\n%s\n", ex.what());
                    } catch (std::out_of_range const& ex) {
                        Rprintf("\n%s\n", ex.what());
                        stop_out_of_range = true;
                    }
                }
            }
            if (spdlog::get_level() == spdlog::level::info && !stop_invalid_argument && !stop_out_of_range) bar.increment();
        }
    }
    if (stop_invalid_argument) Rcpp::stop("Caught 'std::invalid_argument': parsed data columns are not numeric");
    if (stop_out_of_range) Rcpp::stop("Caught 'std::out_of_range': provided column index is too large");
    spdlog::debug("Made matrix in {} s", sw);
    bar.cleanup();

    int mapsize = kh_size(cpg_map);

    std::vector<int> starts_vec(mapsize);

    seqnames = PROTECT(sf_vector(mapsize));
    sf_vec_data& seq_data = sf_vec_data_ref(seqnames);

    spdlog::debug("Created temporary seqnames, samplenames vectors of size {}", kh_size(cpg_map));

    spdlog::info("Creating metadata vectors");
    sw.reset();
    khint_t iter;
    #if defined(_OPENMP)
        #pragma omp parallel for num_threads(nthreads)
    #endif
    for (iter = 0; iter < kh_end(cpg_map); ++iter) {
        if (kh_exist(cpg_map, iter)) {
            CpG cpg = kh_key(cpg_map, iter);
            int row_idx = kh_val(cpg_map, iter) - 1;
            starts_vec[row_idx] = cpg.start;
            seq_data[row_idx] = sfstring(chr_rev_map[cpg.chr], CE_UTF8);
            std::stringstream cpgid_stream;
            cpgid_stream << chr_rev_map[cpg.chr] << ":" << cpg.start + 1;
        }
    }
    start = Rcpp::wrap(starts_vec);

    spdlog::debug("Setting sample names");
    std::string sample_name;
    for (int i = 0; i < sample_names.size(); i++) {
        std::filesystem::path sample_path = bedfile_vec[i];
        sample_name = sample_path.extension() == ".gz" ? sample_path.stem().stem().string() : sample_path.stem().string();
        sample_names[i] = sample_name;
         spdlog::debug("Got {} as sample name from {}", sample_name, bedfile_vec[i]);
    }
    spdlog::debug("Populated seqnames, samplenames vectors in {} s", sw);

    sw.reset();
    int n_rows = bitmat.n_rows;
    if (bitmat.n_rows > mapsize) {
        int diff_rows = bitmat.n_rows - mapsize;
        spdlog::info("{} loci found - {} extra rows allocated with {} resizes", mapsize, diff_rows, resize_count);
        bitmat.resize(mapsize, bedfile_vec.size());
        spdlog::debug("Corrected matrix size in {} s", sw);
    }

    sw.reset();

    switch(type) {
        case BISCUIT: case BISMARK:
            if (sparse) {
                spdlog::info("Creating sparse matrix");
                Rcpp::S4 bitrmat = Rcpp::wrap(bitmat);
                bitrmat.slot("Dimnames") = Rcpp::List::create(R_NilValue, sample_names);
                assays = Rcpp::List::create(
                    Rcpp::_["M"] = bitrmat,
                    Rcpp::_["Cov"] = Rcpp::clone(bitrmat)
                );
                spdlog::debug("Took {}", sw);
            } else {
                spdlog::info("Creating dense matrix");
                Rcpp::NumericMatrix bitrmat = Rcpp::wrap(bitmat);
                Rcpp::colnames(bitrmat) = sample_names;
                assays = Rcpp::List::create(
                    Rcpp::_["M"] = bitrmat,
                    Rcpp::_["Cov"] = Rcpp::clone(bitrmat)
                );
                spdlog::debug("Took {}", sw);
            }
            break;
        default:
            if (sparse) {
                spdlog::info("Creating sparse matrix");
                Rcpp::S4 bitrmat = Rcpp::wrap(bitmat);
                bitrmat.slot("Dimnames") = Rcpp::List::create(R_NilValue, sample_names);
                assays = Rcpp::List::create(
                    Rcpp::_["M"] = bitrmat
                );
                spdlog::debug("Took {}", sw);
            } else {
                spdlog::info("Creating dense matrix");
                Rcpp::NumericMatrix bitrmat = Rcpp::wrap(bitmat);
                Rcpp::colnames(bitrmat) = sample_names;
                assays = Rcpp::List::create(
                    Rcpp::_["M"] = bitrmat
                );
                spdlog::debug("Took {}", sw);
            }
    }


}
template <class Mat>
void QueryAll<Mat>::populate_matrix(RegionQuery& query, int& col_n, const BSType type, const int valInd) {

    int cpg_count = query.cpgs_in_interval.size();
    std::vector<BedRecord> lines;
    std::vector<CpG> ids;
    khmap_m_resize(cpg_map, cpg_count);
    BedRecord parsed_bedline;
    for (std::string cpg_string : query.cpgs_in_interval) {
        switch(type) {
            case BISCUIT:
                parsed_bedline = parseBiscuitRecord(cpg_string);
                break;
            case BISMARK:
                parsed_bedline = parseCovRecord(cpg_string);
                break;
            default:
                parsed_bedline = parseBedRecord(cpg_string, valInd - 1);
        }

        lines.push_back(parsed_bedline);
        spdlog::trace(
            "Parsed {} into chr: {}, start: {}, end: {}",
            cpg_string,
            parsed_bedline.chr,
            parsed_bedline.start,
            parsed_bedline.end,
            parsed_bedline.data[1],
            parsed_bedline.data[2]
        );
        if (!chr_map.count(parsed_bedline.chr)) {
            chr_map.insert({parsed_bedline.chr, ++chr_id});
            chr_rev_map.insert({chr_id, parsed_bedline.chr});
        }

        CpG cpg = CpG{chr_map[parsed_bedline.chr], parsed_bedline.start};

        ids.push_back(cpg);

        khint_t insert_b;
        int absent;
        insert_b = khmap_put(cpg_map, cpg, &absent);
        if (absent) {
            n_cpgs++;
            kh_val(cpg_map, insert_b) = n_cpgs;
        }
    }

    int mapsize = kh_size(cpg_map);
    int cur_nrow = bitmat.n_rows;
    if (cur_nrow < mapsize) {
        resize_mat(cur_nrow, mapsize);
    }

    spdlog::debug("Inserting {} CpGs into matrix", cpg_count);
    for (size_t i = 0; i < lines.size(); i++) {
        khint_t retrieve_b = khmap_get(cpg_map, ids[i]);
        int idx = kh_val(cpg_map, retrieve_b);
        if (lines[i].data[2] == -1) {
            bitmat(idx - 1, col_n) = lines[i].data[0];
        } else {
            bitmat(idx - 1, col_n) = bitpack(lines[i].data[0], lines[i].data[1]);
        }
    }
}

template <class Mat>
void QueryAll<Mat>::resize_mat(int cur_nrow, int mapsize) {
    int diff = mapsize - bitmat.n_rows;
    spdlog::debug("Need {} more rows", diff);
    int extra_rows = diff;
    if (diff < 10) {
        extra_rows = diff * 10000;
    } else if (diff < 10000) {
        extra_rows = diff * 1000;
    } else if (diff < 50000) {
        extra_rows = diff * 500;
    } else {
        extra_rows = diff * 100;
    }
    bitmat.resize(cur_nrow + extra_rows, bitmat.n_cols);
    resize_count++;
    spdlog::debug("Added {} rows to existing {}", extra_rows, cur_nrow);
}


template <class Mat>
int QueryAll<Mat>::bitpack(const float beta_val, const float cov_val) {
    uint16_t cov = cov_val > INT16_MAX ? INT16_MAX : cov_val;
    uint16_t betap = std::round(beta_val * 100);
    betap = betap > INT16_MAX ? INT16_MAX : betap;
    return(betap + (cov << 16));
}

//' Query all CpG info into M and coverage matrices
//' @param bedfiles A vector of bedfiles
//' @param regions A vector of regions
//' @param aligner The aligner used to make the WGBS BED files, only for
//' `make_bsseq_mat`
//' @param valInd The index of the data column needed for the matrix, for `make_mat`
//' @param merged Whether the input strands have been merged/collapsed
//' @param prealloc The number of rows to initialize the matrices with
//' @param nthreads Set number of threads to use overriding the
//' `"iscream.threads"` option. See `?set_threads` for more information.
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List Cpp_query_all(
    std::vector<std::string>& bedfiles,
    std::vector<std::string>& regions,
    const std::string aligner,
    const int valInd,
    const bool merged,
    const bool sparse,
    const int prealloc,
    const int nthreads) {

    BSType type;
    if (aligner == "biscuit") {
        type = BISCUIT;
    } else if (aligner == "bismark" || aligner == "bsbolt") {
        type = BISMARK;
    } else {
        type = GENERAL;
        if (sparse) {
            QueryAll query = QueryAll<arma::sp_fmat>(bedfiles, regions, type, valInd, merged, sparse, prealloc, nthreads);
            return query.ret();
        } else {
            QueryAll query = QueryAll<arma::fmat>(bedfiles, regions, type, valInd, merged, sparse, prealloc, nthreads);
            return query.ret();
        }
    }

    if (sparse) {
        QueryAll query = QueryAll<arma::sp_umat>(bedfiles, regions, type, valInd, merged, sparse, prealloc, nthreads);
        return query.wrap();
    } else {
        QueryAll query = QueryAll<arma::umat>(bedfiles, regions, type, valInd, merged, sparse, prealloc, nthreads);
        return query.wrap();
    }
}
