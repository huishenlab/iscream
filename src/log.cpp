#include "log.hpp"

// [[Rcpp::export]]
void setup_logger(std::string logname = "iscream") {
    auto sp = spdlog::get(logname);
    if (sp == nullptr) sp = spdlog::r_sink_mt(logname);     // or create new one if needed
    spdlog::set_default_logger(sp);                         // and set as default
    spdlog::set_pattern("[%H:%M:%S.%f] [%n] [%l] %v");
}

//' spdlog Logging Lever Setter
//'
//' A helper function to turn a logging level given as string
//' into the current logging level
//'
//' @param name A string with the logging level. Value understood are,
//' in decreasing verbosity \sQuote{trace}, \sQuote{debug}, \sQuote{info},
//' \sQuote{warning}, \sQuote{error}, \sQuote{critical}, and \sQuote{off}.
//' Unrecognised names are equivalent to \sQuote{off}.
//' @return Nothing is returned.
//' @keywords internal
// [[Rcpp::export]]
void Cpp_set_log_level(const std::string &name) {
    spdlog::set_level(spdlog::level::from_str(name));
}

//' @rdname set_log_level
//' @order 2
//' @export
// [[Rcpp::export]]
std::string get_log_level() {
    spdlog::level::level_enum lvl = spdlog::get_level();
    switch (lvl) {
        case spdlog::level::trace:
            return "trace";
        case spdlog::level::info:
            return "info";
        case spdlog::level::debug:
            return "debug";
        case spdlog::level::warn:
            return "warn";
        case spdlog::level::err:
            return "error";
        default:
            return "off";
    }
}
