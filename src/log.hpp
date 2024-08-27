// [[Rcpp::depends(RcppSpdlog)]]
#include <RcppSpdlog>

#ifndef LOG_H
#define LOG_H

# if defined __cplusplus

#include <string>

void set_log_level(const std::string& name);
std::string get_log_level();
void setup_logger(std::string logname);

#endif /* __cplusplus */

#endif /* ifndef LOG_H */
