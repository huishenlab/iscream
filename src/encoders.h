#ifndef ENCODERS_H
#define ENCODERS_H
#include <Rcpp.h>

int encode_meth(uint32_t meth, uint32_t unmeth);
int get_m(uint32_t encoded);
int get_u(uint32_t encoded);

#endif /* ifndef ENCODERS_H */
