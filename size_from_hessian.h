#ifndef SIZE_FROM_HESSIAN_H
#define SIZE_FROM_HESSIAN_H

double* size_from_hessian(
    unsigned nverts,
    unsigned nhess_comps,
    double const* hessians,
    double const* sol_comp_weights,
    double min_h,
    double max_h);

#endif
