#pragma once

#include "rpc.hpp"

double // return fit RMS>0 (pixels) if OK, else <0
fit_dlt(RPC_NS::RPC<double>& rpc,
        double* coeffs11);

