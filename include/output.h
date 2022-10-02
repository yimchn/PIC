#pragma once

#include <fstream>
#include <iostream>
#include <vector>

#include "domain.h"
#include "species.h"

namespace Output {
void fields(Domain &domain, std::vector<Species> &species);
void fields(Domain &domain);
void fields(Domain &domain, int step);
void screenOutput(Domain &domain, std::vector<Species> &species);
void diagOutput(Domain &domain, std::vector<Species> &species);
void Display(Domain &domain, Field<Vec3d> field);
void ProgressBar(double cur, double total);
}  // namespace Output
