#pragma once

#include <fstream>
#include <vector>

#include "domain.h"
#include "species.h"

namespace Output {
void fields(Domain &domain, std::vector<Species> &species);
void fields(Domain &domain);
void screenOutput(Domain &domain, std::vector<Species> &species);
void diagOutput(Domain &domain, std::vector<Species> &species);
void Display(Domain &domain, Field<Vec3d> field);
}  // namespace Output
