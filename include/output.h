#pragma once

#include <fstream>
#include <iostream>
#include <vector>

#include "domain.h"
#include "species.h"

namespace Output {
void fields(Domain &dm, std::vector<Species> &species);
void fields(Domain &dm);
void fields(Domain &dm, int step);
void B(Domain &dm);
void E(Domain &dm);
void J(Domain &dm);
void screenOutput(Domain &dm, std::vector<Species> &species);
void diagOutput(Domain &dm, std::vector<Species> &species);
void Display(Domain &dm, Field<Vec3d> field);
void ProgressBar(Domain dm);
}  // namespace Output
