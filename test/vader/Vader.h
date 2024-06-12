/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <netcdf.h>

#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/meshgenerator.h"

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Variables.h"
#include "oops/runs/Test.h"
#include "oops/util/Logger.h"
#include "test/TestEnvironment.h"

#include "vader/vader.h"
#include "Utils.h"

namespace vader {
namespace test {

/// \brief Top-level options taken by the Vader test.
class VaderTestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VaderTestParameters, Parameters)

 public:
  oops::RequiredParameter<vader::VaderParameters> vader{
        "vader", "vader parameters", this};
  oops::RequiredParameter<eckit::LocalConfiguration> grid{
        "trajectory grid", "trajectory grid description", this};
  oops::RequiredParameter<std::string> filename{"trajectory filename",
        "filename of existing netcdf file with trajectory ingredients", this};
  oops::RequiredParameter<std::vector<std::string>> ingredients{
        "ingredients",
        "ingredient variables that must be available in the trajectory file",
        this};
  oops::RequiredParameter<std::vector<std::string>> products{"products",
        "target variables", this};
  oops::RequiredParameter<double> tolerance{"adjoint test tolerance",
        "adjoint test tolerance", this};
};

// -----------------------------------------------------------------------------
/// \brief Tests adjoint of the recipe.
void testVaderAdjoint() {
  VaderTestParameters params;
  params.validateAndDeserialize(::test::TestEnvironment::config());

  eckit::LocalConfiguration modelVarsConfig;
  modelVarsConfig.set("gas_constant_of_dry_air", 2.8705e2);
  eckit::LocalConfiguration vaderConfig;
  vaderConfig.set(vader::configModelVarsKey, modelVarsConfig);

  vader::Vader vader(params.vader, vaderConfig);
  const oops::Variables ingredientVars(params.ingredients);
  const oops::Variables productVars(params.products);
  oops::Log::info() << "Testing vader" << std::endl;
  oops::Log::info() << " Ingredients (dx): " << ingredientVars << std::endl;
  oops::Log::info() << " Products (dy): " << productVars << std::endl;

  // set up grid and functionspace based on description in the yaml
  const atlas::StructuredGrid grid(params.grid.value());
  const atlas::functionspace::StructuredColumns fs(grid);

  // Open NetCDF file and read all the ingredients for the trajectory
  atlas::FieldSet traj;
  int ncid, retval;
  const std::string & filename = params.filename;
  oops::Log::info() << "Reading trajectory from file: " << filename << std::endl;
  if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid))) ERR(retval);
  std::vector<size_t> ingredientLevels(ingredientVars.size(), 0);
  for (size_t jvar = 0; jvar < ingredientVars.size(); ++jvar) {
    addFieldFromFile(traj, ingredientVars[jvar].name(), fs, grid,
                     ingredientLevels[jvar], ncid);
  }
  // run NL to set trajectory
  oops::Variables vars = productVars;
  vader.changeVarTraj(traj, vars);

  // The test should be setup to do variable change for all variables, exit
  // if it's not
  if (vars.size() > 0) {
    oops::Log::info() << "Not all variables can be converted; no recipes found for "
                      << vars << std::endl;
  }
  EXPECT_EQUAL(vars.size(), 0);

  // Testing whether (dx, K^T dy) == (K dx, dy)
  // Allocating dxin to contain randomized dx (for the ingredient variables)
  atlas::FieldSet dxin;
  for (size_t jvar = 0; jvar < ingredientVars.size(); ++jvar) {
    addRandomField(dxin, ingredientVars[jvar].name(), fs, ingredientLevels[jvar]);
  }
  // fill the product variable with Kdx
  vars = productVars;
  vader.changeVarTL(dxin, vars);

  // Allocating dxout to contain randomized dy (for the product variable).
  // After calling recipe->executeAD, dy (the product variable) will be
  // zeroed out, so the copy is needed.
  atlas::FieldSet dxout, dy;
  for (const auto & productVar : productVars) {
    const std::string & name = productVar.name();
    addRandomField(dxout, name, fs, traj.field(name).shape(1));
    atlas::Field field(name, dxout.field(name).datatype(),
                  dxout.field(name).shape());
    auto from_view = atlas::array::make_view<double, 2>(dxout.field(name));
    auto to_view = atlas::array::make_view<double, 2>(field);
    to_view.assign(from_view);
    dy.add(field);
  }
  // fill the ingredients variables with K^T dxout
  vars = productVars;
  vader.changeVarAD(dxout, vars);

  double zz1 = 0;
  // Compute (dx, K^T dy)
  for (const auto & ingredientVar : ingredientVars) {
    zz1 += dotProduct(dxin[ingredientVar.name()], dxout[ingredientVar.name()]);
  }
  // Compute (Kdx, dy)
  double zz2 = 0;
  for (const auto & productVar : productVars) {
    zz2 += dotProduct(dxin[productVar.name()], dy[productVar.name()]);
  }
  oops::Log::info() << "<dx,KTdy>=" << zz1 << std::endl;
  oops::Log::info() << "<Kdx,dy>=" << zz2 << std::endl;
  oops::Log::info() << "<dx,KTdy>-<Kdx,dy>/<dx,KTdy>="
                    << (zz1-zz2)/zz1 << std::endl;
  oops::Log::info() << "<dx,KTdy>-<Kdx,dy>/<Kdx,dy>="
                    << (zz1-zz2)/zz2 << std::endl;
  const double tol = params.tolerance;
  EXPECT(oops::is_close(zz1, zz2, tol));
}

// -----------------------------------------------------------------------------

class Vader : public oops::Test {
 public:
  Vader() {}
  virtual ~Vader() {}
 private:
  std::string testid() const override {return "vader::test::Vader";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("vader/Vader/testVaderAdjoint")
      { testVaderAdjoint(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace vader
