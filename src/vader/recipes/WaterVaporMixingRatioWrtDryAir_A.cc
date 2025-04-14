/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <iostream>
#include <vector>

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/for_each.h"
#include "atlas/util/Metadata.h"
#include "oops/util/Logger.h"
#include "vader/recipes/WaterVaporMixingRatioWrtDryAir.h"

namespace vader
{
// ------------------------------------------------------------------------------------------------

// Static attribute initialization
const char WaterVaporMixingRatioWrtDryAir_A::Name[] = "WaterVaporMixingRatioWrtDryAir_A";
const oops::Variables WaterVaporMixingRatioWrtDryAir_A::Ingredients{
                      {"water_vapor_mixing_ratio_wrt_moist_air"}};

// Register the maker
static RecipeMaker<WaterVaporMixingRatioWrtDryAir_A> makerWaterVaporMixingRatioWrtDryAir_(
                   WaterVaporMixingRatioWrtDryAir_A::Name);

WaterVaporMixingRatioWrtDryAir_A::WaterVaporMixingRatioWrtDryAir_A(const Parameters_ & params,
                                       const VaderConfigVars & configVariables) :
    configVariables_{configVariables}
{
    oops::Log::trace()
          << "WaterVaporMixingRatioWrtDryAir_A::WaterVaporMixingRatioWrtDryAir_A(params)"
          << std::endl;
}

std::string WaterVaporMixingRatioWrtDryAir_A::name() const
{
    return WaterVaporMixingRatioWrtDryAir_A::Name;
}

oops::Variable WaterVaporMixingRatioWrtDryAir_A::product() const
{
    return oops::Variable{"water_vapor_mixing_ratio_wrt_dry_air"};
}

oops::Variables WaterVaporMixingRatioWrtDryAir_A::ingredients() const
{
    return WaterVaporMixingRatioWrtDryAir_A::Ingredients;
}

size_t WaterVaporMixingRatioWrtDryAir_A::productLevels(const atlas::FieldSet & afieldset) const
{
    return afieldset.field("water_vapor_mixing_ratio_wrt_moist_air").shape(1);
}

atlas::FunctionSpace WaterVaporMixingRatioWrtDryAir_A::productFunctionSpace
                                              (const atlas::FieldSet & afieldset) const
{
    return afieldset.field("water_vapor_mixing_ratio_wrt_moist_air").functionspace();
}

void WaterVaporMixingRatioWrtDryAir_A::executeNL(atlas::FieldSet & afieldset)
{
    oops::Log::trace()
          << "entering WaterVaporMixingRatioWrtDryAir_A::executeNL function"
          << std::endl;

    atlas::field::for_each_value(afieldset["water_vapor_mixing_ratio_wrt_moist_air"],
                                 afieldset["water_vapor_mixing_ratio_wrt_dry_air"],
                                 [&](const double q, double& mixr) {
        mixr = q / (1. - q);
    });

    oops::Log::trace()
          << "leaving WaterVaporMixingRatioWrtDryAir_A::executeNL function"
          << std::endl;
}

}  // namespace vader
