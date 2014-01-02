function result = TDimensionlessQuantity(dimensionalQuantity, lDimensionalUnit, tDimensionalUnit)
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here
    result = dimensionalQuantity * tDimensionalUnit / lDimensionalUnit;