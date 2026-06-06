#pragma once

#include "TypeTraits.hpp"

#include <muParser.h>
#include <string>
#include <iostream>
#include <cmath>

using namespace Traits;

/**
 * @class muParserWrapper
 * @brief This class provides a simple interface for a function wrapper, using
 * the muParser library
 */
class muParserWrapper
{
    private:
    mu::Parser parser;
    Real var_x;
    Real var_y;

    public:
    muParserWrapper(const std::string &expression)
    {
        try{
            parser.DefineVar("x", &var_x);
            parser.DefineVar("y", &var_y);
            parser.DefineConst("pi", M_PI);
            parser.SetExpr(expression);
        } catch (mu::Parser::exception_type &e){
            std::cerr << "muParser Error: " << e.GetMsg() << std::endl;
            exit(1);
        }
    }
    
    Real operator() (Real x, Real y)
    {
        var_x = x;
        var_y = y;
        return parser.Eval();
    }
};
