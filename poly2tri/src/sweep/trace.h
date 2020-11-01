#pragma once

#include <iostream>

namespace p2t
{
constexpr bool ENABLE_TRACE = false;
}

// From Stack Overflow: https://stackoverflow.com/questions/11832960/short-circuit-operator-output-in-c
#define TRACE_OUT if (!p2t::ENABLE_TRACE) {} else std::cerr
