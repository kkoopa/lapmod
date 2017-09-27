/*
 * A test for an implementation of LAPMOD.
 *
 * The test problems are taken from J. E. Beasley,
 * Linear programming on Cray supercomputers
 * Journal of the Operational Research Society 41 (1990) 133-139.
 *
 * MIT License
 *
 * Copyright (c) 2017 Benjamin Byholm
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

#include <algorithm>
#include <chrono>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "lapmod.h"

lapmod::matrix read_data(const std::string &path) {
  std::ifstream ifs(path);
  auto n = 0;
  if (ifs >> n) {
    lapmod::matrix cc(n);
    for (auto c = 0, i = 0, j = 0; ifs >> c; ++j) {
      cc[i][j] = c;
      if (j == n - 1) {
        j = -1;
        ++i;
      }
    }
    return cc;
  } else {
    throw std::logic_error("Bad data file.");
  }
}

int main() {
  static const long answer[] = {305, 475, 626, 804, 991, 1176, 1362, 1552};
  for (auto i = 0; i < static_cast<int>(sizeof(answer) / sizeof(answer[0]));
       ++i) {
    const auto m =
        read_data(std::string("problems/assign") +
                  std::to_string(100 * (i + 1)) + std::string(".txt"));

    auto start = std::chrono::high_resolution_clock::now();

    const lapmod::problem p(&m);
    const auto sol = p.solve();

    std::chrono::duration<double> elapsed_seconds =
        std::chrono::high_resolution_clock::now() - start;

    std::cout << "Solution: " << sol << "\nCost: " << sol.value()
              << "\nDuration: " << elapsed_seconds.count() << " s\n";

    if (sol.value() != answer[i]) {
      throw std::logic_error("The solver is broken.");
    }
  }
  return 0;
}

