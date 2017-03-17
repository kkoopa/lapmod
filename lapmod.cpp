/*
 * An implementation of LAPMOD based on the original Pascal code in
 * A. Volgenant, Linear and Semi Assignment Problems: a core oriented approach,
 * Computers & Operations Research 23 917-932, 1996.
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
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
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
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

class matrix {
public:
  matrix() : height_{}, width_{}, data_{nullptr} {}
  matrix(const matrix &other) = delete;
  matrix(matrix &&other)
      : height_{other.height_}, width_{other.width_}, data_{other.data_} {
    other.data_ = nullptr;
  }
  matrix(int n) : height_{n}, width_{n}, data_{new int[height_ * width_]} {}
  matrix(int height, int width)
      : height_{height}, width_{width}, data_{new int[height_ * width_]} {}
  matrix(std::initializer_list<std::initializer_list<int>> l)
      : height_{static_cast<int>(l.size())},
        width_{static_cast<int>(l.begin()->size())},
        data_{new int[height_ * width_]} {
    if (l.size() > std::numeric_limits<int>::max() ||
        l.begin()->size() > std::numeric_limits<int>::max()) {
      throw std::range_error("Matrix dimensions too large.");
    }
    std::accumulate(l.begin(), l.end(), data_, [](auto *l, const auto &r) {
      return std::copy(r.begin(), r.end(), l);
    });
  }

  matrix &operator=(const matrix &other) = delete;

  matrix &operator=(matrix &&other) {
    height_ = other.height_;
    width_ = other.width_;
    data_ = other.data_;
    other.data_ = nullptr;
    return *this;
  }

  ~matrix() { delete[] data_; }

  int height() const { return height_; }
  int width() const { return width_; }
  int *operator[](std::ptrdiff_t index) { return &data_[index * width_]; }
  const int *operator[](std::ptrdiff_t index) const {
    return &data_[index * width_];
  }

private:
  int height_;
  int width_;
  int *data_;
};

const int inf = std::numeric_limits<int>::max();
const int pp = 39;
const int maxpp = 61;

class problem {
public:
  problem()
      : cost_matrix_{nullptr}, cc_(), kk_(), d_{}, unused_{}, lab_{}, number_{},
        todo_{}, u_{}, v_{}, x_{}, y_{}, ok_{}, z_{-1} {}
  problem(const matrix *cost_matrix)
      : cost_matrix_(cost_matrix), cc_(cost_matrix_->height(), maxpp),
        kk_(cost_matrix_->width(), maxpp), d_{new int[cost_matrix_->height()]},
        unused_{new int[cost_matrix_->height()]},
        lab_{new int[cost_matrix_->height()]},
        number_{new int[cost_matrix_->height()]},
        todo_{new int[cost_matrix_->height()]},
        u_{new int[cost_matrix_->height()]}, v_{new int[cost_matrix_->width()]},
        x_{new int[cost_matrix_->height()]}, y_{new int[cost_matrix_->width()]},
        ok_{new int[cost_matrix_->height()]}, z_{-1} {
    std::fill(v_, v_ + cost_matrix_->width(), inf);
  }
  ~problem() {
    delete[] d_;
    delete[] unused_;
    delete[] lab_;
    delete[] number_;
    delete[] todo_;
    delete[] u_;
    delete[] v_;
    delete[] x_;
    delete[] y_;
    delete[] ok_;
  }

  void selpp_cr() {
    auto t = -1;
    auto end = std::min(cost_matrix_->width(), pp);

    for (auto i = 0; i != cost_matrix_->height(); ++i) {
      auto s = 0ul;
      for (auto j = 0; j != end; ++j) {
        cc_[i][j] = (*cost_matrix_)[i][j];
        kk_[i][j] = j + 1;
        s += (*cost_matrix_)[i][j];
      }

      auto cr = s / end;

      for (auto j = pp; j < cost_matrix_->width(); ++j) {
        if ((*cost_matrix_)[i][j] < cr) {
          auto h = 0;
          do {
            t = t >= pp - 1 ? 0 : t + 1;
            h = cc_[i][t];
          } while (h < cr);
          cc_[i][t] = (*cost_matrix_)[i][j];
          kk_[i][t] = j + 1;
          s = s - h + (*cost_matrix_)[i][j];
          cr = s / pp;
        }
      }

      x_[i] = 0;
      auto diag = (*cost_matrix_)[i][i];

      for (t = 0; t != end; ++t) {
        const auto j = kk_[i][t] - 1;
        if (cc_[i][t] < v_[j]) {
          v_[j] = cc_[i][t];
          y_[j] = i + 1;
        }
        if (j == i) {
          diag = inf;
        }
      }

      if (diag < inf) {
        t = end;
        cc_[i][t] = diag;
        kk_[i][t] = i + 1;
        if (diag < v_[i]) {
          v_[i] = diag;
          y_[i] = i + 1;
        }
      } else {
        t = end - 1;
      }
      number_[i] = t + 1;
    }

    for (auto j = cost_matrix_->width(); j > 0; --j) {
      const auto i = y_[j - 1] - 1;
      if (x_[i] == 0) {
        x_[i] = j;
      } else {
        x_[i] = -std::abs(x_[i]);
        y_[j - 1] = 0;
      }
    }
  }

  int transfer() {
    auto l = -1;

    for (auto i = 0; i != cost_matrix_->height(); ++i) {
      if (x_[i] < 0) {
        x_[i] = -x_[i];
      } else if (x_[i] == 0) {
        unused_[++l] = i;
      } else {
        auto min = inf;
        const auto j1 = x_[i] - 1;
        for (auto t = 0; t != number_[i]; ++t) {
          const auto j = kk_[i][t] - 1;
          if (j != j1 && cc_[i][t] - v_[j] < min) {
            min = cc_[i][t] - v_[j];
          }
        }
        auto t = 0;
        for (; kk_[i][t] != j1 + 1; ++t)
          ;
        v_[j1] = cc_[i][t] - min;
      }
    }

    return l;
  }

  int arr(int l) {
    for (auto cnt = 0; cnt != 2; ++cnt) {
      auto h = 0;
      const auto l0 = l;
      l = -1;

      while (h <= l0) {
        const auto i = unused_[h++];
        auto v0 = inf, vj = inf;
        int j0, j1;

        for (auto t = 0; t != number_[i]; ++t) {
          const auto j = kk_[i][t] - 1, dj = cc_[i][t] - v_[j];

          if (dj < vj) {
            if (dj >= v0) {
              vj = dj;
              j1 = j;
            } else {
              vj = v0;
              v0 = dj;
              j1 = j0;
              j0 = j;
            }
          }
        }

        auto i0 = y_[j0];

        if (v0 < vj) {
          v_[j0] = v_[j0] - vj + v0;
        } else if (i0 > 0) {
          j0 = j1;
          i0 = y_[j0];
        }

        x_[i] = j0 + 1;
        y_[j0] = i + 1;

        if (i0 > 0) {
          if (v0 < vj) {
            unused_[--h] = i0 - 1;
          } else {
            unused_[++l] = i0 - 1;
          }
        }
      }
    }
    return l;
  }

  void augmentation(int l) {
    int i, j;
    do {
      for (auto l0 = l + 1, l = 0; l < l0; ++l) {
        std::fill_n(d_, cost_matrix_->height(), inf);
        std::fill_n(ok_, cost_matrix_->height(), false);
        auto min = inf;
        auto i0 = unused_[l];
        int td1;

        for (auto t = 0; t != number_[i0]; ++t) {
          j = kk_[i0][t] - 1;
          const auto dj = cc_[i0][t] - v_[j];
          d_[j] = dj;
          lab_[j] = i0;

          if (dj <= min) {
            if (dj < min) {
              td1 = -1;
              min = dj;
            }
            todo_[++td1] = j;
          }
        }

        auto last = cost_matrix_->height();
        auto td2 = last - 1;

        for (auto h = 0; h != td1 + 1; ++h) {
          j = todo_[h];
          if (y_[j] == 0) {
            goto augment_2;
          }
          ok_[j] = true;
        }

        do {
          const auto j0 = todo_[td1--], i = y_[j0] - 1;
          todo_[td2--] = j0;

          auto t = 0;
          for (; kk_[i][t] != j0 + 1; ++t)
            ;
          const auto h = cc_[i][t] - v_[j0] - min;
          for (auto t = 0; t != number_[i]; ++t) {
            j = kk_[i][t] - 1;
            if (!ok_[j]) {
              const auto vj = cc_[i][t] - v_[j] - h;
              if (vj < d_[j]) {
                d_[j] = vj;
                lab_[j] = i;
                if (vj == min) {
                  if (y_[j] == 0) {
                    goto augment_1;
                  }
                  todo_[++td1] = j;
                  ok_[j] = true;
                }
              }
            }
          }

          if (td1 == -1) {
            min = inf - 1;
            last = td2 + 1;
            for (j = 0; j != cost_matrix_->width(); ++j) {
              if (d_[j] <= min && !ok_[j]) {
                if (d_[j] < min) {
                  td1 = -1;
                  min = d_[j];
                }
                todo_[++td1] = j;
              }
            }
            for (auto h = 0; h != td1 + 1; ++h) {
              j = todo_[h];
              if (y_[j] == 0) {
                goto augment_1;
              }
              ok_[j] = true;
            }
          }
        } while (true);
      augment_1:
        for (auto k = last; k < cost_matrix_->width(); ++k) {
          const auto j0 = todo_[k];
          v_[j0] = v_[j0] + d_[j0] - min;
        }
      augment_2:
        do {
          i = lab_[j];
          y_[j] = i + 1;
          auto k = j + 1;
          j = x_[i] - 1;
          x_[i] = k;
        } while (i != i0);
      }

      for (i = 0; i != cost_matrix_->height(); ++i) {
        j = x_[i] - 1;
        auto t = 0;
        for (; kk_[i][t] != j + 1; ++t)
          ;
        u_[i] = cc_[i][t] - v_[j];
      }

      l = optcheck();
    } while (l >= 0);
  }

  int optcheck() {
    auto l = -1;
    auto newfree = false;
    for (auto i = 0; i != cost_matrix_->height(); ++i) {
      newfree = false;

      for (auto j = 0; j != cost_matrix_->width(); ++j) {
        if ((*cost_matrix_)[i][j] < u_[i] + v_[j]) {
          number_[i] = number_[i] < cost_matrix_->height() ? number_[i] + 1 : 0;
          newfree = true;
          cc_[i][number_[i]] = (*cost_matrix_)[i][j];
          kk_[i][number_[i]] = j + 1;
        }
      }

      if (newfree) {
        y_[x_[i] - 1] = 0;
        x_[i] = 0;
        unused_[++l] = i;
      }
    }

    return l;
  }

  const int *solution() const { return x_; }
  int size() const { return cost_matrix_->height(); }
  int value() const {
    return z_ < 0
               ? z_ = std::accumulate(
                     v_, v_ + cost_matrix_->width(),
                     std::accumulate(u_, u_ + cost_matrix_->height(), 0))
               : z_;
  }

private:
  const matrix *cost_matrix_;
  matrix cc_;
  matrix kk_;
  int *d_;
  int *unused_;
  int *lab_;
  int *number_;
  int *todo_;
  int *u_;
  int *v_;
  int *x_;
  int *y_;
  int *ok_;
  mutable int z_;
};

matrix read_data(std::string path) {
  std::ifstream ifs(path);
  int n;
  if (ifs >> n) {
    matrix cc(n);
    for (int c, x = 0; ifs >> c; ++x) {
      const auto i = x / n;
      const auto j = x % n;
      cc[i][j] = c;
    }
    return cc;
  } else {
    throw std::logic_error("Bad data file.");
  }
}

int main() {
  int answer[] = {305, 475, 626, 804, 991, 1176, 1362, 1552};
  for (auto i = 0; i < sizeof(answer) / sizeof(answer[0]); ++i) {
    auto m = read_data(std::string("problems/assign") +
                       std::to_string(100 * (i + 1)) + std::string(".txt"));
    problem p(&m);
    p.selpp_cr();
    p.augmentation(p.arr(p.transfer()));

    auto begin = p.solution();
    std::cout << '[' << *begin++ - 1;
    for (auto end = begin + p.size() - 1; begin != end; ++begin) {
      std::cout << ", " << *begin - 1;
    }
    std::cout << "]\n";

    std::cout << p.value() << '\n';

    if (p.value() != answer[i]) {
      throw std::logic_error("The solver is broken.");
    }
  }
  return 0;
}

