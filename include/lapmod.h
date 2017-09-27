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
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace lapmod {
class matrix {
public:
  matrix() noexcept : height_{}, width_{}, data_{} {}
  matrix(const matrix &) = delete;
  matrix(matrix &&other) noexcept
      : height_{other.height_}, width_{other.width_}, data_{std::move(
                                                          other.data_)} {}
  explicit matrix(int n) : matrix(n, n) {}
  matrix(int height, int width)
      : height_{height}, width_{width}, data_{std::make_unique<int[]>(
                                            height_ * width_)} {}
  explicit matrix(std::initializer_list<std::initializer_list<int>> l)
      : height_{static_cast<int>(l.size())}, width_{static_cast<int>(
                                                 l.begin()->size())},
        data_{std::make_unique<int[]>(height_ * width_)} {
    if (l.size() > std::numeric_limits<int>::max() ||
        l.begin()->size() > std::numeric_limits<int>::max()) {
      throw std::range_error("Matrix dimensions too large.");
    }
    std::accumulate(l.begin(), l.end(), data_.get(),
                    [](auto *lhs, const auto &rhs) {
                      return std::copy(rhs.begin(), rhs.end(), lhs);
                    });
  }

  matrix &operator=(const matrix &) = delete;
  matrix &operator=(matrix &&other) = delete;

  int height() const noexcept { return height_; }
  int width() const noexcept { return width_; }
  int *operator[](int index) noexcept { return &data_[index * width_]; }
  const int *operator[](int index) const noexcept {
    return &data_[index * width_];
  }

private:
  friend std::ostream &operator<<(std::ostream &os, const matrix &m) noexcept {
    os << '[' << m[0][0];
    for (auto j = 1; j != m.width_; ++j) {
      os << ", " << m[0][j];
    }
    for (auto i = 1; i != m.height_; ++i) {
      os << "]\n[" << m[i][0];
      for (auto j = 1; j != m.width_; ++j) {
        os << ", " << m[i][j];
      }
    }
    os << ']';
    return os;
  }
  const int height_;
  const int width_;
  std::unique_ptr<int[]> data_;
};

class problem {
public:
  problem() = delete;
  problem(const problem &) = delete;
  explicit problem(const matrix *cost_matrix)
      : cost_matrix_(cost_matrix),
        kk_(cost_matrix_->height()), data_{std::make_unique<int[]>(
                                         7 * cost_matrix_->width())} {
    const auto n = kk_.size() >> 5 < 2 ? kk_.size() : (kk_.size() >> 5) + 1;
    std::for_each(kk_.begin(), kk_.end(), [n](auto &v) { v.reserve(n); });
  }

  problem &operator=(const problem &) = delete;

  class solution {
  public:
    solution() = delete;
    solution(const solution &) = delete;
    solution(solution &&other) = default;

    solution &operator=(const solution &) = delete;
    solution &operator=(solution &&other) = default;

    const int *data() const noexcept { return s_.get(); }
    int size() const noexcept { return l_; }
    long value() const noexcept { return c_; }

  private:
    friend class problem;
    solution(std::unique_ptr<int[]> &&x, const matrix *c, const int *v,
             const int *u) noexcept
        : s_{std::move(x)}, l_{c->height()},
          c_{std::accumulate(v, v + c->width(),
                             std::accumulate(u, u + c->height(), 0l))} {
      std::for_each(s_.get(), s_.get() + c->height(), [](auto &n) { --n; });
    }

    friend std::ostream &operator<<(std::ostream &os,
                                    const solution &s) noexcept {
      auto begin = s.data();
      os << '[' << *begin++;
      for (const auto end = begin + s.size() - 1; begin != end; ++begin) {
        os << ", " << *begin;
      }
      os << "]^T";
      return os;
    }

    std::unique_ptr<int[]> s_;
    int l_;
    long c_;
  };

  solution solve() const {
    const auto u = data_.get() + 4 * cost_matrix_->width();
    const auto v = data_.get() + 5 * cost_matrix_->width();
    selpp_cr();
    augmentation(arr(transfer()));
    return solution(std::move(x_), cost_matrix_, v, u);
  }

private:
  void selpp_cr() const {
    const auto v = data_.get() + 5 * cost_matrix_->width();
    const auto y = data_.get() + 6 * cost_matrix_->width();
    const auto end = cost_matrix_->width() >> 5 < 2
                         ? cost_matrix_->width()
                         : cost_matrix_->width() >> 5;

    std::fill_n(v, cost_matrix_->width(), inf);
    x_ = std::make_unique<int[]>(cost_matrix_->height());

    for (auto i = 0; i != cost_matrix_->height(); ++i) {
      auto s = 0l;
      for (auto j = 0; j != end; ++j) {
        kk_[i].push_back(j);
        s += (*cost_matrix_)[i][j];
      }

      auto cr = s / end;

      for (auto j = end; j != cost_matrix_->width(); ++j) {
        if ((*cost_matrix_)[i][j] < cr) {
          auto h = 0, t = -1;
          do {
            t = t >= end - 1 ? 0 : t + 1;
            h = (*cost_matrix_)[i][kk_[i][t]];
          } while (h < cr);
          kk_[i][t] = j;
          s = s - h + (*cost_matrix_)[i][j];
          cr = s / end;
        }
      }

      x_[i] = 0;
      auto diag = (*cost_matrix_)[i][i];

      for (auto t = 0; t != end; ++t) {
        const auto j = kk_[i][t];
        if ((*cost_matrix_)[i][j] < v[j]) {
          v[j] = (*cost_matrix_)[i][j];
          y[j] = i + 1;
        }
        if (j == i) {
          diag = inf;
        }
      }

      if (diag < inf) {
        kk_[i].push_back(i);
        if (diag < v[i]) {
          v[i] = diag;
          y[i] = i + 1;
        }
      }
    }

    for (auto j = cost_matrix_->width(); j > 0; --j) {
      const auto i = y[j - 1] - 1;
      if (x_[i] == 0) {
        x_[i] = j;
      } else {
        x_[i] = -std::abs(x_[i]);
        y[j - 1] = 0;
      }
    }
  }

  int transfer() const noexcept {
    const auto unused = data_.get() + cost_matrix_->width();
    const auto v = data_.get() + 5 * cost_matrix_->width();
    auto l = -1;

    for (auto i = 0; i != cost_matrix_->height(); ++i) {
      if (x_[i] < 0) {
        x_[i] = -x_[i];
      } else if (x_[i] == 0) {
        unused[++l] = i;
      } else {
        auto min = inf;
        const auto j1 = x_[i] - 1;
        for (auto t = 0u; t != kk_[i].size(); ++t) {
          const auto j = kk_[i][t];
          if (j != j1 && (*cost_matrix_)[i][j] - v[j] < min) {
            min = (*cost_matrix_)[i][j] - v[j];
          }
        }
        v[j1] = (*cost_matrix_)[i][j1] - min;
      }
    }

    return l;
  }

  int arr(int l) const noexcept {
    const auto unused = data_.get() + cost_matrix_->width();
    const auto v = data_.get() + 5 * cost_matrix_->width();
    const auto y = data_.get() + 6 * cost_matrix_->width();

    for (auto cnt = 0; cnt != 2; ++cnt) {
      auto h = 0;
      const auto l0 = l;
      l = -1;

      while (h <= l0) {
        const auto i = unused[h++];
        auto v0 = inf, vj = inf;
        auto j0 = -1, j1 = -1;

        for (auto t = 0u; t != kk_[i].size(); ++t) {
          const auto j = kk_[i][t], dj = (*cost_matrix_)[i][j] - v[j];

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

        auto i0 = y[j0];

        if (v0 < vj) {
          v[j0] = v[j0] - vj + v0;
        } else if (i0 > 0) {
          j0 = j1;
          i0 = y[j0];
        }

        x_[i] = j0 + 1;
        y[j0] = i + 1;

        if (i0 > 0) {
          if (v0 < vj) {
            unused[--h] = i0 - 1;
          } else {
            unused[++l] = i0 - 1;
          }
        }
      }
    }
    return l;
  }

  void augmentation(int l) const {
    const auto d = data_.get();
    const auto unused = data_.get() + cost_matrix_->width();
    const auto lab = data_.get() + 2 * cost_matrix_->width();
    const auto todo = data_.get() + 3 * cost_matrix_->width();
    const auto u = data_.get() + 4 * cost_matrix_->width();
    const auto v = data_.get() + 5 * cost_matrix_->width();
    const auto y = data_.get() + 6 * cost_matrix_->width();

    std::vector<bool> ok(cost_matrix_->height());

    do {
      const auto l0 = l;
      for (l = 0; l <= l0; ++l) {
        auto j = -1;
        std::fill_n(d, cost_matrix_->height(), inf);
        std::fill(ok.begin(), ok.end(), false);
        auto min = inf;
        const auto i0 = unused[l];
        auto td1 = -1;

        for (auto t = 0u; t != kk_[i0].size(); ++t) {
          j = kk_[i0][t];
          const auto dj = (*cost_matrix_)[i0][j] - v[j];
          d[j] = dj;
          lab[j] = i0;

          if (dj <= min) {
            if (dj < min) {
              td1 = -1;
              min = dj;
            }
            todo[++td1] = j;
          }
        }

        auto last = cost_matrix_->height();
        auto td2 = last - 1;

        for (auto h = 0; h != td1 + 1; ++h) {
          j = todo[h];
          if (y[j] == 0) {
            goto augment;
          }
          ok[j] = true;
        }

        do {
          const auto j0 = todo[td1--];
          const auto i = y[j0] - 1;
          todo[td2--] = j0;

          const auto hh = (*cost_matrix_)[i][j0] - v[j0] - min;
          for (auto t = 0; t != static_cast<int>(kk_[i].size()); ++t) {
            j = kk_[i][t];
            if (!ok[j]) {
              const auto vj = (*cost_matrix_)[i][j] - v[j] - hh;
              if (vj < d[j]) {
                d[j] = vj;
                lab[j] = i;
                if (vj == min) {
                  if (y[j] == 0) {
                    goto price_update;
                  }
                  todo[++td1] = j;
                  ok[j] = true;
                }
              }
            }
          }

          if (td1 == -1) {
            min = inf - 1;
            last = td2 + 1;
            for (j = 0; j != cost_matrix_->width(); ++j) {
              if (d[j] <= min && !ok[j]) {
                if (d[j] < min) {
                  td1 = -1;
                  min = d[j];
                }
                todo[++td1] = j;
              }
            }
            for (auto h = 0; h != td1 + 1; ++h) {
              j = todo[h];
              if (y[j] == 0) {
                goto price_update;
              }
              ok[j] = true;
            }
          }
        } while (true);
      price_update:
        for (auto k = last; k < cost_matrix_->width(); ++k) {
          const auto j0 = todo[k];
          v[j0] = v[j0] + d[j0] - min;
        }
      augment:
        auto i = -1;
        do {
          i = lab[j];
          y[j] = i + 1;
          const auto k = j + 1;
          j = x_[i] - 1;
          x_[i] = k;
        } while (i != i0);
      }

      for (auto i = 0; i != cost_matrix_->height(); ++i) {
        const auto j = x_[i] - 1;
        u[i] = (*cost_matrix_)[i][j] - v[j];
      }

      l = optcheck();
    } while (l >= 0);
  }

  int optcheck() const noexcept {
    const auto unused = data_.get() + cost_matrix_->width();
    const auto u = data_.get() + 4 * cost_matrix_->width();
    const auto v = data_.get() + 5 * cost_matrix_->width();
    const auto y = data_.get() + 6 * cost_matrix_->width();

    auto l = -1;

    for (auto i = 0; i != cost_matrix_->height(); ++i) {
      auto newfree = false;

      for (auto j = 0; j != cost_matrix_->width(); ++j) {
        if ((*cost_matrix_)[i][j] < u[i] + v[j]) {
          kk_[i].push_back(j);
          newfree = true;
        }
      }

      if (newfree) {
        y[x_[i] - 1] = 0;
        x_[i] = 0;
        unused[++l] = i;
      }
    }

    return l;
  }

  int size() const noexcept { return cost_matrix_->height(); }

  const matrix *const cost_matrix_;
  mutable std::vector<std::vector<int>> kk_;
  const std::unique_ptr<int[]> data_;
  mutable std::unique_ptr<int[]> x_;
  static const int inf = std::numeric_limits<int>::max();
};
} // namespace lapmod
