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

#ifndef LAPMOD_H_
#define LAPMOD_H_

#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace lapmod {
class matrix {
 public:
    matrix() noexcept : height_{}, width_{}, data_{} {}
    matrix(const matrix &other)
        : height_{other.height_}, width_{other.width_},
          data_{std::make_unique<int[]>(height_ * width_)} {
        std::copy_n(other.data_.get(), height_ * width_, data_.get());
    }
    matrix(matrix &&other) = default;
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
        std::accumulate(l.begin(), l.end(),
                        data_.get(), [](auto *lhs, const auto &rhs) noexcept {
                            return std::copy(rhs.begin(), rhs.end(), lhs);
                        });
    }

    matrix &operator=(matrix other) noexcept {
        using std::swap;
        swap(*this, other);
        return *this;
    }

    int height() const noexcept { return height_; }
    int width() const noexcept { return width_; }
    int *operator[](int index) noexcept { return &data_[index * width_]; }
    const int *operator[](int index) const noexcept {
        return &data_[index * width_];
    }

 private:
    friend std::ostream &operator<<(std::ostream &os,
                                    const matrix &m) noexcept {
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

    friend void swap(matrix &l, matrix &r) noexcept {
        using std::swap;
        swap(l.data_, r.data_);
        swap(l.height_, r.height_);
        swap(l.width_, r.width_);
    }

    int height_;
    int width_;
    std::unique_ptr<int[]> data_;
};

class problem {
    static constexpr int k_d_idx() noexcept { return 0; }
    static constexpr int k_unused_idx() noexcept { return 1; }
    static constexpr int k_lab_idx() noexcept { return 2; }
    static constexpr int k_todo_idx() noexcept { return 3; }
    static constexpr int k_u_idx() noexcept { return 4; }
    static constexpr int k_v_idx() noexcept { return 5; }
    static constexpr int k_y_idx() noexcept { return 6; }
    static constexpr int k_ok_idx() noexcept { return 7; }
    static constexpr int k_field_count() noexcept { return 8; }

    static constexpr int k_max_core() noexcept { return 8192; }

    auto get_field(int idx) const noexcept {
        return data_.get() + idx * cost_matrix_->width();
    }

 public:
    problem(const problem &other)
        : cost_matrix_{other.cost_matrix_}, kk_{other.kk_},
          data_{
              std::make_unique<int[]>(k_field_count() * cost_matrix_->width())},
          x_{} {
        std::copy_n(other.data_.get(), k_field_count() * cost_matrix_->width(),
                    data_.get());
        if (other.x_) {
            x_ = std::make_unique<int[]>(cost_matrix_->height());
            std::copy_n(other.x_.get(), cost_matrix_->height(), x_.get());
        }
    }

    problem(problem &&) = default;

    explicit problem(const matrix *cost_matrix)
        : cost_matrix_(cost_matrix),
          kk_(cost_matrix_->height(),
              std::vector<int>(
                  std::max(std::min(cost_matrix_->width(), k_max_core()),
                           cost_matrix_->width() >> 2u))),
          data_{
              std::make_unique<int[]>(k_field_count() * cost_matrix_->width())},
          x_{} {}

    problem &operator=(problem other) noexcept {
        using std::swap;
        swap(*this, other);
        return *this;
    }

    int size() const noexcept { return cost_matrix_->height(); }

    class solution {
     public:
        solution(const solution &other)
            : s_{std::make_unique<int[]>(other.l_)}, l_{other.l_},
              c_{other.c_} {
            std::copy_n(other.s_.get(), l_, s_.get());
        }
        solution(solution &&other) = default;

        solution &operator=(solution other) noexcept {
            using std::swap;
            swap(*this, other);
            return *this;
        }

        const int *data() const noexcept { return s_.get(); }
        int size() const noexcept { return l_; }
        long value() const noexcept { return c_; }

     private:
        solution() = default;
        friend class problem;

        friend void swap(solution &l, solution &r) noexcept {
            using std::swap;
            swap(l.s_, r.s_);
            swap(l.l_, r.l_);
            swap(l.c_, r.c_);
        }

        template <typename T>
        auto move_helper(T &&x, const matrix *m) noexcept {
            std::for_each(x.get(),
                          x.get() + m->height(), [](auto &n) noexcept { --n; });
            return std::forward<T>(x);
        }

        solution(std::unique_ptr<int[]> &&x, const matrix *c, const int *v,
                 const int *u) noexcept
            : s_{move_helper(std::move(x), c)}, l_{c->height()},
              c_{std::accumulate(v, v + c->width(),
                                 std::accumulate(u, u + c->height(), 0l))} {}

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

    solution solve(bool approximate = false) const {
        const auto u = get_field(k_u_idx());
        const auto v = get_field(k_v_idx());
        selpp_cr();
        augmentation(arr(transfer<decltype(get_field(k_unused_idx()))>()),
                     approximate);
        return solution(std::move(x_), cost_matrix_, v, u);
    }

 private:
    problem() = default;

    void selpp_cr() const {
        const auto v = get_field(k_v_idx());
        const auto y = get_field(k_y_idx());
        const auto end = kk_.empty() ? 0 : static_cast<int>(kk_.front().size());

        std::fill_n(v, cost_matrix_->width(), inf());
        x_ = std::make_unique<int[]>(cost_matrix_->height());

        for (auto i = 0; i != cost_matrix_->height(); ++i) {
            auto s = std::accumulate((*cost_matrix_)[i],
                                     &(*cost_matrix_)[i][end], 0l);
            std::iota(kk_[i].begin(), kk_[i].end(), 0);
            auto cr = static_cast<decltype(s)>(
                static_cast<std::make_unsigned_t<decltype(s)>>(s) /
                static_cast<std::make_unsigned_t<decltype(end)>>(end));

            for (auto j = end; j != cost_matrix_->width(); ++j) {
                if ((*cost_matrix_)[i][j] < cr) {
                    auto h = 0, t = -1;
                    do {
                        t = t >= end - 1 ? 0 : t + 1;
                        h = (*cost_matrix_)[i][kk_[i][t]];
                    } while (h < cr);
                    kk_[i][t] = j;
                    s = s - h + (*cost_matrix_)[i][j];
                    cr = static_cast<decltype(s)>(
                        static_cast<std::make_unsigned_t<decltype(s)>>(s) /
                        static_cast<std::make_unsigned_t<decltype(end)>>(end));
                }
            }

            x_[i] = 0;
            auto diag = (*cost_matrix_)[i][i];

            for (const auto &j : kk_[i]) {
                if ((*cost_matrix_)[i][j] < v[j]) {
                    v[j] = (*cost_matrix_)[i][j];
                    y[j] = i + 1;
                }
                if (j == i) {
                    diag = inf();
                }
            }

            if (diag < inf()) {
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

    template <class BidiIt> BidiIt transfer() const noexcept {
        BidiIt unused = get_field(k_unused_idx());
        const auto v = get_field(k_v_idx());

        for (auto i = 0; i != cost_matrix_->height(); ++i) {
            if (x_[i] < 0) {
                x_[i] = -x_[i];
            } else if (x_[i] == 0) {
                *unused++ = i;
            } else {
                auto min = inf();
                const auto j1 = x_[i] - 1;
                for (const auto &j : kk_[i]) {
                    if (j != j1 && (*cost_matrix_)[i][j] - v[j] < min) {
                        min = (*cost_matrix_)[i][j] - v[j];
                    }
                }
                v[j1] = (*cost_matrix_)[i][j1] - min;
            }
        }

        return unused;
    }

    template <class BidiIt> BidiIt arr(BidiIt uit) const noexcept {
        const BidiIt unused = get_field(k_unused_idx());
        const auto v = get_field(k_v_idx());
        const auto y = get_field(k_y_idx());

        for (auto cnt = 0; cnt != 2; ++cnt) {
            BidiIt unused_begin = unused;
            const BidiIt unused_end = uit;
            uit = unused_begin;

            while (unused_begin != unused_end) {
                const auto i = *unused_begin++;
                auto v0 = inf(), vj = inf();
                auto j0 = -1, j1 = -1;

                for (const auto &j : kk_[i]) {
                    const auto dj = (*cost_matrix_)[i][j] - v[j];

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
                        *--unused_begin = i0 - 1;
                    } else {
                        *uit++ = i0 - 1;
                    }
                }
            }
        }
        return uit;
    }

    template <class BidiIt>
    void augmentation(BidiIt unused_end, bool approximate) const {
        const auto d = get_field(k_d_idx());
        const BidiIt unused = get_field(k_unused_idx());
        const auto lab = get_field(k_lab_idx());
        const auto todo = get_field(k_todo_idx());
        const auto u = get_field(k_u_idx());
        const auto v = get_field(k_v_idx());
        const auto y = get_field(k_y_idx());
        const auto ok = get_field(k_ok_idx());

        do {
            for (auto uit = unused; uit != unused_end; ++uit) {
                auto j = -1;
                std::fill_n(d, cost_matrix_->height(), inf());
                std::fill_n(ok, cost_matrix_->height(), false);
                auto min = inf();
                const auto i0 = *uit;
                auto td1 = -1;

                for (const auto &idx : kk_[i0]) {
                    j = idx;
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

                for (auto it = todo; it != &todo[td1 + 1]; ++it) {
                    j = *it;
                    if (y[j] == 0) {
                        goto augment;
                    }
                    ok[j] = true;
                }

                for (;;) {
                    const auto j0 = todo[td1--];
                    const auto i = y[j0] - 1;
                    todo[td2--] = j0;

                    const auto hh = (*cost_matrix_)[i][j0] - v[j0] - min;
                    for (const auto &idx : kk_[i]) {
                        j = idx;
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
                        min = inf() - 1;
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
                        for (auto it = todo; it != &todo[td1 + 1]; ++it) {
                            j = *it;
                            if (y[j] == 0) {
                                goto price_update;
                            }
                            ok[j] = true;
                        }
                    }
                }
            price_update:
                for (auto it = &todo[last]; it != &todo[cost_matrix_->width()];
                     ++it) {
                    const auto j0 = *it;
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

            if (approximate) {
                break;
            }

            unused_end = optcheck<BidiIt>();
        } while (unused_end != unused);
    }

    template <class BidiIt> BidiIt optcheck() const {
        BidiIt unused = get_field(k_unused_idx());
        const auto u = get_field(k_u_idx());
        const auto v = get_field(k_v_idx());
        const auto y = get_field(k_y_idx());

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
                *unused++ = i;
            }
        }

        return unused;
    }

    friend void swap(problem &l, problem &r) noexcept {
        using std::swap;
        swap(l.cost_matrix_, r.cost_matrix_);
        swap(l.kk_, r.kk_);
        swap(l.data_, r.data_);
        swap(l.x_, r.x_);
    }

    const matrix *cost_matrix_;
    mutable std::vector<std::vector<int>> kk_;
    std::unique_ptr<int[]> data_;
    mutable std::unique_ptr<int[]> x_;

    static constexpr int inf() noexcept {
        return std::numeric_limits<int>::max();
    }
};
} // namespace lapmod

#endif // LAPMOD_H_

