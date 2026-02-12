#pragma once

inline void lpp::detail::generate_combination(const int start, const int n, const int k, std::vector<int>& current,
                                              std::vector<std::vector<int>>& res) {
    if (current.size() == k) {
        res.push_back(current);
        return;
    }
    for (int i = start; i < n; i++) {
        current.push_back(i);
        generate_combination(i + 1, n, k, current, res);
        current.pop_back();
    }
}
